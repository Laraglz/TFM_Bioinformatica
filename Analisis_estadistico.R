
# Estudio de la expresión génica inducida por el ejercicio físico mediante RNA-Seq
# Lara González Pérez
# 2025
# ------------------------------------------------------------------------------

setwd("/Users/tu_usuario/Documents/mi_proyecto")
metadata <- read.csv(file = "metadata.csv")
View(metadata)

# 1. Importar datos y crear matriz de recuentos
# Instalar paquetes necesarios para el análisis
if(!require("BiocManager")) {
  install.packages("BiocManager")
  BiocManager::install("edgeR")
}else{
  library(edgeR)
}

ruta_archivos <- "/Users/tu_usuario/Documents/mi_proyecto"
archivos <- list.files(ruta_archivos, pattern = "\\.tsv$", full.names = FALSE)

# Asignar nombres de muestra (sustituir Run por Sample.Name)
runs <- sub("_counts.tsv", "", basename(archivos))
run_to_sample <- setNames(metadata$Sample.Name, metadata$Run)

# Renombrar y leer los archivos
nombres_muestras <- run_to_sample[runs]
conteos <- readDGE(archivos, columns = c(1, 2), labels = nombres_muestras)

# Ver la matriz de conteo
head(conteos$counts)
nrow(conteos) # número de genes del experimento

# Extraer la matriz de conteos y guardar como CSV
matriz_conteos <- conteos$counts
write.csv(matriz_conteos, file = "matriz_conteos.csv", row.names = TRUE)

# Instalar la librería de genes del genoma humano
if(!require("org.Hs.eg.db")) {
  install.packages("org.Hs.eg.dbr")
}else{
  library(org.Hs.eg.db)
}

# Asociar genes a su nombre
ann <- select(org.Hs.eg.db, keys = row.names(conteos),
              columns = c("ENTREZID", "SYMBOL", "GENENAME"),
              keytype = "ENSEMBL")
head(ann)

# Crear la variable grupo para usarla después teniendo en cuenta que los nombres
# de muestra coincidan con la muestra real
grupo <- metadata$training[match(nombres_muestras, metadata$Sample.Name)]
grupo <- factor(grupo, levels = c("Untrained", "Trained"))


# 2. Crear el objeto DGElist con la matriz de recuentos, los nombres de los genes y
# los grupos y nombres de muestras
DGE_objeto <- DGEList(conteos)
DGE_objeto$genes <- ann
DGE_objeto$samples$group <- grupo
head(DGE_objeto)


# 3. Eliminación de ruido, es decir, genes de expresión muy baja o nula
# la propia función ya hace el cambio a cpm
keep <- filterByExpr(DGE_objeto, group = grupo)
summary(keep)

# No conservar los tamaños de biblioteca originales ni los genes: FALSE
DGE_filtrado <- DGE_objeto[keep, keep.lib.sizes = FALSE]
head(DGE_filtrado)
# Se conservan 14183 genes después del filtrado

# Ordenar el objeto ann para asignar correctamente los nombres de los genes 
ann_ordenado <- ann[match(rownames(DGE_filtrado), ann$ENSEMBL), ]
DGE_filtrado$genes <- ann_ordenado


# 4. Normalizar los conteos
DGE_filtrado <- calcNormFactors(DGE_filtrado)
DGE_filtrado$samples


# 5. MDS plot
colores <- c("Untrained" = "orange", "Trained" = "purple")
par(pty = "s", mar = c(5, 4, 4, 8))
MDS <- plotMDS(DGE_filtrado, col = colores[grupo], pch = 16, top = 1000,
               main = "MDS plot", cex = 1.5)
legend("right",inset = c(-0.8, 0), legend = levels(grupo),
       col = colores[levels(grupo)], pch = 16, xpd = TRUE, bty = "n")
par(pty = "m", mar = c(5, 4, 4, 2)) # restaurar configuración

# Añadir etiquetas a los puntos para ver qué muestra es cada uno
text(MDS, labels = colnames(DGE_filtrado), cex = 0.6, pos = 3)


# 6. PCA plot
logCPM <- cpm(DGE_filtrado, log=TRUE, prior.count=1)
pca <- prcomp(t(logCPM))

# Calcular los porcentajes de varianza explicada
varianza <- pca$sdev^2 / sum(pca$sdev^2) * 100
pc1_var <- round(varianza[1], 1)
pc2_var <- round(varianza[2], 1)

par(pty = "s", mar = c(5, 4, 4, 8))
plot(pca$x[,1:2], col=colores[grupo], pch=16, main="PCA plot", cex = 1.5,
     xlab = paste0("PC1 (", pc1_var, "%)"),
     ylab = paste0("PC2 (", pc2_var, "%)"))
legend("right", inset = c(-0.8, 0), legend = levels(grupo), col = colores[levels(grupo)], pch=16,
       xpd = TRUE, bty = "n")

text(pca$x[,1:2], labels = colnames(DGE_filtrado), cex = 0.6, pos = 3)


# 7. Matriz de diseño
if(!require("statmod")) {
  install.packages("statmod")
}else{
  library(statmod)
}

# Establecer la matriz de diseño
participante <- factor(metadata$Participant)
training <- factor(metadata$training, levels = c("Untrained", "Trained"))

design <- model.matrix(~ participante + training)
design

# 8. Estimar la dispersión
DGE_filtrado <- estimateDisp(DGE_filtrado, design = design, robust = TRUE)

# Gráfico de dispersión
par(pty = "s", mar = c(5, 4, 4, 8), cex = 0.8)
plotBCV(DGE_filtrado)


# Corregir la dispersión y hacer otro gráfico
fit <- glmQLFit(DGE_filtrado, design = design, robust = TRUE)
head(fit$counts)
head(fit$fitted.values)
par(pty = "s", mar = c(5, 4, 4, 8), cex = 0.8)
plotQLDisp(fit)


# 9. Test de expresión diferencial (QLF)
# Definir el contraste
contraste <- makeContrasts(trainingTrained, levels = design)
contraste

res <- glmQLFTest(fit, contrast = contraste)
head(res$table)
dim(res$table)

# 10. Corrección de datos positivos, clasifica los genes de más a menos significativos
res_corrected <- topTags(res, n = Inf) # infinito para que reporte todos los genes
dim(res_corrected)
head(res_corrected)

# 11. Expresión diferencial de genes
# Volcano plot inicial
library(ggplot2)
head(res_corrected$table)
data <- res_corrected$table # facilitar la escritura
ggplot(data, aes(x=logFC, y=-log10(FDR))) + geom_point()

# Añadir una columna nueva a data y modificarla para que aparezca up y down segun 
# los valores de logFC y FDR
data$DE <- "NO"
data$DE[data$logFC > 1 & data$FDR < 0.05] <- "UP"
data$DE[data$logFC < -1 & data$FDR < 0.05] <- "DOWN"
head(data)

# Volcano plot final
if(!require("ggrepel")) {
  install.packages("ggrepel")
}else{
  library(ggrepel)
}

ggplot(data, aes(x = logFC, y = -log10(FDR), col = DE)) +
  geom_point(size = 1) + 
  geom_text_repel(data = subset(data, FDR < 0.05 & abs(logFC) > 1.5), 
                  aes(label = SYMBOL),
                  size = 3,
                  max.overlaps = 20) +
  theme_classic()

# 12. Ontología génica
if(!require("GO.db")) {
  install.packages("GO.db")
}else{
  library(GO.db)
}

# Extraer genes para GO con FDR <0.05 y logFC grande
genes_go <- res_corrected$table
genes_go$DEG <- genes_go$FDR < 0.05 & (genes_go$logFC > 1 | genes_go$logFC < -1)
table(genes_go$DEG)

# Filtrar solo con ENTREZ válidos
genes_go$ENTREZID <- DGE_filtrado$genes$ENTREZID[match(rownames(genes_go), DGE_filtrado$genes$ENSEMBL)]
genes_go <- genes_go[!is.na(genes_go$ENTREZID), ]
head(genes_go)

# Usar goana con esos
go <- goana(de = genes_go$ENTREZID[genes_go$DEG], universe = genes_go$ENTREZID, species = "Hs")

topGO(go, "BP")
topGO(go, "CC")
topGO(go, "MF")

# Tabla con logCPM
logCPM <- cpm(DGE_filtrado, log = TRUE)
head(logCPM)

# Cambiar el nombre de las filas por el nombre de los genes que aparece en el objeto
rownames(logCPM) <- DGE_filtrado$genes$SYMBOL
head(logCPM)

# Cambiar el nombre de las columnas por el tipo de muestra
# Obtener nombres actuales de las muestras (columnas de logCPM)
nombres_muestras <- colnames(DGE_filtrado)

# Obtener participantes en el orden correcto
participantes <- metadata$Participant[match(nombres_muestras, metadata$Sample.Name)]

# Generar los nombres nuevos
colnames(logCPM) <- paste(DGE_filtrado$samples$group, participantes, sep = "-")

head(logCPM)
dim(logCPM)

genes_deg_ensembl <- rownames(genes_go)[genes_go$DEG]  # Solo ENSEMBL IDs de DEGs
symbols <- DGE_filtrado$genes$SYMBOL[match(genes_deg_ensembl, DGE_filtrado$genes$ENSEMBL)]

#Tabla final
logCPM_deg <- logCPM[rownames(logCPM) %in% symbols, ]
logCPM_deg 

# 13. Heatmap
library(pheatmap)

pheatmap(logCPM_deg, 
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "ward.D2",  
         fontsize = 10,                 
         fontsize_row = 10,
         fontsize_col = 10,
         display_numbers = TRUE,
         cellwidth = 25,        
         cellheight = 25)

# 14. Barplots con los términos GO
# Biological Process
bp <- topGO(go, "BP", number = 10)
bp$log10p <- -log10(bp$P.DE)

p1 <- ggplot(bp, aes(x = reorder(Term, log10p), y = log10p)) +
  geom_col(fill = "#6A994E", width = 0.6) +
  coord_flip() +
  labs(title = "Biological Process",
       x = "",
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank())

print(p1)

# Cellular Component
cc <- topGO(go, "CC", number = 10)
cc$log10p <- -log10(cc$P.DE)

p2 <- ggplot(cc, aes(x = reorder(Term, log10p), y = log10p)) +
  geom_col(fill = "#F18F01", width = 0.6) +
  coord_flip() +
  labs(title = "Cellular Component",
       x = "",
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank())

print(p2)

# Molecular Function
mf <- topGO(go, "MF", number = 10)
mf$log10p <- -log10(mf$P.DE)

p3 <- ggplot(mf, aes(x = reorder(Term, log10p), y = log10p)) +
  geom_col(fill = "#CC78BC", width = 0.6) +
  coord_flip() +
  labs(title = "Molecular Function",
       x = "",
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(panel.grid = element_blank())

print(p3)










