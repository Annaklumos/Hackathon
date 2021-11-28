library(DESeq2)
library(factoextra)
library(FactoMineR)

# Chargement des donnees
files <- (Sys.glob("gene_strand_0/*.counts"))
list_files <- lapply(files, function(x) read.table(x, header = T)) 


nb_ech = length(list_files) 
nb_gene = nrow(list_files[[1]])

id_gene <- list_files[[1]]["Geneid"]

df <- rep(0, nb_ech)

for (i in 1:nb_ech) {
   df[i] <- data.frame(list_files[[i]][,7])
}


df <- cbind(id_gene, df)
colnames(df) <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
df_t <- t(as.matrix(df[,-1]))


# Specifier les types 
Type = data.frame("Mutant","Mutant","WT","WT","WT","WT","WT","Mutant")
Type = t(Type)
row.names(Type) <- c(2, 3, 4, 5, 6, 7, 8, 9)
colnames(Type) <- "Type"

# PCA
pdf("PCA.pdf")
fviz_pca_ind(PCA(df_t,graph=F),col.ind = Type)
dev.off()


Des <- DESeqDataSetFromMatrix(countData = df[-1],colData = Type, design = ~Type)
# counts(Des, normalized = TRUE)
analysis = DESeq(Des)
res = results(analysis)
write.table(res, "Deseq2_results_table.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#Normaliser 
vst = getVarianceStabilizedData(analysis)
vst_t <- t(vst)
#PCA(vst_t)


res_mat <- as.matrix(res)
res_mat <- cbind(id_gene, res_mat)

# Lignes dont la pvalue ajustee est < 0.05
ind_padj <- res_mat[,"padj"]
ind_padj <- which(ind_padj<0.05)


# Volcano plot
logp <- -log10(res_mat[,"padj"])
pdf("VolcanoPlot.pdf")
plot(logp~res_mat[,"log2FoldChange"])
dev.off()
