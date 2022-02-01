options(stringsAsFactors = FALSE)
library(sva)
library(DESeq2)
library(pheatmap)
library(MASS)
library(dendextend)

set.seed(10111)

# Récupérer les paramètre de clustring
source("0_Cluster.R")

# Récupérer les fonction necessaire au représentaion graphique et la mise en forme des données
source("3_Visualisation_des_donnees_fonction.R")

condition = names(rnai_list)[3]

analyseName = paste0("Clustering_groupe")
analyseName = paste0(Sys.Date(),"_", analyseName)

path_dir = paste0("./Analyse/",analyseName,"/")
dir.create(path_dir,recursive=T,showWarnings=F)

path = paste0(path_dir,condition ,"/")
dir.create(path,recursive=T,showWarnings=F)
dir.create(paste0(path_dir,condition ,"/Visualisation/Cluster/"),recursive=T,showWarnings=F)

for (condition in names(rnai_list) ){
#### Sans condition de groupe ####
# Ouverture des fichiers avec correction de l'effet Batch sans la condition de groupe
countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_pour_DESeq_v2.tab"), sep="\t",row.names=1,header =  T)

# Ouverture des fichiers avec correction de l'effet Batch avec la condition de groupe
# countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_pour_DESeq_v2.tab"), sep="\t",row.names=1,header =  T)

##### Normalisation des données par rpkm #####
seqlength=read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.transcript.fa.seqlength",h=T)
rownames(seqlength)=sub("PTET.51.1.T","PTET.51.1.G",seqlength$ID)

rpkm = matrix(data = NA,nrow = nrow(countdata),ncol = ncol(countdata))
colnames(rpkm) = colnames(countdata)
rownames(rpkm)=rownames(countdata)

for (i in 1:(ncol(countdata))){
  mapped_reads=sum(countdata[,i])
  rpkm[,i] = (countdata[,i] *1e3) / (seqlength[rownames(countdata),]$LENGTH * (mapped_reads/1e6) )
  
}
write.table(rpkm,paste0(path,condition ,"_expression_table_RPKM.tab"), sep="\t",row.names=T,quote=F)
rpkm = read.table(paste0(path,condition ,"_expression_table_RPKM.tab"), h = T, sep="\t")
# ACP sur les données normalisées
infodata = CreatInfoData3(rpkm, conditions = condition , rnai_list, cluster)

## Couleur par cluster
color = c()
for (r in colnames(rpkm)){
  clust = infodata[r, "Cluster"]
  if (clust == "VEG"){
    color = c(color, clust_color["veg"])
  } else if(clust == "EARLY" ){
    color = c(color, clust_color["early"])
  } else if(clust == "INTER" ){
    color = c(color, clust_color["inter"])
  } else if(clust == "LATE" ){
    color = c(color, clust_color["late"])
  }
}

PCA_plot_generator(rpkm,
                   colors = color,
                   save_path = paste0(path,"Visualisation/ACP_ssGRP/RPKM/"),
                   main = paste0("ACP ", condition," (rpkm)"),
                   sortie = "png")

matDist = as.dist(1-cor(rpkm, method="pearson"))
res = hclust(matDist)
res = as.dendrogram(res)
labels_colors(res)= as.character(color)[order.dendrogram(res)]

png(paste0(path,"/Visualisation/Cluster/", condition,"_Cluster_pearson_rpkm.png"),  width = 800, height = 600)
plot(res, main = "pearson_rpkm")
dev.off()

PCA_plot_generator(log(rpkm+1),
                   colors = color,
                   save_path = paste0(path,"Visualisation/ACP_ssGRP/logRPKM/"),
                   main = paste0("ACP ", condition," (logrpkm)"),
                   sortie = "png")

matDist = as.dist(1-cor(log2(rpkm+1), method="pearson"))
res = hclust(matDist)
res = as.dendrogram(res)
labels_colors(res)= as.character(color)[order.dendrogram(res)]

png(paste0(path,"/Visualisation/Cluster/", condition,"_Cluster_pearson_logrpkm.png"),  width = 800, height = 600)
plot(res, main = "pearson_logrpkm")
dev.off()


##### Analyse DESeq2 ####
# Création du tableau avec les info des colonnes
infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)

# Créataion de l'objet DESeq2
countdata = as.matrix(countdata)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Condition)

# Definition des réplicats techniques
deseq = collapseReplicates(deseq,
                           groupby = deseq$Samples,
                           run     = deseq$Names)

color = c()
for (r in colnames(deseq)){
  clust = infodata[r, "Cluster"]
  if (clust == "VEG"){
    color = c(color, clust_color["veg"])
  } else if(clust == "EARLY" ){
    color = c(color, clust_color["early"])
  } else if(clust == "INTER" ){
    color = c(color, clust_color["inter"])
  } else if(clust == "LATE" ){
    color = c(color, clust_color["late"])
  }
}



# Analyse DESeq2
deseq = DESeq(deseq)

### Recupération des donnée normalisée ####
# Normalisation DESeq2 
norma = as.data.frame(counts(deseq,normalized=TRUE))
write.table(norma,paste0(path,condition ,"_expression_table_DESeq2.tab"), sep="\t",row.names=T,quote=F)

PCA_plot_generator(norma,
                   colors = color,
                   save_path = paste0(path,"Visualisation/ACP_ssGRP/DESeq2/"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")

matDist = as.dist(1-cor(log2(norma+1), method="pearson"))
res = hclust(matDist)
res = as.dendrogram(res)
labels_colors(res)= as.character(color)[order.dendrogram(res)]

png(paste0(path,"/Visualisation/Cluster/", condition,"_Cluster_pearson_DESeq2.png"),  width = 800, height = 600)
plot(res, main = "pearson_DESeq2")
dev.off()


# Variance stabilizing transformation
vsd = assay(vst(deseq, blind = T))
write.table(vsd,paste0(path,condition ,"_expression_table_vst.tab"), sep="\t",row.names=T,quote=F)

PCA_plot_generator(vsd,
                   colors = color,
                   save_path = paste0(path,"Visualisation/ACP_ssGRP/vst/"),
                   main = paste0("ACP ", condition," (vst)"),
                   sortie = "png")

matDist = as.dist(1-cor(log2(vsd+1), method="pearson"))
res = hclust(matDist)
res = as.dendrogram(res)
labels_colors(res)= as.character(color)[order.dendrogram(res)]

png(paste0(path,"/Visualisation/Cluster/", condition,"_Cluster_pearson_vst.png"),  width = 800, height = 600)
plot(res, main = "pearson_vst")
dev.off()

# Regularized log transformation
rlg = assay(rlog(deseq, blind = T))
write.table(rlg,paste0(path,condition ,"_expression_table_rlog.tab"), sep="\t",row.names=T,quote=F)

PCA_plot_generator(rlg,
                   colors = color,
                   save_path = paste0(path,"Visualisation/ACP_ssGRP/rlog/"),
                   main = paste0("ACP ", condition," (rlog)"),
                   sortie = "png")

matDist = as.dist(1-cor(log2(rlg+1), method="pearson"))
res = hclust(matDist)
res = as.dendrogram(res)
labels_colors(res)= as.character(color)[order.dendrogram(res)]

png(paste0(path,"/Visualisation/Cluster/", condition,"_Cluster_pearson_rlog.png"),  width = 800, height = 600)
plot(res, main = "pearson_rlog")
dev.off()
}