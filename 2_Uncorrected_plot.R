###################
#### Visualization of the data befor correction
#### -> On the data used for DESeq2 analysis
###################

options(stringsAsFactors = FALSE)
library(sva)
library(DESeq2)
library(pheatmap)
library(MASS)
library(dendextend)

set.seed(10111)

source("0_Cluster.R") # Groups parameters
source("0_Functions.R") # Library of homemade function

condition = names(rnai_list)[2]

analyseName = paste0(Sys.Date(),"_Uncorrected_", condition)

path = paste0("./Analyse/",analyseName,"/")
dir.create(path,recursive=T,showWarnings=F)


# Open files without ComBat-seq correction
countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_uncorrected.tab"), sep="\t",row.names=1,header =  T)


# Generation the appropriate infodata
infodata = CreatInfoData(countdata, conditions = condition , rnai_list, cluster)

##### DESeq2 analysis ####
# Formatting the data for DESeq2 analysis
countdata = as.matrix(countdata)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Condition)

# Identification of technical replicates
deseq = collapseReplicates(deseq,
                           groupby = deseq$Samples,
                           run     = deseq$Names)

# Saving of the infodata of the collapsed table
infodata_collapse = as.data.frame(colData(deseq))
write.table(infodata_collapse,paste0(path,condition ,"_infodata_collapse.tab"), sep="\t",row.names=T,quote=F)

# DESeq2 analysis
deseq = DESeq(deseq)

# Saving the vst(Variance stabilizing transformation) normalized data
data_tab = assay(vst(deseq, blind = T))
write.table(data_tab,paste0(path,condition ,"_expression_table_vst.tab"), sep="\t",row.names=T,quote=F)

#### Plotting data ####
for (color_type in c("methods","replicates")){
  save_path = paste0(path,color_type,"/")
  dir.create(save_path,recursive=T,showWarnings=F)
  
  # PCA analysis
  PCA_ggplot_generator(data_tab,
                       infodata_collapse,
                       color_type = color_type,
                       collapse = T,
                       save_path = path,
                       main = paste0("PCA ", condition," (DESeq2)"),
                       sortie = "pdf",
                       police_seize = 4,
                       point_seize = 2.5,
                       rename = T,
                       w = 6.5,
                       h = 6)
  
  # Pearson correlation matrix & hierarchical clustering
  matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
  res = hclust(matDist)
  res = as.dendrogram(res)
  
  color = Color_type(data_tab, infodata_collapse, type = color_type)
  labels_colors(res)= as.character(color)[order.dendrogram(res)]
  
  pdf(paste0(save_path, condition ,"_hclust_pearson_vst.pdf"),  width = 12, height = 2.5)
  plot(res, main = "pearson_vst")
  dev.off()
}

sink(paste0(path,"/sessionInfo.txt"))
print(sessionInfo())
sink()