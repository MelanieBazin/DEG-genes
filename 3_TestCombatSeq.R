####
# Validation of ComBat-seq on Paramecium data
# -> Comparison of the data with or without correction
# -> On the data sequenced twice (IC7 and EZL1 time courses)
####

options(stringsAsFactors = FALSE)
library(sva)
library(DESeq2)
library(pheatmap)
library(MASS)
library(dendextend)

set.seed(10111)

source("0_Cluster.R") # Groups parameters
source("0_Functions.R") # Library of homemade function

analyseName = "Test_Combatseq"
analyseName = paste0(Sys.Date(),"_", analyseName)
path_dir = paste0("./Analyse/",analyseName,"/")

condition = names(rnai_list["HiSeqvsNextSeq"])

for (correction in c("corrected","uncorrected")){
  print(paste("Test of conditions",correction, "by ComBatSeq"))
  
  path = paste0(path_dir,"/",condition,"_", correction,"/")
  dir.create(path,recursive=T,showWarnings=F)
  
  # Open files with or without ComBat-seq correction
  countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_",correction,".tab"), sep="\t",row.names=1,header =  T)
  
  # Generation the appropriate infodata
  infodata = CreatInfoData(countdata, conditions = condition , rnai_list, cluster)
  
  # Limit the data to the sample sequenced twice
  time = unique(infodata$Timing[which(infodata$Seq_method == "NextSeq")])
  countdata = countdata[,is.element(infodata$Timing, time)]
  infodata = infodata[is.element(infodata$Timing, time),]
  
  cluster2 = list(
    ICL7 = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
    ICL7bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
    EZL1 = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
    EZL1bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)))
  
  #### DESeq2 analysis ####
  # Formatting the data for DESeq2 analysis
  countdata = as.matrix(countdata)
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ Samples)
  
  # DESeq2 analysis
  print(paste(condition , "-----> DESeq2 analysis"))
  deseq = DESeq(deseq)
  
  # Extracting the vst normalized data
  data_tab = assay(vst(deseq, blind = T))
  
  # Saving the vst normalized data and the correspondinf infodata
  write.table(data_tab,paste0(path,condition,"_expression_table_",correction,".tab"), sep="\t",row.names=T,quote=F)
  write.table(infodata,paste0(path,condition,"_infodata_",correction,".tab"), sep="\t",row.names=T,quote=F)
  
  #### Re-openining of the files ####
  # data_tab = read.table(paste0(path,condition,"_expression_table_",correction,".tab"), sep="\t", header = T)
  # infodata = read.table(paste0(path,condition,"_infodata_",correction,".tab"), sep="\t", header = T)
  
  #### Data visualization  ####
  for (color_type in c("methods","replicates")){
    # PCA plot
    PCA_ggplot_generator(data_tab,
                         infodata,
                         police_seize = 4,
                         point_seize = 3,
                         save_path = path,
                         main = paste0("PCA ", condition," (DESeq2)"),
                         sortie = "pdf",
                         rename = F,
                         color_type = color_type,
                         w = 4,
                         h = 4)
    
    # Pearson correlation matrix & hierarchical clustering
    matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
    res = hclust(matDist)
    res = as.dendrogram(res)
    
    color = Color_type(data_tab, infodata, type = color_type)
    labels_colors(res)= as.character(color)[order.dendrogram(res)]
    
    pdf(paste0(path, color_type, "/", condition ,"_hclust_pearson_vst.pdf"),  width = 12, height = 2.5)
    plot(res, main = "pearson_vst")
    dev.off()
  }
  
}

#### Print R status ####
sink(paste0("./Analyse/",analyseName,"/sessionInfo.txt"))
print(sessionInfo())
sink()
