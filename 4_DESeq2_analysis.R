###################
#### DESeq2 data analysis
#### -> Comparison of the data with or without correction
#### -> On the data sequenced twice (IC7 and EZL1 time courses)
###################
options(stringsAsFactors = FALSE)
library(sva)
library(DESeq2)
library(pheatmap)
library(MASS)

set.seed(10111)

source("0_Cluster.R") # Groups parameters
source("0_Functions.R") # Library of homemade function

##### To reopen files ###
# path_dir = "./Analyse/2022-02-21_Analyse_DESeq2_FC-1.5_pval-0.05/"
# condition = names(rnai_list)[2]
# data_tab = read.table(paste0(path_dir,condition ,"/analyseDE_expression_table_vst.tab"), row.names = 1, sep="\t", header = T)
# infodata = read.table(paste0(path_dir,condition ,"/analyseDE_infodata_collapse.tab"), row.names = 1, sep="\t", header = T)
# path = paste0(path_dir,condition ,"/Visualisation/")
###################################

analyseName = "DESeq2_analysis"
analyseName = paste0(Sys.Date(),"_", analyseName,"/")
path_dir = paste0("./Analyse/",analyseName)

#### Definition of significant threshold for DESeq2 ####
FC = 1.5 #Fold-change (mini: 1.5)
pvalue = 0.05 #p-value (maxi 0.05)


for (condition in names(rnai_list)){
  print(paste("Analysis of", condition , "data -->", paste(rnai_list[[condition]], collapse = ", "), "time courses" ))
  
  path = paste0(path_dir,condition,"_FC-", FC, "_pval-", pvalue ,"/")
  graph_path = paste0(path ,"/Visualization/")
  dir.create(graph_path,recursive=T,showWarnings=F)
  
  # Boxplot visualization of raw data (no ComBat-seq correction nor normalisation)
  countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_uncorrected.tab"), sep="\t",row.names=1,header =  T)
  pdf(paste0(graph_path,"Comptage_bolxplot_uncorrected.pdf"))
    CountBoxplot(countdata, "row", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
  dev.off()
  
  # Open files with or without ComBat-seq correction
  countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_corrected.tab"), sep="\t",row.names=1,header =  T)
  
  # Generation the appropriate infodata
  infodata = CreatInfoData(countdata, conditions = condition , rnai_list, cluster)
  write.table(infodata,paste0(path,condition ,"_infodata.tab"), sep="\t",row.names=T,quote=F)
  
  #### DESeq2 analysis ####
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
  
  # Dispersion estimation
  pdf(paste0(graph_path,"Dispersion_DESeq2.pdf"))
  plotDispEsts(deseq, ylim = c(1e-6, 1e1))
  dev.off()
  
  # Extraction of the normalized data
  data_tab = assay(vst(deseq, blind = T))
  write.table(data_tab,paste0(path,condition ,"_expression_table_vst.tab"), sep="\t",row.names=T,quote=F)
  
  #### Data visualization ####
  print(paste(condition,": Visualization"))
  source("4-1_Data_visualization.R")
  
  
  #### Definition of the deregulated genes ####
  
  # Definition of the time course that will be analysed
  RNAi_list = unique(infodata_collapse$KnockDown)
  RNAi_list = RNAi_list[RNAi_list != "ctrl"]
  
  # Definition of the autogamy stages that will be compared
  stages = infodata$Condition
  stages = list(
    "VEG" = unique(stages[grep("VEG",stages)]),
    "EARLY" = unique(stages[grep("EARLY",stages)]),
    "INTER" = unique(stages[grep("INTER",stages)]),
    "LATE" = unique(stages[grep("LATE",stages)])
  )

  # Definition of deregulated genes at each autogamy stage for each time course
  for (RNAi in RNAi_list){
    print(paste(condition, ": Definition of deregulated genes in", RNAi, "time course" ))
    
    # Saving folder creation
    img_dir=paste0(path,"/DESeq/",RNAi,"/Images/")
    dir.create(img_dir,recursive=T,showWarnings=F)
    
    res_dir=paste0(path,"/DESeq/",RNAi,"/")
    dir.create(res_dir, recursive=T,showWarnings=F)
    
    # For each time course
    source("4-2_DESeq2_results.R") # Extracting the results from DESeq2 analysis
    
  }
  
  #### Comparison with the published data of EZL1 deregulated genes ####
  
  if (condition == "HiSeqvsNextSeq"){
    print("Comparison with the published datas ")
    
    # Using the significant threshold published in Frapporti et al. 2019
    frapp_FC = 2
    frapp_pvalue = 0.05
    
    save_path = paste0(path,"/Frapporti_Comparison_FC-",frapp_FC, "_pval-",frapp_pvalue,"/")
    dir.create(save_path,recursive=T,showWarnings=F)
    
    source("4-3_Comparison_Frapporti.R")
    
  }
}

sink(paste0("./Analyse/",analyseName,"/sessionInfo.txt"))
print(sessionInfo())
sink()