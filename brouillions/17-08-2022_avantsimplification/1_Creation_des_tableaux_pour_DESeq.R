###################
#### Put all data that need will be analysed together in a table
#### + corrected the expression count by ComBat-seq
#### + creation of the infodata table containing all information about the samples
###################

options(stringsAsFactors = FALSE)
library(sva)

set.seed(10111)

source("0_Cluster.R") #Groups parameters
source("0_Visualisation_fonction.R") #Library of homemade function



dir.create("./DATA/Pour_DESeq/",recursive=T,showWarnings=F)

for (i in names(rnai_list)){
  # Opening files
  countdata = ConcatTab(type = "EXPRESSION", conditions = rnai_list[[i]])
  
  # Save the table required for DESeq2 analysis
  infodata = CreatInfoData(countdata, conditions = i , rnai_list, cluster)
  write.table(countdata,paste0("./DATA/Pour_DESeq/",i,"_infodata.tab"), sep="\t",row.names=T,quote=F)
  
  # Before ComBat-seq correction
  countdata = as.matrix(countdata)
  write.table(countdata,paste0("./DATA/Pour_DESeq/",i,"_expression_table_uncorrected.tab"), sep="\t",row.names=T,quote=F)
  
  # After ComBat-seq correction
  batch = paste(infodata$Seq_method,infodata$Labo, sep = "_")
  countdata = ComBat_seq(countdata, batch = batch)
 
  write.table(countdata,paste0("./DATA/Pour_DESeq/",i,"_expression_table_corrected.tab"), sep="\t",row.names=T,quote=F)
   
  print(paste("Table done for",i))
}
