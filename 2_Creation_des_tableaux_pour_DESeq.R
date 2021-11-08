### A Executer sur R version 4

options(stringsAsFactors = FALSE)
source("0_Cluster.R")
source("3_Visualisation_des_donnees_fonction.R")
library(sva)

set.seed(10111)

dir.create("./DATA/Pour_DESeq/",recursive=T,showWarnings=F)
dir.create("./DATA/Pour_DESeq_SansCorrectionBatch/",recursive=T,showWarnings=F)


i = names(rnai_list)[2]
for (i in names(rnai_list)[1:2]){
  ##### Création du tableau de donnée à analyser ensemble ####
  #Ouverture des fichiers et création de l'objet countdata
  countdata = ConcatTab(type = "EXPRESSION", conditions = rnai_list[[i]])

  # Enregistrement des tableau qui seront utiliser pour DESeq
  infodata = CreatInfoData3(countdata, conditions = i , rnai_list, cluster)
  countdata = as.matrix(countdata)
  # Correction de l'effet batch avec ComBat
  if (i =="tout"){
    write.table(countdata,paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",i,"_expression_table_ROW.tab"), sep="\t",row.names=T,quote=F)
    batch = paste(infodata$Batch,infodata$Labo, sep = "_")
    countdata2 = ComBat_seq(countdata, batch = infodata$Batch)
    
    countdata = ComBat_seq(countdata, batch = batch, group = infodata$Cluster)
  } else if(i == "seq2014vs2020"){
    write.table(countdata,paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",i,"_expression_table_ROW.tab"), sep="\t",row.names=T,quote=F)
    timing = unique(infodata$Timing[which(infodata$Batch == "seq_2020")])
    countdata = countdata[,is.element(infodata$Timing, timing)]
    infodata = infodata[is.element(infodata$Timing, timing),]
    countdata2 = ComBat_seq(countdata, batch = infodata$Batch)
    
    
    countdata = ComBat_seq(countdata, batch = infodata$Batch, group = infodata$Timing)
  } else if (length(grep("ICL7", colnames(countdata)))>0){
    write.table(countdata,paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",i,"_expression_table_ROW.tab"), sep="\t",row.names=T,quote=F)
    countdata = ComBat_seq(countdata, batch = infodata$Labo, group = infodata$Cluster)
    countdata2 = ComBat_seq(countdata, batch = infodata$Labo)
  }
  write.table(countdata2,paste0("./DATA/Pour_DESeq/",i,"_expression_table_pour_DESeq_v1.tab"), sep="\t",row.names=T,quote=F)
  write.table(countdata,paste0("./DATA/Pour_DESeq/",i,"_expression_table_pour_DESeq_v2.tab"), sep="\t",row.names=T,quote=F)
  
  print(paste("Tableau pour la condition",i, "termine"))
}

# Pour ouvrir les tableau avec correction BATch
# countdata = read.table(paste0("./DATA/Pour_DESeq/",i,"_expression_table_pour_DESeq.tab"), sep="\t",row.names=1,header =  T)