### A Executer sur R version 4

options(stringsAsFactors = FALSE)
source("0_Cluster.R")
source("3_Visualisation_des_donnees_fonction.R")
library(sva)
set.seed(10111)

dir.create("./DATA/Pour_DESeq/",recursive=T,showWarnings=F)
dir.create("./DATA/Pour_DESeq_SansCorrectionBatch/",recursive=T,showWarnings=F)

rnai_list = rnai_list[2:4]

i = names(rnai_list)[1]
# for (i in names(rnai_list)){
  ##### Création du tableau de donnée à analyser ensemble ####
  #Ouverture des fichiers et création de l'objet countdata
  countdata = ConcatTab(type = "EXPRESSION", conditions = rnai_list[[i]])
  write.table(countdata,paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",i,"_expression_table_ROW.tab"), sep="\t",row.names=T,quote=F)
  
  # Enregistrement des tableau qui seront utiliser pour DESeq
  infodata = CreatInfoData3(countdata, conditions = i , rnai_list, cluster)
  countdata = as.matrix(countdata)
  # Correction de l'effet batch avec ComBat
  if (i =="tout"){
    batch = paste(infodata$Batch,infodata$Labo, sep = "_")
    countdata = ComBat_seq(countdata, batch = batch, group = infodata$Cluster)
  } else if(i == "seq2014vs2020"){
    countdata = ComBat_seq(countdata, batch = infodata$Batch, group = infodata$Timing)
  } else if (length(grep("ICL7", colnames(countdata)))>0){
    countdata = ComBat_seq(countdata, batch = infodata$Labo, group = infodata$Cluster)
  }
  write.table(countdata,paste0("./DATA/Pour_DESeq/",i,"_expression_table_pour_DESeq_v2.tab"), sep="\t",row.names=T,quote=F)
  
  print(paste("Tableau pour la condition",i, "termine"))
# }

# Pour ouvrir les tableau avec correction BATch
# countdata = read.table(paste0("./DATA/Pour_DESeq/",i,"_expression_table_pour_DESeq.tab"), sep="\t",row.names=1,header =  T)