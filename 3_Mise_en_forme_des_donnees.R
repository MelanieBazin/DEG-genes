source("2_Fonction.R")

annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")

#### Definition des paramètres selon l'ARNi étudier ########
if (condition == "CTIP") {
  # Sélectionner les points à regrouper comme faisant parti du même timing
  timing_ctrl = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",2))
  timing_rnai = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",1))

}else if (condition == "XRCC4"){
  # Sélectionner les points à regrouper comme faisant parti du même timing
  timing_ctrl = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",2))
  timing_rnai = timing_ctrl

}else if (condition == "PGM"){
  # Sélectionner les points à regrouper comme faisant parti du même timing
  timing_ctrl_nd7 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3))
  timing_ctrl_icl7 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3))
  timing_ctrl = c(timing_ctrl_nd7, timing_ctrl_icl7)
  timing_rnai = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3))
  
}else if (condition == "KU80c"){
  # Sélectionner les points à regrouper comme faisant parti du même timing
  timing_ctrl_nd7 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3))
  timing_ctrl_icl7 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3))
  timing_ctrl = c(timing_ctrl_nd7, timing_ctrl_icl7)
  timing_rnai = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3))

}

#Ouverture des fichier et création de l'objet countdata
countdata = OpenDataCount(path, condition)

####### Sélection de gènes #####
# selection = c("Ku","ku","PGM","NOWA","PTIWI","mt","TFIIS4","Spo11","Mre11","CER","Rad51", "Lig", "EZL", "SPT", "DCL", "CtIP", "XRCC4", "PDSG2", "PolX", "CAF1")
# 
# selection_ID =c()
# 
# for(i in selection){
#   selection_ID = c(selection_ID,annotation$ID[grep(i,annotation$SYNONYMS)])
#   
# }
# countdata = countdata[is.element(rownames(countdata), selection_ID),]


####### Mise en forme des données pour DESeq2 ##############
countdata =  as.matrix(countdata)

# Création du tableau avec les info des colonnes
infodata = matrix(NA,nrow = ncol(countdata), ncol = 4)
row.names(infodata) = colnames(countdata)
colnames(infodata) = c("Noms", "Feeding", "Timing", "Conditions")
infodata[,"Noms"] = colnames(countdata)

infodata[,"Feeding"] = c(rep("controle",length(grep("CTRL", colnames(countdata)))),
                         rep(paste("RNAi",condition),length(grep("RNAi", colnames(countdata)))))
infodata[,"Timing"] = c(timing_ctrl,timing_rnai)
infodata[,"Conditions"] = c(paste("CTRL", timing_ctrl, sep = "_" ),paste(condition, timing_rnai, sep = "_" ))

infodata = as.data.frame(infodata)
