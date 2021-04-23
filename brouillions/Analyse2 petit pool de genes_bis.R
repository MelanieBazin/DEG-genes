options(stringsAsFactors = FALSE)

annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")

type = "EXPRESSION" #3 possibiltées : "RPM","RPKM", "EXPRESSION"

condition = "CTIP"

#### Ouverture des fichiers ####
files = list.files(path = paste0("./DATA/",type))
files = files[grep(condition,files)]

control = read.table(paste0("./DATA/", type,"/",files[grep("CTRL",files)]),h = T,sep = "\t")
rnai = read.table(paste0("./DATA/", type,"/",files[setdiff(1:length(files),grep("CTRL",files))]),h = T,sep = "\t")

####### Sélection de gènes #####
selection = c("Ku","ku","PGM","NOWA","PTIWI","mt","TFIIS4","Spo11","Mre11","CER","Rad51", "Lig", "EZL", "SPT", "DCL", "CtIP", "XRCC4", "PDSG2", "PolX", "CAF1")
selection = sort(selection)

selection_ID =c()

for(i in selection){
  selection_ID = c(selection_ID,annotation$ID[grep(i,annotation$SYNONYMS)])
  
}
control = control[which(is.element(selection_ID, control$ID)),]
rnai = rnai[which(is.element(selection_ID, rnai$ID)),]

#### Definition des paramètres selon l'ARNi étudier ########

if (condition == "CTIP") {
  # Sélectionner les points à regrouper comme faisant parti du même timing
  timing_ctrl = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",2))
  timing_rnai = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",1))
  
  # Réordonner les colonne du tableau pour les mettre dans l'ordre de la cinétique
  control = control[,c(1,7,2:6)]
  rnai = rnai[,c(1,6,2:5)]
  
  # Definiton des variables DESeq2
  FC = 2
  pvalue = 0.2
  
}
if (condition == "XRCC4"){
  # Sélectionner les points à regrouper comme faisant parti du même timing
  timing_ctrl = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",2))
  timing_rnai = timing_ctrl
  
  # Réordonner les colonne du tableau pour les mettre dans l'ordre de la cinétique
  control = control[,c(1,6,2:5)]
  rnai = rnai[,c(1,6,2:5)]
  
  # Definiton des variables DESeq2
  FC = 2
  pvalue = 0.2
}

####### Mise en forme des données pour DESeq2 ##############
# Fusion des tableau de comptage controle et ARNi
colnames(control) = paste(colnames(control),"CTRL", sep = "_" )
colnames(rnai) = paste(colnames(rnai),"RNAi", sep = "_" )
countdata=merge(control, rnai, by.x = "ID_CTRL", by.y = "ID_RNAi")

row.names(countdata) = countdata$ID_CTRL
countdata = as.matrix(countdata[,-1])

# Création du tableau avec les info des colonnes
infodata = matrix(NA,nrow = ncol(countdata), ncol = 3)
row.names(infodata) = colnames(countdata)
colnames(infodata) = c("Noms", "Conditions", "Timing")
infodata[,"Noms"] = colnames(countdata)
infodata[,"Conditions"] = c(rep("controle",ncol(control)-1),
                            rep(paste("RNAi",condition),ncol(rnai)-1))
infodata[,"Timining"] = c(timing_ctrl,timing_rnai)

infodata = as.data.frame(infodata)

##### Chargement des librairies pour l'analyse #############
library(FactoMineR)
library(factoextra)
library(DESeq2)

##### Analyse DESeq2 ###########
  
# Filtration des gènes avec trop peu de compatage (seuil arbitraire)
# countdata = countdata[rowSums(countdata) > 50,]

# Mise en forme des données
deseq = DESeqDataSetFromMatrix(countData = countdata,
                                      colData  = infodata,
                                      design   = ~ Conditions + Timing)

# Analyse DESeq2
deseq = DESeq(deseq)

###############################################
# Reprise des variables et analyses d'Olivier #
###############################################
notAllZero = (rowSums(counts(deseq)) > 0 )
labels=colnames(countdata)
time_points = c(paste("CTRL", timing_ctrl, sep = "_" ),paste(condition, timing_rnai, sep = "_" ))


# mean for each group of time points
geneNormCountsTable=counts(deseq,normalized=T)
meanGeneNormCountsTable =data.frame(ID=rownames(geneNormCountsTable))


for(p in unique(time_points)) {
  
  if(length(labels[time_points==p])==1)  { 
    meanGeneNormCountsTable[,p] =geneNormCountsTable[, labels[time_points==p]]
  } else { 
    meanGeneNormCountsTable[,p] =apply(geneNormCountsTable[, labels[time_points==p]],1,mean) 
  }
  
}
rownames(meanGeneNormCountsTable)=meanGeneNormCountsTable$ID
meanGeneNormCountsTable=meanGeneNormCountsTable[,-1]
# head(meanGeneNormCountsTable)

meanGeneNormCountsTableInfo=meanGeneNormCountsTable
meanGeneNormCountsTableInfo$ID=rownames(meanGeneNormCountsTableInfo)


