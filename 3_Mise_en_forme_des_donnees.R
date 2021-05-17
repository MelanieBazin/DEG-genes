source("2_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")

# Partie non necessaire si on lance 4_Analyse_DESeq2.R
# 
# analyseName = "DESeq2_test01"
# 
# # Definition des paramètres selon l'ARNi étudier ########
# source("0_Cluster.R")
# 
# # Choix du groupes de donnée à analysé
# i = names(rnai_list)[1]
# print(paste("On analyse le jeu de donnée :", i, "-->", paste(rnai_list[[i]], collapse = ", ") ))



### Création des dossier pour ranger les données ###
path =paste0("./Analyse/",analyseName,"/",i,"/Images/")
dir.create(path,recursive=T,showWarnings=F)


###### Ouverture des fichiers et création de l'objet countdata #####
annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")
countdata = ConcatTab("EXPRESSION", conditions = rnai_list[[i]])

# Passage de la colonne des ID en rowname
if (colnames(countdata)[1]=="ID"){
  rownames(countdata)=countdata$ID
  countdata = countdata[,-1]
}

# Changer le nom des colonnes controles
library("stringr")  
colnames(countdata) = str_replace_all(colnames(countdata),"ND7","ND7_K")
colnames(countdata) = str_replace_all(colnames(countdata),"CTIP_CTRL","ND7_C")
colnames(countdata) = str_replace_all(colnames(countdata),"XRCC4_CTRL","ND7_X")


####### Sélection de gènes #####
selection = c("Ku","ku","PGM","NOWA","PTIWI","mt","TFIIS4","Spo11","Mre11","CER","Rad51", "Lig", "EZL", "SPT", "DCL", "CtIP", "XRCC4", "PDSG2", "PolX", "CAF1")

selection_ID =c()

for(s in selection){
  selection_ID = c(selection_ID,annotation$ID[grep(s,annotation$SYNONYMS)])

}
# countdata = countdata[is.element(rownames(countdata), selection_ID),]
#############################

# Boxplot des comptages non-normalisés
png(paste0(path,i,"_Row_Boxplot.png"))
  CountBoxplot(countdata, "DESeq2")
dev.off()

# Création du tableau avec les info des colonnes
infodata=CreatInfoData1(countdata, conditions = i, rnai_list, cluster)

####### Mise en forme des données pour DESeq2 ##############
countdata =  as.matrix(countdata)
library(DESeq2)

# Deux paramètres possible pour l'analyse
if (parametre_DESeq2 == "Conditions"){
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ Conditions)
}else if (parametre_DESeq2 == "Feeding_Cluster"){
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ Feeding + Cluster)
}



# Analyse DESeq2
deseq = DESeq(deseq)

# Graphique du paramètre de dispersion
png(paste0(path,i,"_dipression_DESeq2.png"))
  plotDispEsts(deseq, ylim = c(1e-6, 1e1))
dev.off()

# Récupération des données de comptage normalisées
tab=counts(deseq,normalized=T)

# Boxplot des comptages normalisés
png(paste0(path,i,"_DESeq2_Boxplot.png"))
  CountBoxplot(tab, "DESeq2")
dev.off()
