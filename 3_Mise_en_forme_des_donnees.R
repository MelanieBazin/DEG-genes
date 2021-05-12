source("0_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")

path = "./Analyse/DESeq2_test01"
annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")

# Definition des paramètres selon l'ARNi étudier ########
source("0_Cluster.R")

# Choix du groupes de donnée à analysé
i= names(rnai_list)[1]

# Ouverture des fichiers et création de l'objet countdata
countdata = ConcatTab("EXPRESSION", conditions = rnai_list[[i]])
row.names(countdata)=countdata$ID
countdata=countdata[,-1]

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
#############################

# Boxplot des comptages non-normalisés
png(paste0(path,"Graph/",i,"_Row_Boxplot.png"))
  CountBoxplot(countdata, "DESeq2")
dev.off()

# Création du tableau avec les info des colonnes
infodata=CreatInfoData2(conditions = rnai_list[[i]])

####### Mise en forme des données pour DESeq2 ##############
countdata =  as.matrix(countdata)
library(DESeq2)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ RNAi + Timing)


# Analyse DESeq2
deseq = DESeq(deseq)

# Graphique du paramètre de dispersion
png(paste0(path,"Graph/",i,"_dipression_DESeq2.png"))
  plotDispEsts(deseq, ylim = c(1e-6, 1e1))
dev.off()

# Récupération des données de comptage normalisées
tab=counts(deseq,normalized=T)

# Boxplot des comptages normalisés
png(paste0(path,"Graph/",i,"_DESeq2_Boxplot.png"))
  CountBoxplot(tab, "DESeq2")
dev.off()
