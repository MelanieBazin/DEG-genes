source("0_Visualisatin  des données_fonction.R", encoding = "UTF-8")

analyseName = "Test_normalisation"
condition = c("CTIP", "XRCC4", "PGM", "KU80c")
condition = "PGM" #4 possibilitées : "CTIP", "XRCC4", "PGM", "KU80c"


#######################################################
# Partie 1 : Comparaison des méthodes de normalisaion #
#######################################################
# Comptage
Type = c("EXPRESSION","RPM", "RPKM")
pdf(paste(analyseName,"Boxplot_Normalsation1.pdf", sep="_"))
for (i in Type){
  tab = ConcatTab(i)
  rownames(tab)=tab$ID
  tab = tab[,-1]
    if (i == "EXPRESSION"){type = "raw data"}else{type = i}
  CountBoxplot(tab, type)
}

############ Normalisation DESeq2 ############
source("2_Fonction.R")


#Ouverture des fichiers et création de l'objet countdata
countdata = ConcatTab("EXPRESSION")
row.names(countdata)=countdata$ID
countdata=countdata[,-1]
countdata =  as.matrix(countdata)

# Création du tableau avec les info des colonnes
infodata=matrix(NA,nrow = ncol(countdata), ncol=3)
row.names(infodata) = colnames(countdata)
colnames(infodata) = c("Name","RNAi","Timing")
infodata[,"Name"] = row.names(infodata)



CTIP = c("T0", "T5.5", "T12.5", "T25", "Veg")
CTIP_CTRL = c("T0", "T5", "T10", "T20", "T30", "Veg")
ICL7 = c("T0", "T5", "T10", "T20", "T35", "T50", "Veg")
KU80c = c( "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
ND7 = c( "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
PGM = c( "T2", "T5", "T10", "T20", "T30", "T40", "Veg")
XRCC4 = c( "T2", "T7", "T22", "T32","Veg")
XRCC4_CTRL = c( "T2", "T7", "T22", "T32","Veg")

infodata[,"Timing"] = c(CTIP, CTIP_CTRL , ICL7, KU80c , ND7, PGM, XRCC4, XRCC4_CTRL)

condi = sub(".tab","",list.files("./DATA/EXPRESSION"))
l = list(CTIP, CTIP_CTRL , ICL7, KU80c , ND7, PGM, XRCC4, XRCC4_CTRL)
rnai = c()
for (i in 1:length(l)){
  rnai = c(rnai, rep(condi[i],length(l[[i]])))
}
infodata[,"RNAi"] = rnai
infodata = as.data.frame(infodata)

# Mise en forme des données
library(DESeq2)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ RNAi + Timing)


# Analyse DESeq2
deseq = DESeq(deseq)

# Récupération des données de comptage normalisées
tab=counts(deseq,normalized=T)
CountBoxplot(tab, "DESeq2")

# Graphique du paramètre de dispersion
plotDispEsts(deseq, ylim = c(1e-6, 1e1))

dev.off()



########################################################
# Partie 2 : Visualisation des données pour clustering #
########################################################
distance_methode =  c("Euclidean", "Correlation")
method_utilisee = c("kmeans", "HCL")


nb_cluster = 4 #VEG, EARLY, INTERMEDIATE, LATE
graph_type = c("heatmap","profils")

#### Visualisation donnée de clustering

pdf(paste0("Lvl3_profil_",sub(".txt","",Nom_de_fichier),"_",nb_cluster,"clusters_",method_utilisee,"_",distance_methode,".pdf"))

Clustering(Nom_de_fichier = Mon_fichier_S,
           distance = distance_methode,
           nb_cluster = nb_cluster,
           method = method_utilisee,
           graph_type = graph_type)

dev.off()
