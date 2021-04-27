source("0_Visualisatin  des données_fonction.R")
#######################################################
# Partie 1 : Comparaison des méthodes de normalisaion #
#######################################################
analyseName = "Test_normalisation"

# RPM

# RPKM

############ Normalisation DESeq2 ############
source("2_Fonction.R")

condition = c("CTIP", "XRCC4", "PGM", "KU80c")
condition = "PGM" #4 possibilitées : "CTIP", "XRCC4", "PGM", "KU80c"

#Ouverture des fichiers et création de l'objet countdata
countdata = OpenDataCount("./DATA/EXPRESSION", condition)
countdata =  as.matrix(countdata)

# Création du tableau avec les info des colonnes
infodata = CreationInfoData(countdata)
infodata = as.data.frame(infodata)

# Mise en forme des données
library(DESeq2)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Feeding + Timing)


# Analyse DESeq2
deseq = DESeq(deseq)

# Récupération des données de comptage normalisées
geneNormCountsTable=counts(deseq,normalized=T)


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
