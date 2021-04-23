####
# Gaëlle LELANDAIS <gaelle.lelandais@universite-paris-saclay.fr>
##
# This content is licensed under CC BY 4.0
# http://creativecommons.org/licenses/by-sa/4.0/ 
####

# -----------------------------------------------------------------------------
# Les données utilisées pour cette mise en application sont présentées 
# dans l'article suivant :
#
# Ostreococcus tauri is a new model green alga for studying iron metabolism 
# in eukaryotic phytoplankton.
#
# Lelandais G, Scheiber I, Paz-Yepes J, Lozano JC, Botebol H, Pilátová J, 
# Žárský V, Léger T, Blaiseau PL, Bowler C, Bouget FY, Camadro JM, Sutak R, 
# Lesuisse E.
#
# BMC Genomics. 2016 May 3;17:319. doi: 10.1186/s12864-016-2666-6.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Lecture des données
# -----------------------------------------------------------------------------

# Deux fichiers sont à lire ici. Ils contiennent des informations complémentaires, 
# les données de comptage (counts) pour l'ensemble des gènes dans l'ensemble 
# des expériences réalisées d'une part, et la description des conditions 
# associées à chacune des expériences.

# Table des comptages (counts)
countData = as.matrix(read.table("mapping_rawdata_allGenes.txt", 
                      row.names = 1, header = T))
# Table de description des conditions
dataInfo  = read.table("sample_info.txt", header = T)

# La première colonne de la table de comptage est retirée, il s'agit de la 
# longueur des gènes (non utilisée dans cette analyse).
countData = countData[,-1]

# -----------------------------------------------------------------------------
# Premières représentations des données
# -----------------------------------------------------------------------------

# Histogramme des valeurs de comptage
hist(countData, xlab = "count values", 
     main = "Histogram of counts\n (raw data - all experiments)")

# Transformation log des données, pour une meilleure lisibilité de 
# l'échelle des valeurs
logCounts = log(countData + 1) # Ajout de la valeur +1 car il y a des 0
hist(logCounts, xlab = "count values (log scale)",
     main = "Histogram of counts - log transformed\n (raw data - all experiments)")
# --> Graphique avec une meilleure lisibilité des données
# --> Des gènes avec très peu de counts sont observés. Ces gènes sont
#     filtrés.

# Filtrage des gènes pour lesquels les valeurs de comptage sont très faibles
# (la valeur 50 est choisie ici arbitrairement)
countData = countData[rowSums(countData) > 50,]
hist(log(countData + 1), main = "Histogram of counts (after filtering)")

# Taille du tableau de données obtenu
dim(countData)
# --> 7610 gènes et 47 expériences

# NOTE : La problématique du filtrage des gènes très faiblement exprimés est 
# importante. Un point de vigilance qui peut être réfléchi/discuté.

# -----------------------------------------------------------------------------
# Normalisation des valeurs de comptage pour tenir compte du biais de la 
# taille de librairie (par expérience)
# -----------------------------------------------------------------------------

# Valeurs de comptage (avant normalisation)
boxplot(log(countData + 1), ylab = "count values (log scale)",
        main = "Count data (raw values)")
# --> Ici des fluctuations sont observées, mais peu importantes. 
#     La correction statistique ne sera pas très forte.

# Mise en application de la normalisation proposée dans DESeq
library(DESeq2)

# Mise en forme des données pour DESeq2, afin d'utiliser ses fonctionnalités
# de normalisation des counts
colData = data.frame(dataInfo)
dds     = DESeqDataSetFromMatrix(countData = countData,
                                  colData  = colData,
                                  design   = ~ Light + Time + Condition + Iron)
# Note : l'argument "design" sera changé par la suite 
# (ne pas en tenir compte ici)

# Calcul des facteurs de normalisation
normParam = sizeFactors(estimateSizeFactors(dds))
# --> Ces paramètres indiquent à quel point la correction des valeurs de
#     comptage doit être plus ou moins forte. Certaines valeurs de counts vont
#     être diminuées (size factor > 1) et d'autres vont être augmentées 
#     (size factor < 1).

# Réalisation de la normalisation (correction) des valeurs de counts. Cela se
# fait par expérience, en utilisant les valeurs de "size factors" calculées
# ci dessus
normDataCount = countData
for(i in 1:ncol(countData)){
  normDataCount[,i] = countData[,i] / normParam[i]
}

# Valeurs de comptage (après normalisation)
boxplot(log(normDataCount + 1), 
        main = "Count data \n(after librairy size correction)")
# --> La normalisation est terminée. L'effet de la taille des librairies a
#     été pris en compte.

# NOTE : La problématique de la normalisation par la taille des librairies est 
# importante. Un point de vigilance qui peut être réfléchi/discuté.

# -----------------------------------------------------------------------------
# Analyse exploratoire des données de comptage (normalisées), quelle
# cohérence entre les réplicats biologiques ?
# -----------------------------------------------------------------------------

# Les expériences RNAseq ont été répétées 3 fois. Ainsi nous disposons de 
# réplicats biologiques. Pour répondre à la question (cohérence entre les
# réplicats), nous allons utiliser un méthode graphique qui est l'ACP 
# (Analyse en Composante Principale).
library(FactoMineR)

# POINT de VIGILANCE pour l'ACP : 
# Les individus doivent positionnés en lignes dans le tableau de données,
# utilisé pour l'ACP. Ici, le tableau est la table de comptage et les individus
# que l'on souhaite considérer sont les expériences (les 47 colonnes).
# Il est donc nécessaire de transposer le tableau des counts normalisés,
# cela est fait avec l'utilisation de la fonction "t()".
resPCA = PCA(t(normDataCount), graph = F)

# La représentation graphique des résultats est réalisée,
# en utilisant la librairie "factoextra", celle ci utilise "ggplot2".
library(factoextra)

# Combien de conditions différentes sont associées à l'ensemble des
# expériences ?
nbCond = length(unique(colData$Name)) 
# --> 16 conditions différentes (avec trois réplicats pour chaque)

# Visualisation du plan principal (axes 1 et 2)
fviz_pca_ind(resPCA,  # Représentation graphique des expériences (47 ici)
             col.ind = colData$Name,  # Les réplicats sont colorés en fonctions des conditions
             palette = rainbow(nbCond),  # Définition d'une palette de couleurs
             addEllipses = F,            # Pas d'ellipses tracées
             label = "none",             # Pas de texte
             legend.title = "Condition", # Titre de la légende
             mean.point = F,             # Pas d'affichage du centre de gravité
             pointsize = 3,              # Taille des points
             pointshape = 19)            # Forme des points
# --> Les réplicats des différentes conditions sont bien proches entre eux.
#     Pas de souci particulier ici.

# -----------------------------------------------------------------------------
# Séparation des données, d'une part la cinétique "Short term" et
# d'autre part la cinétique "Long term".
# -----------------------------------------------------------------------------

# Le séquençage des échantillons biologiques de l'une et l'autre des cinétiques 
# ayant été réalisé simultanément, il était important de réaliser l'ACP 
# précédente sans les séparer. Toutefois, les cinétiques ST et LT 
# correspondent à des conditions de croissance des cellules très différentes. 
# Cela est très bien visualisé en considérant les axes 2 et 3 de l'ACP 
# précédemment réalisée

# Visualisation du plan principal (axes 2 et 3)
fviz_pca_ind(resPCA,  # Représentation graphique des expériences (47 ici)
             axes = c(2,3),
             col.ind = colData$Condition,  # Les réplicats sont colorés en fonctions des conditions
             palette = c("orange", "pink"),  # Définition d'une palette de couleurs
             addEllipses = F,            # Pas d'ellipses tracées
             label = "none",             # Pas de texte
             legend.title = "Time course", # Titre de la légende
             mean.point = F,             # Pas d'affichage du centre de gravité
             pointsize = 3,              # Taille des points
             pointshape = 19)            # Forme des points

# --> Pour la suite des analyses, il est raisonnable de séparer les analyses de
# chacune des deux cinétiques. 

# ----------------------------
# A) Cinétique ST
# ----------------------------

# -- Tables de comptage
countDataST     = countData[,dataInfo[dataInfo[,"Condition"] == "ST", "SampleID"]]
normDataCountST = normDataCount[,dataInfo[dataInfo[,"Condition"] == "ST", "SampleID"]]

# -- Information des échantillons
dataInfoST  = dataInfo[dataInfo[,"Condition"] == "ST",]
dataInfoST  = dataInfoST[,-5] # La colonne "Condition" ne sert plus à rien.

# -- Mise en forme des données pour DESeq2
colDataST = data.frame(dataInfoST)
ddsST     = DESeqDataSetFromMatrix(countData = countDataST,
                                   colData   = colDataST,
                                   design   = ~ Light + Time + Iron)
# ----------------------------
# B) Cinétique LT
# ----------------------------

# -- Tables de comptage
countDataLT = countData[,dataInfo[dataInfo[,"Condition"] == "LT", "SampleID"]]
normDataCountLT = normDataCount[,dataInfo[dataInfo[,"Condition"] == "LT", "SampleID"]]

# -- Information des échantillons
dataInfoLT  = dataInfo[dataInfo[,"Condition"] == "LT",]
dataInfoLT  = dataInfoLT[,-5] # La colonne "Condition" ne sert plus à rien.

# -- Mise en forme des données pour DESeq2
colDataLT = data.frame(dataInfoLT)
ddsLT     = DESeqDataSetFromMatrix(countData = countDataLT,
                                   colData   = colDataLT,
                                   design   = ~ Light + Time + Iron)

# -----------------------------------------------------------------------------
# Analyse exploratoire des données de comptage, quels effets des
# facteurs sur la variabilité des expressions des gènes ?
# -----------------------------------------------------------------------------

# Trois facteurs peuvent avoir un effet sur les valeurs de comptage observées 
# dans les tables ST et LT. Ces facteurs sont : 1) La luminosité, 
# 2) La présence/absence de Fer et 3) Le temps de la cinétique.

# Ici notre objectif est d'évaluer l'impact respectif de ces facteurs sur
# les mesures d'expression des gènes. A nouveau l'ACP sera utilisée.

# ----------------------------
# A) Cinétique ST
# ----------------------------

# Réalisation de l'ACP
resPCA1 = PCA(t(normDataCountST), graph = F) # pas d'affichage graphique

# Représentation des résultats (conditions S1 à S8)
nbCondST = length(unique(colDataST$Name)) # 8 conditions différentes
fviz_pca_ind(resPCA1,                  
             col.ind = colDataST$Name, 
             palette = rainbow(nbCondST), 
             addEllipses = F, 
             label = "none",
             legend.title = "Condition",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 

# Coloration des individus en fonction de la luminosité
fviz_pca_ind(resPCA1,                  
             col.ind = colDataST$Light, 
             palette = c("blue", "yellow"), 
             addEllipses = F, 
             label = "none",
             legend.title = "Light",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 
# --> La luminosité est ici un facteur discriment important des expériences !

# Coloration des individus en fonction du Fer
fviz_pca_ind(resPCA1,                  
             col.ind = colDataST$Iron, 
             palette = c("red", "black"), 
             addEllipses = F, 
             label = "none",
             legend.title = "Light",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 

# Regardons les composantes principales suivantes ...
fviz_pca_ind(resPCA1,      
             axes = c(2,3),
             col.ind = colDataST$Iron, 
             palette = c("red", "black"), 
             addEllipses = F, 
             label = "none",
             legend.title = "Light",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 
# --> Le Fer est également un facteur discriment, mais moins important que
#     la luminosité.

# Coloration des individus en fonction du temps
fviz_pca_ind(resPCA1,  
             axes = c(3,4),
             col.ind = colDataST$Time, 
             palette = c("pink", "purple"), 
             addEllipses = F, 
             label = "none",
             legend.title = "Time",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 
# --> Le Temps est également un facteur discriment, mais moins important
# que les deux autres.

# Ces trois facteurs sont ainsi à prendre en compte, pour que la modélisation
# DESeq2 puisse être la plus pertinente possible.

# ----------------------------
# B) Cinétique LT
# ----------------------------

# Réalisation de l'ACP
resPCA2 = PCA(t(normDataCountLT), graph = F)

# Représentation des résultats
nbCondLT = length(unique(colDataLT$Name)) # 8 conditions différentes
fviz_pca_ind(resPCA2,                  
             col.ind = colDataLT$Name, 
             palette = rainbow(nbCondLT), 
             addEllipses = F, 
             label = "none",
             legend.title = "Condition",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 

# Coloration des individus en fonction de la luminosité
fviz_pca_ind(resPCA2,                  
             col.ind = colDataLT$Light, 
             palette = c("blue", "yellow"), 
             addEllipses = F, 
             label = "none",
             legend.title = "Light",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 
# --> La luminosité est ici un facteur discriment important.

# Coloration des individus en fonction du Fer
fviz_pca_ind(resPCA2,      
             axes = c(2,3),
             col.ind = colDataLT$Iron, 
             palette = c("red", "black"), 
             addEllipses = F, 
             label = "none",
             legend.title = "Light",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 
# --> Le Fer est également un facteur discriment, mais moins important.

# Coloration des individus en fonction du temps
fviz_pca_ind(resPCA2,                  
             col.ind = colDataLT$Time, 
             palette = c("blue", "dark blue", "pink", "purple"), 
             addEllipses = F, 
             label = "none",
             legend.title = "Light",  
             mean.point = F,  
             pointsize = 3, 
             pointshape = 19) 
# --> Le Temps est également un facteur discriment, plus fort ici que dans
# le cas de la ST.

# -----------------------------------------------------------------------------
# Conclusion 
# -----------------------------------------------------------------------------

# Cette première étude exploratoire des données nous a permi :
# 1) De vérifier la cohérence globale des réplicats biologiques entre eux,
# 2) De séparer les données des cinétiques "Short Term" (ST) et "Long Term" (LT),
# 3) D'évaluer les rôles des 3 facteurs (Luminosité, Fer et Temps) sur 
#    les mesures d'expression des gènes. C'est important pour la bonne mise
#    en application de DESeq2.