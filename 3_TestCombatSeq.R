options(stringsAsFactors = FALSE)
library(sva)
library(DESeq2)
library(pheatmap)
library(MASS)

set.seed(10111)

# Récupérer les paramètre de clustring
source("0_Cluster.R")

# Récupérer les fonction necessaire au représentaion graphique et la mise en forme des données
source("3_Visualisation_des_donnees_fonction.R")


analyseName = paste0("Test_Combatseq")

path_dir = paste0("./Analyse/",analyseName,"/")
dir.create(path_dir,recursive=T,showWarnings=F)

condition = names(rnai_list)[2]

print(paste("On analyse le jeu de donnee :", condition , "-->", paste(rnai_list[[condition]], collapse = ", ") ))

path = paste0(path_dir,condition ,"/")
dir.create(path,recursive=T,showWarnings=F)

##### Analyse DESeq2 ####
##### SANS correction de l'effet Batch ####
countdata = read.table(paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",condition ,"_expression_table_ROW.tab"), sep="\t",row.names=1,header =  T)

# Création du tableau avec les info des colonnes
infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)

# Creataion de l'objet DESeq2
countdata = as.matrix(countdata)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Condition)


# Analyse DESeq2
print(paste(condition , "-----> Analyse DESeq2"))
deseq = DESeq(deseq)

# Récupération des données de comptage normalisées
data_tab = counts(deseq,normalized=T)

##### Analyse multi-variée des données pour clustering  #####
# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(paste0(j,"_"), color)]=seq_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"/SANS_ComBat"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")

# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(paste0(j,"_"), color)]=cluster_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"/SANS_ComBat_cluster"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")


#####  AVEC correction de l'effet Batch ####
countdata = read.table(paste0("./DATA/Pour_DESeq/",condition,"_expression_table_pour_DESeq_v2.tab"), sep="\t",row.names=1,header =  T)

# Création du tableau avec les info des colonnes
infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)

# Creataion de l'objet DESeq2
countdata = as.matrix(countdata)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Condition)


# Analyse DESeq2
print(paste(condition , "-----> Analyse DESeq2"))
deseq = DESeq(deseq)

# Récupération des données de comptage normalisées
data_tab = counts(deseq,normalized=T)

##### Analyse multi-variée des données pour clustering  #####
# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(paste0(j,"_"), color)]=seq_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"/AVEC_ComBat"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")

# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(paste0(j,"_"), color)]=cluster_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"/AVEC_ComBat_cluster"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")


#### SEquencer en 2014 ####

condition = names(rnai_list)[3]

print(paste("On analyse le jeu de donnee :", condition , "-->", paste(rnai_list[[condition]], collapse = ", ") ))

path = paste0(path_dir,condition ,"/")
dir.create(path,recursive=T,showWarnings=F)

##### Analyse DESeq2 ####
# Ouverture des fichiers countdata SANS correction de l'effet Batch
countdata = read.table(paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",condition ,"_expression_table_ROW.tab"), sep="\t",row.names=1,header =  T)

# Création du tableau avec les info des colonnes
infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)

# Creataion de l'objet DESeq2
countdata = as.matrix(countdata)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Condition)


# Analyse DESeq2
print(paste(condition , "-----> Analyse DESeq2"))
deseq = DESeq(deseq)

# Récupération des données de comptage normalisées
data_tab = counts(deseq,normalized=T)

##### Analyse multi-variée des données pour clustering  #####
# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(paste0(j,"_"), color)]=seq_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"/SANS_ComBat"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")

# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(paste0(j,"_"), color)]=cluster_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"/SANS_ComBat_cluster"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")


#### SEquencer en 2020 ####

condition = names(rnai_list)[4]

print(paste("On analyse le jeu de donnee :", condition , "-->", paste(rnai_list[[condition]], collapse = ", ") ))

path = paste0(path_dir,condition ,"/")
dir.create(path,recursive=T,showWarnings=F)

##### Analyse DESeq2 ####
# Ouverture des fichiers countdata SANS correction de l'effet Batch
countdata = read.table(paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",condition ,"_expression_table_ROW.tab"), sep="\t",row.names=1,header =  T)

# Création du tableau avec les info des colonnes
infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)

# Creataion de l'objet DESeq2
countdata = as.matrix(countdata)
deseq = DESeqDataSetFromMatrix(countData = countdata,
                               colData  = infodata,
                               design   = ~ Condition)


# Analyse DESeq2
print(paste(condition , "-----> Analyse DESeq2"))
deseq = DESeq(deseq)

# Récupération des données de comptage normalisées
data_tab = counts(deseq,normalized=T)

##### Analyse multi-variée des données pour clustering  #####
# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(paste0(j,"_"), color)]=seq_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"/SANS_ComBat"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")

# Créaction du vecteur de couleur par cluster
color = colnames(data_tab)
for (j in rnai_list[[condition]]){
  color[grep(paste0(j,"_"), color)]=cluster_color[[j]]
}

# Analyse en composante principale
print(paste( condition, "-----> Analyse ACP"))
PCA_plot_generator(data_tab,colors = color,
                   save_path = paste0(path,"/SANS_ComBat_cluster"),
                   main = paste0("ACP ", condition," (DESeq2)"),
                   sortie = "png")
