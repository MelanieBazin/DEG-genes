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

# Récupéreation des fonction d'Olivier pour les analyse des gènes dérégulés
source("3_Functions_AnalyeDESeq2.R")


analyseName = paste0("Analyse_DESeq2_test03_tout_batch")

path_dir = paste0("./Analyse/",analyseName,"/")
dir.create(path_dir,recursive=T,showWarnings=F)

##### Pour analyse DESeq2 uniquement #####
# Definiton des variables DESeq2
FC = 1.5 #Mini 1.5 -> XRCC4 = 2
pvalue = 0.05 #Maxi 0.05 -> XRCC4 = 0.01


### Création  de liste de gènes filtrés (retirés de l'analyse) ###
Filtering= list()
Filtering= NULL
# par exemple : retirer les gènes qui sont DE pendant une manip de silencing,
# ou entre plusieurs manip control ==> des faux positifs


### Vecteur de couleur pour les heatmap
hmcol = colorRampPalette(brewer.pal(10,"RdBu"))(255)
#hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
hmcol = rev(hmcol)

### Limiter le fichier annotation aux gènes avec synonyme ###
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
annotation_synonyms = annotation[-grep("PTET",annotation$NAME),]
# annotation_synonyms = annotation[annotation$SYNONYMS != "",]
rownames(annotation)=annotation$ID

condition = names(rnai_list)[1]

# for (condition in names(rnai_list)){
  print(paste("On analyse le jeu de donnee :", condition , "-->", paste(rnai_list[[condition]], collapse = ", ") ))
  
  path = paste0(path_dir,condition ,"/")
  dir.create(path,recursive=T,showWarnings=F)
  
  ##### Analyse DESeq2 ####
  # Ouverture des fichiers countdata sans correction de l'effet Batch
  countdata = read.table(paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",condition ,"_expression_table_ROW.tab"), sep="\t",row.names=1,header =  T)
  
  # Boxplot des comptages avant normalisation #
  print(paste(condition , "-----> Creation BoxPlot non-normalise"))
  # pdf(paste0(path,condition ,"_row.pdf"))
  png(paste0(path,condition ,"_row.png"))
    CountBoxplot(countdata, "row", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
  dev.off()
  
  # Ouverture des fichiers countdata avec correction de l'effet Batch
  countdata = read.table(paste0("./DATA/Pour_DESeq/",condition ,"_expression_table_pour_DESeq_v2.tab"), sep="\t",row.names=1,header =  T)
  
  # Création du tableau avec les info des colonnes
  infodata = CreatInfoData3(countdata, conditions = condition , rnai_list, cluster)

  # Créataion de l'objet DESeq2
  countdata = as.matrix(countdata)
  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata,
                                 design   = ~ Condition)
  
  # Definition des réplicats techniques
  ddsCollapsed <- collapseReplicates(deseq, 
                                     groupby = deseq$Echantillion, 
                                     run     = deseq$Noms)
  
  
  # Analyse DESeq2
  print(paste(condition , "-----> Analyse DESeq2"))
  deseq = DESeq(deseq)
  
  # Graphique du paramètre de dispersion
  # pdf(paste0(path,condition ,"_dipression_DESeq2.pdf"))
  png(paste0(path,condition ,"_dipression_DESeq2.png"))
  plotDispEsts(deseq, ylim = c(1e-6, 1e1))
  dev.off()

  
  # Récupération des données de comptage normalisées
  data_tab = counts(deseq,normalized=T)
  
  write.table(data_tab,paste0(path,condition ,"_expression_table_normaliserDESeq2.tab"), sep="\t",row.names=T,quote=F)

  mean_data_tab = MeanTabCalculation(data_tab, rnai_list, cluster,condition ) #a supprimer si source réactivé
  
  #### Lancer l'analyse de gènes dérégulés ####

  # Definir les condition à analyser
  RNAi_list = unique(rnai_list[[condition ]][-grep("ND7",rnai_list[[condition ]])])
  if (is.element("ICL7", RNAi_list)){
    RNAi_list = RNAi_list[-grep("ICL7",RNAi_list)]
  }
  RNAi = RNAi_list[1]
  # Regarder la derégulation dnas chaque condition
  for (RNAi in RNAi_list){
    print(paste(condition, ": Analyse des donnees pour ", RNAi ))
    
    #### Création des dossier pour ranger les données ####
    base_img_dir=paste0("./Analyse/",analyseName,"/",condition,"/DESeq/",RNAi,"/Images/")
    dir.create(base_img_dir,recursive=T,showWarnings=F)
    
    base_res_dir=paste0("./Analyse/",analyseName,"/",condition,"/DESeq/",RNAi,"/")
    dir.create(base_res_dir, recursive=T,showWarnings=F)
    
    
    source("3_Analyse_DESeq2.R")
    
  }
  
  #### Lancer les visulalisation des données ####
  data_tab = read.table(paste0("./Analyse/Data_DESeq2_toutBatch/tout/tout_expression_table_normaliserDESeq2.tab"), sep="\t", header = T, row.names = 1)
  data_tab = as.matrix(data_tab)
  path = "./Analyse/Data_DESeq2_toutBatch/tout_Rstudio/"
  dir.create(path,recursive=T,showWarnings=F)
  source("3_Visualisation_des_donnees_new.R")
  # print(paste("Visualisation des donnee fini pour", condition ))
  
  

# }