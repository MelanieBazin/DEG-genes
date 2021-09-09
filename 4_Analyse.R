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


analyseName = paste0("Analyse_DESeq2_tout_CombatON")

#### Definiton des variables DESeq2 ####
FC = 1.5 #Mini 1.5 -> XRCC4 = 2
pvalue = 0.05 #Maxi 0.05 -> XRCC4 = 0.01

analyseName = paste0(Sys.Date(),"_", analyseName, "_FC-", FC, "_pval-", pvalue)

path_dir = paste0("./Analyse/",analyseName,"/")
dir.create(path_dir,recursive=T,showWarnings=F)

#### Limiter le fichier annotation aux gènes avec synonyme ####
annotation_basic = read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.tab",header=T,sep="\t",quote='')
my_annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
annotation = merge(annotation_basic[,c(1,4:7)], my_annotation[,c(1,2,3:7)], by = "ID")[,c(1,7,8,2:6,9:11)]
rm(annotation_basic,my_annotation)
annotation_synonyms = annotation[-grep("PTET",annotation$NAME),]
rownames(annotation)=annotation$ID

### Création  de liste de gènes filtrés (retirés de l'analyse) ###
# par exemple : retirer les gènes qui sont DE pendant une manip de silencing,
# ou entre plusieurs manip control ==> des faux positifs
file_name = "2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05"
condition = "tout"
Ctip_Ku_Pgm_Xrcc4 = read.table(paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection.tab"), sep = "\t", h=T )$ID
Ku_Pgm_Xrcc4 = read.table(paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection_UP.tab"), sep = "\t", h=T )$ID
Ctip_UP = read.table(paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection_CTIP_DOWN-et-UP.tab"), sep = "\t", h=T )$ID
TurboPGM = read.table("./DATA/TurboID/2114003-Pgm-ProteinMeasurements.txt",header=T,sep="\t")$PROTEIN_NAME
TurboPGML4 = read.table("./DATA/TurboID/2114003-PgmL4-ProteinMeasurements.txt",header=T,sep="\t",quote='')$PROTEIN_NAME
# 
# Filtering= list(
#   connu = annotation$ID[which(!is.element(annotation$ID,select_ID))],
#   Ctip_Ku_Pgm_Xrcc4 = annotation$ID[which(!is.element(annotation$ID,Ctip_Ku_Pgm_Xrcc4))],
#   Ku_Pgm_Xrcc4 = annotation$ID[which(!is.element(annotation$ID,Ku_Pgm_Xrcc4))],
#   Ctip_UP = annotation$ID[which(!is.element(annotation$ID,Ctip_UP))],
#   TurboPGM = annotation$ID[which(!is.element(annotation$PROTEIN_NAME,TurboPGM))],
#   TurboPGML4 = annotation$ID[which(!is.element(annotation$PROTEIN_NAME,TurboPGML4))]
# )
Filtering= NULL


#### Vecteur de couleur pour les heatmap ####
hmcol = colorRampPalette(brewer.pal(10,"RdBu"))(255)
#hmcol = colorRampPalette(brewer.pal(9,"GnBu"))(100)
hmcol = rev(hmcol)

condition = names(rnai_list)[1]
# for (condition in names(rnai_list)){
  print(paste("On analyse le jeu de donnee :", condition , "-->", paste(rnai_list[[condition]], collapse = ", ") ))
  
  path = paste0(path_dir,condition ,"/")
  dir.create(path,recursive=T,showWarnings=F)
  dir.create(paste0(path_dir,condition ,"/Visualisation/"),recursive=T,showWarnings=F)
  
  ##### Analyse DESeq2 ####
  # Ouverture des fichiers countdata sans correction de l'effet Batch
  countdata = read.table(paste0("./DATA/Pour_DESeq_SansCorrectionBatch/",condition ,"_expression_table_ROW.tab"), sep="\t",row.names=1,header =  T)
  
  # Boxplot des comptages avant normalisation #
  print(paste(condition , "-----> Creation BoxPlot non-normalise"))
  # pdf(paste0(path,condition ,"_row.pdf"))
  png(paste0(path_dir,condition ,"/Visualisation/Comptage_bolxplot_row.png"))
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
  deseq = collapseReplicates(deseq, 
                             groupby = deseq$Echantillion, 
                             run     = deseq$Noms)

  write.table(as.data.frame(colData(deseq)),paste0(path,condition ,"_infodata_collapse.tab"), sep="\t",row.names=T,quote=F)
  
  # Analyse DESeq2
  print(paste(condition , "-----> Analyse DESeq2"))
  deseq = DESeq(deseq)
  
  
  # Graphique du paramètre de dispersion
  # pdf(paste0(path,condition ,"_dipression_DESeq2.pdf"))
  png(paste0(path_dir,condition ,"/Visualisation/Dipression_DESeq2.png"))
  plotDispEsts(deseq, ylim = c(1e-6, 1e1))
  dev.off()

  
  # Récupération des données de comptage normalisées
  data_tab = counts(deseq,normalized=T)
  
  write.table(data_tab,paste0(path,condition ,"_expression_table_normaliserDESeq2.tab"), sep="\t",row.names=T,quote=F)
  # data_tab = read.table(paste0(path,condition ,"_expression_table_normaliserDESeq2.tab"), row.names = 1, sep="\t", header = T)

  mean_data_tab = MeanTabCalculation(data_tab, rnai_list, cluster,condition ) #necessaire pour les heatmap
  
  #### Lancer l'analyse de gènes dérégulés ####

  # Definir les condition à analyser
  RNAi_list = unique(rnai_list[[condition ]][-grep("ND7",rnai_list[[condition ]])])
  if (is.element("ICL7", RNAi_list)){
    RNAi_list = RNAi_list[-grep("ICL7",RNAi_list)]
  }
  RNAi = RNAi_list[1]
  # Regarder la derégulation dnas chaque condition
  for (RNAi in RNAi_list){
    print(paste(condition, ": Analyse des donnees pour  ----->", RNAi ))
    
    #### Création des dossier pour ranger les données ####
    base_img_dir=paste0("./Analyse/",analyseName,"/",condition,"/DESeq/",RNAi,"/Images/")
    dir.create(base_img_dir,recursive=T,showWarnings=F)
    
    base_res_dir=paste0("./Analyse/",analyseName,"/",condition,"/DESeq/",RNAi,"/")
    dir.create(base_res_dir, recursive=T,showWarnings=F)
    
    
    source("3_Analyse_DESeq2.R")
    
    print(paste(condition, ": Analyse des donnee fini pour ----->",  RNAi))
  }
  
  #### Lancer les visulalisation des données ####
  # data_tab = read.table(paste0("./Analyse/Data_DESeq2_toutBatch/tout/tout_expression_table_normaliserDESeq2.tab"), sep="\t", header = T, row.names = 1)
  data_tab = counts(deseq,normalized=T)
  data_tab = as.matrix(data_tab)
  path = paste0(path_dir,condition ,"/Visualisation/")
  
  
  source("3_Visualisation_des_donnees_new.R")

  
  print(paste("Visualisation des donnee fini pour", condition ))
# }