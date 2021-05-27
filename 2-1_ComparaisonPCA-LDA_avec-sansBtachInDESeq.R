source("0_Cluster.R")
source("2_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")
library(stringr)
library(sva)
library(DESeq2)
library(limma)
library(MASS)
library(ggplot2)
library(klaR)

path = "./Analyse/DESeq2_test02/"
# Utilisation uniquement des normalisations DESeq2
type = "DESeq2"
# RNAi à analyser ensemble
tout = sub("_expression_table_RPKM.tab","",list.files("./DATA/RPKM/"))
rnai_list = list(
  tout = tout,
  controles = tout[which(is.element(tout,c("ND7","ICL7","CTIP_CTRL","XRCC4_CTRL" )))]
  # ,
  # sequencage_2014 = tout[which(is.element(tout,c("ICL7","KU80c","ND7","PGM" )))],
  # controles_2014 = tout[which(is.element(tout,c("ND7", "ICL7")))],
  # sequencage_2014bis = tout[which(is.element(tout,c("KU80c","ND7","PGM" )))],
  # sequencage_2020 = tout[which(is.element(tout,c("CTIP","CTIP_CTRL","XRCC4","XRCC4_CTRL")))],
  # #controles_2020 = tout[which(is.element(tout,c( "CTIP_CTRL","XRCC4_CTRL")))],
  # XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
  # CTIPseul = tout[which(is.element(tout, c("CTIP","CTIP_CTRL")))],
  # CTIPseulctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL", "XRCC4_CTRL")))]
)

i = "tout"

# for (i in names(rnai_list)){

  ##### Création du tableau de donnée à analyser ensemble ####
  #Ouverture des fichiers et création de l'objet countdata
  countdata = ConcatTab("EXPRESSION", conditions = rnai_list[[i]])
  row.names(countdata)=countdata$ID
  countdata=countdata[,-1]
  
  # Boxplot des comptages non-normalisés
  png(paste0(path,i,"_Row_Boxplot.png"))
    CountBoxplot(countdata, "DESeq2", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
  dev.off()
  
  # Création du tableau avec les info des colonnes
  infodata=CreatInfoData3(countdata, conditions = i, rnai_list, cluster)  # La verison 3 ajoute une colonne batch avec les année de séquancage
  
  # Mise en forme des données
  countdata =  as.matrix(countdata)
  countdata = countdata[rowSums(countdata) > 50,]
  
  # Correction de l'effet batch avec ComBat
  countdata = ComBat_seq(countdata, batch = infodata$Batch, group = infodata$Cluster)

  deseq = DESeqDataSetFromMatrix(countData = countdata,
                                 colData  = infodata[,1:(ncol(infodata)-1)],
                                 design   = ~ Feeding + Cluster)
  
  
  # Analyse DESeq2
  deseq = DESeq(deseq)

  
  # Graphique du paramètre de dispersion
  png(paste0(path,i,"_dipression_DESeq2.png"))
    plotDispEsts(deseq, ylim = c(1e-6, 1e1))
  dev.off()
  
  
  # Récupération des données de comptage normalisées
  data_tab=counts(deseq,normalized=T)
  
  # Passage de la colonne des ID en rowname
  if (colnames(data_tab)[1]=="ID"){
    row.names(data_tab)=data_tab$ID
    data_tab = data_tab[,-1]
  }

  # Créaction du vecteur de couleur par cluster
  color = colnames(data_tab)
  for (j in rnai_list[[i]]){
    color[grep(j, color)]=cluster_color[[j]]
  }


  # Changer le nom des colonnes controles
  colnames(data_tab) = str_replace_all(colnames(data_tab),"ND7","ND7_K")
  colnames(data_tab) = str_replace_all(colnames(data_tab),"CTIP_CTRL","ND7_C")
  colnames(data_tab) = str_replace_all(colnames(data_tab),"XRCC4_CTRL","ND7_X")

  # Analyse en composante principale
  print(paste(i, "-----> début ACP"))
  PCA_plot_generator(data_tab,colors = color,
                            save_path = paste0(path,"ACP/DESeq2_",i,"_"),
                            main = paste0("ACP ",i," (",type,")"))
  
  # Analyse de discrimination linéaire (LDA)
  #methode : https://www.towardsdatasciences.com/linear-discriminant_analysis-lda-101-using-r-6a97217a55a6
  #methode : https://www.statology.org/linera-discriminant-analysis-in-r/
  print(paste(i, "-----> début LDA"))
  set.seed(101)
  lda_data_tab=scale(t(data_tab)) #s'assurer que la sd est de 1 et la moyenne à 0 (prédicat des lda)
  summary(apply(lda_data_tab,2,mean)) #verification que la moyenne est à 0 ou très proche
  summary(apply(lda_data_tab,2,sd)) #verification que la sd est à 1
  lda_model = lda(lda_data_tab, grouping = infodata$Cluster)
  lda_model$prior
  summary(lda_model$scaling)

 
  # Evaluer les prédiction du modèle
  lda_pred_train = predict(lda_model)
  table(lda_pred_train$class, infodata$Cluster) #les chiffres sur la diagonal correspondent au rédiction correcte
  
  names(color) = rownames(infodata$Cluster)
  
  png("test.png",width = 480, height = 1000)
    plot(lda_model, dimen = 1, type = "b")
  dev.off()
  plot(lda_model, col = color, dimen = 2)

  ### Installer le placake klaR
  # partimat(paste("Cluster ~", paste(unique(infodata$Cluster), collapse = "+")), data = as.data.frame(lda_train), method = "lda")
  
  
  
  gg_data_tab = cbind(as.data.frame(lda_data_tab), predict(lda_model)$x)
  gp = ggplot(gg_data_tab, aes(LD1, LD2))+
    geom_point(size = 1.5, aes(color = infodata$Cluster)) + 
    # scale_color_manual(labels = names(cluster)) +
    geom_text(size = 3, check_overlap = F, nudge_y = 0.15 , aes(label = row.names(lda_data_tab)))+
    labs(color = "Groupe")
    
  ggsave(paste0("test2.png"), device = "png", plot = gp, width = 20, height = 20, units = "cm")

  #### QDA
  qda_data_tab = t(data_tab)
  summary(apply(qda_data_tab,2,mean)) #verification que la moyenne est à 0 ou très proche
  summary(apply(qda_data_tab,2,sd)) #verification que la sd est à 1
  qda_model = qda(qda_data_tab, grouping = infodata$Cluster)
  
  
  
  # Analyse de clusering
  # for (distance in c("Pearson", "Spearman")){
  #   png(paste0(path,"4Cluster/",i,"_Matrice_",distance,".png"))
  #   # Choisir le mode de calcule des distances
  #   if (distance == "Pearson"){
  #     matDist = as.matrix(cor(data_tab))
  #     pheatmap(matDist, main = paste("Pheatmap Pearson", type, i))
  #     matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
  # 
  #   }else if (distance == "Spearman"){
  #     matDist = as.matrix(cor(data_tab,method="spearman"))
  #     pheatmap(matDist, main = paste("Pheatmap Spearman", type, i))
  #     matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
  #   }
  #   dev.off()
  # 
  # 
  #   for (method in c("kmeans", "HCL")){
  #     print(paste(type,i, "----->",distance, method))
  #     png(paste0(path,"4Cluster/",i,"_Cluster_",method,"_",distance,".png"))
  #     Clustering(matDist = matDist,
  #                nb_cluster = 5,
  #                method = method,
  #                titre = paste(type,i))
  #     dev.off()
    # }
  
  # Boxplot des comptages normalisés divisé par la taille des gènes
  count.seize = DivideByGeneSeize(tab) #NEcessaire que la 1er fois
  write.table(count.seize,paste0("./DATA/DESeq2-seize/",i,"_expression_table_DESEQsurseize.tab"), sep="\t",row.names=F,quote=F)
  
  count.seize = read.table(patse0("./DATA/DESeq2-seize/",i,"_expression_table_DESEQsurseize.tab"), header = T, sep = "\t")
  png(paste0(path,i,"_DESeq-seize_Boxplot.png"))
  CountBoxplot(count.seize, "DESeq2", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
  dev.off()
  
  print("Boxplot terminé")
  }

}