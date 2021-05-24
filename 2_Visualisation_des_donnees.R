# install.packages("stringr") 
# install.packages("pheatmap")
# library(devtools)
# install_github("jokergoo/ComplexHeatmap") #https://jokergoo.github.io/ComplexHeatmap-reference/book/index.html
# install.packages("magick")
# install.packages("viridis")

source("2_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")
source("0_Cluster.R")
library("stringr")  
library("pheatmap")
library("ComplexHeatmap")
library("magick")
library("RColorBrewer")
library(circlize)
# library("viridis")

#######################################################
# Partie 1 : Comparaison des méthodes de normalisaion #
#######################################################
# path = "./Graph/Normalisations/"
# Comptage

# Type = c("EXPRESSION","RPM", "RPKM")
# for (i in Type){
#   png(paste0(path,i,"_Boxplot.png"), width = 500, height = 500)
#     tab = ConcatTab(i)
#     rownames(tab)=tab$ID
#     tab = tab[,-1]
#     tab = tab[,c(12:ncol(tab),1:11)]
#     if (i == "EXPRESSION"){type = "raw data"}else{type = i}
#     
#   CountBoxplot(tab, type, color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
#   dev.off()
# }
# 
# # Récupération des données de comptage normalisées DESeq2
# tab = list.files("./DATA/DESeq2/")[grep("tout",list.files("./DATA/DESeq2/"))]
# tab = tab[grep("normalisation", tab)]
# tab = read.table(paste0("./DATA/DESeq2/", tab), header = T, sep = "\t")
# tab = tab[,c(12:ncol(tab),1:11)]
# 
# png(paste0(path,"DESeq2_Boxplot.png"),  width = 500, height = 500)
#   CountBoxplot(tab, "DESeq2", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21)))
# dev.off()

########################################################
# Partie 2 : Visualisation des données pour clustering #
########################################################
# path = "./Graph/Normalisations/"
# Type = c("RPKM", "RPM","DESeq2")
# # RNAi à analyser ensemble
# tout = sub(paste0("_expression_table_RPKM.tab"),"",list.files("./DATA/RPKM/"))
# rnai_list = list(
#   sequencage_2014 = tout[which(is.element(tout,c("ICL7","KU80c","ND7","PGM" )))],
#   sequencage_2020 = tout[which(is.element(tout,c("CTIP","CTIP_CTRL","XRCC4","XRCC4_CTRL")))],
#   controles = tout[which(is.element(tout,c("ND7", "ICL7", "CTIP_CTRL","XRCC4_CTRL")))],
#   XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
#   XRCC4ctrl2020 = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL","CTIP_CTRL")))],
#   XRCC4tousctrl = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL","ND7", "ICL7", "CTIP_CTRL")))],
#   XRCC4xseq2014 = tout[which(is.element(tout,c("XRCC4","XRCC4_CTRL","ICL7","KU80c","ND7","PGM" )))],
#   CTIPseul = tout[which(is.element(tout, c("CTIP","CTIP_CTRL")))],
#   CTIPctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL","XRCC4_CTRL")))],
#   CTIPtousctrl = tout[which(is.element(tout,c("CTIP","ND7", "ICL7", "CTIP_CTRL","XRCC4_CTRL")))],
#   tout = tout
# )


# pdf(paste0(path,"4Cluster/Tout.pdf"))
# for (type in Type){
#   for (i in names(rnai_list)){
#     # Création du tableau de donnée à analyser ensemble
#     if (type == "DESeq2"){
#       data_tab = read.table(paste0("./DATA/DESeq2/",i,"_",paste(rnai_list[[i]],collapse = "-"),"_normalisation_DESeq2.tab"), header = T, sep="\t")
#     }else{
#       data_tab = ConcatTab("EXPRESSION", conditions = rnai_list[[i]])
#     }
#     
#     # Passage de la colonne des ID en rowname
#     if (colnames(data_tab)[1]=="ID"){
#       row.names(data_tab)=data_tab$ID
#       data_tab = data_tab[,-1]
#     }
#     
#     print(paste(type,i, "- Commencé"))
#     
#     # Analyse en composante principale
#     # print(paste(type,i, "-----> début ACP"))
#     # resExp=PCA_plot_generator(data_tab,colors = NULL, save_path = paste0(path,"/ACP/",type,"_",i,"_"), main = paste0("ACP ",i," (",type,")"))
#     # 
#     # resHCPC = HCPC(resExp, graph = FALSE)
#     # 
#     # png(paste0(path,"/ACP/",type,"_",i,"_PCA_ClusteringHCPC_dendog.png"))
#     # p = fviz_dend(resHCPC, 
#     #               cex = 0.5,                     # Taille du text
#     #               palette = "jco",               # Palette de couleur ?ggpubr::ggpar
#     #               rect = TRUE, rect_fill = TRUE, # Rectangle autour des groupes
#     #               rect_border = "jco",           # Couleur du rectangle
#     #               main = paste0("ACP Dendrogram ",i," (",type,")"), 
#     #               labels_track_height = 2000)   # Augment l'espace pour le texte
#     # print(p)
#     # dev.off()
#     # 
#     # png(paste0(path,"/ACP/",type,"_",i,"_PCA_ClusteringHCPC_clust.png"))
#     # p=fviz_cluster(resHCPC,
#     #                repel = TRUE,            # Evite le chevauchement des textes
#     #                show.clust.cent = F,     # Montre le centre des clusters
#     #                palette = "jco",         # Palette de couleurs, voir ?ggpubr::ggpar
#     #                ggtheme = theme_bw(),
#     #                main = paste0("ACP Factor map ",i," (",type,")"))
#     # print(p)
#     # dev.off()
#     # print(paste(type,i, "-----> fin ACP"))
#     
#     
#     # Analyse de clusering
#     # Choisir le mode de calcule des distances
# if (distance == "Pearson"){
#   matDist = as.matrix(cor(data_tab))
#   pheatmap(matDist, main = paste("Pheatmap Pearson", type, i))
#   matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
#   
# }else if (distance == "Spearman"){
#   matDist = as.matrix(cor(data_tab,method="spearman"))
#   pheatmap(matDist, main = paste("Pheatmap Spearman", type, i))
#   matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
# }
# 
# 
# 
# for (method in c("kmeans", "HCL")){
#   print(paste(type,i, "----->",distance, method))
#   Clustering(matDist = matDist,
#              nb_cluster = 5,
#              method = method,
#              titre = paste(type,i),
#              colors = color)
#        
#       }
#     }
#     
#     print(paste(type,i," - Fait"))
#   }
# }
# dev.off()
# print("-FIN-")

#######################################################
# Partie 3 : Analyse ciblé et coloration de graphique #
#######################################################
# path = "./Graph/Normalisations/Color/"
# # Utilisation uniquement des normalisations DESeq2
# type = "DESeq2"
# # RNAi à analyser ensemble
# tout = sub("_expression_table_RPKM.tab","",list.files("./DATA/RPKM/"))
# rnai_list = list(
#   tout = tout,
#   sequencage_2014 = tout[which(is.element(tout,c("ICL7","KU80c","ND7","PGM" )))],
#   controles_2014 = tout[which(is.element(tout,c("ND7", "ICL7")))],
#   sequencage_2014bis = tout[which(is.element(tout,c("KU80c","ND7","PGM" )))],
#   sequencage_2020 = tout[which(is.element(tout,c("CTIP","CTIP_CTRL","XRCC4","XRCC4_CTRL")))],
#   #controles_2020 = tout[which(is.element(tout,c( "CTIP_CTRL","XRCC4_CTRL")))],
#   XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
#   CTIPseul = tout[which(is.element(tout, c("CTIP","CTIP_CTRL")))],
#   CTIPseulctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL", "XRCC4_CTRL")))]
# )
# 
# for (i in names(rnai_list)){
#   
#   # Création du tableau de donnée à analyser ensemble
#   data_tab = read.table(paste0("./DATA/DESeq2/",i,"_normalisation_DESeq2.tab"), header = T, sep="\t")
# 
#   # Passage de la colonne des ID en rowname
#   if (colnames(data_tab)[1]=="ID"){
#     row.names(data_tab)=data_tab$ID
#     data_tab = data_tab[,-1]
#   }
# 
#   # Créaction du vecteur de couleur par cluster
#   color = colnames(data_tab)
#   for (j in rnai_list[[i]]){
#     color[grep(j, color)]=cluster_color[[j]]
#   }
#   
#   print(paste(i, "- Commencé"))
#   
#   # Changer le nom des colonnes controles
#   colnames(data_tab) = str_replace_all(colnames(data_tab),"ND7","ND7_K")
#   colnames(data_tab) = str_replace_all(colnames(data_tab),"CTIP_CTRL","ND7_C")
#   colnames(data_tab) = str_replace_all(colnames(data_tab),"XRCC4_CTRL","ND7_X")
#   
#   # Analyse en composante principale
#   print(paste(i, "-----> début ACP"))
#   PCA_plot_generator(data_tab,colors = color,
#                             save_path = paste0(path,"/ACP/DESeq2_",i,"_"),
#                             main = paste0("ACP ",i," (",type,")"))
# 
# 
#   # Analyse de clusering
#   for (distance in c("Pearson", "Spearman")){
#     png(paste0(path,"4Cluster/",i,"_Matrice_",distance,".png"))
#     # Choisir le mode de calcule des distances
#     if (distance == "Pearson"){
#       matDist = as.matrix(cor(data_tab))
#       pheatmap(matDist, main = paste("Pheatmap Pearson", type, i))
#       matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
# 
#     }else if (distance == "Spearman"){
#       matDist = as.matrix(cor(data_tab,method="spearman"))
#       pheatmap(matDist, main = paste("Pheatmap Spearman", type, i))
#       matDist = as.dist(1-cor(log2(data_tab+1), method="spearman"))
#     }
#     dev.off()
# 
# 
#     for (method in c("kmeans", "HCL")){
#       print(paste(type,i, "----->",distance, method))
#       png(paste0(path,"4Cluster/",i,"_Cluster_",method,"_",distance,".png"))
#       Clustering(matDist = matDist,
#                  nb_cluster = 5,
#                  method = method,
#                  titre = paste(type,i))
#       dev.off()
#     }
#   }
# 
# }

#######################################################
# Partie 4 : Heatmap après moyenne de points groupés  #
#######################################################
path = "./Graph/Heatmap/"
annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")
annotation = annotation[,c(1,4)]


# RNAi à analyser ensemble
tout = sub("_expression_table_RPKM.tab","",list.files("./DATA/RPKM/"))
rnai_list = list(
  sequencage_2014bis = tout[which(is.element(tout,c("KU80c","ND7","PGM" )))],
  XRCC4seul = tout[which(is.element(tout, c("XRCC4","XRCC4_CTRL")))],
  CTIPseulctrl2020 = tout[which(is.element(tout, c("CTIP","CTIP_CTRL", "XRCC4_CTRL")))]
)

 i = "XRCC4seul"
#for (i in names(rnai_list)){
  print(i)
  # Création du tableau de donnée à analyser ensemble
  data_tab = read.table(paste0("./DATA/DESeq2/",i,"_normalisation_DESeq2.tab"), header = T, sep="\t")
  
  # Calcule des moyennes
  mean_data_tab =data.frame(ID=rownames(data_tab))
  c_split = c()
  
  for (a in rnai_list[[i]]){
    tab = data_tab[,grep(a, colnames(data_tab))] # Récupération des colonnes correspodant a une cinétique
    colnames(tab) = cluster[[a]] # Donner le nom des groupes de points aux colonnes
    
    # Calculer la moyenne pour les points que l'on veux grouper dans la cinétique
    mean_tab =data.frame(ID=rownames(data_tab))
    for (b in unique(cluster[[a]])){
      temp = grep(b, colnames(tab))
      if(length(temp) == 1){
        mean_tab[,b] = tab[,temp]
      }else{
        mean_tab[,b] = apply(tab[,temp], 1, mean)
      }
    }
    
    c_split = c(c_split, rep(a,ncol(mean_tab)-1))
    
    mean_tab = mean_tab[,c(1,ncol(mean_tab),2:(ncol(mean_tab)-1))]
    colnames(mean_tab)[2:ncol(mean_tab)]=paste0(a,"_",colnames(mean_tab)[2:ncol(mean_tab)])
    mean_data_tab = merge.data.frame(mean_data_tab, mean_tab, by = "ID")
    rm(mean_tab, tab)
  }
  
  # Passage de la colonne des ID en rowname
  if (colnames(mean_data_tab)[1]=="ID"){
    rownames(mean_data_tab)=mean_data_tab$ID
    mean_data_tab = mean_data_tab[,-1]
  }
  
  
  # Changer le nom des colonnes controles
  colnames(mean_data_tab) = str_replace_all(colnames(mean_data_tab),"ND7","ND7_K")
  colnames(mean_data_tab) = str_replace_all(colnames(mean_data_tab),"CTIP_CTRL","ND7_C")
  colnames(mean_data_tab) = str_replace_all(colnames(mean_data_tab),"XRCC4_CTRL","ND7_X")
  
  save = mean_data_tab
  mean_data_tab = save
  mean_data_tab = as.matrix(log(mean_data_tab+1))
  
 
  
  color_vec = c(min(mean_data_tab),quantile(mean_data_tab)[[2]],quantile(mean_data_tab)[[3]],
                quantile(mean_data_tab)[[4]], max(mean_data_tab))
  c("white","white", "red","darkred","black")
  color_vec = c(min(mean_data_tab),quantile(mean_data_tab)[[2]],median(mean_data_tab), quantile(mean_data_tab)[[4]],max(mean_data_tab))
  
  # Ordonner les lignes
  if(!is.element(F,annotation$ID == rownames(mean_data_tab))) {
    print("ordre ok")
  }else{print("ATTENTION, ordre des ligne incoherant entre annotation et mean_data_tab")}
  

  h1 = Heatmap(mean_data_tab,
          name = "log(expres)",
          col = colorRamp2(color_vec, c("white","#FEE0D2","#FB6A4A","#BD0026","#67000D")),
          cluster_rows = F, # turn off row clustering
          cluster_columns = F, # turn off column clustering
          column_title = i,
          show_row_names = F,
          row_order = order(annotation$EXPRESSION_PROFIL),
          row_split = annotation$EXPRESSION_PROFIL,
          row_title = "%s", row_title_rot = 0,
          column_split = c_split,
          use_raster = T)
  png(paste0(path,i, "3_heatmap_no_order.png"),width = 400, height = 600)
  print(h1)
  dev.off()
  
  

  
  # h1bis = Heatmap(mean_data_tab,
  #              name = "log(expres)",
  #              col = colorRamp2(color_vec, c("white","#FEE0D2","#FB6A4A","#BD0026","#67000D")),
  #              cluster_rows = F, # turn off row clustering
  #              cluster_columns = F, # turn off column clustering
  #              column_title = i,
  #              show_row_names = F,
  #              row_order = order(annotation$EXPRESSION_PROFIL),
  #              row_split = annotation$EXPRESSION_PROFIL,
  #              row_title = "%s", row_title_rot = 0,
  #              column_split = c_split,
  #              use_raster = F)
  # png(paste0(path,i, "3_heatmap_no_order_sharp.png"),width = 400, height = 600)
  # print(h1bis)
  # dev.off()
  # 
  # 
  # ctrl_expression_order = mean_data_tab[,grep("ND7",colnames(mean_data_tab))]
  # if(ncol(ctrl_expression_order)==1){
  #   ctrl_expression_order = order(mean_data_tab[,grep("VEG",colnames(ctrl_expression_order))], decreasing = T)
  # }else{
  #   ctrl_expression_order = order(apply(mean_data_tab[,grep("ND7",colnames(mean_data_tab))],1, max), decreasing = T)
  # }
  # 
  # 
  # h2 = Heatmap(mean_data_tab,
  #              name = "log(expres)",
  #              col = colorRamp2(color_vec, c("white","#FEE0D2","#FB6A4A","#BD0026","#67000D")),
  #              cluster_rows = F, # turn off row clustering
  #              cluster_columns = F, # turn off column clustering
  #              column_title = i,
  #              show_row_names = F,
  #              row_order = ctrl_expression_order,
  #              row_split = annotation$EXPRESSION_PROFIL[ctrl_expression_order],
  #              row_title = "%s", row_title_rot = 0,
  #              column_split = c_split,
  #              use_raster = F)
  # 
  # png(paste0(path,i, "3_heatmap_order.png"),width = 400, height = 600)
  # print(h2)
  # dev.off()
  


#}

path = "./Graph/Heatmap/" 
for (i in names(rnai_list)){
  print(i)
  # Création du tableau de donnée à analyser ensemble
  data_tab = read.table(paste0("./DATA/DESeq2/",i,"_normalisation_DESeq2.tab"), header = T, sep="\t")
  
  colnames(data_tab)
  data_tab = as.matrix(data_tab)
  # Heatmap sans moyenne ##############
  color_vec = c(min(log(data_tab+1)),quantile(log(data_tab+1))[[2]],median(log(data_tab+1)), 
                quantile(log(data_tab+1))[[4]],max(log(data_tab+1)))
  
  if (i == "CTIPseulctrl2020"){
    c_split = c(rep("RNAi",5),rep("ND7_C",6),rep("ND7_X",5))
    c_order = colnames(data_tab)[c(5, 1:4, 11, 6:10, 16, 12:15)]
  } else if (i == "sequencage_2014bis"){
    c_split = c(rep("Ku80c",7),rep("ND7_K",7),rep("PGM",7))
    c_order = colnames(data_tab)[c(7, 1:6, 14, 8:13, 21, 15:20)]
  }else if (i == "XRCC4seul"){
    c_split = c(rep("RNAi",5),rep("ND7_X",5))
    c_order = colnames(data_tab)[c(5, 1:4, 10, 6:9)]
  }
    
  h = Heatmap(log(data_tab+1),
              name = "log(expres)",
              col = colorRamp2(color_vec, c("white","#FEE0D2","#FB6A4A","#BD0026","#67000D")),
              cluster_rows = F, # turn off row clustering
              cluster_columns = F, # turn off column clustering
              column_title = i,
              show_row_names = F,
              row_order = order(annotation$EXPRESSION_PROFIL),
              row_split = annotation$EXPRESSION_PROFIL,
              row_title = "%s", row_title_rot = 0,
              column_split = c_split,
              column_order = c_order,
              use_raster = T)
    
  png(paste0(path,i, "_noMean_heatmap_no_order.png"),width = 400, height = 600)
    print(h)
  dev.off()
}##############
  