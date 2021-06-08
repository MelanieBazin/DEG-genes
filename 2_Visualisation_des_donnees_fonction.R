options(stringsAsFactors = FALSE)
source("./0_Cluster.R")
library("stringr") 

###### Création/Modification tableau de comptage ####

ConcatTab <- function(type, conditions = NULL){
  annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
  path = paste0("./DATA/", type, "/")
  
  if (type == "EXPRESSION"){
    extention = ".tab"
  }else{
    extention = paste0("_expression_table_",type,".tab")
  }
  
  count = gsub(extention,"",list.files(path))
               
  if (!is.null(conditions)){
    count =  count[which(is.element(count, conditions))]
  }
  
  tab_count = matrix(annotation$ID)
  colnames(tab_count)="ID"
  for (j in count){
    tab = read.table(paste0(path,j,extention), sep = "\t", header = T, quote = "", row.names = 1)
    colnames(tab)=paste(j, colnames(tab), sep = "_")
    tab$ID = rownames(tab)
    tab_count = merge(tab_count, tab, by = "ID")
  } 
  
  if (colnames(tab_count)[1]=="ID"){
    rownames(tab_count) = tab_count$ID
    tab_count = tab_count[,-1]
  }
  return(tab_count)
}

DivideByGeneSeize <- function(countdata){
  annotation = read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.tab",header=T,sep="\t",quote='')
  taille = abs(annotation$NT_START - annotation$NT_END)
  names(taille)=annotation$ID
  
  data_tab_seize  = countdata
  for (j in rownames(data_tab_seize)){
    data_tab_seize[j,] = countdata[j,]/taille[j]
  }
  return(data_tab_seize)
}

MeanTabCalculation <- function(data_tab, rnai_list, cluster,i){
  
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
  mean_data_tab = as.matrix(mean_data_tab)
  return(mean_data_tab)
  
}

##### Création tabeau info avec les métadonnées #####
CreatInfoData1 <- function(countdata, conditions, rnai_list, cluster){
  infodata = matrix(NA,nrow = ncol(countdata), ncol = 5)
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Noms", "Feeding", "Timing","Cluster", "Conditions")
  
  infodata[,"Noms"] = colnames(countdata)

  rnai = names(timing_list)
  timing = c()
  clust = c()
  feeding = colnames(countdata)
  for(r in rnai){
    timing = c(timing,timing_list[[r]])
    clust = c(clust, cluster[[r]])

    for (g in timing_list[r]){
      feeding = str_replace_all(feeding, timing_list[r], "")
    }
  }
  
  infodata[,"Feeding"] = str_sub(feeding, end= -2)
  infodata[,"Timing"] = timing
  infodata[,"Cluster"] = clust
  
  # Colonne condition
  cond = c()
  for (c in 1:nrow(infodata)){
    cond = c(cond, paste(infodata[c,"Feeding"],infodata[c,"Cluster"], sep = "_"))
  }
  
  infodata[,"Conditions"] = cond
  
  infodata = as.data.frame(infodata)
}

CreatInfoData2 <- function(conditions=NULL){
  countdata= ConcatTab("EXPRESSION")
  if (colnames(countdata)[1]=="ID"){
    countdata=countdata[,-1]
  }
  infodata=matrix(NA,nrow = ncol(countdata), ncol=3)
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Name","RNAi","Timing")
  infodata[,"Name"] = row.names(infodata)
  
  infodata[,"Timing"] = unlist(timing_list, use.names = F)
  
  condi = names(timing_list)
  rnai = c()
  for (i in 1:length(timing_list)){
    rnai = c(rnai, rep(condi[i],length(timing_list[[i]])))
  }
  infodata[,"RNAi"] = rnai
  infodata = as.data.frame(infodata)



  if (!is.null(conditions)){
    infodata =  infodata[which(is.element(infodata$RNAi, conditions)),]
  }
  return(infodata)
}

CreatInfoData3 <- function(countdata, conditions, rnai_list, cluster){
  infodata = matrix(NA,nrow = ncol(countdata), ncol = 7)
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Noms", "Feeding", "Timing", "Cluster", "Condition","Batch","Labo")
  
  infodata[,"Noms"] = colnames(countdata)
  
  rnai = rnai_list[[conditions]]
  rnai = rnai[order(rnai)]
  
  timing = c()
  clust = c()
  batch = c()
  feeding = c()
  condition = c()
  labo = c()
  for(r in rnai){
  
    timing = c(timing,timing_list[[r]])
    clust = c(clust, cluster[[r]])
    
    if (length(grep("CTRL",r))>0 | length(grep("ND7",r))>0 | length(grep("ICL7",r))>0 ){
      feeding =c(feeding, rep("ctrl", length(timing_list[[r]])))
      condition = c(condition, paste(cluster[[r]],"ctrl",sep = "_"))
    }else{
      feeding =c(feeding, rep(r, length(timing_list[[r]])))
      condition = c(condition, paste(cluster[[r]],r,sep = "_" ))
    }

    
    if (r == "ND7_K" | r == "PGM"| r == "KU80C" | r == "ICL7"){
      batch = c(batch,rep("seq_2014", length(timing_list[[r]])))
    }else{
      batch = c(batch,rep("seq_2020",length(timing_list[[r]])))
    }
    
    if ( r == "ICL7"){
      labo = c(labo,rep("Duharcourt", length(timing_list[[r]])))
    }else{
      labo = c(labo,rep("Betermier",length(timing_list[[r]])))
    }
  }
  
  infodata[,"Feeding"] = feeding
  infodata[,"Timing"] = timing
  infodata[,"Cluster"] = clust
  infodata[,"Batch"] = batch
  infodata[,"Condition"]= condition
  infodata[,"Labo"]= labo
  
  infodata = as.data.frame(infodata)
  
return(infodata)
}

##### Analyse multi-variée #####
# ACP_plot_generator fait par Olivier
library(FactoMineR)
library("factoextra")
library(ggplot2)
library(gtools)

PCA_plot_generator <- function(Expression_Mat, colors,save_path, main,max_dim=3,barplot_max_dim=3,image_prefix="PCA_",show_barplot=T, vline=0, ...) {
  save_path = paste0(save_path,"/ACP/")
  dir.create(save_path,recursive=T,showWarnings=F)
  
  resExp = PCA(t(Expression_Mat), graph = F)
  
  if(show_barplot) {
    eigenvalues <- resExp$eig
    png(paste0(save_path,image_prefix,"_PCA_Variance.png"))
      barplot(eigenvalues[1:barplot_max_dim, 2], names.arg=1:barplot_max_dim, 
              main = "Variances",
              xlab = "Principal Components",
              ylab = "Percentage of variances",
              col ="gray")
      if(vline!=0) {
        abline(v=vline,lty=2,lwd=2)
      }
    dev.off()
  }    
  
  for (i in 1:dim(combn(1:max_dim,2))[2]) {
    
    gp<-plot.PCA(resExp, axes = combn(1:max_dim,2)[,i], habillage = "ind", col.hab = colors, title = main,
                 ggoptions = list(size=2))
    ggsave(paste0(save_path,image_prefix,i,".png"), device = "png", plot = gp)
    
  }
 return(resExp)
}

library(ggplot2)
library(ggrepel)

LDA_plot_generator <- function(type = "LDA",lda_data_tab,infodata, lda_model, path, condition, color){
  path = paste0(path,"/",type,"/")
  dir.create(path,recursive=T,showWarnings=F)
  
  png(paste0(path,condition,"_",type,"_hist.png"),width = 480, height = 1000)
  plot(lda_model, dimen = 1, type = "b")
  dev.off()
  # plot(lda_model, col = color, dimen = 2)
  prediction = predict(lda_model)$x
  
  
  gg_data_tab = cbind(as.data.frame(lda_data_tab), prediction)
  gp = ggplot(gg_data_tab, aes(LD1, LD2))+
    geom_point(size = 1, aes(color = infodata$Cluster)) +
    geom_text_repel(size = 2, max.overlaps = 30 , aes(label = row.names(lda_data_tab), colour = infodata$Cluster))+
    labs(color = "Groupe")+
    theme_light()+
    scale_color_manual(values = unique(color))
  ggsave(paste0(path,condition,"_",type,"1.png"), device = "png", plot = gp, width = 20, height = 20, units = "cm")
  
  gp = ggplot(gg_data_tab, aes(LD1, LD3))+
    geom_point(size = 1, aes(color = infodata$Cluster)) +
    geom_text_repel(size = 2, max.overlaps = 30 , aes(label = row.names(lda_data_tab), colour = infodata$Cluster))+
    labs(color = "Groupe")+
    theme_light()+
    scale_color_manual(values = unique(color))
  ggsave(paste0(path,condition,"_",type,"2.png"), device = "png", plot = gp, width = 20, height = 20, units = "cm")
  
  gp = ggplot(gg_data_tab, aes(LD3, LD2))+
    geom_point(size = 1, aes(color = infodata$Cluster)) +
    geom_text_repel(size = 2, max.overlaps = 30 , aes(label = row.names(lda_data_tab), colour = infodata$Cluster))+
    labs(color = "Groupe")+
    theme_light()+
    scale_color_manual(values = unique(color))
  ggsave(paste0(path,condition,"_",type,"3.png"), device = "png", plot = gp, width = 20, height = 20, units = "cm")
  
  # gp = ggplot(gg_data_tab, aes(LD1, LD4))+
  #   geom_point(size = 1, aes(color = infodata$Cluster)) +
  #   geom_text_repel(size = 2, max.overlaps = 30 , aes(label = row.names(lda_data_tab), colour = infodata$Cluster))+
  #   labs(color = "Groupe")+
  #   theme_light()+
  #   scale_color_manual(values = unique(color))
  # ggsave(paste0(path,condition,"_",type,"4.png"), device = "png", plot = gp, width = 20, height = 20, units = "cm")
  # 
  # gp = ggplot(gg_data_tab, aes(LD2, LD4))+
  #   geom_point(size = 1, aes(color = infodata$Cluster)) +
  #   geom_text_repel(size = 2, max.overlaps = 30 , aes(label = row.names(lda_data_tab), colour = infodata$Cluster))+
  #   labs(color = "Groupe")+
  #   theme_light()+
  #   scale_color_manual(values = unique(color))
  # ggsave(paste0(path,condition,"_",type,"5.png"), device = "png", plot = gp, width = 20, height = 20, units = "cm")
  # 
  # gp = ggplot(gg_data_tab, aes(LD3, LD4))+
  #   geom_point(size = 1, aes(color = infodata$Cluster)) +
  #   geom_text_repel(size = 2, max.overlaps = 30 , aes(label = row.names(lda_data_tab), colour = infodata$Cluster))+
  #   labs(color = "Groupe")+
  #   theme_light()+
  #   scale_color_manual(values = unique(color))
  # ggsave(paste0(path,condition,"_",type,"6.png"), device = "png", plot = gp, width = 20, height = 20, units = "cm")
  # 
}


EvaluPrediction <- function(type, data_tab, infodata , i, path){
  if (type == "LDA"){
    lda_data_tab=scale(t(data_tab))
    
    keep = c()
    for (l in 1:ncol(lda_data_tab)){
      keep = c(keep, !is.element(T, is.na(lda_data_tab[,l])))
    }
    lda_data_tab = lda_data_tab[,keep]
    
    training_sample = sample(c(T,F), nrow(lda_data_tab), replace = T, prob = c(0.75, 0.25))
    train = lda_data_tab[training_sample,]
    test = lda_data_tab[!training_sample,]
    
    lda_model = lda(train, grouping = infodata$Cluster[training_sample])
    
  }else if (type == "SVM"){
    lda_data_tab=t(data_tab)
    
    training_sample = sample(c(T,F), nrow(lda_data_tab), replace = T, prob = c(0.6, 0.4))
    train = lda_data_tab[training_sample,]
    test = lda_data_tab[!training_sample,]
    
    trctrl = trainControl(method = "repeatedcv", number = 10, repeats = 3)
    lda_model = train(train, y = as.factor(infodata$Cluster[training_sample]), method = "svmLinear",
                      trControl=trctrl,
                      preProcess = c("center", "scale"),
                      tuneLength = 10)
  }
  
  # lda_model$prior
  # summary(lda_model$scaling)
  
  lda_train = predict(lda_model)
  lda_test = predict(lda_model, test)
  
  type = paste0(type,"train")
  
  path = paste0(path,"/",type,"/")
  dir.create(path,recursive=T,showWarnings=F)
  
  LDA_plot_generator(type,train,infodata[training_sample,], lda_model, paste0(path,"Train"), i, color)
  
  info_train = infodata$Cluster[training_sample]
  info_test = infodata$Cluster[!training_sample]
  
  if(type == "LDAtrain"){
    lda_train = as.character(lda_train$class)
    lda_test  = as.character(lda_test$class)
    
  } 
  
  if (length(lda_train)>length(info_train)){
    info_train = c(info_train, rep(NA, length(lda_train)-length(info_train)))
  }else if (length(lda_train)<length(info_train)){
    lda_train = c(lda_train, rep(NA, length(info_train)-length(lda_train)))
  }
  
  if (length(lda_test)>length(info_test)){
    info_test = c(info_test, rep(NA, length(lda_test)-length(info_test)))
  }else if (length(lda_train)<length(info_train)){
    lda_test = c(lda_test, rep(NA, length(info_test)-length(lda_test)))
  }
  
  lda_train = factor(lda_train, levels = unique(infodata$Cluster))
  lda_test  = factor(lda_test, levels = unique(infodata$Cluster))
  info_train = factor(info_train, levels = unique(infodata$Cluster))
  info_test = factor(info_test, levels = unique(infodata$Cluster))
  
  conf_mat_train = confusionMatrix(data = lda_train, reference = info_train)
  conf_mat_test = confusionMatrix(data = lda_test, reference = info_test)
  
  return(print(paste("Précision du modèle : ",signif(conf_mat_test$overall["Accuracy"]*100, digits = 4),"%")))
  
  write(conf_mat_train, paste0(path,i,"_confusionMatrice_train.txt"), sep = "\t")
  write(conf_mat_test, paste0(path,i,"_confusionMatrice_test.txt"), sep = "\t")
  
  
  
  
}


##### Clustering #####

Clustering <- function(matDist, nb_cluster, method, titre, colors = NULL){
  ## Créaction d'un vecteur contenant le clusering calculé a partir de la matrice de distance
  # Choisir le type d'algorithme utilisé pour faire les clusters
  if (method  == "kmeans"){
    res = kmeans(matDist, nb_cluster)
    #Représentataion graphique
    p= fviz_cluster(res, data = matDist, geom = c("point",  "text"), labelsize = 10, repel = T, 
                    show.clust.cent = F, ellipse = T, ggtheme = theme_bw(),
                    main = paste(method, "avec", nb_cluster, "cluster - distance :", distance,"\n", titre), 
                    xlab = "Principal Component 1",
                    ylab = "Principal Component 2")
    print(p)

  }else if(method  == "HCL"){
    res = hclust(matDist)
    #Fait un dendrogramme
    p= plot(res, main = paste(method, "dendrogramme - distance :", distance,"\n", titre))
    print(p)
  }

  
}

###### Création graphique ####

# Crée des boxplot sur toues les gènes par catégories
CountBoxplot <- function (tab, type, color = "lightgray"){
  boxplot(log(tab + 1), ylab = "count values (log scale)",
          main = paste0("Count data (",type,")"), xaxt="n", yaxt="n",
          col = color, outline=FALSE)
  axis(side = 1, labels = FALSE, tick = F)
  axis(side = 2,
       ## Rotate labels perpendicular to y-axis.
       las = 2,
       ## Adjust y-axis label positions.
       mgp = c(3, 0.75, 0))
  text(x = 1:ncol(tab),
       ## Move labels to just below bottom of chart.
       y = par("usr")[3] - 0.1,
       ## Use names from the data list.
       labels = colnames(tab),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 90,
       ## Adjust the labels to almost 100% right-justified.
       adj = 1,
       ## Increase label size.
       cex = 0.5)
}

# Fait les heatmap avec tous les gènes par actégories
library("pheatmap")
library("ComplexHeatmap")
library("magick")
library("RColorBrewer")
library(circlize)
library(gplots)
MyHeatmaps <- function(path, data_tab, moyenne = F, condition, Log = T){
  dir.create(path,recursive=T,showWarnings=F)
  
  annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
  annotation = annotation[which(is.element(annotation$ID, rownames(data_tab))),]
  
  
  if (Log == T){
    data_log =  as.matrix(log(data_tab+1))
    Ylab = "log(expres)"
  }else{
    data_log =  as.matrix(data_tab)
    Ylab = "expres"
  }
  
  color_vec = c(0,quantile(data_log)[[2]],median(data_log), 
                quantile(data_log)[[4]],max(data_log))
  color_vec = colorRamp2(color_vec, 
                         c("white","#FEE0D2","#FB6A4A","#BD0026","#67000D"))
  
  c_split = c()
  c_split_ctr = c()
  c_order = c()
  c_order_ctr = c()
  for (a in rnai_list[[condition]]){
    x = grep(a, colnames(data_log))
    
    
    if (length(grep("ND7",a)) > 0 | length(grep("ICL7",a)) > 0){
      c_order_ctr = c(c_order_ctr,x)
      c_split_ctr = c(c_split_ctr,rep(a,length(x)))
    }else{
      c_order = c(c_order, x)
      c_split = c(c_split,rep(a,length(x)))
    }
  }
  c_split = c(c_split_ctr, c_split)
  c_order = c(c_order_ctr, c_order)
  
  data_log = cbind(apply(data_log, 1, max),data_log)
  data_log = data_log[order(data_log[,1], decreasing = T),]
  data_log = data_log[,is.element(colnames(data_log), colnames(data_tab))]
  data_log = data_log[,c_order]
  
  h = Heatmap(data_log,
              name = Ylab,
              col = color_vec,
              cluster_rows = F, # turn off row clustering
              cluster_columns = F, # turn off column clustering
              column_title = condition,
              show_row_names = F,
              row_order = NULL,
              row_split = order(annotation$EXPRESSION_PROFIL),
              row_title = "%s", row_title_rot = 0,
              column_split = c_split,
              column_order = NULL,
              use_raster = T)
  
  h1 = Heatmap(data_log,
              name = Ylab,
              col = colorRamp2(color_vec, c("white","#FEE0D2","#FB6A4A","#BD0026","#67000D")),
              cluster_rows = F,
              cluster_columns = F, # turn off column clustering
              column_title = condition,
              show_row_names = F,
              row_order = NULL,
              row_split = order(annotation$EXPRESSION_PROFIL),
              row_title = "%s", row_title_rot = 0,
              column_split = c_split,
              column_order = NULL,
              use_raster = F)
  
  h2 = Heatmap(data_log,
               name = Ylab,
               cluster_rows = F, # turn off row clustering
               cluster_columns = F, # turn off column clustering
               column_title = condition,
               show_row_names = F,
               row_order = order(annotation$EXPRESSION_PROFIL),
               row_split = annotation$EXPRESSION_PROFIL,
               row_title = "%s", row_title_rot = 0,
               column_split = c_split,
               column_order = NULL,
               use_raster = T)
  
  if (moyenne == T){
    moyenne = "MOYENNE"
  }else{
    moyenne = ""
  }
  
  png(paste0(path,condition,"_AllPoint_",moyenne,"heatmap_red.png"),width = 500, height = 600)
    print(h)
  dev.off()
  
  png(paste0(path,condition,"_AllPoint_",moyenne,"heatmap_red_unraster.png"),width = 500, height = 600)
    print(h1)
  dev.off()
  
  png(paste0(path,condition, "_AllPoint_",moyenne,"heatmap_blue.png"),width = 500, height = 600)
    print(h2)
  dev.off()
}

# plotGenes fait par Gaëlle
plotGenes <- function(expData, title = "", yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(0, ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(a in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[a,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}

# Dessine les profils des groupes de gènes + trace le profils moyen
ProfilsPDF <- function(save_path = paste0(path,"/profils/"), data_tab, moyenne = F, condition, Log = T){

  annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
  annotation = annotation[which(is.element(annotation$ID, rownames(data_tab))),]
  expr_profil = unique(annotation$EXPRESSION_PROFIL)[-1]
  expr_none = unique(annotation$EXPRESSION_PROFIL)[1]
  
  
  dir.create(save_path,recursive=T,showWarnings=F)
  
  if (moyenne == T){
    moyenne = "MOYENNE"
  }else{
    moyenne = ""
  }
  
  if (Log == T){
    data_log =  log(data_tab)+1
    Ylab = "log(EXPRESSION)"
  }else{
    data_log =  data_tab
    Ylab = "EXPRESSION"
  }
  
  data_log =  log(data_tab+1)
  pdf(paste0(save_path,condition,"_Profils_expression",moyenne,".pdf"))
  for (r in rnai_list[[condition]]){
    rnai = grep(r, colnames(data_tab))
    par(mfrow=c(2,3))
    for( p in expr_profil){
      id = annotation$ID[grep(p, annotation$EXPRESSION_PROFIL)]
      graph = plotGenes(data_tab[id,rnai], title = paste(r,p), yMax = max(data_tab[id,rnai]))
      
    }
    
    for( p in expr_profil){
      id = annotation$ID[grep(p, annotation$EXPRESSION_PROFIL)]
      graph = boxplot(data_log[id,rnai],main = paste(r,p) ,ylab = Ylab,
                      names = str_replace_all(colnames(data_log)[rnai],paste0(r,"_"),""))
      
    }
    
    # Graphiques des "none"
    par(mfrow=c(1,1))
    id = annotation$ID[grep("none", annotation$EXPRESSION_PROFIL)]
    graph = plotGenes(data_tab[id,rnai], title = paste(r,"none"), yMax = max(data_tab[id,rnai]))
    
    graph = boxplot(data_log[id,rnai], main = paste(r,"none"), ylab = Ylab,
                    names = str_replace_all(colnames(data_log)[rnai],paste0(r,"_"),""))
    
  }
  dev.off()
}

ProfilsPNG <- function(save_path = paste0(path,"/profils/"), data_tab, moyenne = F, condition, Log = T){

  annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
  annotation = annotation[which(is.element(annotation$ID, rownames(data_tab))),]
  expr_profil = unique(annotation$EXPRESSION_PROFIL)[-1]
  expr_none = unique(annotation$EXPRESSION_PROFIL)[1]

  dir.create(save_path,recursive=T,showWarnings=F)
  if (moyenne == T){
    moyenne = "MOYENNE"
  }else{
    moyenne = ""
  }
  if (Log == T){
    data_log =  log(data_tab)+1
    Ylab = "log(EXPRESSION)"
  }else{
    data_log =  data_tab
    Ylab = "EXPRESSION"
  }
  
  
  for (r in rnai_list[[condition]]){
    rnai = grep(r, colnames(data_tab))
    
    for( p in expr_profil){
      id = annotation$ID[grep(p, annotation$EXPRESSION_PROFIL)]
      png(paste0(save_path,condition,"_Profil_",moyenne,r,p,".png"))
        graph = plotGenes(data_tab[id,rnai], title = paste(r,p), yMax = max(data_tab[id,rnai]))
        print(graph)
      dev.off()
    }
    
    for( p in expr_profil){
      id = annotation$ID[grep(p, annotation$EXPRESSION_PROFIL)]
      png(paste0(save_path,i,"_Boxplot_",moyenne,r,p,".png"))
        graph = boxplot(data_log[id,rnai],main = paste(r,p) ,ylab = Ylab,
                        xlab = str_replace_all(colnames(data_log)[rnai],paste0(r,"_"),""))
        print(graph)
      dev.off()
    }
    
    # Graphiques des "none"
    par(mfrow=c(1,1))
    id = annotation$ID[grep("none", annotation$EXPRESSION_PROFIL)]
    png(paste0(save_path,i,"_Profils",moyenne,r,"_none.png"))
      graph = plotGenes(data_tab[id,rnai], title = paste(r,"none"), yMax = max(data_tab[id,rnai]))
      print(graph)
    dev.off()
    png(paste0(save_path,i,"_Boxplot",moyenne,r,"_none.png"))
      graph = boxplot(data_log[id,rnai], main = paste(r,"none"), ylab = Ylab,
                      names = str_replace_all(colnames(data_log)[rnai],paste0(r,"_"),""))
      print(graph)
    dev.off()
  }
}

# Dessine les profils d'expression de tous les gènes pour tous les RNai demmandés
ExpressionProfils <- function(type = "RPKM", rnai = NULL){
  dir.create("./Analyse/Profils/",recursive=T,showWarnings=F)
  
  EXPRESSION = ConcatTab(type, conditions = rnai)
  MAX = apply(EXPRESSION,1, max)
  names(MAX)=rownames(EXPRESSION)
  
  if(!is.null(rnai)){
    timing_list = timing_list[rnai]
    rnai_name = paste(rnai, collapse = "_")
  }else{
    rnai_name= "all"
  }
  
  
  pdf(paste0("./Analyse/Profils/Profils",type,"_",rnai_name,".pdf"))
  for (s in 1:length(select_ID)){
    x_axis = unique(unlist(timing_list))
    x_axis_order = c(1,order(as.numeric(gsub("T", "", x_axis[2:length(x_axis)])))+1)
    x_axis = x_axis[x_axis_order]
    
    plot(NULL, xlim = c(1,length(x_axis)),
         ylim=c(0,MAX[select_ID[s]]),
         axes=F,ylab=type,
         xlab="Timing Autogamie",
         main=paste("Profils expression",names(select_ID)[s]))
    
    
    
    axis(1,at=1:length(x_axis),labels=x_axis,las=2)
    axis(2)
    
    legend("topleft",legend=names(timing_list),col=c(1:length(names(timing_list))),lwd=2, cex = 0.75,bty = "n") 
    
    color = 0
    for (i in names(timing_list)){
      color = color+1
      expression = EXPRESSION[select_ID[s],grep(i, colnames(EXPRESSION))]
      positions = c()
      for (g in gsub(paste0(i, "_"),"", colnames(expression))){
        positions = c(positions,which(g == x_axis))
      }
      
      lines(positions,expression,col=color,lwd=2)
      
    } 
    
  }
  dev.off()
}
