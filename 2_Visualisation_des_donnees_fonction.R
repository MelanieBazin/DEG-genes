options(stringsAsFactors = FALSE)

DivideByGeneSeize <- function(countdata){
  annotation = read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.tab",header=T,sep="\t",quote='')
  taille = abs(annotation$NT_START - annotation$NT_END)
  names(taille)=annotation$ID
  
  count.seize  = matrix(NA,ncol = ncol(countdata), nrow = nrow(countdata))
  count.seize = as.data.frame(count.seize)
  colnames(count.seize)= colnames(countdata)
  rownames(count.seize)= annotation$ID
  for (j in rownames(count.seize)){
    count.seize[j,] = countdata[j,]/taille[j]
  }
  return(count.seize)
}


ConcatTab <- function(type, conditions = NULL){
  annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")
  path = paste0("./DATA/", type)
  count = list.files(path)
  if (type == "EXPRESSION"){
    extention = ".tab"
  }else{
    extention = paste0("_expression_table_",type,".tab")
  }
  
  if (!is.null(conditions)){
    count =  count[which(is.element(count, paste0(conditions,extention)))]
  }
  
  count = paste(path,count, sep = "/")
  tab_count = array(annotation$ID, dim = c(nrow(annotation),1))
  colnames(tab_count)="ID"
  for (i in count){
    tab = read.table(i, sep = "\t", header = T)
    name = sub(".tab","",sub(paste0("./DATA/",type,"/"),"",sub(paste0("_expression_table_",type),"",i)))
    colnames(tab)[2:ncol(tab)]=paste(name, colnames(tab)[2:ncol(tab)], sep = "_")
    tab_count = merge(tab_count, tab, by = "ID")
  } 
  return(tab_count)
}

CountBoxplot <- function (tab, type, color = "lightgray"){
  boxplot(log(tab + 1), ylab = "count values (log scale)",
          main = paste0("Count data (",type,")"), xaxt="n", yaxt="n",
          col = color)
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


library(stringr)
CreatInfoData1 <- function(countdata, conditions, rnai_list, cluster){
  infodata = matrix(NA,nrow = ncol(countdata), ncol = 5)
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Noms", "Feeding", "Timing","Cluster", "Conditions")
  
  infodata[,"Noms"] = colnames(countdata)
  
  # Colonne feeding, timing et cluster
  CTIP = c("T0", "T5.5", "T12.5", "T25", "Veg")
  CTIP_CTRL = c("T0", "T5", "T10", "T20", "T30", "Veg")
  ICL7 = c("T0", "T5", "T10", "T20", "T35", "T50", "Veg")
  KU80c = c( "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
  ND7 = c( "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
  PGM = c( "T2", "T5", "T10", "T20", "T30", "T40", "Veg")
  XRCC4 = c( "T2", "T7", "T22", "T32","Veg")
  XRCC4_CTRL = c( "T2", "T7", "T22", "T32","Veg")
  
  rnai = rnai_list[[conditions]]
  rnai = rnai[order(rnai)]
  
  timing = c()
  clust = c()
  feeding = colnames(countdata)
  for(r in rnai){
    t = eval(parse(text = r))
    
    timing = c(timing,t)
    clust = c(clust, cluster[[r]])

    for (g in t){
      feeding = str_replace_all(feeding, t, "")
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



  if (!is.null(conditions)){
    infodata =  infodata[which(is.element(infodata$RNAi, conditions)),]
  }
  return(infodata)
}

CreatInfoData3 <- function(countdata, conditions, rnai_list, cluster){
  infodata = matrix(NA,nrow = ncol(countdata), ncol = 6)
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Noms", "Feeding", "Timing", "Cluster", "Condition","Batch")
  
  infodata[,"Noms"] = colnames(countdata)
  
  # Colonne feeding, timing et cluster
  CTIP = c("T0", "T5.5", "T12.5", "T25", "Veg")
  CTIP_CTRL = c("T0", "T5", "T10", "T20", "T30", "Veg")
  ICL7 = c("T0", "T5", "T10", "T20", "T35", "T50", "Veg")
  KU80c = c( "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
  ND7 = c( "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
  PGM = c( "T2", "T5", "T10", "T20", "T30", "T40", "Veg")
  XRCC4 = c( "T2", "T7", "T22", "T32","Veg")
  XRCC4_CTRL = c( "T2", "T7", "T22", "T32","Veg")
  
  rnai = rnai_list[[conditions]]
  rnai = rnai[order(rnai)]
  
  timing = c()
  clust = c()
  batch = c()
  feeding = c()
  condition = c()
  for(r in rnai){
    t = eval(parse(text = r))
  
    timing = c(timing,t)
    clust = c(clust, cluster[[r]])
    
    if (length(grep("CTRL",r))>0 | length(grep("ND7",r))>0 | length(grep("ICL7",r))>0 ){
      feeding =c(feeding, rep("ctrl", length(t)))
      condition = c(condition, paste(cluster[[r]],"ctrl",sep = "_"))
    }else{
      feeding =c(feeding, rep(r, length(t)))
      condition = c(condition, paste(cluster[[r]],r,sep = "_" ))
    }
    
    
    
    if (r == "ND7" | r == "PGM"| r == "KU80C" | r == "ICL7"){
      batch = c(batch,rep("seq_2014", length(t)))
    }else{
      batch = c(batch,rep("seq_2020",length(t)))
    }
  }
  
  infodata[,"Feeding"] = feeding
  infodata[,"Timing"] = timing
  infodata[,"Cluster"] = clust
  infodata[,"Batch"] = batch
  infodata[,"Condition"]= condition
  
  infodata = as.data.frame(infodata)
  
return(infodata)
}


DA_plot_generator <- function(type,lda_data_tab,infodata, lda_model, path, condition, color){
  path = paste0(path,"/",type,"/")
  dir.create(path,recursive=T,showWarnings=F)
  
  png(paste0(path,condition,"_",type,"_hist.png"),width = 480, height = 1000)
    plot(lda_model, dimen = 1, type = "b")
  dev.off()
  # plot(lda_model, col = color, dimen = 2)
  
  ### Installer le placake klaR
  # partimat(paste("Cluster ~", paste(unique(infodata$Cluster), collapse = "+")), data = as.data.frame(lda_train), method = "lda")
  gg_data_tab = cbind(as.data.frame(lda_data_tab), predict(lda_model)$x)
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
  
}

#############################################
# PCA
#############################################
#Olivier


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


###################
# Clustering
####################


Clustering <- function(matDist, nb_cluster, method, 
                       titre, colors = NULL){
  ## Créaction d'un vecteur contennat le clusering calculé a partir de la matrice de distance
  # Choisir le type d'algorithme utilisé pour faire les clusters
  if (method  == "kmeans"){
    res = kmeans(matDist, nb_cluster)
    #Représentataion graphique
    p= fviz_cluster(res, data = matDist, geom = c("point",  "text"), labelsize = 10, repel = T, 
                    show.clust.cent = F, ellipse = T, ggtheme = theme_bw(),
                    title = paste(method, "avec", nb_cluster, "cluster - distance :", distance,"\n", titre), 
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

##########










