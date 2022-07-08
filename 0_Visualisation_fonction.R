options(stringsAsFactors = FALSE)
library("stringr")
library(dendextend)

###### Modification de tableau ####

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

MeanTabCalculation <- function(data_tab, infodata){
  
  data_tab = OrderColumn(data_tab, infodata)
  infodata = infodata[colnames(data_tab),]
  
  mean_data_tab = as.data.frame(rownames(data_tab))
  rownames(mean_data_tab) = mean_data_tab[,1]
  colnames(mean_data_tab) = "ID"
  
  for (i in unique(infodata$Condition)){
    l = grep(i, infodata$Condition)
    if (length(l) == 1){
      mean_data_tab = cbind(mean_data_tab, data_tab[,l])
    }else {
      tmp = apply(data_tab[,l], 1, mean)
      mean_data_tab = cbind(mean_data_tab, tmp)
    }
  }
  
  mean_data_tab = mean_data_tab[,-grep("ID", colnames(mean_data_tab))]
  colnames(mean_data_tab) = unique(infodata$Condition)
  
  return(mean_data_tab)
  
}

OrderColumn <- function(data_tab, infodata){
  colum_order_ctrl = c()
  colum_order_rnai = c()
  
  infodata$RNAi = str_remove_all(str_split_fixed(infodata$Names, "_T", n=2)[,1], "_Veg")
  
  for (r in unique(infodata$RNAi)){
    rnai_position = str_which(colnames(data_tab),r)
    Timing = infodata$Timing[grep(r, infodata$Names)]
    Timing = sub("Veg","-1", Timing)
    
    cluster_position = rnai_position[order(as.numeric(Timing))]
    
    
    if(infodata$Feeding[rnai_position[1]] == "ctrl"){
      colum_order_ctrl = c(colum_order_ctrl, cluster_position)
    }else{
      colum_order_rnai = c(colum_order_rnai, cluster_position)
    }
  }
  
  colum_order = c(colum_order_ctrl, colum_order_rnai)
  ordered_tab = data_tab[,colum_order]
  
  return(ordered_tab)
  
}

RenameCol_PCA<- function(data_tab){
  name = colnames(data_tab)
  name = str_replace_all(name, "Veg","V")
  name = str_replace_all(name, "ND7_C","NC")
  name = str_replace_all(name, "ND7_X","NX")
  name = str_replace_all(name, "ND7_K","NK")
  name = str_replace_all(name, "ICL7","I")
  name = str_replace_all(name, "CTIP","C")
  name = str_replace_all(name, "XRCC4","X")
  name = str_replace_all(name, "KU80c","K")
  name = str_replace_all(name, "PGM","P")
  colnames(data_tab) = name
  
  return(data_tab)
}


##### Création tabeau info avec les métadonnées #####
CreatInfoData <- function(countdata, conditions, rnai_list, cluster, Timing = NULL){
  infodata = matrix(NA,nrow = ncol(countdata), ncol = 8)
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Names","Samples", "Feeding", "Timing", "Cluster", "Condition","Seq_method","Labo")
  
  infodata[,"Names"] = colnames(countdata)
  
  rnai = rnai_list[[conditions]]
  rnai = rnai[order(rnai)]
  
  timing = c()
  clust = c()
  batch = c()
  feeding = c()
  condition = c()
  labo = c()
  for(r in rnai){
    clust = c(clust, cluster[[r]])
    
    if (length(grep("CTRL",r))>0 | length(grep("ND7",r))>0 | length(grep("ICL7",r))>0 ){
      feeding =c(feeding, rep("ctrl", length(cluster[[r]])))
      condition = c(condition, paste(cluster[[r]],"ctrl",sep = "_"))
    }else if (length(grep("EZL1",r))>0){
      feeding =c(feeding, rep("EZL1", length(cluster[[r]])))
      condition = c(condition, paste(cluster[[r]],"EZL1",sep = "_"))
    }else{
      feeding =c(feeding, rep(r, length(cluster[[r]])))
      condition = c(condition, paste(cluster[[r]],r,sep = "_" ))
    }
    
    
    if (r == "ND7_K" | r == "PGM"| r == "KU80C" | r == "ICL7" | r == "EZL1"){
      batch = c(batch,rep("HiSeq", length(cluster[[r]])))
    }else{
      batch = c(batch,rep("NextSeq",length(cluster[[r]])))
    }
    
    if ( r == "ICL7" | r == "ICL7bis" | r == "EZL1" | r == "EZL1bis"){
      labo = c(labo,rep("Duharcourt", length(cluster[[r]])))
    }else{
      labo = c(labo,rep("Betermier",length(cluster[[r]])))
    }
    
    
    
    if (is.null(Timing)){
      timing = c(timing,timing_list[[r]])
    }else {
      timing = rep(Timing, 4)
    }
  }
  
  
  infodata[,"Feeding"] = feeding
  infodata[,"Timing"] = str_remove_all(timing, "T")
  infodata[,"Cluster"] = clust
  infodata[,"Seq_method"] = batch
  infodata[,"Condition"]= condition
  infodata[,"Labo"]= labo
  infodata[,"Samples"] = str_remove_all(infodata[,"Names"],"bis")
  
  infodata = as.data.frame(infodata)
  
  return(infodata)
}

##### Analyse multi-variée #####
# ACP_plot_generator fait par Olivier
library(FactoMineR)
library("factoextra")
library(ggplot2)
library(gtools)
library(ggrepel)
library(ggbiplot)

PCA_plot_generator <- function(data_tab, colors,save_path, main,max_dim=3,barplot_max_dim=3,
                               image_prefix="PCA_",show_barplot=T, selection = NULL, vline=0, sortie = "pdf", 
                               label = c("all","none","ind","ind.sup","quali","var","quanti.sup"),police_seize = 3, ...) {
  
  dir.create(save_path,recursive=T,showWarnings=F)
  
  resExp = RenameCol_PCA(data_tab)
  resExp = PCA(t(resExp), graph = F)
  
  if(show_barplot) {
    eigenvalues <- resExp$eig
    if (sortie == "pdf") {
      pdf(paste0(save_path,image_prefix,"_PCA_Variance.pdf"))
    }else if (sortie =="pngf"){
      png(paste0(save_path,image_prefix,"_PCA_Variance.png"))
    }
    
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
    
    gp<-plot.PCA(resExp, axes = combn(1:max_dim,2)[,i], habillage = "ind", 
                 label = label, col.hab = colors, 
                 title = main,
                 ggoptions = list(size=police_seize),
                 select = selection,
                 unselect = 0.85)
    ggsave(paste0(save_path,image_prefix,i,".",sortie), device = sortie, plot = gp)
    
  }
  
  return(resExp)
}

PCA_ggplot_generator <- function(data_tab, infodata,save_path, main,max_dim=3,barplot_max_dim=3,
                               image_prefix="PCA_",show_barplot=T, selection = NULL, vline=0, sortie = "pdf", 
                               label = c("all","none","ind","ind.sup","quali","var","quanti.sup"),police_seize = 3, ...) {
  
  dir.create(save_path,recursive=T,showWarnings=F)
  
  data_tab2 = t(RenameCol_PCA(data_tab))
  resExp = PCA(data_tab2, graph = F)
  
  if(show_barplot) {
    eigenvalues <- resExp$eig
    if (sortie == "pdf") {
      pdf(paste0(save_path,image_prefix,"_PCA_Variance.pdf"))
    }else if (sortie =="pngf"){
      png(paste0(save_path,image_prefix,"_PCA_Variance.png"))
    }
    
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
  
  
  for (i in max_dim:1){
    data_tab2 = cbind(resExp$ind$coord[, i],data_tab2 )
    colnames(data_tab2)[1] = paste0("Dim",i)
  }
  
  for (i in 1:dim(combn(1:max_dim,2))[2]){
    gp = ggplot(data = as.data.frame(data_tab2), aes(x = get(paste0("Dim",combn(1:max_dim,2)[1,i])), y = get(paste0("Dim",combn(1:max_dim,2)[2,i])), 
                                                     color = infodata$Cluster, 
                                                     label = rownames(data_tab2))) + 
      geom_point(alpha = 0.8, size = 0.5) +
      scale_color_manual(values = c("VEG" = "darkorange1", "EARLY" = "deepskyblue", "INTER" = "chartreuse3", "LATE"= "red" )) +
      theme_classic() +
      geom_hline(yintercept = 0, lty = 1) +
      geom_vline(xintercept = 0, lty = 1) +
      geom_text_repel(size=police_seize) +
      theme(legend.position="none") +
      labs(title = main,
           x = paste0("Dim",combn(1:max_dim,2)[1,i]),
           y = paste0("Dim",combn(1:max_dim,2)[2,i]))
    
    gp
    
    ggsave(paste0(save_path ,image_prefix,i,".",sortie), device = sortie, plot = gp)
    
  }
  
  return(resExp)
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
library(magick)
library("RColorBrewer")
library(circlize)
library(gplots)
MyHeatmaps <- function(path, data_tab,infodata, condition, moyenne = F,color = "red", Log = T, sortie = "pdf", raster = F){
  dir.create(path,recursive=T,showWarnings=F)
  
  if (moyenne == T){
    moyenne = "MOYENNE"
  }else{
    moyenne = ""
  }
  
  if (raster == T){
    r = ""
  }else{
    r = "_unraster"
  }
  
  if (Log == T){
    data_log =  as.data.frame(log(data_tab+1))
    Ylab = "log(expres)"
  }else{
    data_log =  as.data.frame(data_tab)
    Ylab = "expres"
  }
  
  rnai = sub("Veg","",colnames(data_log)[grep("Veg", colnames(data_log), ignore.case = T)], ignore.case = T)
  if (is.element("_ctrl",rnai)){
    rnai = sub("_","",rnai)
  }else{
    str_sub(rnai,-1) = ""
  }
  
  c_split = c()
  for (a in rnai){
    x = grep(a, colnames(data_log))
    c_split = c(c_split,rep(a,length(x)))
  }
  
  # Ajout des annotation au tableau
  annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t", row.names = 1)
  data_log = merge(data_log, annotation, by = 0)
  rownames(data_log) = data_log$Row.names
  data_log= data_log[,-1]
  
  # Ordonner les lignes par valeur de log du plus grand au plus petit
  data_log = cbind(apply(data_log[,1:ncol(data_tab)], 1, max),data_log)
  colnames(data_log)[1] = "MAX"
  data_log = data_log[order(data_log$MAX, decreasing = T),]
  
  # Ordonner les lignes par profils d'expression
  data_log = data_log[order(data_log$EXPRESSION_PROFIL),]
  
  # Retirer les colonnes superflues et faire un matrice numerique pour les heatmap
  data_log_mat = data_log[,is.element(colnames(data_log), colnames(data_tab))]
  data_log_mat = as.matrix(data_log_mat)
  
  # Definition des couleurs
  if (color == "red"){
    color_vec = c(0,quantile(data_log_mat)[[2]],median(data_log_mat),
                  max(data_log_mat))
    color_vec = colorRamp2(color_vec,
                           c("#FFF7EC","#FDBB84","#EF6548","#990000"))
    
    # color_vec = brewer.pal(8, "OrRd")
  }else if (color == "blue"){
    color_vec = rev(colorRampPalette(brewer.pal(10,"RdBu"))(255))
  }
  
  
  h = Heatmap(data_log_mat,
              name = Ylab,
              column_title = condition,
              col = color_vec,
              cluster_rows = F, cluster_columns = F, # turn off  clustering
              cluster_row_slices = F, cluster_column_slices = F, # turn off the clustering on slice
              show_row_names = F,
              
              row_order = nrow(data_log):1,
              row_split = data_log$EXPRESSION_PROFIL,
              row_title_rot = 0,
              
              column_split = factor(c_split, levels = unique(c_split)),
              column_order = 1:ncol(data_log_mat),
              
              use_raster = raster #reduce the original image seize
  )
  
  
  if (sortie == "pdf"){
    pdf(paste0(path,condition,"_AllPoint_",moyenne,"heatmap_",color,r, ".pdf"),width = 800, height = 800)
    draw(h)
    dev.off()
    
    
  }else if (sortie == "pdf"){
    pdf(paste0(path,condition,"_AllPoint_",moyenne,"heatmap_",color,r, ".pdf"),width = 800, height = 800)
    draw(h)
    dev.off()
  }
  
}

MyHeatmaps.2 <- function(path, data_tab, infodata, condition){
  data_log = log(data_tab+1)
  
  rnai = sub("Veg","",colnames(data_log)[grep("Veg", colnames(data_log), ignore.case = T)], ignore.case = T)
  if (is.element("_ctrl",rnai)){
    rnai = sub("_","",rnai)
  }else{
    str_sub(rnai,-1) = ""
  }
  
  c_split = c()
  for (a in rnai){
    x = grep(a, colnames(data_log))
    c_split = c(c_split,max(x))
  }
  
  # Ajout des annotation au tableau
  annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t", row.names = 1)
  data_log = merge(data_log, annotation, by = 0)
  rownames(data_log) = data_log$Row.names
  data_log= data_log[,-1]
  
  # Ordonner les lignes par valeur de log du plus grand au plus petit
  data_log = cbind(apply(data_log[,1:ncol(data_tab)], 1, max),data_log)
  colnames(data_log)[1] = "MAX"
  data_log = data_log[order(data_log$MAX, decreasing = T),]
  
  # Ordonner les lignes par profils d'expression
  data_log = data_log[order(data_log$EXPRESSION_PROFIL),]
  
  r_split =c()
  for (e in unique(data_log$EXPRESSION_PROFIL)){
    x = grep(e, data_log$EXPRESSION_PROFIL)
    r_split = c(r_split,max(x))
  }
  
  # Retirer les colonnes superflues et faire un matrice numerique pour les heatmap
  data_log_mat = data_log[,is.element(colnames(data_log), colnames(data_tab))]
  data_log_mat = as.matrix(data_log_mat)
  
  pdf("Test.pdf")
  heatmap.2(data_log_mat,
            Rowv=F, Colv=F,
            rowsep = r_split,
            colsep = c_split,
            scale="none",
            labRow="",
            # col=hmcol,
            trace="none",
            dendrogram = "none",
            main=paste("Test"))
  dev.off()
}




# Dessine les profils d'expression de tous les gènes pour tous les RNai demmandés
library("RColorBrewer")
ExpressionProfils <- function(type , condition, file, name = NULL, select_ID = NULL, rnai = NULL, infodata = NULL){
  
  if (type  == "vst"){
    path = paste0(file,condition,"/")
    EXPRESSION = read.table(paste0(path,condition ,"_expression_table_vst.tab"))
    infodata = read.table(paste0(path,condition ,"_infodata_collapse.tab"))
    
    if (is.null(rnai)){
      rnai = rnai_list[[condition]]
      RNAi = NULL
    }else{
      RNAi = paste(rnai, collapse = "-")
    }
    rnai = rnai[-grep("bis",rnai)]
    
  }else{
    path = "./Analyse/Profils/"
    dir.create(path,recursive=T,showWarnings=F)
    EXPRESSION = ConcatTab(type, conditions = rnai)
    
  }
  path = paste0(path,"Visualisation/profils/")
  dir.create(path,recursive=T,showWarnings=F)
  
  if(!is.null(select_ID)){
    EXPRESSION = EXPRESSION[select_ID,]
    
    NAMES = c()
    for (g in select_ID){
      NAMES = c(NAMES,annotation$NAME[grep(g,annotation$ID)])
    }
    names(select_ID) = NAMES
    select = "selectedID"
  
  } else {
    select_ID = annotation$ID
    names(select_ID) = annotation$NAME
    select = "all_genes"
  
  }
  
  EXPRESSION = OrderColumn(EXPRESSION, infodata)
  infodata = infodata[colnames(EXPRESSION),]
  
  MAX = apply(EXPRESSION,1, max)
  names(MAX)=rownames(EXPRESSION)
  
  timing_list = list()
  for (r in rnai){
    timing = infodata$Timing[grep(r, infodata$Names)]
    timing_list = c(timing_list, list(timing))
  }
  names(timing_list)=rnai
  timing_list = timing_list[rnai]
  
  
  pdf(paste0(path,"Profils_",type,"_",RNAi,"_", select,"_", name, ".pdf"))
  for (s in rownames(EXPRESSION)){
    x_axis = unique(unlist(timing_list))
    x_axis_order = c(1,order(as.numeric(gsub("T", "", x_axis[2:length(x_axis)])))+1)
    x_axis = x_axis[x_axis_order]
    
    p = plot(NULL, xlim = c(1,length(x_axis)),
             ylim=c(0,MAX[s]),
             axes=F,ylab=type,
             xlab="Timing Autogamie (h)",
             main=paste(names(select_ID)[grep(s,select_ID)]))
    
    
    
    axis(1,at=1:length(x_axis),labels=x_axis,las=2)
    axis(2)
    
    col_pal = brewer.pal(n = length(names(timing_list)), name = "Set1")
    legend("topleft",legend=names(timing_list),col=col_pal,lwd=2, cex = 0.75,bty = "n") 
    
    color = 0
    for (i in names(timing_list)){
      color = color+1
      expression = EXPRESSION[s,grep(paste0(i,"_"), colnames(EXPRESSION))]
      positions = c()
      for (g in gsub("T","",gsub(paste0(i, "_"),"", colnames(expression)))){
        positions = c(positions,grep(g, x_axis)[1])
      }
      
      lines(positions,expression,col=col_pal[color],lwd=2)
      
    } 
    print(p)
    p = NULL
  }
  dev.off()
}


library(ggplot2)
Mean_expression <- function( condition, file, select_ID = NULL, rnai = NULL){
  
  path = paste0(file,condition,"/")
  EXPRESSION = read.table(paste0(path,condition ,"_expression_table_vst.tab"))
  infodata = read.table(paste0(path,condition ,"_infodata_collapse.tab"))
  
  path = paste0(path,"Visualisation/profils/")
  dir.create(path,recursive=T,showWarnings=F)
  
  if (rnai == "ctrl"){
    rnai = rnai_list[[condition]][c(grep("ND7",rnai_list[[condition]]),grep("ICL7",rnai_list[[condition]]))]
  }else {
    rnai = rnai_list[[condition]]
  }
  rnai = rnai[-grep("bis",rnai)]
  
  if(!is.null(select_ID)){
    EXPRESSION = EXPRESSION[select_ID,]
  } else {
    select_ID = annotation$ID
    names(select_ID) = annotation$NAME
  }
  
  timing_list = list()
  col_order = c()
  for(r in rnai){
    timing = str_remove_all(colnames(EXPRESSION)[grep(r, colnames(EXPRESSION))], paste0(r,"_"))
    timing = timing[c(length(timing), order(as.numeric(gsub("T", "", timing[1:(length(timing)-1)]))))]
    timing_list = c(timing_list, list(timing))
    for (t in timing){
      col_order = c(col_order,
                    grep(paste(r,t,sep="_"), colnames(EXPRESSION))[1])
    } 
  }
  names(timing_list)=rnai
  
  EXPRESSION = OrderColumn(EXPRESSION, infodata)
  
  
  timing_list = timing_list[rnai]
  rnai_name = paste(rnai, collapse = "_")
  
  tab = t(EXPRESSION)
  
  ggplot(EXPRESSION[1,]) 
  geom_point(aes(color = infodata$Feeding))
}

BoxnBarpolt_repartion <- function(LIST, path){
  for(i in 1:ncol(LIST[[1]])){
    pdf(paste0(path, colnames(LIST[[1]][i]),".pdf"))
    boxplot(c(LIST[[1]][i],LIST[[2]][i]), 
            names = c("Avec_Motif","Sans_motif"),
            main =  colnames(LIST[[1]][i]),
            outline = F)
    dev.off()
    
    pdf(paste0(path, colnames(LIST[[1]][i]),"_avecMotif.pdf"))
    hist(LIST[[1]][[i]],
         main =  paste0(colnames(LIST[[1]][i]),"_avecMotif"),
         breaks = 50,
         xlim = c(1,8))
    abline(v = mean(LIST[[1]][[i]]), col = "grey", lty = "dashed")
    abline(v = median(LIST[[1]][[i]]), col = "red", lty = "dashed")
    dev.off()
    
    pdf(paste0(path, colnames(LIST[[2]][i]),"_sansMotif.pdf"))
    hist(LIST[[2]][[i]],
         main =  paste0(colnames(LIST[[2]][i]),"_sansMotif"),
         breaks = 50,
         xlim = c(1,8))
    abline(v = mean(LIST[[2]][[i]]), col = "grey", lty = "dashed")
    abline(v = median(LIST[[2]][[i]]), col = "red", lty = "dashed")
    dev.off()
  }
}

###### Pour études des gènes DEG ######
Profile_Barplot <- function(filtre_list, nom, path){
  row_order = c("Early peak", "Intermediate peak", "Late peak", "Early repression" ,"Late induction", "Late repression", "none" )
  colors = c("purple3","red2","chartreuse4","dodgerblue3","deeppink","darkorange","snow3")
  
  profil = as.data.frame(table(annotation$EXPRESSION_PROFIL))
  for (n in names(filtre_list)){
    tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, filtre_list[[n]]))]))
    profil = merge(profil, tab, by = "Var1", all = T)
    
  }
  
  rownames(profil) = profil$Var1
  profil = profil[,-1]
  colnames(profil) = c("ALL", names(filtre_list))
  
  # Réordonner les lignes
  profil = profil[row_order,]
  
  
  ### Histogramme empilés
  pdf(paste0(path,"Profils_barplot_",nom,".pdf"),width = 550, height = 500)
  barplot(as.matrix(profil),
          col = colors,
          main = "Profil repartition",
          ylab = "gene nb")
  
  legend("topright",
         legend = rownames(profil),
         fill = colors,
         bty = "n")
  dev.off()
  
  # Création d'un tableau avec ses pourcentages
  profil_prct = profil
  for (n in 1:ncol(profil)){
    profil_prct[,n] = profil_prct[,n]/sum(profil[,n])*100
  }
  
  pdf(paste0(path,"Profils_barplot_",nom,"_prct.pdf"),width = 550, height = 500)
  barplot(as.matrix(profil_prct),
          col = colors,
          main = "Profil repartition",
          ylab = "% of genes",
          names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
  dev.off()
}

Profile_EnrichmentBarplot <- function (filtre_list, path, names){
  colors = c("purple3","red2","chartreuse4","dodgerblue3","deeppink","darkorange","snow3")
  
  profil = as.data.frame(table(annotation$EXPRESSION_PROFIL))
  for(up in names(filtre_list)){
    tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, filtre_list[[up]]))]))
    tab = merge(profil, tab, by = "Var1")
    rownames(tab) = tab[,"Var1"]
    tab = tab[,-1]
    colnames(tab) = c("ALL",up)
    
    tab_prct = tab
    tab_prct[,"ALL"] = tab[,"ALL"]/length(annotation$ID)*100
    for (n in rownames(tab)){
      tab_prct[n,2] = tab[n,2]/tab[n,"ALL"]*100
    }
    
    
    pdf(paste0(path,"Enrichment_barplot_",names,".pdf"),width = 800, height = 500)
    barplot(t(as.matrix(tab_prct)),
            beside = T,
            main = paste("Profil repartition of",up ),
            ylab = "% of genes",
            ylim = c(0,60),
            col = c("grey", "indianred2"),
            names.arg = sub(" "," \n ",rownames(tab)))
    
    legend("topleft",
           legend = paste0(colnames(tab)," (", apply(tab, 2, sum)," genes)"),
           fill = c("grey", "indianred2"),
           bty = "n")
    dev.off()
  }
}

IES_Barplot <- function(filtre_list, path, nom){
  profil = as.data.frame(c(sum(annotation$NB_IES != 0), 
                           sum(annotation$NB_IES == 0)), 
                         row.names = c("IES+", "IES-"))
  for (n in names(filtre_list)){
    tab = c(sum(annotation$NB_IES[which(is.element(annotation$ID, filtre_list[[n]]))] != 0),
            sum(annotation$NB_IES[which(is.element(annotation$ID, filtre_list[[n]]))] == 0))
    profil = cbind(profil, tab)
    
  }
  
  colnames(profil) = c("ALL", names(filtre_list))
  
  
  ### Histogramme empilés
  pdf(paste0(path,"Profils_barplot",nom,".pdf"),width = 550, height = 500)
  barplot(as.matrix(profil),
          col = colors,
          main = "Profil repartition",
          ylab = "gene nb")
  
  legend("topright",
         legend = rownames(profil),
         fill = colors,
         bty = "n")
  dev.off()
  
  # Création d'un tableau avec ses pourcentages
  profil_prct = profil
  for (n in 1:ncol(profil)){
    profil_prct[,n] = profil_prct[,n]/sum(profil[,n])*100
  }
  
  pdf(paste0(path,"Profils_barplot_",nom,"_prct.pdf"),width = 550, height = 500)
  barplot(as.matrix(profil_prct),
          col = colors,
          main = "Profil repartition",
          ylab = "% of genes",
          names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
  dev.off()
}

IES_EnrichmentBarplot <- function (filtre_list, path){
  profil = as.data.frame(c(sum(annotation$NB_IES != 0), 
                           sum(annotation$NB_IES == 0)), 
                         row.names = c("IES+", "IES-"))
  for(up in names(filtre_list)){
    tab = c(sum(annotation$NB_IES[which(is.element(annotation$ID, filtre_list[[up]]))] != 0),
            sum(annotation$NB_IES[which(is.element(annotation$ID, filtre_list[[up]]))] == 0))
    tab = cbind(profil, tab)
    colnames(tab) = c("ALL",up)
    
    tab_prct = tab
    tab_prct[,"ALL"] = tab[,"ALL"]/length(annotation$ID)*100
    for (n in rownames(tab)){
      tab_prct[n,2] = tab[n,2]/tab[n,"ALL"]*100
    }
    
    
    pdf(paste0(path,"Profils_barplot_",up,".pdf"),width = 400, height = 500)
    barplot(t(as.matrix(tab_prct)),
            beside = T,
            main = "Profil repartition of UP deregulated genes",
            ylab = "% of genes",
            ylim = c(0,60),
            col = c("grey", "indianred2"),
            names.arg = sub(" "," \n ",rownames(tab)))
    
    legend("topleft",
           legend = paste0(colnames(tab)," (", apply(tab, 2, sum)," genes)"),
           fill = c("grey", "indianred2"),
           bty = "n")
    dev.off()
  }
}

###### Pour études des motifs ######
PositionHistogram <- function (filtre_list, path, name){
  pos_list = list()
  for (n in names(filtre_list)){
    filtre = filtre_list[[n]]
    if (length(filtre) != 0){
      position = prom_motif$START[is.element(prom_motif$ID, filtre)]
      pos_list = c(pos_list, list(position) ) 
      pdf(paste0(path, "Histogramme_STARTposition_", n,".pdf"))
      hist(position, breaks = 75, xlim = c(-150,0), axes = F,
           xlab = paste("Distance from", debut),
           ylab = "Nb of motif",
           main = n)
      axis(2)
      axis(1, at = seq(-150,0,10))
      dev.off()
    }
  }
  pdf(paste0(path, "Boxplot_STARTposition_", name,".pdf"))
  names(pos_list) = names(filtre_list)
  boxplot(pos_list)
  dev.off()
}



PilBarplot <- function(strand_tab,filtre_list, nom, path, colors, row_order){
  
  tab_stand = as.data.frame(table(strand_tab[,2]))
  for (n in names(filtre_list)){
    tab = as.data.frame(table(strand_tab[which(is.element(strand_tab$ID, filtre_list[[n]])),2]))
    tab_stand = merge(tab_stand, tab, by = "Var1", all = T)
    
  }
  
  rownames(tab_stand) = tab_stand$Var1
  tab_stand = tab_stand[,-1]
  colnames(tab_stand) = c("ALL", names(filtre_list))
  tab_stand = tab_stand[row_order,]
  
  ### Histogramme empilés
  pdf(paste0(path,"Pil_barplot_",nom,".pdf"),width = 550, height = 500)
  barplot(as.matrix(tab_stand),
          col = colors)
  legend("topright",
         legend = rev(rownames(tab_stand)),
         fill = rev(colors),
         bty = "n")
  dev.off()
  
  # Création d'un tableau avec ses pourcentages
  profil_prct = tab_stand
  for (n in 1:ncol(tab_stand)){
    profil_prct[,n] = profil_prct[,n]/sum(tab_stand[,n])*100
  }
  
  pdf(paste0(path,"Pil_barplot_",nom,"_prct.pdf"),width = 550, height = 500)
  barplot(as.matrix(profil_prct),
          col = colors,
          main = "Profil repartition",
          ylab = "% of genes",
          names.arg = paste(colnames(profil_prct), apply(tab_stand, 2, sum), sep = "\n"))
  dev.off()
}

EnrichmentBarplot <- function (strand_tab, filtre_list, names, path, colors){
  
  tab_stand = as.data.frame(table(strand_tab[,2]))
  for(up in names(filtre_list)){
    tab = as.data.frame(table(strand_tab[which(is.element(strand_tab$ID, filtre_list[[up]])),2]))
    tab = merge(tab_stand, tab, by = "Var1")
    rownames(tab) = tab[,"Var1"]
    tab = tab[,-1]
    colnames(tab) = c("ALL",up)
    
    tab_prct = tab
    tab_prct[,"ALL"] = tab[,"ALL"]/length(strand_tab$ID)*100
    for (n in rownames(tab)){
      tab_prct[n,2] = tab[n,2]/tab[n,"ALL"]*100
    }
    
    
    pdf(paste0(path,"Enrichment_barplot_",names,"_",up,".pdf"),width = 800, height = 500)
    barplot(t(as.matrix(tab_prct)),
            beside = T,
            main = paste("Profil repartition of", up ),
            ylab = "% of genes",
            ylim = c(0,60),
            col = c("grey", "indianred2"),
            names.arg = sub(" "," \n ",rownames(tab)))
    
    legend("topleft",
           legend = paste0(colnames(tab)," (", apply(tab, 2, sum)," genes)"),
           fill = c("grey", "indianred2"),
           bty = "n")
    dev.off()
  }
}
