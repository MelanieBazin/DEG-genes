####
# Homemade function used in the analysis
####
options(stringsAsFactors = FALSE)
library(stringr)
library(RColorBrewer)
library(dendextend)
library(FactoMineR)
library(factoextra)
library(ggplot2)
library(gtools)
library(ggrepel)
library(ggbiplot)
library(topGO)

###### Table creation/modification ####
# Merge count table in one table to be analysed by DESeq2
ConcatTab <- function(type, conditions = NULL){
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

# Calculate the mean value of the replicates of each autogamy stages
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

# Creation of the table with the information about the sample in the column of a count table
CreatInfoData <- function(countdata, conditions, rnai_list, cluster, Timing = NULL){
  infodata = matrix(NA,nrow = ncol(countdata), ncol = 8)
  row.names(infodata) = colnames(countdata)
  colnames(infodata) = c("Names","Samples", "KnockDown", "Timing", "AutogStage", "Condition","Seq_method","Labo")
  
  infodata[,"Names"] = colnames(countdata)
  
  rnai = rnai_list[[conditions]]
  rnai = rnai[order(rnai)]
  
  timing = c()
  clust = c()
  batch = c()
  KnockDown = c()
  condition = c()
  labo = c()
  for(r in rnai){
    clust = c(clust, cluster[[r]])
    
    if (length(grep("CTRL",r))>0 | length(grep("ND7",r))>0 | length(grep("ICL7",r))>0 ){
      KnockDown =c(KnockDown, rep("ctrl", length(cluster[[r]])))
      condition = c(condition, paste(cluster[[r]],"ctrl",sep = "_"))
    }else if (length(grep("EZL1",r))>0){
      KnockDown =c(KnockDown, rep("EZL1", length(cluster[[r]])))
      condition = c(condition, paste(cluster[[r]],"EZL1",sep = "_"))
    }else{
      KnockDown =c(KnockDown, rep(r, length(cluster[[r]])))
      condition = c(condition, paste(cluster[[r]],r,sep = "_" ))
    }
    
    
    if (r == "ND7_K" | r == "PGM"| r == "KU80c" | r == "ICL7" | r == "EZL1"){
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
  
  
  infodata[,"KnockDown"] = KnockDown
  infodata[,"Timing"] = str_remove_all(timing, "T")
  infodata[,"AutogStage"] = clust
  infodata[,"Seq_method"] = batch
  infodata[,"Condition"]= condition
  infodata[,"Labo"]= labo
  infodata[,"Samples"] = str_remove_all(infodata[,"Names"],"bis")
  
  infodata = as.data.frame(infodata)
  
  return(infodata)
}


# Reordering the table by time and time course
OrderColumn <- function(data_tab, infodata){
  colum_order_ctrl = c()
  colum_order_rnai = c()
  
  infodata$RNAi = str_remove_all(str_split_fixed(infodata$Names, "_T", n=2)[,1], "_Veg")
  
  for (r in unique(infodata$RNAi)){
    rnai_position = str_which(colnames(data_tab),r)
    Timing = infodata$Timing[grep(r, infodata$Names)]
    Timing = sub("Veg","-1", Timing)
    
    cluster_position = rnai_position[order(as.numeric(Timing))]
    
    
    if(infodata$KnockDown[rnai_position[1]] == "ctrl"){
      colum_order_ctrl = c(colum_order_ctrl, cluster_position)
    }else{
      colum_order_rnai = c(colum_order_rnai, cluster_position)
    }
  }
  
  colum_order = c(colum_order_ctrl, colum_order_rnai)
  ordered_tab = data_tab[,colum_order]
  
  return(ordered_tab)
  
}

#### List modification ####
# Cross two list
Crossinglist <- function (list1, list2){
  LIST = list()
  list_name = c()
  for (l in names(list1)){
    for (m in names(list2)){
      LIST = c(LIST, list(intersect(list1[[l]],list2[[m]])))
      list_name = c(list_name, paste(names(list1[l]),names(list2[m]), sep = "x"))
    }
  }
  names(LIST)=list_name
  
  return(LIST)
}


#### Data exploration ####
# Define the color vector that will be used
Color_type <- function(data_tab, infodata, type){
  if(type == "methods"){
    color_list = method_color
    
    method_list = infodata$Seq_method
    names(method_list) = infodata$Names
    
    if(colnames(infodata)[ncol(infodata)] == "runsCollapsed"){
      method_list[grep("ICL7",names(method_list))] = "both"
      method_list[grep("EZL1",names(method_list))] = "both"
    }
    
  }else if(type == "replicates"){
    color_list = cluster_color
    
    method_list = infodata$AutogStage
    names(method_list) = infodata$Names
  }
  
  color = c()
  for (samp in names(method_list)){
    color = c(color, color_list[method_list[samp]])
  }
  
  names(color) = colnames(data_tab)
  return(color)
  
}

# PCA analysis using ggplot representation
PCA_ggplot_generator <- function(data_tab, infodata, save_path, color_type, main,  max_dim=3,barplot_max_dim=3,
                                 image_prefix="PCA_",show_barplot=T, selection = NULL, vline=0, sortie = "pdf", h = 5, w = 5,
                                 label = c("all","none","ind","ind.sup","quali","var","quanti.sup"),police_seize = 3, point_seize = 0.5, 
                                 rename = F, collapse = F, ...) {
  
  save_path = paste0(save_path,"/",color_type,"/ggPCA/")
  dir.create(save_path,recursive=T,showWarnings=F)
  
  # Identification of the different time courses to be ploted
  TimeCourses = str_split_fixed(colnames(data_tab),"_T", n= 2)[,1]
  TimeCourses = str_split_fixed(TimeCourses,"_V", n= 2)[,1]
  # Definition of the shape to apply
  shape_palette = c("CTIP" = 4,
                    "ICL7" = 1,
                    "KU80c" = 15,
                    "ND7_K" = 0,
                    "ND7_X" = 2,
                    "ND7_L" = 3,
                    "PGM" = 18,
                    "XRCC4" = 17,
                    "EZL1" = 16)
  
  # Shorten the name of the points
  shortname = colnames(data_tab)
  if(rename == T){
    shortname = str_replace_all(shortname, "Veg","V")
    shortname = str_remove_all(shortname, "ND7_")
    shortname =str_split_fixed(shortname,"_",n=2)[,2]
  }
  shortname = str_remove_all(shortname, "bis")
  
  # PCA analisis
  data_tab2 = t(data_tab)
  resExp = PCA(data_tab2, graph = F)
  
  # % of varience barplot
  if(show_barplot) {
    eigenvalues <- resExp$eig
    if (sortie == "pdf") {
      pdf(paste0(save_path,image_prefix,"_PCA_Variance.pdf"))
    }else if (sortie =="png"){
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
  
  # Definition of the PCA dimension to represent
  for (i in max_dim:1){
    data_tab2 = cbind(resExp$ind$coord[, i],data_tab2 )
    colnames(data_tab2)[1] = paste0("Dim",i)
  }
  
  # Definition of the color to apply
  if (color_type == "replicates"){
    color_info = infodata$AutogStage
    color_palette = cluster_color
  } else if (color_type == "methods"){
    if (collapse == T){
      for (t in c("T0", "T10"))
        infodata$Seq_method[grep(",", infodata$runsCollapsed)] = "Both"
    }
    color_info = infodata$Seq_method
    color_palette = method_color
  }
  
  # Definition of the shape to apply

  
  # PCA plot
  for (i in 1:dim(combn(1:max_dim,2))[2]){
    gp = ggplot(data = as.data.frame(data_tab2), aes(x = get(paste0("Dim",combn(1:max_dim,2)[1,i])), y = get(paste0("Dim",combn(1:max_dim,2)[2,i])), 
                                                     color = color_info, 
                                                     label = shortname,
                                                     shape = TimeCourses)) + 
      geom_point(alpha = 0.85, size = point_seize) +
      scale_shape_manual(values = shape_palette) +
      scale_color_manual(values = color_palette) +
      theme_classic(base_size = 9, base_line_size = 0.25) +
      geom_hline(yintercept = 0, lty = 2, size = 0.25) +
      geom_vline(xintercept = 0, lty = 2, size = 0.25) +
      geom_text_repel(size=police_seize) +
      guides(color = FALSE) +
      labs(title = main,
           x = paste0("Dim",combn(1:max_dim,2)[1,i]),
           y = paste0("Dim",combn(1:max_dim,2)[2,i]))
    
    if (rename != T){
      gp = gp+theme(legend.position="none")
    }
    gp
    
    ggsave(paste0(save_path,image_prefix,i,".",sortie), plot = gp, height = h, width = w, device = sortie, )
    
  }
  
  return(gp)
}

# Boxplot of the reads count for each time point of each time course
CountBoxplot <- function (tab, type, color = "lightgray"){
  boxplot(log(tab + 1), ylab = "count values (log scale)",
          main = paste0("Count data (",type,")"), xaxt="n", yaxt="n",
          col = color, outline=FALSE)
  axis(side = 1, labels = FALSE, tick = F)
  axis(side = 2, # Rotate labels perpendicular to y-axis.
       las = 2, # Adjust y-axis label positions.
       mgp = c(3, 0.75, 0))
  text(x = 1:ncol(tab), # Move labels to just below bottom of chart.
       y = par("usr")[3] - 0.1, # Use names from the data list.
       labels = colnames(tab), # Change the clipping region.
       xpd = NA, # Rotate the labels by 35 degrees.
       srt = 90, # Adjust the labels to almost 100% right-justified.
       adj = 1, # Increase label size.
       cex = 0.5)
}

# Draw the expression profiles for all or selected genes in all time courses
ExpressionProfils <- function(type , condition, path, name = NULL, select_ID = NULL, rnai = NULL, infodata = NULL){
  
  if (type  == "vst"){
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
    path2 = "./Analyse/Profils/"
    dir.create(path2,recursive=T,showWarnings=F)
    EXPRESSION = ConcatTab(type, conditions = rnai)
    
  }
  path = paste0(path,"Visualization/Profils/")
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

# Creation of pil-chart representing the proportion of each expression cluster for lists of genes
Profile_Barplot <- function(filtre_list, nom, save_path, w= 6.8, h = 7.2){
  
  # Create a table with the number of genes in each category
  profil = as.data.frame(table(annotation$EXPRESSION_PROFIL))
  colnames(profil) = c("Expression_cluster", "ALL")
  for (n in names(filtre_list)){
    tab = as.data.frame(table(annotation$EXPRESSION_PROFIL[which(is.element(annotation$ID, filtre_list[[n]]))]))
    colnames(tab) = c("Expression_cluster",n)
    profil = merge(profil, tab, by = "Expression_cluster", all = T)
    
  }
  
  profil[is.na(profil)] = 0
  rownames(profil) = profil$Expression_cluster
  profil = profil[,-grep("Expression_cluster", colnames(profil))]
  
  write.table(profil, paste0(save_path, "Profils_tab_",nom,".tab"), sep="\t",row.names=T,quote=F)
  
  # Transform the table in percentages
  profil_prct = t(t(profil)/apply(profil, 2, sum)*100)
  profil_prct = profil_prct[names(profile_color),] # reorder the lines
  
  # Creation of a piled-barplot with the different expression profile proportion
  pdf(paste0(save_path,"Profils_barplot_",nom,"_prct.pdf"),width = w, height = h)
  barplot(as.matrix(profil_prct),
          col = profile_color,
          main = "Profils repartition",
          ylab = "% of genes",
          names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
  dev.off()
  
  # Add a file with the legend
  pdf(paste0(save_path,"Profils_barplot_legends.pdf"),width = 4.3, height = 5.3)
  plot.new()
  legend("center",
         x.intersp=0.1,
         legend = rev(rownames(profil)),
         fill = rev(profile_color),
         bty = "n")
  dev.off()
  
}

# Creation of pil-chart representing the proportion of genes with IES in lists of genes
IES_Barplot <- function(filtre_list, nom, save_path, w= 6.8, h = 7.2){
  
  # Create a table with the number of genes in each category
  profil = as.data.frame(c(sum(annotation$NB_IES != 0), 
                           sum(annotation$NB_IES == 0)), 
                         row.names = c("IES+", "IES-"))
  for (n in names(filtre_list)){
    tab = c(sum(annotation$NB_IES[which(is.element(annotation$ID, filtre_list[[n]]))] != 0),
            sum(annotation$NB_IES[which(is.element(annotation$ID, filtre_list[[n]]))] == 0))
    profil = cbind(profil, tab)
    
  }
  colnames(profil) = c("ALL", names(filtre_list))
  write.table(profil, paste0(save_path, "IES_tab_",nom,".tab"), sep="\t",row.names=T,quote=F)
  
  # Transform the table in percentages
  profil_prct = t(t(profil)/apply(profil, 2, sum)*100)
  
  # Creation of a piled-barplot with the different proportion of genes with IES
  pdf(paste0(save_path,"IES_barplot_",nom,"_prct.pdf"),width = w, height = h)
  barplot(as.matrix(profil_prct),
          col = c("navyblue","snow2"),
          main = "Profil repartition",
          ylab = "% of genes",
          names.arg = paste(colnames(profil_prct), apply(profil, 2, sum), sep = "\n"))
  dev.off()
  
  # Add a file with the legend
  pdf(paste0(save_path,"IES_barplot_legends.pdf"),width = 4.3, height = 5.3)
  plot.new()
  legend("center",
         x.intersp=0.1,
         legend = rev(rownames(profil)),
         fill = rev(c("navyblue","snow2")),
         bty = "n")
  dev.off()
}


#### Statistics ####
# Calculate the enrichment in a profile of a gene list with correction for multiple tests
Enrichment_padj <- function(LIST, data_tab, nb_simulation = 1000){
  data_tab = as.data.frame(data_tab)
  colnames(data_tab) = c("ID", "PROFIL")
  
  TAB_list = list()
  
  print("Table of all data")
  print(table(data_tab$PROFIL))
  
  Profils = unique(data_tab$PROFIL)
  all = length(data_tab$ID)
  for (l in names(LIST)){
    ## Calculate the experimental p-value
    pval = c()
    nb_genes = length(LIST[[l]][which(is.element(LIST[[l]], data_tab$ID))])
    for (p in Profils){
      nb_profil = sum(data_tab$PROFIL == p)
      nb_genes_profil = length(LIST[l][which(is.element(LIST[[l]], data_tab$ID[data_tab$PROFIL == p]))])
      
      pval = c(pval, 1- phyper(nb_genes_profil-1, nb_profil, all - nb_profil, nb_genes))
    }
    names(pval) = Profils
    
    ## p-values correction by simulation for multiple-test
    nb_simulation = 1000
    
    # Create random data set and calculate theoretical p-values
    if (length(data_tab$ID)>=length(LIST[[l]])){
      theoric_ID = matrix(NA, nrow = length(LIST[[l]]), ncol = nb_simulation)
      for (c in 1:ncol(theoric_ID)){
        theoric_ID[,c] = sample(data_tab$ID, length(LIST[[l]]))
      }
    }else{
      theoric_ID = matrix(NA, nrow = length(data_tab$ID)*2/3, ncol = nb_simulation)
      for (c in 1:ncol(theoric_ID)){
        theoric_ID[,c] = sample(data_tab$ID, length(data_tab$ID)*2/3)
      }
    }
    theoric_pval = matrix(NA, nrow = length(Profils), ncol = nb_simulation)
    row.names(theoric_pval) = Profils
    
    for (s in 1:nb_simulation){
      genes = theoric_ID[,s]
      pvals = c()
      for (p in Profils){
        nb_profil = sum(data_tab$PROFIL == p)
        nb_genes = length(genes)
        nb_genes_profil = length(genes[which(is.element(genes, data_tab$ID[data_tab$PROFIL == p]))])
        pvals = c(pvals, 1- phyper(nb_genes_profil-1, nb_profil, all - nb_profil, nb_genes))
      }
      theoric_pval[,s] = pvals
    }
    
    # Calculate the adjusted p-value based on the simulated data
    pval_adj = c()
    for (p in Profils){
      th_pvals = theoric_pval[p,]
      pval_adj = c(pval_adj, sum(th_pvals< pval[p])/length(th_pvals))
    }
    
    names(pval_adj) = Profils
    pval_adj = format(pval_adj, scientific = T, digit = 5)
    
    print(l)
    freq = as.data.frame(table(data_tab$PROFIL[which(is.element(data_tab$ID, LIST[[l]]))]))
    tab = cbind(as.data.frame(pval_adj),as.numeric(pval_adj) < 0.01)
    tab = merge(freq, tab, by.x = "Var1", by.y = 0)
    
    rownames(tab) = tab$Var1
    tab = tab[,-1]
    colnames(tab) = c("nb_event", "pval_adj", "pval_adj < 0.01")
    print(tab)
    
    TAB_list = c(TAB_list,list(tab))
  }
  
  names(TAB_list) = names(LIST)
  
  return(TAB_list)
}

# Calculate the enrichment in a profile of a gene list
Chi2_pvalue <- function(LIST, data_tab){
  data_tab = as.data.frame(data_tab)
  colnames(data_tab) = c("ID", "PROFIL")
  
  print("Table of all data")
  print(table(data_tab$PROFIL))
  
  Profils = unique(data_tab$PROFIL)
  all = length(data_tab$ID)
  for (l in names(LIST)){
    ## Calculate the p-value
    pval = c()
    signif = c()
    nb_genes = length(LIST[[l]][which(is.element(LIST[[l]], data_tab$ID))])
    for (p in Profils){
      nb_profil = sum(data_tab$PROFIL == p)
      nb_genes_profil = length(LIST[l][which(is.element(LIST[[l]], data_tab$ID[data_tab$PROFIL == p]))])
      
      # Format the data for chi2
      data_mat = matrix(data=NA, ncol=2, nrow=2)
      data_mat[ , 1] = c(nb_genes_profil, nb_profil)
      data_mat[ , 2] = c(nb_genes-nb_genes_profil , all-nb_profil)
      rownames(data_mat) = c("IN_LIST", "ALL")
      colnames(data_mat) = c("With_profile","Without_profile")
      
      # Calculate the p-value
      chi2 = chisq.test(data_mat, correct = F)
      
      pval=c(pval, chi2$p.value)
      
      # Add stars to the significant value
      if(chi2$p.value < 1e-200) {sign="****"
      }else if(chi2$p.value < 1e-100) {sign="***"
      }else if(chi2$p.value < 1e-20) {sign="**"
      }else if(chi2$p.value < 1e-10) {sign="*"
      }else {sign="ns"}
      signif = c(signif, sign)
      
    }
    names(pval) = Profils
    pval = format(pval, scientific = T, digit = 5)
    nb = table(data_tab$PROFIL[which(is.element(data_tab$ID, LIST[[l]]))])
    tab = t(rbind(nb,rbind(pval,signif)))
    
    print(l)
    print(tab)
    
  }
}

# GO enrichment calculation

top_go_enrichment <- function(map_file,myInterestingGenes,vocabulary="MF", 
                              random=FALSE, pvalue=0.05, adjust_pvalue=TRUE,prefix="",
                              geneUniverse=c()) {

  geneID2GO = map_file
  GO2geneID = inverseList(geneID2GO)
  
  geneNames = names(geneID2GO)
  if(random==TRUE) {
    print("Random")
    myInterestingGenes = sample(geneNames, length(myInterestingGenes))
  }
  
  if(length(geneUniverse) != 0) {
    myInterestingGenes = intersect(geneNames,geneUniverse)
  }
  
  geneList = factor(as.integer(geneNames %in% myInterestingGenes))
  names(geneList) = geneNames

  GOdata = new("topGOdata",
               ontology = vocabulary,
               allGenes = geneList,
               annot = annFUN.gene2GO,
               gene2GO = geneID2GO)
  
  resultFisher = runTest(GOdata, algorithm = "classic", statistic = "fisher")
  results.table = GenTable(GOdata, classicFisher = resultFisher , topNodes = length(resultFisher@score))
  
  results.table.bh = results.table[results.table$classicFisher<=pvalue,]
  if(adjust_pvalue==TRUE) {
    results.table.bh = results.table[which(p.adjust(results.table[,"classicFisher"],method="BH")<=pvalue),] 
  }
  
  results.table.bh
}

