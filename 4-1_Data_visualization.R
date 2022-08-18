###################
#### DESeq2 data visualization
#### -> Boxplot of the normalized count
#### -> PCA and hierarchical clustering
###################
options(stringsAsFactors = FALSE)

#### Boxplot of normalized count ####
print(paste("-----> Normalized counts boxplot"))
pdf(paste0(graph_path, "Vst_Counts_bolxplot.pdf"))
CountBoxplot(data_tab, "vst counts", color = c(rep("darkolivegreen2",28), rep("chartreuse4",21))) 
dev.off()

#### Cluster visualization (PCA and hclust)  ####
for (color_type in c("methods","replicates")){
  
  # PCA analysis
  print(paste("-----> PCA colored by :", color_type))
  PCA_ggplot_generator(data_tab,
                       infodata_collapse,
                       color_type,
                       police_seize = 4,
                       point_seize = 2.5,
                       save_path = graph_path,
                       main = paste0("PCA ", condition," (DESeq2)"),
                       sortie = "pdf",
                       rename = T,
                       w = 6.5,
                       h = 6,
                       collapse = T)
  
  # Pearson correlation matrix & hierarchical clustering
  print(paste("-----> Hierarchical clustering colored by :", color_type))
  matDist = as.dist(1-cor(log2(data_tab+1), method="pearson"))
  res = hclust(matDist)
  res = as.dendrogram(res)
  
  color = Color_type(data_tab, infodata_collapse, type = color_type)
  labels_colors(res)= as.character(color)[order.dendrogram(res)]
  
  pdf(paste0(graph_path, color_type,"/", condition ,"_hclust_pearson_vst.pdf"),  width = 12, height = 2.5)
  plot(res, main = "pearson_vst")
  dev.off()
}


#### Boxplot of mean expression level ####
print(paste("-----> Boxplot of expression level"))

# Order the time course to plot
rnai = unique(infodata_collapse$KnockDown)
rnai = c(rnai[grep("ctrl", rnai)],rnai[-grep("ctrl", rnai)])

profils = unique(annotation$EXPRESSION_PROFIL)

# Calculate the mean expression of replicates
mean_data_tab = MeanTabCalculation(data_tab, infodata_collapse)
write.table(mean_data_tab,paste0(path,condition ,"_MEANexpression_table_vst.tab"), sep="\t",row.names=T,quote=F)

pdf(paste0(graph_path,condition , "_Boxplot_ExpressionVST", ".pdf"))
par(mfrow = c(1,length(rnai)))

# For each expression profile
for (p in profils){
  data = mean_data_tab[annotation$ID[annotation$EXPRESSION_PROFIL == p],]
  
  # For each time course
  for (r in rnai){
    data_box = data[,grep(r,colnames(data))]
    
    mediane = apply(data_box, 2, median)
    
    # Create a boxplot + median of the expression profile
    plot(NULL, main = paste0(p, "\n", "n = ", nrow(data_box)), sub = r,
         ylim = c(0, 20),
         ylab = "Expression level (vst)",
         xlim = c(0.5,4.5),
         xlab = "",
         xaxt = "n")
    
    if(r == "PGM"){ # Exception for PGM since there is no EARLY stage
      boxplot(data_box, add = T, boxwex =0.5, col = "white", outline = F,
              names = str_remove_all(colnames(data_box), paste0("_",r)),
              at = c(1, 3, 4))
      lines(c(mediane[1], mean(c(mediane[1], mediane[2])), mediane[2:3]), col = profile_color[p], lwd = 4)
    }else{
      boxplot(data_box, add = T, boxwex =0.5, col = "white", outline = F,
              names = str_remove_all(colnames(data_box), paste0("_",r)))
      lines(mediane, col = profile_color[p], lwd = 4)
    }
  }
}

dev.off()






