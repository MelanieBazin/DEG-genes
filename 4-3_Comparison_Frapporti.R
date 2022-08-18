####
# Comparison of the deregulated genes obtain #### 
# with this procedure 
# and those obtain in Frapporti et al. 2019
####
options(stringsAsFactors = FALSE)
library(ggvenn)
library("VennDiagram")

#### To re-open the files ####
# save_path = "./Analyse/2022-02-21_Analyse_DESeq2/HiSeqvsNextSeq_FC-1.5_pval-0.05/Frapporti_comparison_FC-2_pval-0.05/"
# data_tab = read.table(paste0(save_path, "HiSeqvsNextSeq_expression_table_vst.tab"), sep = "\t", header = T)
# infodata_collapse = read.table(paste0(save_path, "HiSeqvsNextSeq_infodata_collapse.tab"), sep = "\t", header = T)
###################################

for (deg in c("DOWN", "UP")){
  # Open published table of deregulated genes in EZL1 KD
  Frapporti = read.table(paste0("./DATA/",deg,"_EZL1-Frapporti_2019.csv"), sep = ";", header = T, quote = "")
  Frapporti = Frapporti[!apply(Frapporti == "", 1, all), ] # Delet empty lines 
  
  for(timing in c("EARLY", "LATE")){    
    #### Definition of significant deregulated genes ####
    
    # In Frapporti's data
    sign_Frapporti = Frapporti[ !is.na(Frapporti[paste0("Pvalue",timing)]) & Frapporti[paste0("Pvalue",timing)] < frapp_pvalue ,]
    sign_Frapporti = sign_Frapporti[sign_Frapporti[paste0("LogFC_",timing)] > log2(frapp_FC) | sign_Frapporti[paste0("LogFC_",timing)] < log2(1/frapp_FC),]
    
    # In my analysis
    MyData = read.table(paste0(path, "DESeq/EZL1/DEgenes_HiSeqvsNextSeq_",timing,".tab"), sep = "\t", header = T, quote = "")
    
    sign_MyData = MyData[ !is.na(MyData$padj) & MyData$padj < frapp_pvalue ,]
    sign_MyData = sign_MyData[sign_MyData$log2FoldChange > log2(frapp_FC) | sign_MyData$log2FoldChange < log2(1/frapp_FC),]
    
    # Save the merge
    CrossData = merge(sign_MyData, sign_Frapporti, by = "ID")
    write.table(CrossData, paste0(save_path, "Common_",deg,"_genes_", timing,".tab"), sep = "\t")
    
    #### Splitting UP and DOWN deregualted genes ####
    if (deg == "UP"){
      sign_MyData = sign_MyData[sign_MyData$REGULATION == "Up-regulated",]
      sign_Frapporti = sign_Frapporti[sign_Frapporti[paste0("LogFC_",timing)] > 0,]
      My_Up = sign_MyData
      Frap_Up = sign_Frapporti
    }else if (deg == "DOWN"){
      sign_MyData = sign_MyData[sign_MyData$REGULATION == "Down-regulated",]
      sign_Frapporti = sign_Frapporti[sign_Frapporti[paste0("LogFC_",timing)] < 0,]
    }
    
    #### Graphical data comparison ####
    # Venn diagram
    pdf(paste0(save_path,"Venn_",timing,"_sig",deg,".pdf"))
    grid.newpage() 
    draw.pairwise.venn(area1 = length(sign_MyData$ID),
                       area2 = length(sign_Frapporti$ID),
                       cross.area = length(intersect(sign_MyData$ID, sign_Frapporti$ID)),
                       category = c("My_data", "Frapporti"),
                       lwd = rep(0.5, 2),
                       fill = c("dodgerblue", "gold1"),
                       alpha = rep(0.5, 2),
                       scaled = TRUE)
    dev.off()

  }
  
}

### Comparison of the 
data_tab = OrderColumn(data_tab, infodata_collapse)

# Split the deregulated gene list :
diff_ID = list( not_Frapp = setdiff(My_Up$ID, Frap_Up$ID), # Genes deregulated in Frapporti but not in this study
                not_My = setdiff(Frap_Up$ID, My_Up$ID), # Genes deregulated in this study but not in Frapporti
                both = intersect(My_Up$ID, Frap_Up$ID)) # Genes deregulated in both study

for (d in names(diff_ID)){
  # Restrict the normalized expression table to one part of the Venn diagramm
  diff_My_Frapp = diff_ID[[d]]
  diff_My_Frapp = data_tab[diff_My_Frapp,]
  
  # Add the profile information to the normalized expression table
  rownames(annotation) = annotation$ID
  profil = annotation[, c("NAME" , "EXPRESSION_PROFIL")]
  profil = profil[rownames(diff_My_Frapp),]
  
  diff_My_Frapp = merge(profil, diff_My_Frapp, by= 0)
  
  
  profils = unique(diff_My_Frapp$EXPRESSION_PROFIL)
  
  pdf(paste0(save_path, "ExpressionVST_",d, ".pdf"))
  par(mfrow = c(1,2))
  
  # For each expression profile
  for (p in profils){
    data = diff_My_Frapp[diff_My_Frapp$EXPRESSION_PROFIL == p,]
    data = data[,which(is.element(colnames(data),colnames(data_tab)))]
    data = as.matrix(data)

    # For each time course
    for (r in c("ICL7", "EZL1")){
      data_box = data[,grep(r,colnames(data))]
      colnames(data_box) = str_remove_all(colnames(data_box), paste0(r,"_"))
      
      mediane = apply(data_box, 2, median)
      
      # Create a boxplot + median of the expression profile
      plot(NULL, main = p, sub = r,
           ylim = c(0, round(max(data))+1),
           ylab = "Expression level (vst)",
           xlim = c(0,ncol(data_box)+1),
           xlab = "",
           xaxt = "n")
      
      boxplot(data_box, add = T, boxwex =0.5, col = "white", outline = F)
      lines(mediane, col = profile_color[p], lwd = 4)
    }
  }
  dev.off()
}

