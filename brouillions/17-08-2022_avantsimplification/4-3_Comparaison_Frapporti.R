options(stringsAsFactors = FALSE)
library(ggvenn)
library("VennDiagram")

# path = "./Analyse/2022-02-21_Analyse_DESeq2_FC-1.5_pval-0.05/HiSeqvsNextSeq/"

save_path = paste0(path,"Compared_Frapporti/")
dir.create(save_path,recursive=T,showWarnings=F)

for (deg in c("DOWN", "UP")){
  Frapporti = read.table(paste0("./DATA/",deg,"_EZL1-Frapporti_2019.csv"), sep = ";", header = T, quote = "")
  # Supprimer les lignes vides
  Frapporti = Frapporti[!apply(Frapporti == "", 1, all), ]  
  
  for(timing in c("EARLY", "LATE")){    
    # Definir les genes significatif d'après les donnée Frapporti
    sign_Frapporti = Frapporti[ !is.na(Frapporti[paste0("Pvalue",timing)]) & Frapporti[paste0("Pvalue",timing)] < pvalue ,]
    sign_Frapporti = sign_Frapporti[sign_Frapporti[paste0("LogFC_",timing)] > log2(FC) | sign_Frapporti[paste0("LogFC_",timing)] < log2(1/FC),]
    
    # Ouverture de mon fichier de donné correspondant
    MyData = read.table(paste0(path, "DESeq/EZL1/DEgenes_HiSeqvsNextSeq_",timing,"_NoFilter.tab"), sep = "\t", header = T, quote = "")
    
    sign_MyData = MyData[ !is.na(MyData$padj) & MyData$padj < pvalue ,]
    sign_MyData = sign_MyData[sign_MyData$log2FoldChange > log2(FC) | sign_MyData$log2FoldChange < log2(1/FC),]
    
    if (deg == "UP"){
      sign_MyData = sign_MyData[sign_MyData$REGULATION == "Up-regulated",]
      sign_Frapporti = sign_Frapporti[sign_Frapporti[paste0("LogFC_",timing)] > 0,]
      My_Up = sign_MyData
      Frap_Up = sign_Frapporti
    }else if (deg == "DOWN"){
      sign_MyData = sign_MyData[sign_MyData$REGULATION == "Down-regulated",]
      sign_Frapporti = sign_Frapporti[sign_Frapporti[paste0("LogFC_",timing)] < 0,]
    }
    
    selection = list(sign_MyData$ID, sign_Frapporti$ID)
    names(selection)=c("MyData", "Frapporti")
    p = ggvenn(selection,
               fill_color = c("dodgerblue", "gold1"),
               stroke_size = 0.5,
               set_name_size = 6,
               show_percentage = F,
               text_size = 6,
               set_name_color = c("dodgerblue", "gold1"))
    
    pdf(paste0(save_path,"Venn_",timing,"_sig",deg,".pdf"))
    print(p)
    dev.off()
    
    pdf(paste0(save_path,"VennProp_",timing,"_sig",deg,".pdf"))
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
    
    
    
    CrossData = merge(sign_MyData, sign_Frapporti, by = "ID")
    
    write.table(CrossData, paste0(save_path, "Common_",deg,"_genes_", timing,".tab"), sep = "\t")
  }
  
}

vst = read.table(paste0(path, "HiSeqvsNextSeq_expression_table_vst.tab"), sep = "\t", header = T)
infodata = read.table(paste0(path, "HiSeqvsNextSeq_infodata_collapse.tab"), sep = "\t", header = T)
vst = OrderColumn(vst, infodata)

diff_ID = list( not_Frapp = setdiff(My_Up$ID, Frap_Up$ID),
                not_My = setdiff(Frap_Up$ID, My_Up$ID),
                both = intersect(My_Up$ID, Frap_Up$ID))

for (d in names(diff_ID)){
  diff_My_Frapp = diff_ID[[d]]
  diff_My_Frapp = vst[setdiff(My_Up$ID, Frap_Up$ID),]

  rownames(annotation) = annotation$ID
  profil = annotation[, c("NAME" , "EXPRESSION_PROFIL")]
  profil = profil[rownames(diff_My_Frapp),]
  
  diff_My_Frapp = merge(profil, diff_My_Frapp, by= 0)
  
  profils = unique(diff_My_Frapp$EXPRESSION_PROFIL)
  
  pdf(paste0(save_path, "ExpressionVST_",d, ".pdf"))
  par(mfrow = c(1,2))
  for (p in profils){
    data = diff_My_Frapp[diff_My_Frapp$EXPRESSION_PROFIL == p,]
    data = data[,which(is.element(colnames(data),colnames(vst)))]
    
    max = apply(data[,grep("ICL7", colnames(data))], 1, max)
    data = data[order(max, decreasing = T),]


    data = as.matrix(data)

    # heatmap.2(data,
    #           Rowv=F, Colv=F,dendrogram = "none", key=T,
    #           colsep = as.integer(length(grep("ICL7",  colnames(vst)))),
    #           labRow="",
    #           trace="none",
    #           col=colorRampPalette(c("white","goldenrod1","red2","darkred","black"))(1000),
    #           main=p)
    

    for (r in c("ICL7", "EZL1")){
      
      data_box = data[,grep(r,colnames(data))]
      colnames(data_box) = str_remove_all(colnames(data_box), paste0(r,"_"))
      
      
      mediane = apply(data_box, 2, median)
      
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

