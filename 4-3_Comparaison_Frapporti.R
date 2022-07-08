options(stringsAsFactors = FALSE)
library(ggvenn)

save_path = paste0(path,"Compared_Frapporti/")
dir.create(save_path,recursive=T,showWarnings=F)

for (deg in c("UP", "DOWN")){
  Frapporti = read.table(paste0("./DATA/",deg,"_EZL1-Frapporti_2019.csv"), sep = ";", header = T, quote = "")
  # Supprimer les lignes vides
  Frapporti = Frapporti[!apply(Frapporti == "", 1, all), ]  
  
  for(timing in c("EARLY", "LATE")){    
    # Definir les genes significatif d'après les donnée Frapporti
    sign_Frapporti = Frapporti[ !is.na(Frapporti[paste0("Pvalue",timing)]) & Frapporti[paste0("Pvalue",timing)] < pvalue ,]
    sign_Frapporti = sign_Frapporti[sign_Frapporti[paste0("LogFC_",timing)] >= log2(FC) | sign_Frapporti[paste0("LogFC_",timing)] <= log2(1/FC),]
    
    # Ouverture de mon fichier de donné correspondant
    MyData = read.table(paste0(path, "/DESeq/EZL1/DEgenes_HiSeqvsNextSeq_",timing,"_NoFilter.tab"), sep = "\t", header = T, quote = "")
    
    sign_MyData = MyData[ !is.na(MyData$padj) & MyData$padj < pvalue ,]
    sign_MyData = sign_MyData[sign_MyData$log2FoldChange > log2(FC) | sign_MyData$log2FoldChange < log2(1/FC),]
    
    if (deg == "UP"){
      sign_MyData = sign_MyData[sign_MyData$REGULATION == "Up-regulated",]
    }else if (deg == "DOWN"){
      sign_MyData = sign_MyData[sign_MyData$REGULATION == "Down-regulated",]
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
    
    CrossData = merge(sign_MyData, sign_Frapporti, by = "ID")
    
    write.table(CrossData, paste0(save_path, "Common_",deg,"_genes_", timing,".tab"), sep = "\t")
  }
  
}
