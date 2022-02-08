options(stringsAsFactors = FALSE)
library(ggvenn)


date = "02-08"


file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
FC = 1.5
pvalue = 0.05

path = paste0("./Analyse/",file_name, "/HiSeqvsNextSeq/DESeq/EZL1/")


Frapporti = read.table("./DATA/DE_EZL1-Frapporti_2019.csv", sep = ";", header = T)
# Supprimer les lignes vides
Frapporti = Frapporti[!apply(Frapporti == "", 1, all), ]  


for(timing in c("EARLY", "LATE")){
  # Definir les genes significatif d'après les donnée Frapporti
  sign_Frapporti = Frapporti[ !is.na(Frapporti[paste0("Pvalue",timing)]) & Frapporti[paste0("Pvalue",timing)] < pvalue ,]
  sign_Frapporti = sign_Frapporti[sign_Frapporti[paste0("LogFC_",timing)] >= log2(FC) | sign_Frapporti[paste0("LogFC_",timing)] <= log2(1/FC),]
  
  # Ouverture de mon fichier de donné correspondant
  MyData = read.table(paste0(path, "DEgenes_HiSeqvsNextSeq_",timing,"_NoFilter.tab"), sep = "\t", header = T, quote = "")
  
  selection = list(MyData$ID, sign_Frapporti$ID)
  names(selection)=c("MyData", "Frapporti")
  p = ggvenn(selection,
         fill_color = c("dodgerblue", "gold1"),
         stroke_size = 0.5,
         set_name_size = 6,
         show_percentage = F,
         text_size = 6,
         set_name_color = c("dodgerblue", "gold1"))
  
  png(paste0(path,"Venn_",timing,"sig.png"))
  print(p)
  dev.off()
  
  CrossData = merge(MyData, sign_Frapporti, by = "ID")
  
  write.table(CrossData, paste0(path, "Common_genes_", timing,".tab"), sep = "\t")
  
}
