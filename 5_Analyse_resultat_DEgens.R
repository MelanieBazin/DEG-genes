options(stringsAsFactors = FALSE)
library("stringr") 
library(ggvenn)

source("0_Cluster.R")

file_name = "2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05"

annotation_basic = read.table("./DATA/ptetraurelia_mac_51_annotation_v2.0.tab",header=T,sep="\t",quote='')
my_annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
annotation = merge(annotation_basic[,c(1,4:7)], my_annotation[,c(1,2,3:7)], by = "ID")[,c(1,7,8,2:6,9:11)]
rm(annotation_basic,my_annotation)

TurboPGM = read.table("./DATA/TurboID/2114003-Pgm-ProteinMeasurements.txt",header=T,sep="\t")
TurboPGML4 = read.table("./DATA/TurboID/2114003-PgmL4-ProteinMeasurements.txt",header=T,sep="\t",quote='')

ID_tab = annotation
condition =  names(rnai_list)[1]

# for (condition in names(rnai_list)){
  RNAi_list = unique(rnai_list[[condition ]][-grep("ND7",rnai_list[[condition ]])])
  RNAi_list = RNAi_list[-grep("ICL7",RNAi_list)]
  RNAi_list = RNAi_list[-grep("bis",RNAi_list)]

  
  # RNAi = RNAi_list[1]
  # Regarder la derégulation dans chaque condition
  for (RNAi in RNAi_list){
    if (RNAi == "CTIP"){
      timing = c("EARLY","INTER")
    }else{
      timing = "LATE"
    }
    for (t in timing){
    tab = read.table(paste0("Analyse/",file_name,"/",condition,"/DESeq/",RNAi,"/NoFilter/DEgenes_",condition,"_",t,"_NoFilter.tab"), header = T, sep = "\t")
    colnames(tab)[3:5]= paste(RNAi,t,colnames(tab)[3:5], sep = "_")
    ID_tab = merge(ID_tab, tab[,c(1,3:5)], by = "ID", all = T)
    colnames(ID_tab) = str_remove_all(str_remove_all(colnames(ID_tab), ".x"),".y")
    }
  }
   colnames(TurboPGM)[2:ncol(TurboPGM)]= paste("TurboPGM",colnames(TurboPGM)[2:ncol(TurboPGM)], sep="_")
   colnames(TurboPGML4)[2:ncol(TurboPGML4)]= paste("TurboPGML4",colnames(TurboPGML4)[2:ncol(TurboPGML4)], sep="_")
   ID_tab = merge(ID_tab, TurboPGM, by = "PROTEIN_NAME", all = T)
   ID_tab = merge(ID_tab, TurboPGML4, by = "PROTEIN_NAME", all = T)
   ID_tab = ID_tab[,c(2:3,1,4:ncol(ID_tab))]
  
  write.table(as.matrix(ID_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes.tab"), sep = "\t", row.names = T)  

  #### Croiser les RNAi ####
  selection = list(
    CTIP_down = ID_tab$ID[unique(
      grep("Down-regulated",ID_tab$CTIP_EARLY_REGULATION),
      grep("Down-regulated",ID_tab$CTIP_INTER_REGULATION))],
    PGM_up = ID_tab$ID[grep("Up-regulated",ID_tab$PGM_LATE_REGULATION)],
    KU80c_up = ID_tab$ID[grep("Up-regulated",ID_tab$KU80c_LATE_REGULATION)],
    XRCC4_up = ID_tab$ID[grep("Up-regulated",ID_tab$XRCC4_LATE_REGULATION)]
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_selection.png"))
  ggvenn(selection,
         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
         stroke_size = 0.5,
         set_name_size = 5,
         set_name_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"))
  dev.off()
  
  select_ID = Reduce(intersect, selection)
  select_tab = ID_tab[is.element(ID_tab$ID, select_ID),]
  
  write.table(as.matrix(select_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection.tab"), sep = "\t", row.names = F )
  
  selection1 = list(
    PGM_up = ID_tab$ID[grep("Up-regulated",ID_tab$PGM_LATE_REGULATION)],
    KU80c_up = ID_tab$ID[grep("Up-regulated",ID_tab$KU80c_LATE_REGULATION)],
    XRCC4_up = ID_tab$ID[grep("Up-regulated",ID_tab$XRCC4_LATE_REGULATION)]
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_selection_UP.png"))
  ggvenn(selection1,
         fill_color = c( "#EFC000FF", "#868686FF", "#CD534CFF"),
         stroke_size = 0.5,
         set_name_size = 5,
         set_name_color = c("#EFC000FF", "#868686FF", "#CD534CFF"))
  dev.off()
  
  select1_ID = Reduce(intersect, selection1)
  select1_tab = ID_tab[is.element(ID_tab$ID, select1_ID),]
  
  write.table(as.matrix(select1_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection_UP.tab"), sep = "\t", row.names = F )
  
  selection2 = list(
    CTIP_down = ID_tab$ID[unique(
      grep("Down-regulated",ID_tab$CTIP_EARLY_REGULATION),
      grep("Down-regulated",ID_tab$CTIP_INTER_REGULATION))],
    UP = select1_ID
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_selection_CTIP_DOWN-et-UP.png"))
  ggvenn(selection2,
         fill_color = c( "#0073C2FF", "coral2"),
         stroke_size = 0.5,
         set_name_size = 5,
         set_name_color = c("#0073C2FF", "coral2"))
  dev.off()
  
  select2_ID = Reduce(intersect, selection2)
  select2_tab = ID_tab[is.element(ID_tab$ID, select2_ID),]
  
  write.table(as.matrix(select2_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection_CTIP_DOWN-et-UP.tab"), sep = "\t", row.names = F )
  
  
  
  #### Croiser avec les résultat de BioID ####
  turbo = list(
    CTIP_down = ID_tab$PROTEIN_NAME[unique(
      grep("Down-regulated",ID_tab$CTIP_EARLY_REGULATION),
      grep("Down-regulated",ID_tab$CTIP_INTER_REGULATION))],
    TurboPGM = TurboPGM$PROTEIN_NAME,
    TurboPGML4 = TurboPGML4$PROTEIN_NAME,
    UP = ID_tab$PROTEIN_NAME[which(is.element(ID_tab$ID,select1_ID))]
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_turbo_selection.png"))
  ggvenn(turbo,
         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "coral2"),
         stroke_size = 0.5,
         set_name_size = 5,
         set_name_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "coral2"))
  dev.off()
  
  turbo_ID = Reduce(intersect, turbo)
  turbo_tab = ID_tab[is.element(ID_tab$PROTEIN_NAME, turbo_ID),]
  
  write.table(as.matrix(turbo_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_turbo_selection.tab"), sep = "\t", row.names = F )
  
  
  
  turbo1 = list(
    CTIP_down = ID_tab$PROTEIN_NAME[unique(
      grep("Down-regulated",ID_tab$CTIP_EARLY_REGULATION),
      grep("Down-regulated",ID_tab$CTIP_INTER_REGULATION))],
    TurboPGM = TurboPGM$PROTEIN_NAME,
    TurboPGML4 = TurboPGML4$PROTEIN_NAME
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_turbo_ctip.png"))
  ggvenn(turbo1,
         fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
         stroke_size = 0.5,
         set_name_size = 5,
         set_name_color = c("#0073C2FF", "#EFC000FF", "#868686FF"))
  dev.off()
  
  turbo1_ID = Reduce(intersect, turbo1)
  turbo1_tab = ID_tab[is.element(ID_tab$PROTEIN_NAME, turbo1_ID),]
  
  write.table(as.matrix(turbo1_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_turbo_ctip.tab"), sep = "\t", row.names = F )
  
  
  turbo2 = list(
    TurboPGM = TurboPGM$PROTEIN_NAME,
    TurboPGML4 = TurboPGML4$PROTEIN_NAME,
    UP = ID_tab$PROTEIN_NAME[which(is.element(ID_tab$ID,select1_ID))]
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_turbo_up.png"))
  ggvenn(turbo2,
         fill_color = c("#0073C2FF", "#EFC000FF","coral2"),
         stroke_size = 0.5,
         set_name_size = 5,
         set_name_color = c("#0073C2FF", "#EFC000FF", "coral2"))
  dev.off()
  
  turbo2_ID = Reduce(intersect, turbo2)
  turbo2_tab = ID_tab[is.element(ID_tab$PROTEIN_NAME, turbo2_ID),]
  
  write.table(as.matrix(turbo2_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_turbo_UP.tab"), sep = "\t", row.names = F )
  
  
  
  
# }








