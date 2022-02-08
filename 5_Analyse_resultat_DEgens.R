options(stringsAsFactors = FALSE)
library("stringr") 
library(ggvenn)

source("0_Cluster.R")

date = "02-08"

file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
annotation = annotation[,c(1,3:5,13,6:11,2)]
rownames(annotation)=annotation$ID

TurboPGM = read.table("./DATA/TurboID/2114003-Pgm-ProteinMeasurements.txt",header=T,sep="\t")
TurboPGML4 = read.table("./DATA/TurboID/2114003-PgmL4-ProteinMeasurements.txt",header=T,sep="\t",quote='')

ID_tab = annotation
condition =  names(rnai_list)[2]

# for (condition in names(rnai_list)){
  RNAi_list = unique(rnai_list[[condition ]][-grep("ND7",rnai_list[[condition ]])])
  RNAi_list = RNAi_list[-grep("bis",RNAi_list)]
  RNAi_list = RNAi_list[-grep("ICL7",RNAi_list)]
  

  
  # RNAi = RNAi_list[1]
  # Regarder la derégulation dans chaque condition
  for (RNAi in RNAi_list){
    if (RNAi == "CTIP"){
      timing = c("EARLY","INTER")
    }else{
      timing = "LATE"
    }
    for (t in timing){
    tab = read.table(paste0("./Analyse/",file_name,"/",condition,"/DESeq/",RNAi,"/DEgenes_",condition,"_",t,"_NoFilter.tab"), header = T, sep = "\t")
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
         fill_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"),
         stroke_size = 0.5,
         set_name_size = 7,
         show_percentage = F,
         text_size = 7,
         set_name_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"))
  dev.off()
  
  select_ID = Reduce(intersect, selection)
  select_tab = ID_tab[is.element(ID_tab$ID, select_ID),]
  
  write.table(as.matrix(select_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection.tab"), sep = "\t", row.names = F )
  
  #### Selection des gènes UP DEG en PGM, KU & XRCC4 ####
  selection1 = list(
    PGM_up = ID_tab$ID[grep("Up-regulated",ID_tab$PGM_LATE_REGULATION)],
    KU80c_up = ID_tab$ID[grep("Up-regulated",ID_tab$KU80c_LATE_REGULATION)],
    XRCC4_up = ID_tab$ID[grep("Up-regulated",ID_tab$XRCC4_LATE_REGULATION)]
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_selection_UP.png"))
  ggvenn(selection1,
         fill_color = c( "chocolate", "#868686FF", "darkolivegreen3"),
         stroke_size = 0.5,
         set_name_size = 7,
         show_percentage = F,
         text_size = 7,
         set_name_color = c("chocolate", "#868686FF", "darkolivegreen3"))
  dev.off()
  
  select1_ID = Reduce(intersect, selection1)
  select1_tab = ID_tab[is.element(ID_tab$ID, select1_ID),]
  
  write.table(as.matrix(select1_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection_UP.tab"), sep = "\t", row.names = F )
  
  #### Selection des gènes DOWN DEG en PGM, KU & XRCC4 ####
  selection1bis = list(
    PGM_up = ID_tab$ID[grep("Down-regulated",ID_tab$PGM_LATE_REGULATION)],
    KU80c_up = ID_tab$ID[grep("Down-regulated",ID_tab$KU80c_LATE_REGULATION)],
    XRCC4_up = ID_tab$ID[grep("Down-regulated",ID_tab$XRCC4_LATE_REGULATION)]
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_DOWN.png"))
  ggvenn(selection1bis,
         fill_color = c( "chocolate", "#868686FF", "darkolivegreen3"),
         stroke_size = 0.5,
         set_name_size = 7,
         show_percentage = F,
         text_size = 7,
         set_name_color = c("chocolate", "#868686FF", "darkolivegreen3"))
  dev.off()
  
  select1bis_ID = Reduce(intersect, selection1bis)
  select1bis_tab = ID_tab[is.element(ID_tab$ID, select1bis_ID),]
  
  write.table(as.matrix(select1bis_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_DOWN.tab"), sep = "\t", row.names = F )
  
  #### Selection des gènes UP DEG en PGM, KU & XRCC4 et DOWN DEG en CTIP ####
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
         set_name_size = 7,
         show_percentage = F,
         text_size = 10,
         set_name_color = c("#0073C2FF", "coral2"))
  dev.off()
  
  select2_ID = Reduce(intersect, selection2)
  select2_tab = ID_tab[is.element(ID_tab$ID, select2_ID),]
  
  write.table(as.matrix(select2_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_selection_CTIP_DOWN-et-UP.tab"), sep = "\t", row.names = F )
  
  select2_expression = read.table("./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/tout_expression_table_normaliserDESeq2.tab",sep = "\t", h=T, row.names = 1)
  col = c(grep("ICL7", colnames(select2_expression)),grep("ND7", colnames(select2_expression)),grep("CTIP", colnames(select2_expression)))
  select2_expression = select2_expression[is.element(rownames(select2_expression),select2_ID),col]
  col = c(7,1,5,2:4,6,13,8,12,9:11,20,14,19,15:18,25,21,24,22:23,30,26,29,27:28)
  select2_expression = select2_expression[,col]
  col = grep("Veg",colnames(select2_expression))
  select2_expression = select2_expression[,-col]
  
  MyHeatmaps(path = paste0("Analyse/",file_name,"/",condition,"/"),
             data_tab = select2_expression,
             condition = condition)
  
  low = select2_expression[select2_expression$CTIP_T12.5 < 200,]
  select2_low = ID_tab[is.element(ID_tab$ID, rownames(low)),]
  MyHeatmaps(path = paste0("Analyse/",file_name,"/",condition,"/low_FC2"),
             data_tab = low,
             condition = condition)
  
  write.table(as.matrix(select2_low),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_low_CTIP_DOWN-et-UP.tab"), sep = "\t", row.names = F )
  
  summary(select2_expression)
  
  #### Croiser avec les résultat de BioID ####
  #### Selection des gènes UP DEG en PGM, KU & XRCC4 et DOWN DEG en CTIP et présent dans les bioID de Marc ####
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
         fill_color = c("#0073C2FF", "chocolate", "plum", "coral2"),
         stroke_size = 0.5,
         set_name_size = 7,
         show_percentage = F,
         text_size = 7,
         set_name_color = c("#0073C2FF", "chocolate", "plum", "coral2"))
  dev.off()
  
  turbo_ID = Reduce(intersect, turbo)
  turbo_tab = ID_tab[is.element(ID_tab$PROTEIN_NAME, turbo_ID),]
  
  write.table(as.matrix(turbo_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_turbo_selection.tab"), sep = "\t", row.names = F )
  
  
  #### Selection des gènes DOWN DEG en CTIP et preésent dans les bioID de Marc ####
  turbo1 = list(
    CTIP_down = ID_tab$PROTEIN_NAME[unique(
      grep("Down-regulated",ID_tab$CTIP_EARLY_REGULATION),
      grep("Down-regulated",ID_tab$CTIP_INTER_REGULATION))],
    TurboPGM = TurboPGM$PROTEIN_NAME,
    TurboPGML4 = TurboPGML4$PROTEIN_NAME
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_turbo_ctip.png"))
  ggvenn(turbo1,
         fill_color = c("#0073C2FF", "chocolate", "plum"),
         stroke_size = 0.5,
         set_name_size = 7,
         show_percentage = F,
         text_size = 7,
         set_name_color = c("#0073C2FF", "chocolate", "plum"))
  dev.off()
  
  turbo1_ID = Reduce(intersect, turbo1)
  turbo1_tab = ID_tab[is.element(ID_tab$PROTEIN_NAME, turbo1_ID),]
  
  write.table(as.matrix(turbo1_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_turbo_ctip.tab"), sep = "\t", row.names = F )
  
  #### Selection des gènes UP DEG en PGM, KU & XRCC4 et présent ans les BioID de Marc ####
  turbo2 = list(
    TurboPGM = TurboPGM$PROTEIN_NAME,
    TurboPGML4 = TurboPGML4$PROTEIN_NAME,
    UP = ID_tab$PROTEIN_NAME[which(is.element(ID_tab$ID,select1_ID))]
  )
  png(paste0("Analyse/",file_name,"/",condition,"/Venn_turbo_up.png"))
  ggvenn(turbo2,
         fill_color = c("chocolate", "plum","coral2"),
         stroke_size = 0.5,
         set_name_size = 7,
         show_percentage = F,
         text_size = 7,
         set_name_color = c("chocolate", "plum", "coral2"))
  dev.off()
  
  turbo2_ID = Reduce(intersect, turbo2)
  turbo2_tab = ID_tab[is.element(ID_tab$PROTEIN_NAME, turbo2_ID),]
  
  write.table(as.matrix(turbo2_tab),paste0("Analyse/",file_name,"/",condition,"/Resumer_DEgenes_turbo_UP.tab"), sep = "\t", row.names = F )
  
  
# }








