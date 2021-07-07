options(stringsAsFactors = FALSE)
library("stringr") 

source("0_Cluster.R")

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
for (i in grep("PTIWI",annotation$SYNONYMS)){
  annotation$NAME[i]= str_split(annotation$SYNONYMS[i],",")[1]
}


ID_tab = annotation[,1:5]
condition =  names(rnai_list)[1]

# for (condition in names(rnai_list)){
  RNAi_list = unique(rnai_list[[condition ]][-grep("ND7",rnai_list[[condition ]])])
  if (is.element("ICL7", RNAi_list)){
    RNAi_list = RNAi_list[-grep("ICL7",RNAi_list)]
  }
  
  # RNAi = RNAi_list[1]
  # Regarder la derégulation dnas chaque condtion
  for (RNAi in RNAi_list){
    if (RNAi == "CTIP"){
      timing = c("EARLY","INTER")
    }else{
      timing = "LATE"
    }
    for (t in timing){
    tab = read.table(paste0("Analyse/Analyse_DESeq2_test03_tout_batch/",condition,"/DESeq/",RNAi,"/NoFilter/DEgenes_",condition,"_",t,"_NoFilter.tab"), header = T, sep = "\t")
    colnames(tab)[3:5]= paste(RNAi,t,colnames(tab)[3:5], sep = "_")
    ID_tab = merge(ID_tab, tab[,c(1,3:5)], by = "ID", all = T)
    }
  }
# }
ID_tab = as.matrix(ID_tab)
write.table(ID_tab,"Analyse/Analyse_DESeq2_test03_tout_batch/tout/Resumer_DEgenes.tab", sep = "\t", row.names = T)


ID_ctip = ID_tab[which(
  ID_tab$CTIP_EARLY_REGULATION == "Down-regulated" |
  ID_tab$CTIP_INTER_REGULATION == "Down-regulated"
  ),c(1,grep("CTIP", colnames(ID_tab)))] 
ID_up = ID_tab[which(
  ID_tab$PGM_LATE_REGULATION =="Up-regulated" & 
  ID_tab$KU80c_LATE_REGULATION == "Up-regulated" &
  ID_tab$XRCC4_LATE_REGULATION == "Up-regulated"
  ),-grep("CTIP", colnames(ID_tab))] 

select_ET = merge(ID_up, ID_ctip, by = "ID")
select_ET = as.matrix(select_ETselect_ET)
write.table(select_ET,"Analyse/Analyse_DESeq2_test03_tout_batch/tout/Resumer_DEgenes_select_ET.tab", sep = "\t", row.names = F)

select_OU = ID_tab[which(
  ID_tab$CTIP_EARLY_REGULATION == "Down-regulated" |
  ID_tab$CTIP_INTER_REGULATION == "Down-regulated"|
  ID_tab$PGM_LATE_REGULATION =="Up-regulated" | 
  ID_tab$KU80c_LATE_REGULATION == "Up-regulated" |
  ID_tab$XRCC4_LATE_REGULATION == "Up-regulated"
),] 

write.table(select_OU,"Analyse/Analyse_DESeq2_test03_tout_batch/Resumer_DEgenes_select_OU.tab", sep = "\t", row.names = F )



rownames(select_OU)= 1:nrow(select_OU)
nm = grep("PTET", select_OU$NAME)
syn = c(grep("PTET", select_OU$SYNONYMS),grep("PTMB", select_OU$SYNONYMS),grep("rab", select_OU$SYNONYMS))
a = intersect(nm, syn)

abs =as.numeric(rownames(select_OU)[which(select_OU$SYNONYMS == "")])
b = intersect(nm, abs)
c = c(a, b)

select_OUbis = select_OU[-c,]
write.table(select_OUbis,"Analyse/Analyse_DESeq2_test03_tout_batch/Resumer_DEgenes_select_OUbis.tab", sep = "\t", row.names = F )
