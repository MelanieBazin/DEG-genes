
annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")

CTIP = read.table("Analyse_DESeq2_101genes/CTIP_EXPRESSION/NoFilter/DEgenes_CTIP_NoFilter.tab",h=T,sep="\t")
colnames(CTIP)[2]="CTIP_REGULATION"
XRCC4 = read.table("Analyse_DESeq2_101genes/XRCC4_EXPRESSION/NoFilter/DEgenes_XRCC4_NoFilter.tab",h=T,sep="\t")
colnames(XRCC4)[2]="XRCC4_REGULATION"
PGM = read.table("Analyse_DESeq2_101genes/PGM_EXPRESSION/NoFilter/DEgenes_PGM_NoFilter.tab",h=T,sep="\t")
colnames(PGM)[2]="PGM_REGULATION"
KU80c = read.table("Analyse_DESeq2_101genes/KU80c_EXPRESSION/NoFilter/DEgenes_KU80c_NoFilter.tab",h=T,sep="\t")
colnames(KU80c)[2]="KU80c_REGULATION"

CTIP_man = read.table("DATA/Manuel/DEG_manuel_CTIP.txt",h=T,sep="\t")
XRCC4_man = read.table("DATA/Manuel/DEG_manuel_XRCC4.tab",h=T,sep="\t")

selection = c("Ku","ku","PGM","NOWA","PTIWI","mt","TFIIS4","Spo11","Mre11","CER","Rad51", "Lig", "EZL", "SPT", "DCL", "CtIP", "XRCC4", "PDSG2", "PolX", "CAF1")

selection_ID =c()

for(i in selection){
  selection_ID = c(selection_ID,annotation$ID[grep(i,annotation$SYNONYMS)])
  
}


tab=annotation[is.element(selection_ID, annotation$ID)],c(1,3,4,2,5)]
tab = merge(tab, PGM[,1:2], by = "ID", all = T)

tab = merge(PGM[,1:2], KU80c[,1:2], by = "ID", all = T)
tab = merge(tab, CTIP_man, by = "ID", all = T)
tab = merge(tab, CTIP[,1:2], by = "ID", all = T)
tab = merge(tab, XRCC4_man, by = "ID", all = T)
tab = merge(tab, XRCC4[,1:2], by = "ID", all = T)
tab = merge(annotation[,c(1,3,4,2,5,6)], tab, by= "ID")
tab = tab[order(tab$ID),c(1:5,7,6,8:12)]


write.table(tab, "Comparaison données.tab",sep="\t", row.names = F, quote = F)
write.csv(tab, "Comparaison données.csv",sep="\t", row.names = F, quote = F)
