options(stringsAsFactors = FALSE)
source("6_Filtres.R")

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"
path_DESeq = paste0(path,"DESeq/")
path_motif = paste0(path, "MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/")

colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")


RNAi = list.files(path_DESeq)
tab = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")[,c(1,3,13)]

for (r in RNAi){
  path_tmp = paste0(path_DESeq,r,"/NoFilter/")
  data = list.files(path_tmp, pattern = "DEgenes")
  

  for (d in data){
    data_name = sub("DEgenes_tout","",sub("_NoFilter.tab","",d))
    if (data_name != ""){
      
      tab_tmp = read.table(paste0(path_tmp,d), header = T, sep = "\t")
      
      tab = merge(tab,tab_tmp[,c(1,5)], by = "ID", all = T)
      colnames(tab)[ncol(tab)] = paste0(r,data_name)
      print(colnames(tab)[ncol(tab)])
    }
  }
}

tab_select = tab[which(is.element(tab$ID,UP_PKXE)),]

write.table(tab_select, paste0(path, "CinetiqueDEGUP_PKXE.tab"),sep = "\t", row.names = F)



file = list.files(path_motif, pattern = "IN_MAC_CDS.gff")
prom_with_motif = read.table(paste0(path_motif, file), header=F, sep="\t")
colnames(prom_with_motif)=colnamesgff3

suspect = c()
for (r in c("PGM","KU80c","XRCC4","EZL1")){

  s = is.element(tab_select[,grep(paste0(r,"_INTER"), colnames(tab_select))], "Down-regulated")
  s = tab_select$ID[s]
  suspect = c( suspect,s)
}

suspect = tab_select[is.element(tab_select$ID,suspect),]

suspect_motif = intersect(suspect$ID, prom_with_motif$ID)

source("0_Cluster.R")
source("3_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")

selectID = suspect
names(selectID)= tab_select$NAME[is.element(tab_select$ID,suspect)]

ExpressionProfils(type = "DESeq2", 
                  condition = "tout", 
                  file = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/",
                  select_ID = selectID)


