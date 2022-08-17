path = "./Analyse/2022-02-21_Analyse_DESeq2_FC-1.5_pval-0.05/analyseDE/Motif/From_TSS_IN_MAC/UP_CTIP_inter/STREME/"
file_name = list.files(path)[grep("sequences",list.files(path))]

enrichment_STREME = mean(112/187, 112/187, 86/187, 85/187, 90/187)
STREME = mean(112, 112, 86, 85, 90)

genes_motif = NULL

for (f in file_name){
  liste = read.table(paste0(path, f), sep = "\t", header = T)
  liste = unique(liste$seq_ID)
  genes_motif = c(genes_motif,list(liste))
}

names(genes_motif) = file_name

ID = intersect(intersect(intersect(intersect(genes_motif[[1]], genes_motif[[2]]), genes_motif[[3]]), genes_motif[[4]]), genes_motif[[5]])

enrichissement = length(ID)/187

STREME
enrichment_STREME
length(ID)
enrichissement
