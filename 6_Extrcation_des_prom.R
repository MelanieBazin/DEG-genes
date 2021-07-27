options(stringsAsFactors = FALSE)

library(seqinr)
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")

IES = NULL
debut = "TSS"
promoteur = read.fasta(paste0("./DATA/Promoteur/IN_MAC_", IES, "upstream_150nt_",debut,".fa"))

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"
UP = read.table(paste0(path,"Resumer_DEgenes_selection_UP.tab"), header = T)

prom_UP = promoteur[which(is.element(names(promoteur),UP$ID))]
# write.fasta(sequences = prom_UP, names = names(prom_UP), file.out = paste0(path, "PromUP.fa") )

Inter = annotation[which(annotation$EXPRESSION_PROFIL == "Intermediate peak"),]

prom_UP_int = prom_UP[which(is.element(names(prom_UP),Inter$ID))]
write.fasta(sequences = prom_UP_int, names = names(prom_UP_int), file.out = paste0(path, "PromUP_inter.fa") )

for (i in 1:5){
  prom_rand = sample(names(promoteur),length(names(prom_UP_int)), replace = F)
  prom_rand = promoteur[which(is.element(names(promoteur),prom_rand))]
  
  write.fasta(sequences = prom_rand, names = names(prom_rand), file.out = paste0(path, "PromRand_",i,".fa") )
  
}

# prom_UP_rand500 =   sample(names(prom_UP_int),500, replace = F)
# prom_UP_rand500 = promoteur[which(is.element(names(prom_UP_int),prom_UP_rand500))]
# write.fasta(sequences = prom_UP_rand500, names = names(prom_UP_rand500), file.out = paste0(path, "PromUP500.fa") )
