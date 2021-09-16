source("6_Filtres.R")
library(seqinr)

#### Génération de promoteurs aléatoir avec les meme effectifs que ls gènes d'intéret ####
promoteur = read.fasta(paste0("./DATA/Promoteur/IN_MAC_upstream_150nt_TSS.fa"))

prom_UP_int = promoteur[which(is.element(names(promoteur),UP_PKXE_inter))]
write.fasta(sequences = prom_UP_int, names = names(prom_UP_int), file.out = paste0(path_motif, "PromUP_PKXE_int.fa") )

for (i in 1:5){
  prom_rand = sample(names(promoteur),length(names(prom_UP_int)), replace = F)
  prom_rand = promoteur[which(is.element(names(promoteur),prom_rand))]

  write.fasta(sequences = prom_rand, names = names(prom_rand), file.out = paste0(path_motif, "PromRand_PKXE_",i,".fa") )

}

prom_UP_DOWN_int = promoteur[which(is.element(names(promoteur),UP_DOWN_inter))]
write.fasta(sequences = prom_UP_DOWN_int, names = names(prom_UP_DOWN_int), file.out = paste0(path_motif, "PromUP_DOWN_int.fa") )
for (i in 1:5){
  prom_rand = sample(names(promoteur),length(names(prom_UP_DOWN_int)), replace = F)
  prom_rand = promoteur[which(is.element(names(promoteur),prom_rand))]

  write.fasta(sequences = prom_rand, names = names(prom_rand), file.out = paste0(path_motif, "PromRand_PKXE_CTIP_",i,".fa") )

}