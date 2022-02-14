options(stringsAsFactors = FALSE)
library(seqinr)
library(memes)
source("0_Cluster.R")

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")

# Definition des fichier promoteur à ouvrir
IES = NULL
debut = "TSS"
promoteur = read.fasta(paste0("./DATA/Promoteur/IN_MAC", IES, "_upstream_150nt_",debut,".fa"))

# Definitir les fichiers d'analyse à ouvrir
date = "02-08"
condition =  names(rnai_list)[2]

# Localiser les donner
file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
save_path = paste0("./Analyse/",file_name, "/", condition, "/Motif/From_",debut, "_IN_MAC",IES,"/")


# Ouvrir les filtres sur les dérégulation
RNAi = rnai_list[[condition]]
RNAi = RNAi[-grep("bis", RNAi)]
RNAi = RNAi[-grep("ICL7", RNAi)]
RNAi = RNAi[-grep("ND7", RNAi)]
source("5-1_Filtres.R")

##### Partie 1 - définir les promoteur pour analyse STREME ####
path = paste0(save_path, "Promotors/")
dir.create(path,recursive=T,showWarnings=F)

prom_UP = promoteur[which(is.element(names(promoteur),up_pkx))]
# write.fasta(sequences = prom_UP, names = names(prom_UP), file.out = paste0(path, "PromUP.fa") )

prom_UP_int = prom_UP[which(is.element(names(prom_UP),inter_genes))]
write.fasta(sequences = prom_UP_int, names = names(prom_UP_int), file.out = paste0(path, "PromUP_inter.fa") )

for (i in 1:5){
  prom_rand = sample(names(promoteur),length(names(prom_UP_int)), replace = F)
  prom_rand = promoteur[which(is.element(names(promoteur),prom_rand))]
  
  write.fasta(sequences = prom_rand, names = names(prom_rand), file.out = paste0(path, "PromRand_",i,".fa") )
  
}

##### Partie 2 - Déintion logo par STREME ####
# Definition lieu sauvegarde
path = paste0(save_path, "STREME/")
dir.create(save_path,recursive=T,showWarnings=F)

runStreme(input = paste0(save_path, "Promotors/PromUP_inter.fa"),
          control = paste0(save_path, "Promotors/PromRand_",i,".fa"),
          outdir = path,
          meme_path)




print("Faire les recherche de motif avec STREME avant de lancer 6-1_MeanMotif.R /n 
      Parametres : /n 
      Minimum width = 5 ; Minimum width = 20 /n 
      p-value threshold = 0.05 /n 
      Markov order -> default /n 
      Align sequences on their Right Ends")

# Définition d'un motif moyen
# source("6-1_MeanMotif.R")

##### Partie 4 - définition des prom avec motif par FIMO ####
path2 = paste0(save_path,"FIMO/")
dir.create(path2,recursive=T,showWarnings=F)

print("Faire les recherche de motif avec STREME avant de lancer 6-1_MeanMotif.R /n 
      Parametres : /n 
      Minimum width = 5 ; Minimum width = 20 /n 
      p-value threshold = 0.05 /n 
      Markov order -> default /n 
      Align sequences on their Right Ends")




