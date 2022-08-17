options(stringsAsFactors = FALSE)
library(seqinr)
library(universalmotif)
library(ggplot2)
library(ggseqlogo)
source("0_Cluster.R")

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")
colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")

# Definition des fichier promoteur à ouvrir
IES = NULL
debut = "TSS"
promoteur = read.fasta(paste0("./DATA/Promoteur/IN_MAC", IES, "_upstream_150nt_",debut,".fa"))

# Definitir les fichiers d'analyse à ouvrir
date = Sys.Date()
date = "2022-02-21"
condition =  names(rnai_list)[2]

# Localiser les donner
file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
save_path = paste0("./Analyse/",file_name, "/", condition, "/Motif/From_",debut, "_IN_MAC",IES,"/UP_CTIP_inter/")


# Ouvrir les filtres sur les dérégulation
RNAi = rnai_list[[condition]]
RNAi = RNAi[-grep("bis", RNAi)]
RNAi = RNAi[-grep("ICL7", RNAi)]
RNAi = RNAi[-grep("ND7", RNAi)]
source("5-1_Filtres.R")

##### Partie 1 - définir les promoteur pour analyse STREME ####
print("List of promotors creation") 
path = paste0(save_path, "Promotors/")
dir.create(path,recursive=T,showWarnings=F)


prom_UP = promoteur[which(is.element(names(promoteur),stdCTIP$DOWN_UP))]
# write.fasta(sequences = prom_UP, names = names(prom_UP), file.out = paste0(path, "PromUP.fa") )

# prom_UP_not_int = prom_UP[which(!is.element(names(prom_UP),AUTOGAMY$inter_peak))]
# write.fasta(sequences = prom_UP_not_int, names = names(prom_UP), file.out = paste0(path, "PromUP_not_inter.fa") )

prom_UP_int = prom_UP[which(is.element(names(prom_UP),AUTOGAMY$inter_peak))]
write.fasta(sequences = prom_UP_int, names = names(prom_UP_int), file.out = paste0(path, "PromUP_CTIP_inter.fa") )

for (i in 1:5){
  prom_rand = sample(names(promoteur),length(names(prom_UP_int)), replace = F)
  prom_rand = promoteur[which(is.element(names(promoteur),prom_rand))]
  
  write.fasta(sequences = prom_rand, names = names(prom_rand), file.out = paste0(path, "PromRand_",i,".fa") )
  
}

##### Partie 2 - Définition du motif par STREME ####
# Definition lieu sauvegarde
path = paste0(save_path, "STREME/")
dir.create(path,recursive=T,showWarnings=F)

# Attendre que ce soit fait avant de passer à la suite
readline(prompt=paste("1- Search enriched motif on STREME",
                      "Parametres :",
                      "  Minimum width = 5 ; Minimum width = 20",
                      "  p-value threshold = 0.05",
                      "  Markov order -> default" ,
                      "  Align sequences on their Right Ends",
                      "2- Save the 5 Minimal MEME motifs in : ", 
                      path,
                      "Press [enter] to continue", sep ="\n"))

print("3- Definition of a mean motif") 

# Ouverture des fichiers MimimalMEME
files = list.files(path, pattern = ".meme")

motifs=c()
for (f in files){
  motif = read_meme(paste0(path,f))
  motif["name"] = paste("Motif",str_sub(f,1,1))
  motifs = c(motifs,motif)
}


# Compraisoan des différentes motifs obtenus
comp_logo = compare_motifs(motifs)
write.table(comp_logo, paste0(path, "PreasonCorrelation.tab"), sep = "\t")

#Création d'un motif moyen
MergeMotif = merge_motifs(motifs)
MergeMotif["name"] = "Mean_motif"
write_meme(MergeMotif, paste0(path, "0-MeanMotif.meme"))

motifs = c(motifs, MergeMotif)

# Création de logo pour les motifs
png(paste0(path,"All_Motifs_logo.png"),width = 800, height = 1000, bg = "transparent")
view_motifs(motifs)
dev.off()

for (m in motifs){
  v=view_motifs(m, show.positions = F)
  png(paste0(path,m["name"], "_logo.png"),width = 1000, height = 400, bg = "transparent")
  print(v)
  dev.off()
  
  #Image du reverse complement
  logo_mean_rev = motif_rc(m)
  v = view_motifs(logo_mean_rev, show.positions = F)
  png(paste0(path,m["name"], "_logo-rev.png"),width = 1000, height = 400, bg = "transparent")
  print(v)
  dev.off()
}

dev.off()


##### Partie 4 - définition des prom avec motif par FIMO ####
path2 = paste0(save_path,"FIMO_1E-4/")
dir.create(path2,recursive=T,showWarnings=F)

print("5- Search for motif in promotors on FIMO") 
print("Parametres :")
print(paste("  Uploade the input motif :", path, "MeanMotif.meme"))  
print(paste0("  Uploade the input sequence : ./DATA/Promoteur/IN_MAC", IES, "_upstream_150nt_",debut,".fa"))
print("  Match p-value < 1E-4"  )
print(paste0("6- Save the identified motif GFF table in : ", path2))


sink(paste0(save_path,"/sessionInfo_defintionmotif.txt"))
print(sessionInfo())
sink()
