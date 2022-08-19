####
# Motif definition step
####
options(stringsAsFactors = FALSE)
library(seqinr)
library(universalmotif)
library(ggplot2)
library(ggseqlogo)

source("0_Cluster.R") # Groups parameters

IES = NULL
debut = "TSS"
# Open the promoter files
promoteur = read.fasta(paste0("./DATA/Promoteur/IN_MAC", IES, "_upstream_150nt_",debut,".fa"))

# Define the analyse to open
date = Sys.Date()
date = "2022-08-17"
condition =  names(rnai_list)[2]

source("5-1_Filters_candidats.R")

save_path = paste0(path, "Motif/From_",debut, "_IN_MAC",IES,"/")

#### Part 1 - Motif definition by STREME ####
save_path2 = paste0(save_path, "Promotors/")
dir.create(save_path2,recursive=T,showWarnings=F)

# File for "Input the sequences"
prom_int = promoteur[which(is.element(names(promoteur),candidats))]
write.fasta(sequences = prom_int, names = names(prom_int), file.out = paste0(path, "Prom_interest.fa") )

# Files for "Input the control sequences" for the 5 round of STREME
for (i in 1:5){
  prom_rand = sample(names(promoteur),length(names(prom_int)), replace = F)
  prom_rand = promoteur[which(is.element(names(promoteur),prom_rand))]
  
  write.fasta(sequences = prom_rand, names = names(prom_rand), file.out = paste0(path, "PromRand_",i,".fa") )
  
}

# Files to save the STREME data
save_path2 = paste0(save_path, "STREME/")
dir.create(save_path2,recursive=T,showWarnings=F)

readline(prompt=paste("STRME Parametres :",
                      "  Minimum width = 5 ; Minimum width = 20",
                      "  p-value threshold = 0.05",
                      "  Markov order -> default" ,
                      "  Align sequences on their Right Ends",
                      "2- Save the 5 Minimal MEME motifs and the sequences.tsv table in : ", 
                      path,
                      "Press [enter] to continue", sep ="\n"))

#### Part 2 - Identification of the genes present in all STREME analysis ####
file_name = list.files(save_path2)[grep("sequences",list.files(save_path2))]
genes_motif = NULL
for (f in file_name){
  liste = read.table(paste0(save_path2, f), sep = "\t", header = T)
  liste = unique(liste$seq_ID)
  genes_motif = c(genes_motif,list(liste))
}
names(genes_motif) = file_name

liste = genes_motif[[1]]
for (f in 2:length(genes_motif)){
  liste = intersect(liste, genes_motif[[f]])
}

write(liste, paste0(save_path2,"Genes_with_all_motif.txt"))

#### Part 3 - Formation of a mean motif ####
# Open the meme files
files = list.files(save_path2, pattern = ".meme")

motifs=list()
for (f in files){
  motif = read_meme(paste0(save_path2,f))
  motif["name"] = paste("Motif",str_sub(f,1,1))
  motifs = c(motifs, list(motif))
}
names(motifs) = files

# Correlation of the motif obtain in the different round of STREME
comp_logo = compare_motifs(motifs)
write.table(comp_logo, paste0(save_path2, "PreasonCorrelation.tab"), sep = "\t")

# Trimmed and merged motif
motifs_t = trim_motifs(motifs, min.ic = 1,trim.from = "both")

MergeMotif = merge_motifs(motifs_t)
MergeMotif["name"] = "Merge_motif"
write_meme(MergeMotif, paste0(save_path2, "0-MergeMotif.meme"), overwrite = T)

motifs = c(motifs, Merge = list(MergeMotif))

# Logos view of the motifs
pdf(paste0(save_path2,MergeMotif["name"], "_logo.pdf"),width = 6, height = 2)
view_motifs(MergeMotif, show.positions = F)
dev.off()

Merge_rev = motif_rc(MergeMotif)
pdf(paste0(save_path2,MergeMotif["name"], "_logo-rev.pdf"),width = 6, height = 2)
view_motifs(Merge_rev, show.positions = F)
dev.off()

pdf(paste0(save_path2,"All_Motifs_logo.pdf"),width = 7, height = 6)
view_motifs(motifs)
dev.off()


#### Part 4 - Detection of the motif in all promoters using FIMO ####
save_path3 = paste0(save_path,"FIMO/")
dir.create(save_path3,recursive=T,showWarnings=F)

print("5- Search for motif in promoters on FIMO") 
print("Parametres :")
print(paste("  Uploade the input motif :", save_path2, "MeanMotif.meme"))  
print(paste0("  Uploade the input sequence : ./DATA/Promoteur/IN_MAC", IES, "_upstream_150nt_",debut,".fa"))
print("  Match p-value < 1E-4"  )
print(paste0("6- Save the .tsv table of identified motif in : ", save_path3))

#### Print R status ####
sink(paste0(save_path2,"/sessionInfo_defintionmotif.txt"))
print(sessionInfo())
sink()
