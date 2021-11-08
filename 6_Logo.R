# install.packages("ggseqlogo")
options(stringsAsFactors = FALSE)

library(universalmotif)

library(ggplot2)
library(ggseqlogo)

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"
path_logo = paste0(path, "MOTIF/Motif_dans_prom/Motif_intermed/STREME/")


files = list.files(path_logo, pattern = ".txt")


logo1 = read.table(paste0(path_logo,files[4]), header=F, sep=" ")
logo2 = read.table(paste0(path_logo,files[5]), header=F, sep=" ")
logo3 = read.table(paste0(path_logo,files[2]), header=F, sep=" ")
logo4 = read.table(paste0(path_logo,files[1]), header=F, sep=" ")
logo5 = read.table(paste0(path_logo,files[3]), header=F, sep=" ")

logo_mean = matrix(nrow = 18, ncol = 4)

for (l in 1:nrow(logo_mean)){
  for(c in 1:ncol(logo_mean)){
    logo_mean[l,c] = mean(logo1[l,c],logo2[l,c],logo3[l,c],logo4[l,c],logo5[l,c]) 
  }
}

write.table(logo_mean, paste0(path_logo,"Mean_logo.txt"), sep = " ", row.names = F, col.names = F)


logo_mean = create_motif(t(logo_mean), alphabet = "DNA", strand = "+")




png(paste0(path_logo,"Mean_logo.png"),width = 1200, height = 400)
view_motifs(logo_mean)
dev.off()


logo_mean_rev = motif_rc(logo_mean)


png(paste0(path_logo,"Mean_logo_rev.png"),width = 1200, height = 400)
view_motifs(logo_mean_rev)
dev.off()
