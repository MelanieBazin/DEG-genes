# install.packages("ggseqlogo")
options(stringsAsFactors = FALSE)

library(ggplot2)
library(ggseqlogo)

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"
path_logo = paste0(path, "MOTIF/Motif_dans_prom/Motif_intermed/STREME_CTIP/")


files = list.files(path_logo, pattern = ".txt")


logo1 = read.table(paste0(path_logo,files[1]), header=T, sep=" ")
logo2 = read.table(paste0(path_logo,files[2]), header=T, sep=" ")
logo3 = read.table(paste0(path_logo,files[3]), header=T, sep=" ")
logo4 = read.table(paste0(path_logo,files[4]), header=T, sep=" ")
logo5 = read.table(paste0(path_logo,files[5]), header=T, sep=" ")

logo_mean = matrix(nrow = 18, ncol = 4)

for (l in 1:18){
  for(c in 1:4){
    logo_mean[l,c] = mean(logo1[l,c],logo2[l,c],logo3[l,c],logo4[l,c],logo5[l,c]) 
  }
}

write.table(logo_mean, paste0(path_logo,"Mean_logo.txt"), sep = " ", row.names = F, col.names = F)


logo_mean = t(logo_mean)
row.names(logo_mean)= c("A","C", "G", "T")

cs1 = make_col_scheme(chars=c('A', 'C', 'G', 'T'), 
                      cols=c('red2', 'blue3', 'goldenrod2', 'green4'))

png(paste0(path_logo,"Mean_logo.png"),width = 1200, height = 400)
ggseqlogo(logo_mean, col_scheme= cs1)
dev.off()
