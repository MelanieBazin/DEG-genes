# install.packages("ggseqlogo")
options(stringsAsFactors = FALSE)

library(universalmotif)

library(ggplot2)
library(ggseqlogo)

files = list.files(path, pattern = ".meme")


logo1 = read_meme(paste0(path,files[1]))
logo2 = read.table(paste0(path,files[2]), header=F, sep=" ")
logo3 = read.table(paste0(path,files[3]), header=F, sep=" ")
logo4 = read.table(paste0(path,files[4]), header=F, sep=" ")
logo5 = read.table(paste0(path,files[5]), header=F, sep=" ")





files = list.files(path, pattern = "_freqs.txt")


logo1 = read.table(paste0(path,files[1]), header=F, sep=" ")
logo2 = read.table(paste0(path,files[2]), header=F, sep=" ")
logo3 = read.table(paste0(path,files[3]), header=F, sep=" ")
logo4 = read.table(paste0(path,files[4]), header=F, sep=" ")
logo5 = read.table(paste0(path,files[5]), header=F, sep=" ")

min = min(nrow(logo1), nrow(logo2),nrow(logo3),nrow(logo4),nrow(logo5))

logo_mean = matrix(nrow = min, ncol = 4)

for (l in 1:nrow(logo_mean)){
  for(c in 1:ncol(logo_mean)){
    logo_mean[l,c] = mean(logo1[l,c],logo2[l,c],logo3[l,c],logo4[l,c],logo5[l,c]) 
  }
}

write.table(logo_mean, paste0(path,"Mean_logo_freqs.txt"), sep = " ", row.names = F, col.names = F)


logo_mean = create_motif(t(logo_mean), alphabet = "DNA", strand = "+")

png(paste0(path,"Mean_logo.png"),width = 1200, height = 400)
view_motifs(logo_mean)
dev.off()


logo_mean_rev = motif_rc(logo_mean)


png(paste0(path,"Mean_logo_rev.png"),width = 1200, height = 400)
view_motifs(logo_mean_rev)
dev.off()
