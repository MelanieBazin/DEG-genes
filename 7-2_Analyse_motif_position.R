#### Position motif ####
path = paste0(save_path,"Positions/")
dir.create(path ,recursive=T,showWarnings=F)


PositionHistogram(MOTIF, path, "MOTIF")
PositionHistogram(MOTIFxUP_PKX, path, "MOTIF_UP")
PositionHistogram(MOTIFxCTIP, path, "MOTIF_CTIP")
PositionHistogram(MOTIFxAUTOG, path, "MOTIF_AUTOG")
PositionHistogram(MOTIFxUP_PKXxAUTOG, path, "MOTIF_UP_AUTO")
PositionHistogram(MOTIFxCTIPxAUTOG, path, "MOTIF_CTIP_AUTO")
PositionHistogram(SUPP[1:3], path, "SUPP")

PositionHistogram(MOTIFxnotUP_PKX, path, "MOTIF_notUP")
PositionHistogram(MOTIFxnotCTIP, path, "MOTIF_notCTIP")
PositionHistogram(MOTIFxnotUP_PKXxAUTOG, path, "MOTIF_notUP_auto")

# Faire une séléaction de gènes et voir les localisation
pos_inter = prom_motif$START[which(is.element(prom_motif$ID,SUPP$Inter_motif))]
pos_not_inter = prom_motif$START[which(is.element(prom_motif$ID,SUPP$not_Inter_motif))]

pdf(paste0(path, "Histogramme_STARTposition_intervs_notinter.pdf"),width = 1000, height = 800)
par(mfrow=c(2,1))
hist(pos_not_inter,breaks = 75, xlim = c(-150,0))
hist(pos_inter,breaks = 75, xlim = c(-150,0))
dev.off()

t_not_up = table(pos_not_inter)
t_up = table(pos_inter)
tab = merge(as.data.frame(t_not_up), as.data.frame(t_up), by = 1, all = T)
tab[is.na(tab)] = 0
colnames(tab) = c("start", "not_UP", "UP")
mat = matrix(c(tab$not_UP, tab$UP),2,length(tab$not_UP),byrow=T)

sink(paste0(path,"/chi2_pos_UPvsNOTUP.txt"))
chisq.test(mat)
sink()

pdf(paste0(path, "Histogramme_STARTposition_not_inter_rand.pdf"),width = 1700, height = 900)
par(mfrow=c(4,5))

rand_pos = list()
rand_mean = c()
for (i in 1:20){
  rand = sample(1:length(pos_not_inter),length(pos_inter), replace = F)
  pos = pos_not_inter[rand]
  rand_pos = c(rand_pos, list(pos))
  rand_mean= c(rand_mean, median(pos))
  
  t_not_up = table(pos)
  tab = merge(as.data.frame(t_not_up), as.data.frame(t_up), by = 1, all = T)
  tab[is.na(tab)] = 0
  colnames(tab) = c("start", "not_UP", "UP")
  mat = matrix(c(tab$not_UP, tab$UP),2,length(tab$not_UP),byrow=T)
  
  chi2 = chisq.test(mat)
  
  hist(pos,breaks = 75, xlim = c(-150,0), axes = F,
       main = paste("pvalue :", round(chi2$p.value, 4)))
  axis(2)
  axis(1, at = seq(-150,0,10))
  
}
dev.off()

pdf(paste0(path, "Histogramme_STARTposition_not_inter_rand_last.pdf"))
par(mfrow=c(1,1))
hist(pos,breaks = 75, xlim = c(-150,0), axes = F,
     main = paste("pvalue :", round(chi2$p.value, 4)))
axis(2)
axis(1, at = seq(-150,0,10))
dev.off()

## Calcul enrichissement
# fractionner les positions
fraction = 1

pos_fraction = c()
for(l in 1:nrow(prom_motif)){
  for(n in -seq(15, 150-fraction, by=fraction)){
    if(prom_motif$START[l] < n & prom_motif$START[l] >= n-fraction){
      pos_fraction = c(pos_fraction, paste(n,n-fraction))
    }
  }
}


my_data = as.data.frame(cbind(prom_motif$ID, pos_fraction))
LIST = c(UP_PKX, stdCTIP)
LIST = c(LIST, list(UP_CTIP_inter = intersect(stdCTIP$DOWN_UP, AUTOGAMY$inter_peak)))

sink(paste0(path,"/Enrichissement_position_autogamy_",fraction,"frac.txt"))
Enrichment_padj(AUTOGAMY, my_data)
sink()
sink(paste0(path,"/Enrichissement_position_motif_",fraction,"frac.txt"))
Enrichment_padj(LIST, my_data)
sink()

all_pos = prom_motif[prom_motif$START <= b & prom_motif$START >= a,]
all_pos_DOWN_UP = all_pos$ID[which(is.element(all_pos$ID, stdCTIP$DOWN_UP))]
all_pos_other = all_pos$ID[which(!is.element(all_pos$ID, stdCTIP$DOWN_UP))]

motif_DOWN_UP = prom_motif$ID[which(is.element(prom_motif$ID, stdCTIP$DOWN_UP))]
motif_other = prom_motif$ID[which(!is.element(prom_motif$ID, stdCTIP$DOWN_UP))]


tab = c(nrow(all_pos), nrow(prom_motif), nrow(all_pos)/nrow(prom_motif)*100)
tab = cbind(tab, c(length(all_pos_DOWN_UP), length(motif_DOWN_UP), length(all_pos_DOWN_UP)/length(motif_DOWN_UP)*100))
tab = cbind(tab, c(length(all_pos_other), length(motif_other), length(all_pos_other)/length(motif_other)*100))
colnames(tab) = c("ALL", "DOWN_UP", "other")
rownames(tab) = c(paste(a,b), "anywhere", "prct")
write.table(tab, paste0(path,"nb_of_motif.tab"), sep = "\t")
