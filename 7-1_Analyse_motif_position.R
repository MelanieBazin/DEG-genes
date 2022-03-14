# source("7_Ouverture_fichier_motif_filtres.R")

####  Répartition des motifs sur les promoteurs ####
print("Histogram and box plot of motif position")
path = paste0(save_path,"Positions/")
dir.create(path ,recursive=T,showWarnings=F)

# Parmis les gènes d'intéret

PositionHistogram(MOTIF, path, "MOTIF")
PositionHistogram(MOTIFxUP_PKX, path, "MOTIF_UP")
PositionHistogram(MOTIFxCTIP, path, "MOTIF_CTIP")

PositionHistogram(MOTIFxAUTOG, path, "MOTIF_AUTOG")
PositionHistogram(MOTIFxUP_PKXxAUTOG, path, "MOTIF_UP_AUTO")
PositionHistogram(MOTIFxCTIPxAUTOG, path, "MOTIF_CTIP_AUTO")

PositionHistogram(SUPP[1:3], path, "SUPP")

# Parmis les autres gènes

PositionHistogram(MOTIFxnotUP_PKX, path, "MOTIF_notUP")
PositionHistogram(MOTIFxnotCTIP, path, "MOTIF_notCTIP")
PositionHistogram(MOTIFxnotUP_PKXxAUTOG, path, "MOTIF_notUP_auto")
PositionHistogram(MOTIFxnotCTIPxAUTOG, path, "MOTIF_not_CTIP_auto")

# Faire une séléaction de gènes et voir les localisation
pos_up = prom_motif$START[which(is.element(prom_motif$ID,MOTIFxUP_PKX$MotifxUP_ALL))]
pos_not_up = prom_motif$START[which(is.element(prom_motif$ID,MOTIFxnotUP_PKX$Motifxnot_up_PKX))]

png(paste0(path, "Histogramme_STARTposition_UPvs_notUP.png"),width = 1000, height = 800)
par(mfrow=c(2,1))
hist(pos_not_up,breaks = 75, xlim = c(-150,0))
hist(pos_up,breaks = 75, xlim = c(-150,0))
dev.off()

t_not_up = table(pos_not_up)
t_up = table(pos_up)
tab = merge(as.data.frame(t_not_up), as.data.frame(t_up), by = 1, all = T)
tab[is.na(tab)] = 0
colnames(tab) = c("start", "not_UP", "UP")
mat = matrix(c(tab$not_UP, tab$UP),2,length(tab$not_UP),byrow=T)

sink(paste0(path,"/chi2_pos_UPvsNOTUP.txt"))
chisq.test(mat)
sink()

png(paste0(path, "Histogramme_STARTposition_not_UP_rand.png"),width = 1700, height = 900)
par(mfrow=c(4,5))

rand_pos = list()
rand_mean = c()
for (i in 1:20){
  rand = sample(1:length(pos_not_up),length(pos_up), replace = F)
  pos = pos_not_up[rand]
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

png(paste0(path, "Histogramme_STARTposition_not_UP_rand_last.png"))
par(mfrow=c(1,1))
hist(pos,breaks = 75, xlim = c(-150,0), axes = F,
     main = paste("pvalue :", round(chi2$p.value, 4)))
axis(2)
axis(1, at = seq(-150,0,10))
dev.off()
##
pos_not_inter = prom_motif$START[which(is.element(prom_motif$ID,SUPP$not_Inter_motif))]
png(paste0(path, "Histogramme_STARTposition_not_inter_rand.png"),width = 1700, height = 900)
par(mfrow=c(4,5))

rand_pos = list()
rand_mean = c()
for (i in 1:20){
  rand = sample(1:length(pos_not_inter),length(pos_up), replace = F)
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

png(paste0(path, "Histogramme_STARTposition_not_inter_rand_last.png"))
par(mfrow=c(1,1))
hist(pos,breaks = 75, xlim = c(-150,0), axes = F,
     main = paste("pvalue :", round(chi2$p.value, 4)))
axis(2)
axis(1, at = seq(-150,0,10))
dev.off()



## Calcul enrichissement
# fractionner les positions
pos_fraction = c()
for(l in 1:nrow(prom_motif)){
  if(prom_motif$START[l] < 0 & prom_motif$START[l] >=-20){
    pos_fraction = c(pos_fraction, "0-20")
  }else if(prom_motif$START[l] < -20 & prom_motif$START[l] >=-50){
    pos_fraction = c(pos_fraction, "-21-50")
  }else if(prom_motif$START[l] < -50 & prom_motif$START[l] >=-80){
    pos_fraction = c(pos_fraction, "-51-80")
  }else if(prom_motif$START[l] < -80 & prom_motif$START[l] >=-110){
    pos_fraction = c(pos_fraction, "-81-110")
  }else if(prom_motif$START[l] < -110 & prom_motif$START[l] >=-150){
    pos_fraction = c(pos_fraction, "-111-150")
  }else{
    print(paste("ERROR line",l))
  }
}


my_data = as.data.frame(cbind(prom_motif$ID, pos_fraction))
LIST = c(UP_PKX, stdCTIP)
names(LIST)= c(names(UP_PKX), names(stdCTIP))

sink(paste0(path,"/Enrichissement_position_autogamy.txt"))
Enrichment_padj(AUTOGAMY, my_data)
sink()
sink(paste0(path,"/Enrichissement_position_motif.txt"))
Enrichment_padj(LIST, my_data)
sink()


##### Etat de R #####
sink(paste0(path,"/Analyse_sessionInfo.txt"))
print(sessionInfo())
sink()
