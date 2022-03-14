# source("7_Ouverture_fichier_motif_filtres.R")

#### Barplot répartition des motifs +/- ####
print("Repartition of 50.80 motifs")
path = paste0(save_path,"Barplot_Motifs_plusVSmoins/")
dir.create(path ,recursive=T,showWarnings=F)

# Sur UP PGM KU80c & XRCC4
STRAND = list(
  plus = setdiff(MOTIF_uniq$Motif_plus, MOTIF_uniq$motif_moins),
  moins = setdiff(MOTIF_uniq$motif_moins, MOTIF_uniq$Motif_plus),
  both = intersect(MOTIF_uniq$Motif_plus, MOTIF_uniq$motif_moins),
  none = setdiff(annotation$ID, MOTIF_uniq$Motif))

strand_tab = matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("ID", "Motif_STRAND")))
for (S in names(STRAND)){
  plus_tab = rep(S, length(STRAND[[S]]))
  plus_tab = cbind(STRAND[[S]], plus_tab)
  
  strand_tab = rbind(strand_tab,plus_tab)
}
strand_tab =  as.data.frame(strand_tab)
row_order = c("moins", "both", "plus", "none")
color = c("royalblue1", "gold1", "chartreuse3", "grey")

PilBarplot(strand_tab,AUTOGAMY, "Motifs", path, color, row_order)
PilBarplot(strand_tab,UP_PKX, "UP_PKX", path, color, row_order)
PilBarplot(strand_tab,stdCTIP, "CTIP", path, color, row_order)
PilBarplot(strand_tab,not_UP_PKX, "not_UP_PKX", path, color, row_order)
PilBarplot(strand_tab,not_stdCTIP, "not_CTIP", path, color, row_order)

EnrichmentBarplot(strand_tab,AUTOGAMY,"strand", path, color)
EnrichmentBarplot(strand_tab,UP_PKX, "strand", path, color)
EnrichmentBarplot(strand_tab,stdCTIP, "strand", path, color)
EnrichmentBarplot(strand_tab,not_UP_PKX, "strand", path, color)
EnrichmentBarplot(strand_tab,not_stdCTIP, "strand", path, color)

#### Réapartition des gènes à la "bonne" position ####
print("Repartition of +/- motifs")
path = paste0(save_path,"Barplot_Motifs_position/")
dir.create(path ,recursive=T,showWarnings=F)

POS = list(
  in50.80 = MOTIF_uniq$motif_50.80,
  elseware = setdiff(MOTIF_uniq$Motif, MOTIF_uniq$motif_50.80),
  none = setdiff(annotation$ID, MOTIF_uniq$Motif)
)
pos_tab = matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("ID", "Motif_pos")))
for (S in names(POS)){
  plus_tab = rep(S, length(POS[[S]]))
  plus_tab = cbind(POS[[S]], plus_tab)
  
  pos_tab = rbind(pos_tab,plus_tab)
}
pos_tab =  as.data.frame(pos_tab)
row_order = names(POS)
color = c("chartreuse4", "deeppink3", "grey")

PilBarplot(pos_tab,AUTOGAMY, "Motifs_position", path, color, row_order)
PilBarplot(pos_tab,UP_PKX, "UP_PKX_position", path, color, row_order)
PilBarplot(pos_tab,stdCTIP, "CTIP_position", path, color, row_order)
PilBarplot(pos_tab,not_UP_PKX, "not_UP_PKX_position", path, color, row_order)
PilBarplot(pos_tab,not_stdCTIP, "not_CTIP_position", path, color, row_order)


EnrichmentBarplot(pos_tab,AUTOGAMY, "pos", path, color)
EnrichmentBarplot(pos_tab,UP_PKX, "pos", path, color)
EnrichmentBarplot(pos_tab,stdCTIP, "pos", path, color)
EnrichmentBarplot(pos_tab,not_UP_PKX, "pos", path, color)
EnrichmentBarplot(pos_tab,not_stdCTIP, "pos", path, color)

##### Etat de R #####
sink(paste0(path,"/Analyse_sessionInfo.txt"))
print(sessionInfo())
sink()