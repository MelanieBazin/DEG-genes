options(stringsAsFactors = FALSE)
library(ggvenn)
library(ggplot2) 
library(RColorBrewer)
source("0_Cluster.R")
source("0_Visualisation_fonction.R")

# Definition des fichier promoteur à ouvrir
IES = NULL
debut = "TSS"

# Definitir les fichiers d'analyse à ouvrir
date = "02-08"
condition =  names(rnai_list)[2]
p_valueFIMO = "1E-5"

# Localiser les donner
file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
save_path = paste0("./Analyse/",file_name, "/", condition, "/Motif/From_",debut, "_IN_MAC",IES,"/FIMO_",p_valueFIMO)

# Ouvrir les filtres sur les dérégulation
RNAi = rnai_list[[condition]]
RNAi = RNAi[-grep("bis", RNAi)]
RNAi = RNAi[-grep("ICL7", RNAi)]
RNAi = RNAi[-grep("ND7", RNAi)]
source("5-1_Filtres.R")

# Ouverture du fichier FIMO
colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")
path = paste0(save_path,"FIMO/")
prom_motif = read.table(paste0(path,"fimo.gff"), sep = "\t", header = T)
colnames(prom_motif) = colnamesgff3
prom_motif$START = prom_motif$START-150
prom_motif$END= prom_motif$END-150

#### Création de filtre supplémentaire pour les motifs ####
print("Addition of new filter")
      
UP_PKX = c(UP_PKX, UP_ALL = list(up_pkx))

MOTIF = list(
  Motif = prom_motif$ID,
  Motif_plus = prom_motif$ID[prom_motif$STRAND == "+"],
  motif_moins = prom_motif$ID[prom_motif$STRAND == "-"],
  motif_50.80 = prom_motif$ID[prom_motif$START > -80 & prom_motif$START < -50 ]
)

MOTIF_uniq = list(
  Motif = unique(prom_motif$ID),
  Motif_plus = unique(prom_motif$ID[prom_motif$STRAND == "+"]),
  motif_moins = unique(prom_motif$ID[prom_motif$STRAND == "-"]),
  motif_50.80 = unique(prom_motif$ID[prom_motif$START > -80 & prom_motif$START < -50 ])
)

# Croiser les filtres
MOTIFxAUTOG = Crossinglist(MOTIF_uniq, AUTOGAMY)

MOTIFxUP_PKX = Crossinglist(MOTIF_uniq, UP_PKX)
MOTIFxUP_PKXxAUTOG = Crossinglist(MOTIFxUP_PKX, AUTOGAMY)

MOTIFxCTIP = Crossinglist(MOTIF_uniq, stdCTIP)
MOTIFxCTIPxAUTOG = Crossinglist(MOTIFxCTIP, AUTOGAMY)

MOTIFxnotUP_PKX = Crossinglist(MOTIF_uniq, not_UP_PKX)
MOTIFxnotUP_PKXxAUTOG = Crossinglist(MOTIFxnotUP_PKX, AUTOGAMY)

MOTIFxnotCTIP = Crossinglist(MOTIF_uniq, not_stdCTIP)
MOTIFxnotCTIPxAUTOG = Crossinglist(MOTIFxnotCTIP, AUTOGAMY)

#### Répartition des profils parmis les gènes avec et sans motifs ####
print("Profiles repartition barplot")
path = paste0(save_path,"Barplot_profil/")
dir.create(path ,recursive=T,showWarnings=F)

# Sur motif
Profile_Barplot(MOTIF_uniq, "Motif", path)
Profile_EnrichmentBarplot(MOTIF_uniq, path, "Motif")

# Sur Motif + UP PGM KU80c & XRCC4 
Profile_Barplot(MOTIFxUP_PKX, "Motif_UP", path)
Profile_EnrichmentBarplot(MOTIFxUP_PKX, path, "Motif")

# Sur Motif + DOWN CTIP + UP PKX 
Profile_Barplot(MOTIFxCTIP, "Motif_CTIP", path)
Profile_EnrichmentBarplot(MOTIFxCTIP, path, "Motif")

#### Motif parmis les différents filtres ####
print("Histogram of motif position")
path = paste0(save_path,"Hitogramme_positions/")
dir.create(path ,recursive=T,showWarnings=F)

# Répartition des motifs sur les promoteurs

PositionHistogram(MOTIF)
PositionHistogram(MOTIFxUP_PKX)
PositionHistogram(MOTIFxCTIP)
PositionHistogram(MOTIFxnotUP_PKX)
PositionHistogram(MOTIFxnotCTIP)

PositionHistogram(MOTIFxAUTOG)
PositionHistogram(MOTIFxUP_PKXxAUTOG)
PositionHistogram(MOTIFxCTIPxAUTOG)
PositionHistogram(MOTIFxnotUP_PKXxAUTOG)
PositionHistogram(MOTIFxnotCTIPxAUTOG)

#### Diagramme de Venn ####
print("Venn Diagramm in progress")
path = paste0(save_path,"Venn_Diagramm/")
dir.create(path ,recursive=T,showWarnings=F)

# Croiser Motif avec genes de l'autoagmie
LIST = list(Autogamy = AUTOGAMY$autogamy,
            Motif =  MOTIF_uniq$Motif,
            Motif50.80 =  MOTIF_uniq$motif_50.80)
png(paste0(path,"Venn_Motif_autogamy.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

# Croiser Motif avec UP PGM, KU80c ou XRCC4
LIST = c(UP_PKX[-4],   Motif = list(MOTIF_uniq[[1]]))
png(paste0(path,"Venn_Motif_UP.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = c(UP_PKX[-4],   Motif_50.80 = list(MOTIF_uniq[["motif_50.80"]]))
png(paste0(path,"Venn_Motif50-80_UP.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(UP_PKX = UP_PKX$UP_ALL,
            IntermediatePeak = AUTOGAMY$inter_peak,
            Motif =  MOTIF_uniq$Motif)
png(paste0(path,"Venn_Motifs_UP_inter.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(UP_PKX = UP_PKX$UP_ALL,
            EarlyPeak = AUTOGAMY$early_peak,
            Motif =  MOTIF_uniq$Motif)
png(paste0(path,"Venn_Motifs_UP_early.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()


LIST = list(UP_PKX = UP_PKX$UP_ALL,
            IntermediatePeak = AUTOGAMY$inter_peak,
            Motif50.80 =  MOTIF_uniq$motif_50.80)
png(paste0(path,"Venn_Motifs50.80_UP_inter.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(UP_PKX = UP_PKX$UP_ALL,
            EarlyPeak = AUTOGAMY$early_peak,
            Motif50.80 =  MOTIF_uniq$motif_50.80)
png(paste0(path,"Venn_Motifs50.80_UP_early.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()


# Croiser Motif avec UP_PKX + DOWN CTIP 
LIST = list(UP_PKX = UP_PKX$UP_ALL,
         DOWN_CTIP = stdCTIP$DOWN_CTIP,
           Motif =  MOTIF_uniq$Motif)
png(paste0(path,"Venn_Motif_UP-CTIP.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP,
            Motif =  MOTIF_uniq$Motif,
            IntermediatePeak = AUTOGAMY$inter_peak)
png(paste0(path,"Venn_Motif_UP-CTIP_inter.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP,
            Motif =  MOTIF_uniq$Motif,
            EarlyPeak = AUTOGAMY$early_peak)
png(paste0(path,"Venn_Motif_UP-CTIP_early.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP,
            Motif.80.50 =  MOTIF_uniq$motif_50.80,
            EarlyPeak = AUTOGAMY$early_peak)
png(paste0(path,"Venn_Motif80-50_UP-CTIP_early.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP,
            Motif.80.50 =  MOTIF_uniq$motif_50.80,
            IntermediatePeak = AUTOGAMY$inter_peak)
png(paste0(path,"Venn_Motif80-50_UP-CTIP_inter.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP,
            Motif.80.50 =  MOTIF_uniq$motif_50.80,
            Turbo = TURBO$turbo_OU)
png(paste0(path,"Venn_Motif80-50_UP-CTIP_Turbo.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()

LIST = list(InterPeak = AUTOGAMY$inter_peak,
            UP_DOWN = stdCTIP$DOWN_UP,
            Motif.80.50 =  MOTIF_uniq$motif_50.80,
            Turbo = TURBO$turbo_OU)
png(paste0(path,"Venn_Motif80-50_UP-DOWN_Turbo_inter.png"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set2"))
dev.off()


#### Barplot répartition des motifs +/- ####
print("Repartition of +/- motifs")
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
colors = c("royalblue1", "gold1", "chartreuse3", "grey")

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

####Création de tableaux récapitulatifs ####
print("Creation of summary table")
summary_tab = read.table(paste0("./Analyse/",file_name, "/", condition, "/Summary_",condition,".tab"), sep = '\t', header = T)

# Ajout des info su les gènes pour chaque motifs
mini_tab = prom_motif[,c("ID","SCORE", "START", "STRAND")]
colnames(mini_tab)[-1] = paste0("Motif_",colnames(mini_tab)[-1])

summary_tab2 = merge(summary_tab, mini_tab, by = "ID", all = T)
write.table(summary_tab,paste0("Analyse/",file_name,"/",condition,"/Summary_Motif_",condition,".tab"), sep = "\t", row.names = F) 

# Ajout des information du nombre de motif et de la position pour chaque gènes
nb_motif = table(MOTIF$Motif)
nb_motif = as.data.frame(nb_motif)
colnames(nb_motif)= c("ID", "Motif")
summary_tab = merge(summary_tab, nb_motif, by = "ID", all = T)
summary_tab = merge(summary_tab, strand_tab, by = "ID")
summary_tab$Motif_50_80 = is.element(summary_tab$ID, MOTIF_uniq$motif_50.80)
write.table(summary_tab,paste0("Analyse/",file_name,"/",condition,"/Summary2_",condition,".tab"), sep = "\t", row.names = F) 


# Création de séléction de gènes d'intérêt
selection_tab = summary_tab[which(is.element(summary_tab$ID, stdCTIP$DOWN_UP)),]
selection_tab = selection_tab[which(is.element(selection_tab$ID, MOTIF$motif_50.80)),]
write.table(selection_tab,paste0("Analyse/",file_name,"/",condition,"/Selection1_UP_DOWN-",condition,".tab"), sep = "\t", row.names = F)

selection_tab = selection_tab[which(str_detect(selection_tab$NAME, "PTET")),]
selection_tab = selection_tab[which(str_detect(selection_tab$SYNONYMS, "PTMB")| selection_tab$SYNONYMS == ""),]
write.table(selection_tab,paste0("Analyse/",file_name,"/",condition,"/Selection2_unknown_",condition,".tab"), sep = "\t", row.names = F)

selection_tab = selection_tab[which(is.element(selection_tab$ID, TURBO$turbo_OU)),]
write.table(selection_tab,paste0("Analyse/",file_name,"/",condition,"/Selection3_Tubo_",condition,".tab"), sep = "\t", row.names = F)

sink(paste0(save_path,"/Analyse_sessionInfo.txt"))
print(sessionInfo())
sink()
