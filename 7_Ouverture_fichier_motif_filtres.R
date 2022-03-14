options(stringsAsFactors = FALSE)
library(ggvenn)
library(ggplot2) 
library(RColorBrewer)
source("0_Cluster.R")
source("0_Visualisation_fonction.R")
source("0_Stat_function.R")

# Definition des fichier promoteur à ouvrir
IES = NULL
debut = "TSS"

# Definitir les fichiers d'analyse à ouvrir
date = Sys.Date()
date = "2022-02-21"
condition =  names(rnai_list)[2]
p_valueFIMO = "1E-4"
additional_folder = "/UP_inter"

# Localiser les donner
file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
save_path = paste0("./Analyse/",file_name, "/", condition, "/Motif/From_",debut, "_IN_MAC",IES,additional_folder,"/FIMO_",p_valueFIMO, "/")

# Ouvrir les filtres sur les dérégulation
RNAi = rnai_list[[condition]]
RNAi = RNAi[-grep("bis", RNAi)]
RNAi = RNAi[-grep("ICL7", RNAi)]
RNAi = RNAi[-grep("ND7", RNAi)]
source("5-1_Filtres.R")

# Open the fimo.tsv file
prom_motif = read.table(paste0(save_path,"fimo.tsv"), sep = "\t", header = T)
colnames(prom_motif) = c("MOTIF_ID", "MOTIF_alt_ID", "ID","START","END","STRAND","SCORE", "p.value", "q.value", "MOTIF_SEQUENCE")
prom_motif$START = prom_motif$START-150
prom_motif$END= prom_motif$END-150

# Saving of the annotated file
prom_motif2 = merge(prom_motif, annotation[,c("ID","NAME", "SYNONYMS")], by = "ID")
write.table(prom_motif2,paste0(save_path,"/fimoTSV_merge.tab"), sep = "\t", row.names = F)

#### Find the p_value corresponding to the enrichment given by STREME ####
if (additional_folder == "/UP_inter"){
  enrichment = mean(97/215, 103/215, 95/215, 98/215, 95/215)
}else if (additional_folder == "/neg_intermediate_notUP2"){
  enrichment = mean(102/215, 98/215, 80/215, 95/215, 135/215)
}
query = intersect(AUTOGAMY$inter_peak, stdCTIP$UP_ALL)
prom_UP_motif = prom_motif[which(is.element(prom_motif$ID, query)),]
prom_UP_motif = prom_UP_motif[order(prom_UP_motif$p.value),]
prom_UP_motif = prom_UP_motif[!duplicated(prom_UP_motif$ID),]
nb_selected = round(length(query)*enrichment)
p_value = round(prom_UP_motif$p.value[nb_selected], digit = 6)

# Check the enrichmùent obtain
prom_motif2 = prom_motif[prom_motif$p.value <= p_value,]
prom_UP_motif = unique(prom_motif2$ID[which(is.element(prom_motif2$ID, query))])
enrichment = length(prom_UP_motif)/length(query)*100
enrichment

# Restrain the data to the significative ones
save_path = paste0(save_path, "p-value_",p_value, "/")
prom_motif = prom_motif[prom_motif$p.value <= p_value,]

#### Création de filtre supplémentaire pour les motifs ####
print("Addition of new filter")

UP_PKX = c(UP_PKX, UP_ALL = list(up_pkx))

not_UP_PKX = c(not_UP_PKX, not_up_PKX = list(not_up_pkx))

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

SUPP = list(
  Inter_motif = intersect(AUTOGAMY$inter_peak,MOTIF_uniq$Motif),
  not_Inter_motif = setdiff(MOTIF_uniq$Motif,AUTOGAMY$inter_peak),
  Inter_UP_motif = intersect(intersect(AUTOGAMY$inter_peak,UP_PKX$UP_ALL), MOTIF_uniq$Motif),
  Inter_UP_ssmotif = setdiff(intersect(AUTOGAMY$inter_peak,UP_PKX$UP_ALL), MOTIF_uniq$Motif)
)

# Croiser les filtres
MOTIFxAUTOG = Crossinglist(MOTIF_uniq["Motif"], AUTOGAMY)

MOTIFxUP_PKX = Crossinglist(MOTIF_uniq["Motif"], UP_PKX)
MOTIFxUP_PKXxAUTOG = Crossinglist(MOTIFxUP_PKX, AUTOGAMY)

MOTIFxCTIP = Crossinglist(MOTIF_uniq["Motif"], stdCTIP)
MOTIFxCTIPxAUTOG = Crossinglist(MOTIFxCTIP, AUTOGAMY)

MOTIFxnotUP_PKX = Crossinglist(MOTIF_uniq["Motif"], not_UP_PKX)
MOTIFxnotUP_PKXxAUTOG = Crossinglist(MOTIFxnotUP_PKX, AUTOGAMY)

MOTIFxnotCTIP = Crossinglist(MOTIF_uniq["Motif"], not_stdCTIP)
MOTIFxnotCTIPxAUTOG = Crossinglist(MOTIFxnotCTIP, AUTOGAMY)

### Analyses ###
code = list.files("./")[grep("7-",list.files("./"))]

for (c in code){
  print(paste("Loading =====>", c))
  source(c)
}
