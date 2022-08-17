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
additional_folder = "/UP_CTIP_inter"

# Localiser les données
file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
save_path = paste0("./Analyse/",file_name, "/", condition, "/Motif/From_",debut, "_IN_MAC",IES,additional_folder,"/FIMO_",p_valueFIMO, "/")

# prom_motif = read.table(paste0(save_path, "/p-value_1.4e-05/Motifs_1.4e-05.tab"), sep = "\t", header = T)

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
}else if (additional_folder == "/UP_CTIP_inter"){
  enrichment = mean(112/187, 112/187, 86/187, 85/187, 90/187)
}

query = query = intersect(AUTOGAMY$inter_peak, stdCTIP$UP_ALL)
prom_inter = prom_motif[which(is.element(prom_motif$ID, query)),]
prom_inter = prom_inter[order(prom_inter$p.value),]
prom_inter = prom_inter[!duplicated(prom_inter$ID),]
nb_selected = round(length(query)*enrichment)
p_value = round(prom_inter$p.value[nb_selected], digit = 7)

# Check the enrichment obtain
prom_motif2 = prom_motif[prom_motif$p.value <= p_value,]
prom_inter = unique(prom_motif2$ID[which(is.element(prom_motif2$ID, query))])
enrichment = length(prom_inter)/length(query)*100


# Venn Diagram
source("7-0_Venn_digram.R")


# Restrain the data to the significant ones
save_path = paste0(save_path, "p-value_",p_value, "/")
prom_motif = prom_motif2
write.table(prom_motif, paste(save_path, "Motifs_",p_value,".tab"), sep = "\t", row.names = F)


#### Création de filtre supplémentaire pour les motifs ####
print("Addition of new filter")
source("7-1_New_filter.R")

# Barplot enrichissement
tab = c(length(SUPP2$NotCandidats), length(setdiff(AUTOGAMY$all_genes, stdCTIP$DOWN_UP)))
tab = rbind(tab,c(length(SUPP2$Candidats), length(stdCTIP$DOWN_UP)))
tab = rbind(tab,c(length(SUPP2$NotCandidats_inter), length(setdiff(AUTOGAMY$inter_peak, stdCTIP$DOWN_UP))))
colnames(tab) = c("Motif", "Total")
rownames(tab) = c("NotCandidats", "Candidats" ,"NotCandidats_inter")
tab[,1] =  tab[,1]/tab[,"Total"]*100

pdf(paste0(save_path,"Motif_enrichment.pdf"), width = 2.5, height = 6)
barplot(tab[,1],
        ylim = c(0,100),
        ylab = "% of gene with motif",
        col = "darkgoldenrod2",
        las =  2)
dev.off()


#### Position motif ####
source("7-2_Analyse_motif_position.R")

# position enrichment
a = -70
b = -50

prom_motif_pos = prom_motif[prom_motif$START > a & prom_motif$START < b ,]
write.table(prom_motif_pos, paste0(save_path, "Motifs_",p_value,"_",a,"_",b,".tab"), sep = "\t", row.names = F)

MOTIF = c(MOTIF, list(motif_pos = prom_motif$ID[prom_motif$START > a & prom_motif$START < b ]))
MOTIF_uniq = c(MOTIF_uniq, list(motif_pos = unique(prom_motif$ID[prom_motif$START > a & prom_motif$START < b ])))

source("7-3_Motif_among_gene_cathegories.R")
source("7-4_Motif_strand.R")
source("7-5_Summary_table.R")
