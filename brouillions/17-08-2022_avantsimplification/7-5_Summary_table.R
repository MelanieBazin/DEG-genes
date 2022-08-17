# source("7_Ouverture_fichier_motif_filtres.R")

####Création de tableaux récapitulatifs ####
print("Creation of summary table")
summary_tab = read.table(paste0("./Analyse/",file_name, "/", condition, "/Summary_",condition,".tab"), sep = '\t', header = T)

# Ajout des info su les gènes pour chaque motifs
mini_tab = prom_motif[,c("ID","SCORE", "START", "STRAND")]
colnames(mini_tab)[-1] = paste0("Motif_",colnames(mini_tab)[-1])

summary_tab2 = merge(summary_tab, mini_tab, by = "ID", all = T)
write.table(summary_tab,paste0(save_path,"/Summary_Motif_",condition,".tab"), sep = "\t", row.names = F) 

# Ajout des information du nombre de motif et de la position pour chaque gènes
nb_motif = table(MOTIF$Motif)
nb_motif = as.data.frame(nb_motif)
colnames(nb_motif)= c("ID", "Motif")
summary_tab = merge(summary_tab, nb_motif, by = "ID", all = T)
summary_tab = merge(summary_tab, strand_tab, by = "ID")
summary_tab$Motif_50_70 = is.element(summary_tab$ID, MOTIF_uniq$motif_pos)
write.table(summary_tab,paste0(save_path,"/Summary2_",condition,".tab"), sep = "\t", row.names = F) 


# Création de séléction de gènes d'intérêt
selection_tab = summary_tab[which(is.element(summary_tab$ID, stdCTIP$UP_ALL)),]
selection_tab = selection_tab[which(is.element(selection_tab$ID, MOTIF$motif_pos)),]
write.table(selection_tab,paste0(save_path, "/Selection1_UP-",condition,".tab"), sep = "\t", row.names = F)

selection_tab = summary_tab[which(is.element(summary_tab$ID, stdCTIP$DOWN_UP)),]
selection_tab = selection_tab[which(is.element(selection_tab$ID, MOTIF$motif_pos)),]
write.table(selection_tab,paste0(save_path, "/Selection2_UP_DOWN-",condition,".tab"), sep = "\t", row.names = F)

selection_tab = selection_tab[which(str_detect(selection_tab$NAME, "PTET")),]
selection_tab = selection_tab[which(str_detect(selection_tab$SYNONYMS, "PTMB")| selection_tab$SYNONYMS == ""),]
write.table(selection_tab,paste0(save_path,"/Selection3_unknown_",condition,".tab"), sep = "\t", row.names = F)

selection_tab = selection_tab[which(is.element(selection_tab$ID, TURBO$turbo_OU)),]
write.table(selection_tab,paste0(save_path,"/Selection4_Tubo_",condition,".tab"), sep = "\t", row.names = F)



##### Etat de R #####
sink(paste0(path,"/Analyse_sessionInfo.txt"))
print(sessionInfo())
sink()

