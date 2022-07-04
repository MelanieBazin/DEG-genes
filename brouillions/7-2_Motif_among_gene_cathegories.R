# source("7_Ouverture_fichier_motif_filtres.R")

#### Répartition des profils parmis les gènes avec et sans motifs ####
print("Repartition of expression profiles barplots")
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

##### Etat de R #####
sink(paste0(path,"/Analyse_sessionInfo.txt"))
print(sessionInfo())
sink()


