# source("7_Ouverture_fichier_motif_filtres.R")

#### Répartition des profils parmis les gènes avec et sans Motif ####
print("Repartition of expression profiles barplots")
path = paste0(save_path,"Barplot_profil/")
dir.create(path ,recursive=T,showWarnings=F)

# Sur motif
Profile_Barplot(MOTIF_uniq, "Motif", path)
Profile_EnrichmentBarplot(MOTIF_uniq, path, "Motif")

# Sur Motif + UP PGM KU70c & XRCC4 
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
            Motifpos =  MOTIF_uniq$motif_pos)
pdf(paste0(path,"Venn_Motif_autogamy.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

# Croiser Motif avec genes du pic intermediaire
LIST = list(IntermediatePeak = AUTOGAMY$inter_peak,
            Motif =  MOTIF_uniq$Motif,
            Motifpos =  MOTIF_uniq$motif_pos)

pdf(paste0(path,"Venn_Motif_inter.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()


# Croiser Motif avec UP PGM, KU70c ou XRCC4
LIST = c(UP_PKX[-4],   Motif = list(MOTIF_uniq[[1]]))
pdf(paste0(path,"Venn_Motif_UP.pdf"))
ggvenn(LIST,
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7)
dev.off()

LIST = c(UP_PKX[-4],   Motif_pos = list(MOTIF_uniq[["motif_pos"]]))
pdf(paste0(path,"Venn_Motif",a,b,"_UP.pdf"))
ggvenn(LIST,
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7)
dev.off()

LIST = list(IntermediatePeak = AUTOGAMY$inter_peak,
            Motif =  MOTIF_uniq$Motif,
            UP_PKX = UP_PKX$UP_ALL)
pdf(paste0(path,"Venn_Motif_UP_inter.pdf"))
ggvenn(LIST,
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7)
dev.off()

LIST = list(EarlyPeak = AUTOGAMY$early_peak,
            Motif =  MOTIF_uniq$Motif,
            UP_PKX = UP_PKX$UP_ALL)
pdf(paste0(path,"Venn_Motif_UP_early.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()


LIST = list(IntermediatePeak = AUTOGAMY$inter_peak,
            Motifpos =  MOTIF_uniq$motif_pos,
            UP_PKX = UP_PKX$UP_ALL)
pdf(paste0(path,"Venn_Motif",a,b,"_UP_inter.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

LIST = list(EarlyPeak = AUTOGAMY$early_peak,
            Motifpos =  MOTIF_uniq$motif_pos,
            UP_PKX = UP_PKX$UP_ALL)
pdf(paste0(path,"Venn_Motif",a,b,"_UP_early.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()


# Croiser Motif avec UP_PKX + DOWN CTIP 
LIST = list(DOWN_CTIP = stdCTIP$DOWN_CTIP,
            Motif =  MOTIF_uniq$Motif,
            UP_PKX = UP_PKX$UP_ALL)
pdf(paste0(path,"Venn_Motif_UP-CTIP.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

LIST = list(IntermediatePeak = AUTOGAMY$inter_peak,
            Motif =  MOTIF_uniq$Motif,
            UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP)
pdf(paste0(path,"Venn_Motif_UP-CTIP_inter.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

LIST = list(EarlyPeak = AUTOGAMY$early_peak,
            Motif =  MOTIF_uniq$Motif,
            UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP)
pdf(paste0(path,"Venn_Motif_UP-CTIP_early.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

LIST = list(EarlyPeak = AUTOGAMY$early_peak,
            Motifpos =  MOTIF_uniq$motif_pos,
            UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP)
pdf(paste0(path,"Venn_Motif",a,b,"_UP-CTIP_early.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

LIST = list(IntermediatePeak = AUTOGAMY$inter_peak,
            Motifpos =  MOTIF_uniq$motif_pos,
            UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP)
pdf(paste0(path,"Venn_Motif",a,b,"_CTIP_inter.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

LIST = list(Turbo = TURBO$turbo_OU,
            Motifpos =  MOTIF_uniq$motif_pos,
            UP_PKX = UP_PKX$UP_ALL,
            DOWN_CTIP = stdCTIP$DOWN_CTIP)
pdf(paste0(path,"Venn_Motif",a,b,"_UP-CTIP_Turbo.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

LIST = list(InterPeak = AUTOGAMY$inter_peak,
            Motifpos =  MOTIF_uniq$motif_pos,
            UP_DOWN = stdCTIP$DOWN_UP,
            Turbo = TURBO$turbo_OU)
pdf(paste0(path,"Venn_Motif",a,b,"_UP-DOWN_Turbo_inter.pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()

##### Etat de R #####
sink(paste0(path,"/Analyse_sessionInfo.txt"))
print(sessionInfo())
sink()


