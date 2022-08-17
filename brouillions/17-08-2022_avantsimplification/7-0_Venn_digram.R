path = paste0(save_path,"Venn_Diagramm/")
dir.create(path ,recursive=T,showWarnings=F)

col2 = brewer.pal(n = 3, name = "Set1")
col2 = col2[-3]

# Croiser Motif avec genes de l'autoagmie
LIST = list(IntermediatePeak = AUTOGAMY$inter_peak,
            Motif =  prom_motif$ID)
pdf(paste0(path,"Venn_Motif_inter_10-4.pdf"))
ggvenn(LIST,
       fill_color = col2,
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = col2)
dev.off()


# Croiser Motif avec genes de l'autoagmie
LIST = list(IntermediatePeak = AUTOGAMY$inter_peak,
            Motif =  prom_motif2$ID)
pdf(paste0(path,"Venn_Motif_inter_",p_value,".pdf"))
ggvenn(LIST,
       fill_color = col2,
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = col2)
dev.off()

LIST = list(IntermediatePeak = AUTOGAMY$inter_peak,
            Motif =  prom_motif2$ID,
            UP_PKX = up_pkx)
pdf(paste0(path,"Venn_Motif_inter_UP_",p_value,".pdf"))
ggvenn(LIST,
       fill_color = brewer.pal(n = length(LIST), name = "Set1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 7,
       set_name_color = brewer.pal(n = length(LIST), name = "Set1"))
dev.off()
