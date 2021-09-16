source("6_Filtres.R")
library(ggvenn)

#### Faire diagramme de Venne pour les gènes selectionnée comme étant la population de gènes positif pour STREM ####
selection = list(UP_PGM, UP_KU80c, UP_XRCC4, Pic_intermediaire)
names(selection)=c("UP_PGM", "UP_KU80c", "UP_XRCC4", "Pic_intermediaire")
png(paste0(path,"Venn_selection.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"))
dev.off()

selection = list(UP_PGM, UP_KU80c, UP_XRCC4, UP_EZL1)
names(selection)=c("UP_PGM", "UP_KU80c", "UP_XRCC4", "UP_EZL1")
png(paste0(path,"Venn_UP.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"))
dev.off()

selection = list(UP_PKXE, Pic_intermediaire)
names(selection)=c("UP_PKXE","Pic_intermediaire")
png(paste0(path,"Venn_UP_inter.png"))
ggvenn(selection,
       fill_color = c("darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("darkorange","darkolivegreen3"))
dev.off()

selection = list(UP_PKXE, Pic_early)
names(selection)=c("UP_PKXE","Pic_early")
png(paste0(path,"Venn_UP_early.png"))
ggvenn(selection,
       fill_color = c("darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("darkorange","darkolivegreen3"))
dev.off()

selection = list(UP_PKXE, Early_rep)
names(selection)=c("UP_PKXE","Early_rep")
png(paste0(path,"Venn_UP_early_rep.png"))
ggvenn(selection,
       fill_color = c("darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("darkorange","darkolivegreen3"))
dev.off()

selection = list(UP_PKXE, DOWN_CTIP, Pic_intermediaire)
names(selection)=c("UP_PKXE", "DOWN_CTIP", "Pic_intermediaire")
png(paste0(path,"Venn_UP_DOWN_inter.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 7,
       show_percentage = F,
       text_size = 7,
       set_name_color = c("#0073C2FF", "darkorange","darkolivegreen3"))
dev.off()

selection = list(UP_PKXE, DOWN_CTIP, Pic_early)
names(selection)=c("UP_PKXE", "DOWN_CTIP", "Pic_early")
png(paste0(path,"Venn_UP_DOWN_early.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 7,
       show_percentage = F,
       text_size = 7,
       set_name_color = c("#0073C2FF", "darkorange","darkolivegreen3"))
dev.off()

selection = list(UP_PKXE, DOWN_CTIP, Early_rep)
names(selection)=c("UP_PKXE", "DOWN_CTIP", "Early_rep")
png(paste0(path,"Venn_UP_DOWN_early_rep.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 7,
       show_percentage = F,
       text_size = 7,
       set_name_color = c("#0073C2FF", "darkorange","darkolivegreen3"))
dev.off()

selection = list(DOWN_CTIP, Early_rep)
names(selection)=c( "DOWN_CTIP", "Early_rep")
png(paste0(path,"Venn_DOWN_early_rep.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 7,
       show_percentage = F,
       text_size = 7,
       set_name_color = c("#0073C2FF", "darkorange","darkolivegreen3"))
dev.off()



selection = list(UP_PKXE, DOWN_CTIP)
names(selection)=c("UP_PKXE", "DOWN_CTIP")
png(paste0(path,"Venn_UP_DOWN.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("#0073C2FF", "darkorange"))
dev.off()

#### Faire diagramme de Venne pour les gènes non-selectionnée ####
selection = list(DOWN_PGM, DOWN_KU80c, DOWN_XRCC4, Pic_intermediaire)
names(selection)=c("DOWN_PGM", "DOWN_KU80c", "DOWN_XRCC4", "Pic_intermediaire")
png(paste0(path,"Venn_DOWN_selection.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"))
dev.off()

selection = list(DOWN_PGM, DOWN_KU80c, DOWN_XRCC4, DOWN_EZL1)
names(selection)=c("DOWN_PGM", "DOWN_KU80c", "DOWN_XRCC4", "DOWN_EZL1")
png(paste0(path,"Venn_DOWN.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"))
dev.off()

selection = list(DOWN_PKXE, Pic_intermediaire)
names(selection)=c("DOWN_PKXE","Pic_intermediaire")
png(paste0(path,"Venn_DOWN_inter.png"))
ggvenn(selection,
       fill_color = c("darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("darkorange","darkolivegreen3"))
dev.off()

selection = list(DOWN_PKXE, UP_CTIP, Pic_intermediaire)
names(selection)=c("DOWN_PKXE", "UP_CTIP", "Pic_intermediaire")
png(paste0(path,"Venn_DOWN_UP_inter.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 7,
       show_percentage = F,
       text_size = 7,
       set_name_color = c("#0073C2FF", "darkorange","darkolivegreen3"))
dev.off()

selection = list(DOWN_PKXE, UP_CTIP, Pic_early)
names(selection)=c("DOWN_PKXE", "UP_CTIP", "Pic_early")
png(paste0(path,"Venn_DOWN_UP_early.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 7,
       show_percentage = F,
       text_size = 7,
       set_name_color = c("#0073C2FF", "darkorange","darkolivegreen3"))
dev.off()

selection = list(DOWN_PKXE, UP_CTIP, Early_rep)
names(selection)=c("DOWN_PKXE", "UP_CTIP", "Early_rep")
png(paste0(path,"Venn_DOWN_UP_early-rep.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 7,
       show_percentage = F,
       text_size = 7,
       set_name_color = c("#0073C2FF", "darkorange","darkolivegreen3"))
dev.off()


selection = list(DOWN_PKXE, UP_CTIP)
names(selection)=c("DOWN_PKXE", "UP_CTIP")
png(paste0(path,"Venn_DOWN_UP.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("#0073C2FF", "darkorange"))
dev.off()


