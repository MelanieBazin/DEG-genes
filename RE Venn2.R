options(stringsAsFactors = FALSE)
library(ggvenn)
library(ggplot2) 
library(RColorBrewer)

path = "./Analyse/2022-02-21_Analyse_DESeq2_FC-1.5_pval-0.05/analyseDE/NewGraphs_UP_CTIP_inter/"
# path = "./NewGraphs/"
dir.create(path ,recursive=T,showWarnings=F)


file = "./Analyse/2022-02-21_Analyse_DESeq2_FC-1.5_pval-0.05/analyseDE/Motif/From_TSS_IN_MAC/UP_CTIP_inter/FIMO_1E-4/p-value_1.4e-05/Summary2_analyseDE.tab"
# file = "./Summary2_analyseDE.tab"

summary = read.table(file, header = T, sep = "\t")

file2 = "./Analyse/2022-02-21_Analyse_DESeq2_FC-1.5_pval-0.05/analyseDE/Motif/From_TSS_IN_MAC/UP_CTIP_inter/FIMO_1E-4/p-value_1.4e-05/Motifs_1.4e-05_-70_-50.tab"
motif_pos = read.table(file2, header = T, sep = "\t")


KU = summary$ID[which(summary$KU80c_REG == "Up-regulated")]
PGM = summary$ID[which(summary$PGM_REG == "Up-regulated")]
XRCC4 = summary$ID[which(summary$XRCC4_REG == "Up-regulated")]
CTIP_early = summary$ID[which(summary$CTIP_early_REG == "Down-regulated")]
CTIP_inter = summary$ID[which(summary$CTIP_inter_REG == "Down-regulated")]

UP = intersect(intersect(KU, PGM), XRCC4)
CTIP = unique(c(CTIP_early, CTIP_inter))

UP_CTIP = intersect(UP, CTIP)

INTER = summary$ID[which(summary$EXPRESSION_PROFIL == "Intermediate peak")]
not_INTER = summary$ID[which(!is.element(summary$ID, INTER))]

UP_CTIP_inter = intersect(UP_CTIP, INTER)
write(UP_CTIP_inter, paste0(path,"ID_UPpkx_DOWNc_inter.txt"))

notUP_CTIP_inter = setdiff(summary$ID, UP_CTIP_inter)


MOTIF = summary$ID[which(!is.na(summary$Motif))]
MOTIF_pos = summary$ID[which(summary$Motif_50_70 == TRUE)]
MOTIF_notPos = setdiff(MOTIF,MOTIF_pos)

MOTIF_plus = summary$ID[which(summary$Motif_STRAND == "plus")]
MOTIF_minus = summary$ID[which(summary$Motif_STRAND == "moins")]
MOTIF_both = summary$ID[which(summary$Motif_STRAND == "both")]

MOTIF_pos_inter = intersect(MOTIF_pos,UP_CTIP_inter)

MOTIF_pos_plus = unique(motif_pos$ID[which(motif_pos$STRAND == "+")])
MOTIF_pos_moins = unique(motif_pos$ID[which(motif_pos$STRAND == "-")])
MOTIF_pos_both = intersect(MOTIF_pos_plus, MOTIF_pos_moins)
MOTIF_pos_plus = setdiff(MOTIF_pos_plus, MOTIF_pos_both)
MOTIF_pos_moins = setdiff(MOTIF_pos_moins, MOTIF_pos_both)

#Venn with UP
LIST = list(RNAi_KU80c = KU,
            RNAi_PGM = PGM,
            RNAi_XRCC4 = XRCC4)
png(paste0(path,"Venn_UP.png"))
ggvenn(LIST,
       fill_color = c("deeppink1","dodgerblue3","ivory2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()



#Venn with UP& CTIP
LIST = list(UP = UP,
            RNAi_CTIP = CTIP,
            Intermediate = INTER)
png(paste0(path,"Venn_UP_CTIP_inter.png"))
ggvenn(LIST,
       fill_color = c("plum4", "chartreuse3", "red1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

#Venn with UP& CTIP
LIST = list(UP = UP,
            RNAi_CTIP = CTIP)
png(paste0(path,"Venn_UP_CTIP.png"))
ggvenn(LIST,
       fill_color = c("plum4", "chartreuse3"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

#Venn motif intermed
LIST = list(UP_CTIP = UP_CTIP,
            Intermediate = INTER,
            Motif = MOTIF)
png(paste0(path,"Venn_Motif_inter_up_ctip.png"))
ggvenn(LIST,
       fill_color = c("chartreuse4", "red1","darkgoldenrod2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

LIST = list(UP_CTIP_inter = UP_CTIP_inter,
            Motif = MOTIF)
png(paste0(path,"Venn2_Motif_inter_up_ctip_CTIP.png"))
ggvenn(LIST,
       fill_color = c("indianred","darkgoldenrod2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()


#Venn motif intermed with motif ar right position and up deregulated
LIST = list(UP_CTIP = UP_CTIP,
            Inter = INTER,
            Motif = MOTIF_pos)
png(paste0(path,"Venn_Motif_pos_inter_up_ctip_CTIP.png"))
ggvenn(LIST,
       fill_color = c("chartreuse4", "red1","darkgoldenrod2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

LIST = list(UP_CTIP_inter = UP_CTIP_inter,
            Motif = MOTIF_pos)
png(paste0(path,"Venn2_Motif_pos_inter_up_ctip_CTIP.png"))
ggvenn(LIST,
       fill_color = c("indianred","darkgoldenrod2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

A = intersect(intersect(UP_CTIP,INTER),MOTIF_pos)
write.table(summary[which(is.element(summary$ID,A)),],paste0(path,"Venn_Motif_pos_inter_up_ctip_CTIP.tab"), sep = "\t" , row.names = F)

#Venn motif intermed UP
LIST = list(Motif = MOTIF,
            Intermediate = INTER,
            UP = UP)
png(paste0(path,"Venn_Motif_inter_UP.png"))
ggvenn(LIST,
       fill_color = c("darkgoldenrod2", "red1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

#Venn motif intermed with motif ar right position and up deregulated
LIST = list(Motif = MOTIF_pos,
            Inter = INTER,
            UP = UP)
png(paste0(path,"Venn_Motif_pos_inter_up.png"))
ggvenn(LIST,
       fill_color = c("darkgoldenrod2", "red1","plum4"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

#### Barplot Motif position proportion ####
tab = c(length(MOTIF_pos), length(MOTIF_notPos), length(MOTIF))
tab = rbind(tab,c(length(intersect(MOTIF_pos, UP_CTIP_inter)),length(intersect(MOTIF_notPos, UP_CTIP_inter)), length(intersect(MOTIF,UP_CTIP_inter))))
tab = rbind(tab,c(length(intersect(MOTIF_pos, notUP_CTIP_inter)),length(intersect(MOTIF_notPos, notUP_CTIP_inter)),length(intersect(MOTIF,notUP_CTIP_inter))))
colnames(tab)= c("-50,-70","out","sum")
rownames(tab)= c("ALL","UP_CTIP_inter", "notUP_CTIP_inter")
write.table(tab, paste0(path,"Motif_pos.tab"), sep = "\t")

tab2 = tab[,1:2]/tab[,"sum"]*100
write.table(tab2, paste0(path,"Motif_pos_prct.tab"), sep = "\t")

png(paste0(path,"BarPlot_Motif_position.png"))
barplot(t(tab2),
        col = c("goldenrod3","moccasin"))
dev.off()

#Bar plot motif orientation UP+CTIP+INTER
tab = c(length(MOTIF_plus),
        length(MOTIF_both), 
        length(MOTIF_minus),
        length(MOTIF))
# tab =  rbind(tab,c(length(MOTIF_pos_plus)
#                    ,length(MOTIF_pos_both),
#                    length(MOTIF_pos_moins),
#                    length(MOTIF_pos)))
tab = rbind(tab,c(length(intersect(MOTIF_plus, UP_CTIP_inter)),
                  length(intersect(MOTIF_both, UP_CTIP_inter)),
                  length(intersect(MOTIF_minus, UP_CTIP_inter)),
                  length(intersect(MOTIF,UP_CTIP_inter))))
# tab = rbind(tab,c(length(intersect(MOTIF_pos_plus, UP_CTIP_inter)),
#                   length(intersect(MOTIF_pos_both, UP_CTIP_inter)),
#                   length(intersect(MOTIF_pos_moins, UP_CTIP_inter)),
#                   length(intersect(MOTIF_pos,UP_CTIP_inter))))
colnames(tab)= c("plus","both","minus","sum")
rownames(tab)= c("ALL","ALL_pos","UP_CTIP_inter","UP_CTIP_inter_pos")
rownames(tab)= c("ALL","UP_CTIP_inter")
write.table(tab, paste0(path,"Motif_orientation_UP_CTIP_inter.tab"), sep = "\t")

tab2 = tab[,1:3]/tab[,"sum"]*100
write.table(tab2, paste0(path,"Motif_orientation_UP_CTIP_inter_prct.tab"), sep = "\t")

png(paste0(path,"BarPlot_Motif_orientation_UP_CTIP_inter.png"),width = 240, height = 480)
barplot(t(tab2),
        col = c("goldenrod2","olivedrab","cyan4"))
dev.off()

#Bar plot motif orientation INTER
tab = c(length(MOTIF_plus),
        length(MOTIF_both), 
        length(MOTIF_minus),
        length(MOTIF))
# tab =  rbind(tab,c(length(MOTIF_pos_plus)
#                    ,length(MOTIF_pos_both),
#                    length(MOTIF_pos_moins),
#                    length(MOTIF_pos)))
tab = rbind(tab,c(length(intersect(MOTIF_plus, INTER)),
                  length(intersect(MOTIF_both, INTER)),
                  length(intersect(MOTIF_minus, INTER)),
                  length(intersect(MOTIF,INTER))))
# tab = rbind(tab,c(length(intersect(MOTIF_pos_plus, INTER)),
#                   length(intersect(MOTIF_pos_both, INTER)),
#                   length(intersect(MOTIF_pos_moins, INTER)),
#                   length(intersect(MOTIF_pos,INTER))))
colnames(tab)= c("plus","both","minus","sum")
rownames(tab)= c("ALL","ALL_pos","INTER","INTER_pos")
rownames(tab)= c("ALL","INTER")
write.table(tab, paste0(path,"Motif_orientation_INTER.tab"), sep = "\t")

tab2 = tab[,1:3]/tab[,"sum"]*100
write.table(tab2, paste0(path,"Motif_orientation_INTER_prct.tab"), sep = "\t")

png(paste0(path,"BarPlot_Motif_orientation_INTER.png"),width = 240, height = 480)
barplot(t(tab2),
        col = c("goldenrod2","olivedrab","cyan4"))
dev.off()

