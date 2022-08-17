options(stringsAsFactors = FALSE)
library(ggvenn)
library(ggplot2) 
library(RColorBrewer)
library("VennDiagram")

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
MOTIF_uniq = unique(MOTIF)

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

Candidats = intersect(UP_CTIP,MOTIF_uniq)
NotCandidats = setdiff(MOTIF_uniq, UP_CTIP)
NotCandidats_inter = intersect(INTER,setdiff(MOTIF_uniq, UP_CTIP))


#Venn with UP
LIST = list(RNAi_KU80c = KU,
            RNAi_PGM = PGM,
            RNAi_XRCC4 = XRCC4)
pdf(paste0(path,"Venn_UP.pdf"))
ggvenn(LIST,
       fill_color = c("deeppink1","dodgerblue3","ivory2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

#Venn with UP& CTIP
LIST = list(UP = UP,
            RNAi_CTIP = CTIP,
            Intermediate = INTER)
pdf(paste0(path,"Venn_UP_CTIP_inter.pdf"))
ggvenn(LIST,
       fill_color = c("plum4", "chartreuse3", "red1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

#Venn with UP& CTIP
LIST = list(UP = UP,
            RNAi_CTIP = CTIP)
pdf(paste0(path,"Venn_UP_CTIP.pdf"))
ggvenn(LIST,
       fill_color = c("plum4", "chartreuse3"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

#Venn motif intermed
LIST = list(UP_CTIP = UP_CTIP,
            Intermediate = INTER,
            Motif = MOTIF)
pdf(paste0(path,"Venn_Motif_inter_up_ctip.pdf"))
ggvenn(LIST,
       fill_color = c("chartreuse4", "red1","darkgoldenrod2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

LIST = list(UP_CTIP_inter = UP_CTIP_inter,
            Motif = MOTIF)
pdf(paste0(path,"Venn2_Motif_inter_up_ctip_CTIP.pdf"))
ggvenn(LIST,
       fill_color = c("indianred","darkgoldenrod2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()


#Venn motif intermed with motif ar right position and up deregulated
LIST = list(UP_CTIP = UP_CTIP,
            Inter = INTER,
            Motif = MOTIF_pos)
pdf(paste0(path,"Venn_Motif_pos_inter_up_ctip_CTIP.pdf"))
ggvenn(LIST,
       fill_color = c("chartreuse4", "red1","darkgoldenrod2"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

LIST = list(UP_CTIP_inter = UP_CTIP_inter,
            Motif = MOTIF_pos)
pdf(paste0(path,"Venn2_Motif_pos_inter_up_ctip_CTIP.pdf"))
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
pdf(paste0(path,"Venn_Motif_inter_UP.pdf"))
ggvenn(LIST,
       fill_color = c("darkgoldenrod2", "red1"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()

#Venn motif intermed with motif ar right position and up deregulated
LIST = list(Motif = MOTIF_pos,
            Inter = INTER,
            UP = UP)
pdf(paste0(path,"Venn_Motif_pos_inter_up.pdf"))
ggvenn(LIST,
       fill_color = c("darkgoldenrod2", "red1","plum4"),
       stroke_size = 0.5, set_name_size = 6, show_percentage = F, text_size = 10)
dev.off()


#### Venn porportionels ####
# UP in PGM, KU80c, XRCC4
pdf(paste0(path,"VennProp_UP.pdf"))
grid.newpage() 
overrideTriple=T
draw.triple.venn(area1 = length(KU),
                 area2 = length(PGM),
                 area3 = length(XRCC4),
                 n12 = length(intersect(KU, PGM)),
                 n23 = length(intersect(PGM, XRCC4)),
                 n13 = length(intersect(KU, XRCC4)),
                 n123 = length(intersect(intersect(KU, PGM),XRCC4)),
                 category = c("KU80c RNAi", "PGM RNAi", "XRCC4 RNAi"),
                 col = c("deeppink1","dodgerblue3","ivory2"),
                 fill = c("deeppink1","dodgerblue3","ivory2"),
                 alpha = rep(0.5, 3),
                 scaled = TRUE)
dev.off()

# UP in PGM, KU80c, XRCC4 & DOWN CTIP
pdf(paste0(path,"VennProp_UP_CTIP.pdf"))
grid.newpage() 
draw.pairwise.venn(area1 = length(UP),
                   area2 = length(CTIP),
                   cross.area = length(intersect(UP, CTIP)),
                   category = c("UP", "CTIP"),
                   col = c("plum4", "chartreuse3"),
                   fill = c("plum4", "chartreuse3"),
                   alpha = rep(0.5, 2),
                   scaled = TRUE)
dev.off()

# UP in PGM, KU80c, XRCC4 & DOWN CTIP + Motif + Intermediate peak
pdf(paste0(path,"VennProp_Motif_inter_up_ctip.pdf"))
grid.newpage() 
overrideTriple=T
draw.triple.venn(area1 = length(UP_CTIP),
                 area2 = length(INTER),
                 area3 = length(MOTIF),
                 n12 = length(intersect(UP_CTIP, INTER)),
                 n23 = length(intersect(INTER, MOTIF)),
                 n13 = length(intersect(UP_CTIP, MOTIF)),
                 n123 = length(intersect(intersect(UP_CTIP, INTER),MOTIF)),
                 category = c("UP_CTIP", "IntermediatePeak", "Motif"),
                 col = c("chartreuse4", "red1","darkgoldenrod2"),
                 fill = c("chartreuse4", "red1","darkgoldenrod2"),
                 alpha = rep(0.5, 3),
                 scaled = TRUE)
dev.off()

pdf(paste0(path,"VennProp2_Motif_inter_up_ctip_CTIP.pdf"))
grid.newpage() 
draw.pairwise.venn(area1 = length(UP_CTIP_inter),
                   area2 = length(MOTIF),
                   cross.area = length(intersect(UP_CTIP_inter, MOTIF)),
                   category = c("UP_CTIP_inter", "Motif"),
                   col = c("indianred","darkgoldenrod2"),
                   fill = c("indianred","darkgoldenrod2"),
                   alpha = rep(0.5, 2),
                   scaled = TRUE)
dev.off()


pdf(paste0(path,"VennProp_FrapportiTRASH.pdf"))
grid.newpage() 
draw.pairwise.venn(area1 = (1358+1579),
                   area2 = (1358+147),
                   cross.area = 1358,
                   category = c("My_data", "Frapporti"),
                   col = c("dodgerblue", "gold1"),
                   fill = c("dodgerblue", "gold1"),
                   alpha = rep(0.5, 2),
                   scaled = TRUE)
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

pdf(paste0(path,"BarPlot_Motif_position.pdf"))
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

pdf(paste0(path,"BarPlot_Motif_orientation_UP_CTIP_inter.pdf"),width = 240, height = 480)
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

pdf(paste0(path,"BarPlot_Motif_orientation_INTER.pdf"),width = 240, height = 480)
barplot(t(tab2),
        col = c("goldenrod2","olivedrab","cyan4"))
dev.off()

#Bar plot motif orientation all candidats et inter
tab = c(length(MOTIF_plus),
        length(MOTIF_both), 
        length(MOTIF_minus),
        length(MOTIF))
tab = rbind(tab,c(length(intersect(MOTIF_plus, NotCandidats)),
                  length(intersect(MOTIF_both, NotCandidats)),
                  length(intersect(MOTIF_minus, NotCandidats)),
                  length(intersect(MOTIF,NotCandidats))))
tab = rbind(tab,c(length(intersect(MOTIF_plus, Candidats)),
                  length(intersect(MOTIF_both, Candidats)),
                  length(intersect(MOTIF_minus, Candidats)),
                  length(intersect(MOTIF,Candidats))))
tab = rbind(tab,c(length(intersect(MOTIF_plus, NotCandidats_inter)),
                  length(intersect(MOTIF_both, NotCandidats_inter)),
                  length(intersect(MOTIF_minus, NotCandidats_inter)),
                  length(intersect(MOTIF,NotCandidats_inter))))

colnames(tab)= c("plus","both","minus","sum")
rownames(tab)= c("ALL","NotCandidats","Candidats","NotCandidats_inter")
write.table(tab, paste0(path,"Motif_orientation_SUPP2.tab"), sep = "\t")

tab2 = tab[,1:3]/tab[,"sum"]*100
write.table(tab2, paste0(path,"Motif_orientation_SUPP2_prct.tab"), sep = "\t")

pdf(paste0(path,"BarPlot_Motif_orientation_SUPP2.pdf"),width = 5, height = 8)
barplot(t(tab2[2:4,]),
        col = c("goldenrod2","olivedrab","cyan4"),
        ylab = "% of genes with motif")
dev.off()


Candidats = intersect(UP_CTIP,MOTIF_uniq)
NotCandidats = setdiff(MOTIF_uniq, UP_CTIP)
NotCandidats_inter = intersect(INTER,setdiff(MOTIF_uniq, UP_CTIP))
