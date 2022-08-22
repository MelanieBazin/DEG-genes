####
# Analysis of the motifs
####
options(stringsAsFactors = FALSE)
library(ggvenn)

source("0_Cluster.R") # Groups parameters
source("0_Functions.R") # Library of homemade function

# Define the analyse to open
date = Sys.Date()
date = "2022-08-17"
condition =  names(rnai_list)[2]

IES = NULL
debut = "TSS"

source("5-1_Filters_candidats.R")

save_path = paste0(path, "Motif/From_",debut, "_IN_MAC",IES,"/")

# Open the fimo.tsv file
prom_motif_raw = read.table(paste0(save_path,"/FIMO/fimo.tsv"), sep = "\t", header = T)
colnames(prom_motif_raw) = c("MOTIF_ID", "MOTIF_alt_ID", "ID","START","END","STRAND","SCORE", "p.value", "q.value", "MOTIF_SEQUENCE")
prom_motif_raw$START = prom_motif_raw$START-150
prom_motif_raw$END= prom_motif_raw$END-150

#### Find the p_value corresponding to the enrichment given by STREME ####

# Calculation of the enrichment factor based on genes with motif in the 5 STREME analysis
# enrichment = mean(112/187, 112/187, 86/187, 85/187, 90/187)

minimal_genes = read.table(paste0(save_path,"/STREME/Genes_with_all_motif.txt"))[,1]
enrichment_STREME = length(minimal_genes)/length(candidats)

# Define the p_value threshold to obtain th previous enrichment
prom = prom_motif_raw[which(is.element(prom_motif_raw$ID, candidats)),]
prom = prom[order(prom$p.value),]
prom = prom[!duplicated(prom$ID),]
nb_selected = round(length(candidats)*enrichment_STREME)
p_value = round(prom$p.value[nb_selected], digit = 7)

save_path = paste0(save_path,"p-value_",p_value,"/")

prom_motif = prom_motif_raw[prom_motif_raw$p.value <= p_value,]
prom_motif = merge(annotation[,c("ID","Name", "Aliases")], prom_motif[,-(1:2)], by = "ID")
write.table(prom_motif,paste0(save_path,"fimoTSV_merge_pval-",p_value,".tab"), sep = "\t", row.names = F)

#### Venn Diagram ####
save_path2 = paste0(save_path,"Caracteristics/")
dir.create(save_path2 ,recursive=T,showWarnings=F)

# Motif among the genes of interest
a1 = DE_genes$UPpkx_Dc
a2 = AUTOGAMY$`Intermediate peak`
a3 = prom_motif$ID

pdf(paste0(save_path2,"VennProp_UPpkx.pdf"))
grid.newpage() 
overrideTriple=T
draw.triple.venn(area1 = length(a1),
                 area2 = length(a2),
                 area3 = length(a3),
                 n12 = length(intersect(a1, a2)),
                 n23 = length(intersect(a2, a3)),
                 n13 = length(intersect(a1, a3)),
                 n123 = length(intersect(intersect(a1, a2),a3)),
                 category = c(paste("UPpkx_Dc",length(a1), sep = "\n"),
                              paste("IntermediatePeak",length(a2), sep = "\n"),
                              paste("Motif",length(a3), sep = "\n")),
                 fill = c("chartreuse4", "red1","darkgoldenrod2"),
                 alpha = rep(0.5, 3),
                 scaled = TRUE)
dev.off()

#### New filters ####
# Genes with motifs
MOTIF = list(MOTIF = prom_motif$ID,
             uniqMOTIF = unique(prom_motif$ID))

ForAnalysis = list(Allx = annotation$ID[which(!is.element(annotation$ID, candidats))],
                   Candidats = candidats,
                   InterPx = AUTOGAMY$`Intermediate peak`[which(!is.element(AUTOGAMY$`Intermediate peak`, candidats))])

MOTIF.ForAnalysis = Crossinglist(list("MOTIF" = MOTIF$uniqMOTIF), ForAnalysis)

# Motif strands
MOTIF_strand = list(plus = unique(prom_motif$ID[prom_motif$STRAND == "+"]),
                    minus = unique(prom_motif$ID[prom_motif$STRAND == "-"]))

MOTIF_strand = c(MOTIF_strand,
                 both = list(intersect(MOTIF_strand$plus,MOTIF_strand$minus)))

MOTIF_strand$plus = setdiff(MOTIF_strand$plus, MOTIF_strand$both)
MOTIF_strand$minus = setdiff(MOTIF_strand$minus, MOTIF_strand$both)

MOTIF_strand.ForAnalysis = Crossinglist(MOTIF_strand, ForAnalysis)

#### Motif enrichment ####

# Barplot of the enrichment
tab = matrix(NA,length(ForAnalysis),2)
colnames(tab) = c("Motif", "Total")
rownames(tab) = names(ForAnalysis)
for (r in 1:length(ForAnalysis)){
  tab[r,] = c(length(MOTIF.ForAnalysis[[r]]), length(ForAnalysis[[r]]))
}
write.table(tab,paste0(save_path2,"Motif_enrichment_tab.tab"), sep = "\t", row.names = T)

tab[,1] =  tab[,1]/tab[,"Total"]*100

pdf(paste0(save_path2,"Motif_enrichment.pdf"), width = 2.5, height = 4)
barplot(tab[,1],
        ylim = c(0,100),
        ylab = "% of gene with motif",
        col = "darkgoldenrod2",
        las =  2)
dev.off()

# Stat of the enrichment
my_data = cbind(annotation$ID, is.element(annotation$ID, MOTIF$uniqMOTIF))

sink(paste0(save_path2,"Chi2_candidates_vs_ALL.txt"))
Chi2_pvalue(ForAnalysis, my_data)
sink()

my_data = my_data[which(is.element(my_data[,1], AUTOGAMY$`Intermediate peak`)),]

sink(paste0(save_path2,"Chi2_candidates_vs_InterP.txt"))
Chi2_pvalue(ForAnalysis, my_data)
sink()

#### Motif orientation ####
tab = matrix(NA,length(ForAnalysis),4)
colnames(tab) = c("plus","both","minus","all_motif")
rownames(tab) = names(ForAnalysis)
for (r in names(ForAnalysis)){
  tab[r,] = c(length(MOTIF_strand.ForAnalysis[[paste0("plusx",r)]]),
              length(MOTIF_strand.ForAnalysis[[paste0("bothx",r)]]), 
              length(MOTIF_strand.ForAnalysis[[paste0("minusx",r)]]), 
              length(MOTIF.ForAnalysis[[paste0("MOTIFx",r)]]))
}
write.table(tab, paste0(save_path2,"Motif_strand_candidats.tab"), sep = "\t")
tab = tab[,1:3]/tab[,"all_motif"]*100
write.table(tab, paste0(save_path2,"Motif_strand_candidats_prct.tab"), sep = "\t")

# Barplot of the proportion for each orientation
pdf(paste0(save_path2,"Motif_strand_candidats_barplot.pdf"),width = 3.7, height = 5.5)
barplot(t(tab),
        col = c("goldenrod2","olivedrab","cyan4"),
        ylab = "% of genes with motif")
dev.off()


# Add a file with the legend
pdf(paste0(save_path2,"PMotif_strand_legends.pdf"),width = 7, height = 7)
plot.new()
legend("center",
       x.intersp=0.1,
       legend = rev(rownames(tab)),
       fill = c("goldenrod2","olivedrab","cyan4"),
       bty = "n")
dev.off()

source("7-4_Motif_strand.R")

#### Motif position ####
# Plot density position
density_tab = NULL
for (n in names(ForAnalysis)){
  filtre = ForAnalysis[[n]]
  position = prom_motif$START[is.element(prom_motif$ID, filtre)]
  
  h = hist(position, breaks = 150, plot = F)
  h_tab = as.data.frame(cbind(h$mids, h$density*100))
  h_tab = cbind(h_tab, rep(n, nrow(h_tab)))
  
  density_tab = rbind(density_tab, h_tab)
}
colnames(density_tab) = c("mids","density", "group")

col_plot = c("chartreuse3","black","red")
names(col_plot) = names(ForAnalysis)

gp = ggplot(data = density_tab, aes(x = mids , y = density, color = group))+
  theme_classic() +
  scale_color_manual("Gene groups", values = col_plot) +
  labs(x = paste("Distance from", debut), y = "Motif density") +
  geom_smooth(span = 0.5,
              se = FALSE) +
  scale_x_continuous(limits = c(-150, 0), breaks = seq(-150, 0, 10))

ggsave(paste0(save_path2,"/Position_density.pdf"), device = "pdf", plot = gp,
       width = 7, height = 4, dpi = 300)

## Enrichment statistic
# Fractions of the positions
fraction = 5
pos_fraction = c()
for(l in 1:nrow(prom_motif)){
  for(n in -seq(15, 150-fraction, by=fraction)){
    if(prom_motif$START[l] < n & prom_motif$START[l] >= n-fraction){
      pos_fraction = c(pos_fraction, paste(n,n-fraction))
    }
  }
}

my_data = as.data.frame(cbind(prom_motif$ID, pos_fraction))

sink(paste0(save_path2,"/Position_enrichment_autogamy_",fraction,"frac.txt"))
Enrichment_padj(AUTOGAMY, my_data)
sink()
sink(paste0(save_path2,"/Position_enrichment_candidats_",fraction,"frac.txt"))
enrch = Enrichment_padj(ForAnalysis, my_data)
sink()

#### Set the enriched position of the motif ####
enrch_pos = enrch$Candidats
enrch_pos = rownames(enrch_pos)[enrch_pos$`pval_adj < 0.01`]
enrch_pos = paste(enrch_pos, collapse = " ")
enrch_pos = str_split(enrch_pos," ")[[1]]
enrch_pos = as.numeric(enrch_pos)

a = min(enrch_pos)
b = max(enrch_pos)

save_path2 = paste0(save_path,"Position",a,"_",b,"/")
dir.create(save_path2 ,recursive=T,showWarnings=F)

prom_motif_pos = prom_motif[prom_motif$START > a & prom_motif$START < b ,]
write.table(prom_motif_pos, paste0(save_path2, "Motifs_",p_value,"_",a,"_",b,".tab"), sep = "\t", row.names = F)

#### New filters  ####
MOTIF_pos = list(MOTIF = prom_motif_pos$ID,
                 uniqMOTIF = unique(prom_motif_pos$ID))

MOTIF_pos.ForAnalysis = Crossinglist(list("MOTIF_pos" = MOTIF_pos$uniqMOTIF), ForAnalysis)

# Motif strands
MOTIF_pos_strand = list(plus = unique(prom_motif_pos$ID[prom_motif_pos$STRAND == "+"]),
                    minus = unique(prom_motif_pos$ID[prom_motif_pos$STRAND == "-"]))

MOTIF_pos_strand = c(MOTIF_pos_strand,
                 both = list(intersect(MOTIF_pos_strand$plus,MOTIF_pos_strand$minus)))

MOTIF_pos_strand$plus = setdiff(MOTIF_pos_strand$plus, MOTIF_pos_strand$both)
MOTIF_pos_strand$minus = setdiff(MOTIF_pos_strand$minus, MOTIF_pos_strand$both)

MOTIF_str_pos.ForAnalysis = Crossinglist(MOTIF_pos_strand, ForAnalysis)

#### Motif at the position ####
tab = matrix(NA,length(ForAnalysis),3)
colnames(tab) = c(paste("Motif",a,b),"Motif", "Total")
rownames(tab) = names(ForAnalysis)
for (r in 1:length(ForAnalysis)){
  tab[r,] = c(length(MOTIF_pos.ForAnalysis[[r]]),length(MOTIF.ForAnalysis[[r]]), length(ForAnalysis[[r]]))
}
write.table(tab,paste0(save_path2,"Motif_pos_count_tab.tab"), sep = "\t", row.names = T)

#### Motif orientation at the position ####
tab = matrix(NA,length(ForAnalysis),4)
colnames(tab) = c("plus","both","minus","all_motif")
rownames(tab) = names(ForAnalysis)
for (r in names(ForAnalysis)){
  tab[r,] = c(length(MOTIF_str_pos.ForAnalysis[[paste0("plusx",r)]]),
              length(MOTIF_str_pos.ForAnalysis[[paste0("bothx",r)]]), 
              length(MOTIF_str_pos.ForAnalysis[[paste0("minusx",r)]]), 
              length(MOTIF_pos.ForAnalysis[[paste0("MOTIFx",r)]]))
}
write.table(tab, paste0(save_path2,"Motif_strand_pos_candidats.tab"), sep = "\t")
tab = tab[,1:3]/tab[,"all_motif"]*100
write.table(tab, paste0(save_path2,"Motif_strand_pos_candidats_prct.tab"), sep = "\t")

# Barplot of the proportion for each orientation
pdf(paste0(save_path2,"Motif_strand_pos_candidats_barplot.pdf"),width = 3.7, height = 5.5)
barplot(t(tab),
        col = c("goldenrod2","olivedrab","cyan4"),
        ylab = "% of genes with motif")
dev.off()

#### Summary table creation ####
summary_tab = read.table(paste0(path, "/Summary_",condition,".tab"), sep = '\t', header = T)[1:41533,]

# Add information about motifs
mini_tab = prom_motif[,c("ID","SCORE", "START", "STRAND")]
colnames(mini_tab)[-1] = paste0("Motif_",colnames(mini_tab)[-1])

summary_tab2 = merge(summary_tab, mini_tab, by = "ID", all = T)
write.table(summary_tab,paste0(save_path,"/Summary_Motif_",condition,".tab"), sep = "\t", row.names = F) 

# Addition of the number of motif founded
nb_motif = table(MOTIF$MOTIF)
nb_motif = as.data.frame(nb_motif)
colnames(nb_motif)= c("ID", "Motif")
summary_tab = merge(summary_tab, nb_motif, by = "ID", all = T)

# Addition of the motif strand
strand_tab = matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c("ID", "Motif_STRAND")))
for (S in names(MOTIF_strand)){
  tab = rep(S, length(MOTIF_strand[[S]]))
  tab = cbind(MOTIF_strand[[S]], tab)
  
  strand_tab = rbind(strand_tab,tab)
}
strand_tab = as.data.frame(strand_tab)
summary_tab = merge(summary_tab, strand_tab, by = "ID", all = T)

# Addition of a boolen to know if at least on of the motif is in the enriched position
summary_tab$Pos = is.element(summary_tab$ID, MOTIF_pos$uniqMOTIF)
colnames(summary_tab)[ncol(summary_tab)] = paste0("Motif_",a,"_",b)

write.table(summary_tab,paste0(save_path,"/Summary2_",condition,".tab"), sep = "\t", row.names = F) 

summary_tab = summary_tab[which(is.element(summary_tab$ID, candidats)),]
write.table(summary_tab,paste0(save_path,"/Summary2-candidats_",condition,".tab"), sep = "\t", row.names = F) 

#### Print R status ####
sink(paste0(save_path,"/sessionInfo.txt"))
print(sessionInfo())
sink()
