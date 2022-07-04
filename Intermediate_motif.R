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
p_valueFIMO = "1.4e-05"
additional_folder = "/UP_CTIP_inter"

# Localiser les donner
file_name = list.files("./Analyse/")[grep(paste0(date,"_Analyse_DESeq2"),list.files("./Analyse/"))]
save_path = paste0("./Analyse/",file_name, "/", condition, "/Motif/From_",debut, "_IN_MAC",IES,additional_folder,"/FIMO_1E-4/p-value_",p_valueFIMO, "/")


summary_tab = read.table(paste0(save_path,"Summary2_",condition,".tab"), sep = '\t', header = T)
prom_motif = read.table(paste0(save_path,"Motifs_",p_valueFIMO,".tab"), sep = '\t', header = T)

MotifColor <- function(tab){
  motif_color = cbind(is.element(rownames(tab), MOTIF_uniq$Motif),
                      is.element(rownames(tab), MOTIF_uniq$motif_pos))
  motif_color = apply(motif_color, 1, sum)
  motif_color = sub("0","grey",sub("1","deepskyblue",sub("2","red", motif_color)))
}

UPdownColor <- function(tab){
  up_color = cbind(is.element(rownames(tab), UP_PKX$UP_ALL),
                   is.element(rownames(tab), unique(downCTIP$CTIP_early, downCTIP$CTIP_inter)))
  
  up_color = apply(up_color, 1, sum)
  up_color = sub("0","grey",sub("1","red",sub("2","deepskyblue", up_color)))
}

# Ouvrir les filtres sur les dérégulation
RNAi = rnai_list[[condition]]
RNAi = RNAi[-grep("bis", RNAi)]
RNAi = RNAi[-grep("ICL7", RNAi)]
RNAi = RNAi[-grep("ND7", RNAi)]
source("5-1_Filtres.R")
MOTIF_uniq = list(
  Motif = unique(prom_motif$ID),
  Motif_plus = unique(prom_motif$ID[prom_motif$STRAND == "+"]),
  motif_moins = unique(prom_motif$ID[prom_motif$STRAND == "-"]),
  motif_pos = unique(prom_motif$ID[prom_motif$START > -70 & prom_motif$START < -50 ])
)


SUPP = list(
  Inter_motif = intersect(AUTOGAMY$inter_peak,MOTIF_uniq$Motif),
  not_Inter_motif = setdiff(MOTIF_uniq$Motif,AUTOGAMY$inter_peak),
  Inter_UP_motif = intersect(intersect(AUTOGAMY$inter_peak,UP_PKX$UP_ALL), MOTIF_uniq$Motif),
  Inter_UP_ssmotif = setdiff(intersect(AUTOGAMY$inter_peak,UP_PKX$UP_ALL), MOTIF_uniq$Motif)
)

#### Comparison FC and pvalue ####
print("Comparison FC and pvalue")
path = paste0(save_path,"FC_padj_comparison/")
dir.create(path ,recursive=T,showWarnings=F)

Inter_UP_motif = summary_tab[which(is.element(summary_tab$ID, SUPP$Inter_UP_motif)),c(1,19:27)]
Inter_UP_ssmotif = summary_tab[which(is.element(summary_tab$ID, SUPP$Inter_UP_ssmotif)),c(1,19:27)]

Inter_UP_motif_FC = Inter_UP_motif[,grep("FC",colnames(Inter_UP_motif))]
Inter_UP_ssmotif_FC = Inter_UP_ssmotif[,grep("FC",colnames(Inter_UP_ssmotif))]

BoxnBarpolt_repartion(list(Inter_UP_motif_FC = Inter_UP_motif[,grep("FC",colnames(Inter_UP_motif))],
                           Inter_UP_ssmotif_FC = Inter_UP_ssmotif[,grep("FC",colnames(Inter_UP_ssmotif))]),
                      path)

Inter_UP_motif_padj = Inter_UP_motif[,grep("padj",colnames(Inter_UP_motif))]
Inter_UP_ssmotif_padj = Inter_UP_ssmotif[,grep("padj",colnames(Inter_UP_ssmotif))]

BoxnBarpolt_repartion(list(Inter_UP_motif_padj = Inter_UP_motif[,grep("padj",colnames(Inter_UP_motif))],
                           Inter_UP_ssmotif_padj = Inter_UP_ssmotif[,grep("padj",colnames(Inter_UP_ssmotif))]),
                      path)

Inter_UP_motif50.80_FC = Inter_UP_motif[which(is.element(Inter_UP_motif$ID, MOTIF$motif_pos)),grep("FC",colnames(Inter_UP_motif))]
Inter_UP_motif50.80_padj = Inter_UP_motif[which(is.element(Inter_UP_motif$ID, MOTIF$motif_pos)),grep("padj",colnames(Inter_UP_motif))]

# Histograms with the position limitation
for(i in 1:ncol(Inter_UP_motif50.80_FC)){
  
  png(paste0(path, colnames(Inter_UP_motif50.80_FC[i]),"_avecMotif_80-50.png"))
  hist(Inter_UP_motif50.80_FC[[i]],
       main =  paste0(colnames(Inter_UP_motif50.80_FC[i]),"_avecMotif"),
       breaks = 50,
       xlim = c(1,8))
  abline(v = mean(Inter_UP_motif50.80_FC[[i]]), col = "grey", lty = "dashed")
  abline(v = median(Inter_UP_motif50.80_FC[[i]]), col = "red", lty = "dashed")
  dev.off()
  
  
  png(paste0(path, colnames(Inter_UP_motif50.80_padj[i]),"_avecMotif_80-50.png"))
  hist(Inter_UP_motif50.80_padj[[i]],
       main =  paste0(colnames(Inter_UP_motif50.80_padj[i]),"_avecMotif"),
       breaks = 50,
       xlim = c(0,0.05))
  abline(v = mean(Inter_UP_motif50.80_padj[[i]]), col = "grey", lty = "dashed")
  abline(v = median(Inter_UP_motif50.80_padj[[i]]), col = "red", lty = "dashed")
  dev.off()
  
}

### PCA on FC and pvalue
INTER_UP = summary_tab[which(is.element(summary_tab$ID, intersect(AUTOGAMY$inter_peak, UP_PKX$UP_ALL))),]
rownames(INTER_UP)= INTER_UP$ID

# Color defintion by motif
motif_color = MotifColor(INTER_UP)


# Split table between fold-change ans pvalue
INTER_UP_FC = INTER_UP[,grep("FC",colnames(INTER_UP))[3:5]]
PCA_plot_generator(t(INTER_UP_FC),
                   motif_color,
                   main = "Intermediate_FC",
                   save_path = paste0(path,"PCA_FC/"),
                   label = "none")

INTER_UP_padj = INTER_UP[,grep("padj",colnames(INTER_UP))[3:5]]
PCA_plot_generator(t(INTER_UP_padj),
                   motif_color,
                   main = "Intermediate_padj",
                   save_path = paste0(path,"PCA_padj/"),
                   label = "none")

#### Comparison FC and pvalue with higher FC threshold  ####
path = paste0(save_path,"FC_padj_comparison/FC_2.5/")
dir.create(path ,recursive=T,showWarnings=F)

FC_threshold = 2.5

print(paste("Add threshold of",FC_threshold))

Inter_UP_motif_FC = Inter_UP_motif
Inter_UP_ssmotif_FC = Inter_UP_ssmotif
for (i in 1:3){
  Inter_UP_motif_FC = Inter_UP_motif_FC[Inter_UP_motif_FC[grep("FC", colnames(Inter_UP_motif_FC))[i]]>= FC_threshold,]
  Inter_UP_ssmotif_FC = Inter_UP_ssmotif_FC[Inter_UP_ssmotif_FC[grep("FC", colnames(Inter_UP_ssmotif_FC))[i]]>= FC_threshold,]
}

BoxnBarpolt_repartion(list(Inter_UP_motif_FC_FC = Inter_UP_motif_FC[,grep("FC",colnames(Inter_UP_motif_FC))],
                           Inter_UP_ssmotif_FC_FC = Inter_UP_ssmotif_FC[,grep("FC",colnames(Inter_UP_ssmotif_FC))]),
                      path)

BoxnBarpolt_repartion(list(Inter_UP_motif_FC_padj = Inter_UP_motif_FC[,grep("padj",colnames(Inter_UP_motif_FC))],
                           Inter_UP_ssmotif_FC_padj = Inter_UP_ssmotif_FC[,grep("padj",colnames(Inter_UP_ssmotif_FC))]),
                      path)

Inter_UP_motif50.80_FC = Inter_UP_motif_FC[which(is.element(Inter_UP_motif_FC$ID, MOTIF$motif_pos)),grep("FC",colnames(Inter_UP_motif_FC))]
Inter_UP_motif50.80_padj = Inter_UP_motif_FC[which(is.element(Inter_UP_motif_FC$ID, MOTIF$motif_pos)),grep("padj",colnames(Inter_UP_motif_FC))]

# Histograms with the position limitation
for(i in 1:ncol(Inter_UP_motif50.80_FC)){
  
  png(paste0(path, colnames(Inter_UP_motif50.80_FC[i]),"_avecMotif_80-50.png"))
  hist(Inter_UP_motif50.80_FC[[i]],
       main =  paste0(colnames(Inter_UP_motif50.80_FC[i]),"_avecMotif"),
       breaks = 50,
       xlim = c(1,8))
  abline(v = mean(Inter_UP_motif50.80_FC[[i]]), col = "grey", lty = "dashed")
  abline(v = median(Inter_UP_motif50.80_FC[[i]]), col = "red", lty = "dashed")
  dev.off()
  
  
  png(paste0(path, colnames(Inter_UP_motif50.80_padj[i]),"_avecMotif_80-50.png"))
  hist(Inter_UP_motif50.80_padj[[i]],
       main =  paste0(colnames(Inter_UP_motif50.80_padj[i]),"_avecMotif"),
       breaks = 50,
       xlim = c(0,0.05))
  abline(v = mean(Inter_UP_motif50.80_padj[[i]]), col = "grey", lty = "dashed")
  abline(v = median(Inter_UP_motif50.80_padj[[i]]), col = "red", lty = "dashed")
  dev.off()
  
}

INTER_UP = INTER_UP[which(is.element(INTER_UP$ID, unique(c(Inter_UP_motif_FC$ID, Inter_UP_ssmotif_FC$ID)))),]

# Color defintion by motif
motif_color = MotifColor(INTER_UP)

# Split table between fold-change ans pvalue
INTER_UP_FC = INTER_UP[,grep("FC",colnames(INTER_UP))[3:5]]
PCA_plot_generator(t(INTER_UP_FC),
                   motif_color,
                   main = "Intermediate_FC",
                   save_path = paste0(path,"PCA_FC/"),
                   label = "none")

INTER_UP_padj = INTER_UP[,grep("padj",colnames(INTER_UP))[3:5]]
PCA_plot_generator(t(INTER_UP_padj),
                   motif_color,
                   main = "Intermediate_padj",
                   save_path = paste0(path,"PCA_padj/"),
                   label = "none")

#### Comparison of intermediate peak expression profile with and without motif ####
print("Profil drawing")
path = paste0(save_path,"Profils/")
dir.create(path ,recursive=T,showWarnings=F)

ID_motif = SUPP$Inter_motif
names(ID_motif) = annotation$NAME[which(is.element(annotation$ID,ID_motif))]

ID_NOmotif = setdiff(AUTOGAMY$inter_peak,SUPP$Inter_motif)
names(ID_NOmotif) = annotation$NAME[which(is.element(annotation$ID,ID_NOmotif))]

# All conditions
pdf(paste0(path,"Expression_profile_intermed_motif.pdf"))
ExpressionProfils(type = "vst",
                  condition,
                  file = paste0("./Analyse/",file_name, "/"),
                  select_ID = ID_motif)
dev.off()


pdf(paste0(path,"Expression_profile_intermed_NO_motif.pdf"))
ExpressionProfils(type = "vst",
                  condition,
                  file = paste0("./Analyse/",file_name, "/"),
                  select_ID = ID_NOmotif)
dev.off()

# Control conditions only
pdf(paste0(path,"Expression_profile_CTRL_intermed_motif.pdf"))
ExpressionProfils(type = "vst",
                  condition,
                  file = paste0("./Analyse/",file_name, "/"),
                  select_ID = ID_motif,
                  rnai = rnai_list[[condition]][c(grep("ND7",rnai_list[[condition]]),grep("ICL7",rnai_list[[condition]]))] )
dev.off()


pdf(paste0(path,"Expression_profile_CTRL_intermed_NO_motif.pdf"))
ExpressionProfils(type = "vst",
                  condition,
                  file = paste0("./Analyse/",file_name, "/"),
                  select_ID = ID_NOmotif,
                  rnai = rnai_list[[condition]][c(grep("ND7",rnai_list[[condition]]),grep("ICL7",rnai_list[[condition]]))] )
dev.off()


# PCA and volcanoplot
print("PCA on profile and Volcanopolot")
data_path = paste0("./Analyse/",file_name, "/",condition,"/")
EXPRESSION = read.table(paste0(data_path,condition ,"_expression_table_vst.tab"))

DEG_PROFIL = summary_tab
rownames(DEG_PROFIL) = DEG_PROFIL$ID
DEG_PROFIL = DEG_PROFIL[,13:27]

# Limit defintion for volcanoplot
FC_lim = na.omit(DEG_PROFIL[,grep("FC", colnames(DEG_PROFIL))])
FC_lim = c(trunc(min(FC_lim)),ceiling(max(FC_lim)))

padj_lim = -log(na.omit(DEG_PROFIL[,grep("padj", colnames(DEG_PROFIL))]))
padj_lim = c(trunc(min(padj_lim)),ceiling(max(padj_lim)))

tgrey = rgb(t(col2rgb("grey")), max = 255, alpha = 1, names = "t.Gray")

# Profils comparison by PCA
rnai = rnai_list[[condition]]
for (ctrl in c("ND7","ICL7")){
  rnai = rnai[-grep(ctrl, rnai)]
}
cond = rnai_list[[condition]][-grep("bis", rnai_list[[condition]])]

for (r in cond){
  expression = EXPRESSION[AUTOGAMY$inter_peak,grep(r, colnames(EXPRESSION))]
  expression = na.omit(expression)
  motif_color = MotifColor(expression)
  deg_color = UPdownColor(expression)
  
  PCA_plot_generator(as.matrix(t(expression)),
                     motif_color,
                     main = r,
                     save_path = paste0(path,"PCA_profil_interP_",r,"/Motifcolor/"),
                     label = "none")
  
  PCA_plot_generator(as.matrix(t(expression)),
                     motif_color,
                     main = r,
                     save_path = paste0(path,"PCA_profil_interP_",r,"/Motifcolor/motif/"),
                     label = "none",
                     selection = SUPP$Inter_motif)
  
  PCA_plot_generator(as.matrix(t(expression)),
                     deg_color,
                     main = r,
                     save_path = paste0(path,"PCA_profil_interP_",r,"/DEGcolor/"),
                     label = "none")
  
  PCA_plot_generator(as.matrix(t(expression)),
                     deg_color,
                     main = r,
                     save_path = paste0(path,"PCA_profil_interP_",r,"/DEGcolor/motif/"),
                     label = "none",
                     selection = SUPP$Inter_motif)
  
  if (is.element(r, rnai)){
    if(r == "CTIP"){
      r = c("CTIP_early", "CTIP_inter")
    }
    for (s in r){
      deg_profil = DEG_PROFIL[AUTOGAMY$inter_peak, grep(s, colnames(DEG_PROFIL))]
      deg_profil = na.omit(deg_profil)
      motif_color2 = MotifColor(deg_profil)
      log2fc = deg_profil[,grep("log2FC", colnames(deg_profil))]
      padj = deg_profil[,grep("padj", colnames(deg_profil))]
      
      png(paste0(path,"Volcanoplot_interP_",s,".png"), width = 6, height = 6, units = 'in', res = 300)
      plot(log2fc,-log(padj),log="y",col=motif_color2,
           xlab=paste0("log2(ctrl/",s,")"),
           ylab="-log(p-value)",
           xlim = FC_lim,
           ylim = padj_lim,
           pch=20,main=s,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      dev.off()
      
      deg_profil = deg_profil[MOTIF_uniq$Motif,]
      deg_profil = na.omit(deg_profil)
      motif_color = MotifColor(deg_profil)
      log2fc = deg_profil[,grep("log2FC", colnames(deg_profil))]
      padj = deg_profil[,grep("padj", colnames(deg_profil))]
      
      png(paste0(path,"Volcanoplot_interP_",s,"_motif_only.png"), width = 6, height = 6, units = 'in', res = 300)
      plot(log2fc,-log(padj),log="y",col=motif_color,
           xlab=paste0("log2(ctrl/",s,")"),
           ylab="-log(p-value)",
           xlim = FC_lim,
           ylim = padj_lim,
           pch=20,main=s,cex=1.3,cex.axis=1.3,cex.lab=1.3)
      dev.off()
      
    }
  }
}

####



#### Profils of selected genes ####
rownames(summary_tab) = summary_tab$ID
info_data = read.table(paste0(data_path,condition ,"_infodata_collapse.tab"))
EXPRESSION = OrderColumn(EXPRESSION, info_data)
info_data = info_data[colnames(EXPRESSION),]


for (motif in c("motif", "ssmotif")){
  for (type in c("RNAi", "CTRL")){
    # Definition of the selected ID
    selection = SUPP[[paste0("Inter_UP_",motif)]]
    
    # Defintion of the RNAi to be plot
    if (type == "RNAi"){
      cond = unique(info_data$Feeding)[c(-1,-2)]
      
    }else if (type == "CTRL"){
      rnai = rnai_list[[condition]]
      cond = c()
      for (r in c("ND7", "ICL7")){
        cond = c(cond, rnai[grep(r, rnai)])
      }
    }
    
    if(any(grepl("bis",cond))){
      cond = cond[-grep("bis", cond)]
    }
    
    pdf(paste0(path,"Gplot_",motif,"_",type,".pdf"))
    for(id in selection){
      expression = EXPRESSION[id,]
      
      tab = c("Time", "Expression", "RNAi")
      for(c in cond){
        rnai_position = str_which(colnames(expression),c)
        tmp = cbind(info_data$Timing[rnai_position],t(expression[rnai_position]), rep(c, length(rnai_position)))
        tab = rbind(tab,tmp)
      }
      
      colnames(tab) = tab[1,]
      tab = as.data.frame(tab[-1,])
      
      tab$Expression = as.numeric(tab$Expression)
      tab$Time = as.numeric(sub("Veg","-1", tab$Time))
      
      tab = tab[order(tab$Time),]
      
      
      gp = ggplot(tab, aes(x = Time, y = Expression))+ 
        geom_point(aes(color = RNAi))+ 
        geom_smooth(alpha = 0.1, color = "black", size = 0.6)+
        theme_minimal() +
        ggtitle(summary_tab[id, "NAME"])
      
      print(gp)
      
    }
    dev.off()
  }
}


##### R status #####
sink(paste0(path,"/Analyse_sessionInfo.txt"))
print(sessionInfo())
sink()
