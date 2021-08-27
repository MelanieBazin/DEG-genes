options(stringsAsFactors = FALSE)
library(ggvenn)
library(stringr) 


colors=c("purple", "orange", "red3","forestgreen", "royalblue","deeppink", "grey")
colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")


path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"
path_motif = paste0(path, "MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/")

#### Initialisation des filtres ####
infogenes = read.table(paste0(path,"Resumer_DEgenes.tab"), header = T)


Tous_les_genes = unique(na.omit(infogenes$ID))
Genes_codant = unique(na.omit(infogenes$ID[infogenes$TYPE == "Coding Gene"]))
UP_PGM = unique(na.omit(infogenes$ID[infogenes$PGM_LATE_REGULATION == "Up-regulated"]))
UP_KU80c = unique(na.omit(infogenes$ID[infogenes$KU80c_LATE_REGULATION == "Up-regulated"]))
UP_XRCC4 = unique(na.omit(infogenes$ID[infogenes$XRCC4_LATE_REGULATION == "Up-regulated"]))
UP_EZL1 = unique(na.omit(infogenes$ID[infogenes$EZL1_LATE_REGULATION == "Up-regulated"]))
DOWN_CTIP = unique(na.omit(infogenes$ID[infogenes$CTIP_EARLY_REGULATION == "Down-regulated" | infogenes$CTIP_INTER_REGULATION == "Down-regulated"]))
UP_PGM_KU80c = unique(intersect(UP_PGM, UP_KU80c))
UP_PKX = unique(intersect(UP_PGM_KU80c, UP_XRCC4))
UP_PKXE= unique(intersect(UP_PKX, UP_EZL1))
UP_PKX_DOWN_CTIP = unique(intersect(UP_PKX, DOWN_CTIP))
UP_DOWN = unique(intersect(UP_PKXE, DOWN_CTIP))
DE_autogamie = unique(na.omit(infogenes$ID[infogenes$EXPRESSION_PROFIL != "none"]))
Pic_intermediaire = unique(na.omit(infogenes$ID[infogenes$EXPRESSION_PROFIL == "Intermediate peak"]))
UP_PKX_inter = unique(intersect(Pic_intermediaire,UP_PKX))
UP_PKXE_inter = unique(intersect(Pic_intermediaire,UP_PKXE))
UP_PKX_DOWN_inter = unique(intersect(Pic_intermediaire,UP_PKX_DOWN_CTIP))
UP_DOWN_inter = unique(intersect(Pic_intermediaire,UP_DOWN))

FILTRES = list(Tous_les_genes,Genes_codant,
               UP_PGM, UP_KU80c, UP_XRCC4, UP_EZL1, DOWN_CTIP,
               UP_PGM_KU80c, UP_PKX, UP_PKXE, 
               UP_PKX_DOWN_CTIP, UP_DOWN ,
               DE_autogamie, Pic_intermediaire,
               UP_PKX_inter, UP_PKXE_inter,
               UP_PKX_DOWN_inter, UP_DOWN_inter)

names(FILTRES) = c("Tous_les_genes","Genes_codant",
                   "UP_PGM", "UP_KU80c", "UP_XRCC4", "UP_EZL1", "DOWN_CTIP",
                   "UP_PGM_KU80c", "UP_PKX", "UP_PKXE", 
                   "UP_PKX_DOWN_CTIP", "UP_DOWN" ,
                   "DE_autogamie", "Pic_intermediaire",
                   "UP_PKX_inter", "UP_PKXE_inter",
                   "UP_PKX_DOWN_inter", "UP_DOWN_inter")

#### Faire diagramme de Venne pour les gènes selectionnée comme étant la population de gènes positif pour STREM ####
selection = list(UP_PGM, UP_KU80c, UP_XRCC4, Pic_intermediaire)
names(selection)=c("UP_PGM", "UP_KU80c", "UP_XRCC4", "Pic_intermediaire")
png(paste0(path_motif,"Venn_selection.png"))
ggvenn(selection,
       fill_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"),
       stroke_size = 0.5,
       set_name_size = 6,
       show_percentage = F,
       text_size = 6,
       set_name_color = c("#0073C2FF", "darkorange", "#868686FF", "darkolivegreen3"))
dev.off()



### Ouverture des liste de gènes avec motifs ####
for (type in c("fuzznuc", "fimo", "FIMO2")){
  # type = "fimo"
  if(type == "fimo"){
    dossier  = "FIMO/"
  }else if (type == "fuzznuc"){
    dossier  = "FUZZNUC/"
  }else if (type == "FIMO2"){
    dossier  = "FIMO2/"
  }
  
  files = list.files(paste0(path_motif,dossier), pattern = ".gff")
  
  save_path = paste0(path_motif,type,"_Distribution_Motif/")
  dir.create(save_path, recursive=T,showWarnings=F)
  
  for (f in files){
    # f = files[1]
    print(f)
    
    # Création du tableau de donnée pour stocker les effectifs pour chaque filtres
    effectif_tab = data.frame(as.numeric(summary(FILTRES)[1:length(FILTRES)]), row.names = rownames(summary(FILTRES)))
    colnames(effectif_tab) = "Effectifs"
    effectif_tab[,"Sans_Motif"]= NA
    effectif_tab[,"Avec_Motif"]= NA
    effectif_tab[,"Motif_-"]= NA
    effectif_tab[,"Motif_+"]= NA
    
    effectif_tab_Noduplicat = effectif_tab
    effectif_tab_Noduplicat[,"Motif_+/-"] = 0
    
    prom_with_motif = read.table(paste0(path_motif,dossier,f), header=F, sep="\t")
    colnames(prom_with_motif)=colnamesgff3
    
    pdf(paste0(save_path,sub(".gff","",sub(".gff3","",f)),"_graphs.pdf"))
    par(mfrow=c(2,2), cex.axis = 0.75)
    for(filtre in names(FILTRES)){
      selection = FILTRES[[filtre]]
      filter_motif = prom_with_motif[is.element(prom_with_motif$ID, selection),]
      
      #### Distribution du start du motif dans les promoteur ####
      start_motif_plus = filter_motif$START[filter_motif$STRAND == "+"]-150
      hist(start_motif_plus, 
           xlab = "Position début du motif", 
           main = paste("Motif +","\n",filtre),
           sub = paste("Effectif =",length(start_motif_plus)),
           breaks = 30,
           xlim = c(-150,0),
           xaxp = c(-150,0,15))
      
      start_motif_moins = filter_motif$START[filter_motif$STRAND == "-"]-150
      hist(start_motif_moins, 
           xlab = "Position début du motif", 
           main = paste("Motif -","\n",filtre),
           sub = paste("Effectif =",length(start_motif_moins)),
           breaks = 30,
           xlim = c(-150,0),
           xaxp = c(-150,0,15))
      
      #### Calcul des effectifs pour chaque condition ####
      effectif = effectif_tab[filtre, "Effectifs"]
      motifs_nb = table(filter_motif$STRAND)
      effectif_tab[filtre,2:ncol(effectif_tab)] = c(effectif-sum(motifs_nb),sum(motifs_nb), motifs_nb )
      
      ### Faire une liste de genes uniques ###
      # Les gènes avec plus d'un motif sont réparti de la facon suivent :
      # - les gènes avec uniquemnt les motifs plus sont dans "Motif_+"
      # - les gènes avec uniquemnt les motifs moins sont dans "Motif_-"
      # - les gènes avec au moins une forme plus et une forme moins sont dans "Motif_+/-"
      multi = unique(filter_motif$ID[duplicated(filter_motif$ID)])
      if (length(multi)>0){
        ## Retirer tous les genes avec plus d'un motif ##
        filter_motif_NOduplicat = filter_motif[which(!is.element(filter_motif$ID,multi)),]
        motif_uniq_tot = length(unique(filter_motif$ID))
        motifs_nb = table(filter_motif_NOduplicat$STRAND)
        
        ## Compter les genes ayant plus d'un motif ##
        motif_both = 0
        for (id in multi){
          
          motif_tab_temp = filter_motif[grep(id, filter_motif$ID),]
          motif_strands = names(table(motif_tab_temp$STRAND))
          if(length(motif_strands)==1){
            if (motif_strands == "-"){
              motifs_nb["-"] = motifs_nb["-"]+1
            }else if (motif_strands == "+"){
              motifs_nb["+"] = motifs_nb["+"]+1
            }}else{
              motif_both = motif_both+1
            }
        }
        
        ## Editer le tableau des effectifs ##
        if(motif_uniq_tot == sum(motifs_nb, motif_both)){
          effectif_tab_Noduplicat[filtre,2:ncol(effectif_tab_Noduplicat)] = c(effectif-motif_uniq_tot,sum(motifs_nb), motifs_nb,motif_both)
        }
      }
    }
    write.table(effectif_tab, paste0(save_path,sub(".gff","",sub(".gff3","",f)),"_Effectifs_avec_duplicat.tab"), sep = "\t")
    
    
    par(mfrow=c(1,2), cex.axis = 0.5, cex.main = 0.75, cex.lab = 0.75)
    ### Effectif avec les gnes porteur plus d'un motif ###
    barplot(t(as.matrix(effectif_tab[,c("Sans_Motif","Motif_-", "Motif_+")])),
            las = 2,
            ylab = "Nombre de motifs",
            main = paste0("Representation des motifs","\n",f),
            legend.text =T)
    barplot(t(as.matrix(effectif_tab[,c("Sans_Motif","Motif_-", "Motif_+")]/effectif_tab$Effectifs*100)),
            las = 2,
            ylab = "% de motifs",
            main = paste0("Representation des motifs","\n",f))
    
    
    ### Effectif avec les gènes sans duplicat ###
    if (sum(effectif_tab_Noduplicat$`Motif_+/-`)>0){
      write.table(effectif_tab_Noduplicat, paste0(save_path,sub(".gff","",sub(".gff3","",f)),"_sans_duplicat.tab"), sep = "\t")
      ## Comparaison avec le nombre de gènes sans motif de la cathégorie ##
      barplot(t(as.matrix(effectif_tab_Noduplicat[,c("Sans_Motif","Motif_-", "Motif_+","Motif_+/-")])),
              las = 2,
              ylab = "Nombre de promoteur",
              main = paste0("Genes avec motif","\n",f),
              legend.text =T)
      barplot(t(as.matrix(effectif_tab_Noduplicat[,c("Sans_Motif","Motif_-", "Motif_+","Motif_+/-")]/effectif_tab_Noduplicat$Effectifs*100)),
              las = 2,
              ylab = "% de promoteur",
              main = paste0("Genes avec un motif","\n",f))
      
    }
    dev.off()
  }
  
}


#### Faire la réciproque : chercher un enrichissement dans certain groupe parmis les gènes avec motif ####
