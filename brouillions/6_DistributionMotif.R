options(stringsAsFactors = FALSE)
library(seqinr)

colors=c("purple", "orange", "red3","forestgreen", "royalblue","deeppink", "grey")
colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")


path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"
path_motif = paste0(path, "MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/")

#### Initialisation des filtres ####
source("6_Filtres.R")


### Ouverture des liste de gènes avec motifs ####
type = "FIMO2"
dossier  = "FIMO2/"
# for (type in c("fuzznuc", "fimo", "FIMO2")){
#   # type = "fimo"
#   if(type == "fimo"){
#     dossier  = "FIMO/"
#   }else if (type == "fuzznuc"){
#     dossier  = "FUZZNUC/"
#   }else if (type == "FIMO2"){
#     dossier  = "FIMO2/"
#   }
  
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
  
# }


#### Faire la réciproque : chercher un enrichissement dans certain groupe parmis les gènes avec motif ####
type = "FIMO2"
dossier  = "FIMO2/"

files = list.files(paste0(path_motif,dossier), pattern = ".gff")

save_path = paste0(path_motif,type,"_Distribution_Motif/")
dir.create(save_path, recursive=T,showWarnings=F)

for (f in files){
  # f = files[1]
  print(f)
  
  # Création du tableau de donnée pour stocker les effectifs pour chaque filtres
  effectif_tab = data.frame(as.numeric(summary(ANTIFILTRES)[1:length(ANTIFILTRES)]), row.names = rownames(summary(ANTIFILTRES)))
  colnames(effectif_tab) = "Effectifs"
  effectif_tab[,"Sans_Motif"]= NA
  effectif_tab[,"Avec_Motif"]= NA
  effectif_tab[,"Motif_-"]= NA
  effectif_tab[,"Motif_+"]= NA
  
  effectif_tab_Noduplicat = effectif_tab
  effectif_tab_Noduplicat[,"Motif_+/-"] = 0
  
  prom_with_motif = read.table(paste0(path_motif,dossier,f), header=F, sep="\t")
  colnames(prom_with_motif)=colnamesgff3
  
  for(filtre in names(ANTIFILTRES)){
    # filtre = names(ANTIFILTRES)[3]
    
    selection = ANTIFILTRES[[filtre]]
    filter_motif = prom_with_motif[is.element(prom_with_motif$ID, selection),]
    
    #### Calcul des effectifs pour chaque condition ####
    effectif = effectif_tab[filtre, "Effectifs"]
    motifs_nb = table(filter_motif$STRAND)
    if (is.na(motifs_nb["+"])){
      motifs_nb["+"] = 0
    }
    if (is.na(motifs_nb["-"])){
      motifs_nb["-"] = 0
    }
    effectif_tab[filtre,2:ncol(effectif_tab)] = c(effectif-sum(motifs_nb),sum(motifs_nb), motifs_nb["-"], motifs_nb["+"])
    
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
      if (is.na(motifs_nb["+"])){
        motifs_nb["+"] = 0
      }
      if (is.na(motifs_nb["-"])){
        motifs_nb["-"] = 0
      }
      
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
    }else{
      motif_uniq_tot = 0
      motif_both = 0
      motifs_nb = c(0,0)
      effectif_tab_Noduplicat[filtre,2:ncol(effectif_tab_Noduplicat)] = c(effectif-motif_uniq_tot,sum(motifs_nb), motifs_nb,motif_both)
    }
  }
  write.table(effectif_tab, paste0(save_path,sub(".gff","",sub(".gff3","",f)),"_Effectifs_avec_duplicat_DOWN.tab"), sep = "\t")

  
  ### Effectif avec les gènes sans duplicat ###
  if (sum(effectif_tab_Noduplicat$`Motif_+/-`)>0){
    write.table(effectif_tab_Noduplicat, paste0(save_path,sub(".gff","",sub(".gff3","",f)),"_sans_duplicat_DOWN.tab"), sep = "\t")

    
  }
 
}
