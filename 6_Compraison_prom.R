options(stringsAsFactors = FALSE)
library(seqinr)
library(stringr) 

# f = files[1]
# type = "fimo" #fuzznuc ou fimo

path = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/"
save_path = paste0(path, "MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/")

annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")

colnamesgff3=c("ID","SOURCE","TYPE","START","END","SCORE","STRAND","PHASE", "ATTRIBUTES")

UP = read.table(paste0(path,"Resumer_DEgenes_selection_UP.tab"), header = T)

for (type in c("fuzznuc", "fimo")){
  
  if(type == "fimo"){
    dossier  = "FIMO/"
    
  }else if (type == "fuzznuc"){
    dossier  = "FUZZNUC/"
  }
  
  files = list.files(paste0(save_path,dossier), pattern = ".gff")
  
  for (f in files){
    print(f)
    
    prom_with_motif = read.table(paste0(save_path,dossier,f), header=F, sep="\t")
    colnames(prom_with_motif)=colnamesgff3
    
    tab_annot = merge(prom_with_motif, annotation, by = "ID")
    tab_annot$START = tab_annot$START-150
    tab_annot$END = tab_annot$END-150
    
    if(type == "fimo"){
      prom_with_motif2 = read.table(paste0(save_path,dossier,sub(".gff","",f),".tsv"), header=T, sep="\t")
      tab_annot2 = merge(prom_with_motif2, annotation, by.x = "sequence_name", by.y ="ID")
      tab_annot2$start = tab_annot2$start-150
      tab_annot2$stop = tab_annot2$stop-150
      write.table(tab_annot2, paste0(save_path, sub(".gff","",f),"_motif.tab"), sep = "\t", row.names = F) 
      
      tab_annot = tab_annot[order(tab_annot2$p.value),]
      
    }else if (type == "fuzznuc"){
      write.table(tab_annot, paste0(save_path, sub(".gff3","",f),"_motif.tab"), sep = "\t", row.names = F) 
    }
    
    debut = str_sub(sub(".gff","",sub(".gff3","",f)), -3)
    if(debut == "CDS"){
      debut = paste0("CDS_start")
    }
    if (grepl("IES",f)){
      IES = "IES_"
    }else{IES = ""}
    
    promoteur = read.fasta(paste0("./DATA/Promoteur/IN_MAC_", IES, "upstream_150nt_",debut,".fa"))
    
    
    for (strand in c("+","-")){
      tab_annot_select = tab_annot[tab_annot$STRAND == strand,]
      for(j in c("local", "loin","proche", "all", "UP", "UP_int")){
        if (j == "all"){
          tab_annot_select_plus = tab_annot_select
        }else if (j == "local"){
          tab_annot_select_plus = tab_annot_select[tab_annot_select$START >= -70 & tab_annot_select$START<= -30,]
        }else if(j == "loin"){
          tab_annot_select_plus = tab_annot_select[tab_annot_select$START < -70,]
        }else if(j == "proche"){
          tab_annot_select_plus = tab_annot_select[tab_annot_select$START > -30,]
        }else if (j =="UP"){
          tab_annot_select_plus = tab_annot_select[which(is.element(tab_annot_select$ID, UP$ID)),]
        }else if(j == "UP_int"){
          tab_annot_select_plus = tab_annot_select[which(is.element(tab_annot_select$ID, UP$ID)),]
          tab_annot_select_plus = tab_annot_select_plus[tab_annot_select_plus$EXPRESSION_PROFIL == "Intermediate peak",]
        }
        if (nrow(tab_annot_select_plus) > 5){
          if (nrow(tab_annot_select_plus)>=100){
            tab_annot_select_plus = tab_annot_select_plus[1:100,]
          }

          prom_select = promoteur[which(is.element(names(promoteur),tab_annot_select_plus$ID))]
          for (i in na.omit(names(prom_select))){
            repet = max(tab_annot_select_plus$START)-max(tab_annot_select_plus$START[grep(i,tab_annot_select_plus$ID)])
            prom_select[[i]] = c(rep("-",repet ),prom_select[[i]])
          }

          name = c()
          for (n in names(prom_select)){
            name = c(name, tab_annot$NAME[grep(n, tab_annot$ID)])
          }
          
          write.fasta(sequences = prom_select, names = name, file.out = paste0(save_path, sub(".gff","",sub(".gff3","",f)),"_prom",strand,"_",j,".fa") )
        }
      }
    }
  }
}
print ("- FINI ! -")

###### Notes #####
#get_sequence_from_gff.pl -gff ~/DATA/ptetraurelia_mac_51_annotation_v2.0.gff3 -genome /data/PARAMECIUM/REFERENCES/tetraurelia/ptetraurelia_mac_51.fa -upstream 150 -type CDS -out IN_MAC_upstream_150nt_CDS_start.tab > IN_MAC_upstream_150nt_CDS_start.fa
# promotors=read.table("SEARCH MOTIF/IN_MAC_upstream_150nt_CDS_start.tab",header=F, sep="\t")
# colnames(promotors)=colnamesgff3[1:4]

# Motif = AAAATC[AC]TT[AGT][AT]A[AT]T[AT][AT]TT
# exéctué dans : cd PROJECTS/DEG-genes/DATA/Promoteur/

# fuzznuc -sequence IN_MAC_upstream_150nt_CDS_start.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]T[AT][AT]TT -outfile fuzznuc_IN_MAC_CDS.gff3 -complement -rformat2 gff
# fuzznuc -sequence IN_MAC_upstream_150nt_TSS.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]T[AT][AT]TT -outfile fuzznuc_IN_MAC_TSS.gff3 -complement -rformat2 gff
# fuzznuc -sequence IN_MAC_IES_upstream_150nt_CDS_start.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]T[AT][AT]TT -outfile fuzznuc_IN_MAC_IES_CDS.gff3 -complement -rformat2 gff
# fuzznuc -sequence IN_MAC_IES_upstream_150nt_TSS.fa -pattern AAAATC[AC]TT[AGT][AT]A[AT]T[AT][AT]TT -outfile fuzznuc_IN_MAC_IES_TSS.gff3 -complement -rformat2 gff
