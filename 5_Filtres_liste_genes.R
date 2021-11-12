
annotation = read.table("./DATA/My_annotation2.tab",header=T,sep="\t")

# Liste des gènes d'interets
ALL = read.table(paste0("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/Resumer_DEgenes.tab"), sep = "\t", h=T )
ALL = ALL[-41534,]

PGM_UP = na.omit(ALL$ID[ALL$PGM_LATE_REGULATION == "Up-regulated"])
KU80c_UP = na.omit(ALL$ID[ALL$KU80c_LATE_REGULATION == "Up-regulated"])
XRCC4_UP = na.omit(ALL$ID[ALL$XRCC4_LATE_REGULATION == "Up-regulated"])
EZL1_UP = na.omit(ALL$ID[ALL$EZL1_LATE_REGULATION == "Up-regulated"])

UP = Reduce(intersect, list(PGM_UP,KU80c_UP, XRCC4_UP, EZL1_UP))

PGM_DOWN = na.omit(ALL$ID[ALL$PGM_LATE_REGULATION == "Down-regulated"])
KU80c_DOWN = na.omit(ALL$ID[ALL$KU80c_LATE_REGULATION == "Down-regulated"])
XRCC4_DOWN = na.omit(ALL$ID[ALL$XRCC4_LATE_REGULATION == "Down-regulated"])
EZL1_DOWN = na.omit(ALL$ID[ALL$EZL1_LATE_REGULATION == "Down-regulated"])

DOWN = Reduce(intersect, list(PGM_DOWN,KU80c_DOWN, XRCC4_DOWN, EZL1_DOWN))

CTIP_e_DOWN = na.omit(ALL$ID[ALL$CTIP_EARLY_REGULATION == "Down-regulated"])
CTIP_i_DOWN = na.omit(ALL$ID[ALL$CTIP_INTER_REGULATION == "Down-regulated"])

CTIP_DOWN = intersect(CTIP_e_DOWN, CTIP_i_DOWN)

CTIP_e_UP = na.omit(ALL$ID[ALL$CTIP_EARLY_REGULATION == "Down-regulated"])
CTIP_i_UP = na.omit(ALL$ID[ALL$CTIP_INTER_REGULATION == "Down-regulated"])

CTIP_UP = intersect(CTIP_e_UP, CTIP_i_UP)

FILTRE_ARNi = list(
  PGM_UP = PGM_UP,
  KU80c_UP = KU80c_UP,
  XRCC4_UP = XRCC4_UP,
  EZL1_UP = EZL1_UP,
  UP = UP,
  CTIP_DOWN = CTIP_DOWN,
  UP_DOWN = intersect(UP, CTIP_DOWN),
  PGM_DOWN = PGM_DOWN,
  KU80c_DOWN = KU80c_DOWN,
  XRCC4_DOWN = XRCC4_DOWN,
  EZL1_DOWN = EZL1_DOWN,
  DOWN = DOWN,
  CTIP_UP = CTIP_UP,
  DOWN_UP = intersect(DOWN, CTIP_UP)
)

MOTIF = read.table("Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/tout/MOTIF/Motif_dans_prom/Motif_intermed/Parmis_tous_les_genes/FIMO2/FIMO2_IN_MAC_CDS.tsv",sep = "\t", h=T)
MOTIF$start =  MOTIF$start-150
MOTIF$stop = MOTIF$stop-150

ALL2 = ALL
colnames(ALL2)[1] = "sequence_name"
FILTRE_MOTIF = list(
  MOTIF = MOTIF,
  MOTIF_pos= MOTIF[MOTIF$start < -50 & MOTIF$start > -80,],
  MOTIF_plus = MOTIF[MOTIF$strand == "+",],
  MOTIF_moins = MOTIF[MOTIF$strand == "-",],
  MOTIF_pos_plus = MOTIF[MOTIF$strand == "+" & MOTIF$start < -50 & MOTIF$start > -80,],
  MOTIF_pos_moins = MOTIF[MOTIF$strand == "-" & MOTIF$start < -50 & MOTIF$start > -80,]
  
)

rm(ALL2,PGM_UP, KU80c_UP, XRCC4_UP, EZL1_UP, UP, 
   CTIP_DOWN, PGM_DOWN, KU80c_DOWN, XRCC4_DOWN, EZL1_DOWN, DOWN, 
   CTIP_UP, CTIP_e_DOWN, CTIP_e_UP, CTIP_i_DOWN, CTIP_i_UP,
   MOTIF)

#### Fonctions pour analyse enrichissement ####
pVal_hypergeometric <- function(data_base, genes_list){
  tab_base = table(data_base$EXPRESSION_PROFIL)
  
  tab = data_base[genes_list,]
  tab = table(tab$EXPRESSION_PROFIL)
  
  # Regarder l'enrichissement pour chacun des groupes de gènes
  N = nrow(data_base)
  n = length(genes_list)
  
  p_vec = c()
  for (t in names(tab)){
    # Nombre de gènes porteur de la fonction dans la base de donnée et dans la liste
    N1 = tab_base[[t]]
    k = tab[[t]]
    
    # Calcul de la probabilité d'observer cette valeur k (si H0 est vraie)
    pval = 1 - phyper(k - 1, N1, N - N1, n)
    p_vec = c(p_vec, pval)
    
  }
  names(p_vec) = names(tab)
  
  return(p_vec)
}


Enrichissement <- function(data_base, genes_list, nb_sim_corr ){
  rownames(data_base) = data_base$ID
  # Verifier que tout les gènes d'interet sont dans la base
  if(length(setdiff(genes_list, data_base$ID))>0){
    return(print("Interest genes that are not in the database are removed"))
    genes_list = intersect(genes_list, data_base$ID)
  }
  
  p_vec = pVal_hypergeometric(data_base, genes_list)
  
  ### Corection des p-value par simulation
  # Tirages de gènes au hasard
  sim_genes_tab = matrix(NA, nrow = length(genes_list), ncol = nb_sim_corr)
  for(i in 1:ncol(sim_genes_tab)){
    sim_genes_tab[,i] = sample(row.names(data_base), length(genes_list))
  }
  
  # Simulation des p-values
  groups = unique(data_base$EXPRESSION_PROFIL)
  sim_pVal = matrix(NA, nrow = length(groups), ncol = nb_sim_corr)
  row.names(sim_pVal) = groups
  
  for (s in 1:nb_sim_corr){
    p_vec_sim = pVal_hypergeometric(data_base,  sim_genes_tab[,s])
    
    sim_pVal[,s] = p_vec_sim
  }
  
  # Correctioon des p-value
  pVal_adj = c()
  
  for(g in groups){
    
    p_vec_sim = sim_pVal[g,]
    
    padj2 = sum(p_vec_sim < p_vec[g])/length(p_vec_sim)    
    # --> Proportion de pval hasard < pvalue réelle  
    
    pVal_adj = c(pVal_adj, padj2)
  }
  
  names(pVal_adj) = groups
  
  return(pVal_adj)
}