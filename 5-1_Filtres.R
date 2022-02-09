all_genes = annotation$ID
inter_genes = annotation$ID[annotation$EXPRESSION_PROFIL == "Intermediate peak"]
early_genes = annotation$ID[annotation$EXPRESSION_PROFIL == "Early peak"]
auto_genes = annotation$ID[annotation$EXPRESSION_PROFIL != "none"]

# Les gènes UP dérégulés en PGM, KU80c et/OU XRCC4
UP_PKX = list()
for (R in RNAi[-1]){
  tab = TAB[[R]]
  tab = tab[tab$REGULATION == "Up-regulated", "ID"]
  UP_PKX = c(UP_PKX, list(tab))
}
names(UP_PKX) = paste("UP", RNAi[-1], sep = "_")

# Les gènes NOT UP dérégulés en PGM, KU80c et/OU XRCC4
not_UP_PKX = list()
for (l in names(UP_PKX)){
  up_genes = UP_PKX[[l]]
  not_UP_PKX = c(not_UP_PKX, list(setdiff(all_genes, up_genes)))
}
names(not_UP_PKX) = paste("not",names(UP_PKX), sep = "_")

# Les gènes DOWN dérégulés en PGM, KU80c et/OU XRCC4
DOWN_PKX = list()
for (R in RNAi[-1]){
  tab = TAB[[R]]
  tab = tab[tab$REGULATION == "Down-regulated", "ID"]
  DOWN_PKX = c(DOWN_PKX, list(tab))
}
names(DOWN_PKX) = paste("DOWN", RNAi[-1], sep = "_")

# Les gènes NOT dérégulés en PGM, KU80c et/OU XRCC4
not_DE_PKX = list()
for (l in names(UP_PKX)){
  up_genes = UP_PKX[[l]]
  down_genes = DOWN_PKX[[l]]
  not_DE_PKX = c(not_DE_PKX, list(setdiff(all_genes, unique(down_genes,up_genes))))
}
names(not_DE_PKX) = paste("not_DE",RNAi[-1], sep = "_")

# Les gènes DOWN dérégulés en CTIP early ou inter
CTIP = rbind(TAB[["CTIP_early"]],TAB[["CTIP_inter"]])

DOWN_C = unique(CTIP[CTIP$REGULATION == "Down-regulated", "ID"])

# Les gènes UP dérégulés en CTIP
UP_C = unique(CTIP[CTIP$REGULATION == "Up-regulated", "ID"])

# Les gènes NOT DE dérégulés en CTIP
not_DE_C = setdiff(all_genes, unique(DOWN_C,UP_C))



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