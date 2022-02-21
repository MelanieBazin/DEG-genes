
Khi2_intermed <- function(genes_list, profil_list){
  nb_genes = length(genes_list[which(is.element(genes_list, profil_list))])
  other = length(genes) - nb_genes
  # Définition des probabilitée théoriques
  proba = length(profil_list)/length(AUTOGAMY$all_genes) 
  proba = c(proba, 1-proba)
  chi2 = chisq.test(c(nb_genes, other), p  = proba)
  
  pv=chi2$p.value
  signif="ns"
  if(pv < 1e-200) {
    signif="****"
  } else {
    if(pv < 1e-100) {
      signif="***"
    } else {   
      
      if(pv < 1e-20) {
        signif="**"
      } else {
        if(pv < 1e-10) {
          signif="*"
        }
      }
    }
  }
  
  return(paste(format(pv,digits=3), signif))
}

Enrichment_padj <- function(LIST, data_tab, nb_simulation = 1000){

  data_tab = as.data.frame(data_tab)
  colnames(data_tab) = c("ID", "PROFIL")
  
  print("Table of all data")
  print(table(data_tab$PROFIL))

  Profils = unique(data_tab$PROFIL)
  all = length(data_tab$ID)
  for (l in names(LIST)){
    ## Calcul des p-value
    pval = c()
    for (p in Profils){
      nb_profil = sum(data_tab$PROFIL == p)
      nb_genes = length(LIST[[l]])
      nb_genes_profil = length(LIST[l][which(is.element(LIST[[l]], data_tab$ID[data_tab$PROFIL == p]))])
      
      
      pval = c(pval, 1- phyper(nb_genes_profil-1, nb_profil, all - nb_profil, nb_genes))
    }
    names(pval) = Profils
    
    ##  Correction des p-value
    nb_simulation = 1000
    
    # Création de données séléction aléatoirement
    theoric_ID = matrix(NA, nrow = length(LIST[[l]]), ncol = nb_simulation)
    for (c in 1:ncol(theoric_ID)){
      theoric_ID[,c] = sample(data_tab$ID, length(LIST[[l]]))
    }
    theoric_pval = matrix(NA, nrow = length(Profils), ncol = nb_simulation)
    row.names(theoric_pval) = Profils
    
    for (s in 1:nb_simulation){
      genes = theoric_ID[,s]
      pvals = c()
      for (p in Profils){
        nb_profil = sum(data_tab$PROFIL == p)
        nb_genes = length(genes)
        nb_genes_profil = length(genes[which(is.element(genes, data_tab$ID[data_tab$PROFIL == p]))])
        pvals = c(pvals, 1- phyper(nb_genes_profil-1, nb_profil, all - nb_profil, nb_genes))
      }
      theoric_pval[,s] = pvals
    }
    
    # Ajuster les pvalue
    pval_adj = c()
    for (p in Profils){
      th_pvals = theoric_pval[p,]
      pval_adj = c(pval_adj, sum(th_pvals< pval[p])/length(th_pvals))
    }
    
    names(pval_adj) = Profils
    
    print(l)
    print(table(data_tab$PROFIL[which(is.element(data_tab$ID, LIST[[l]]))]))
    print(cbind(as.data.frame(pval_adj),pval_adj < 0.05))
    
  }
}
