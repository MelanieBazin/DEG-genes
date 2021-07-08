source("0_Cluster.R")
source("3_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")


ExpressionProfils("RPKM", select_ID = select_ID)


for (i in names(rnai_list)){
  rnai = rnai_list[[i]]
  ExpressionProfils("RPKM", rnai, select_ID = select_ID)
}



ExpressionProfils(type = "DESeq2", 
                  condition = "tout", 
                  file = "./Analyse/2021-07-07_Analyse_DESeq2_tout_CombatON_FC-1.5_pval-0.05/",
                  select_ID = select_ID)

