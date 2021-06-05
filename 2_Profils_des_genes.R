source("2_Visualisation_des_donnees_fonction.R", encoding = "UTF-8")


ExpressionProfils("RPKM")


for (i in names(rnai_list)){
  rnai = rnai_list[[i]]
  ExpressionProfils("RPKM", rnai)
}



