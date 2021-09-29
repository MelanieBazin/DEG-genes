# Effectif observé
UP = c(564,
       516,
       716,
       652,
       157,
       910,
       4427
)

# Effectif théorique
ALL = c(1974,
        4536,
        2037,
        3741,
        468,
        4434,
        24343)
pALL = ALL/sum(ALL)


tab = matrix(c(ALL,UP),2,length(ALL), byrow = T)
colnames(tab) = c("Early peak",
                  "Early repression",
                  "Intermediate peak",
                  "Late induction",
                  "Late peak",
                  "Late repression",
                  "none")
rownames(tab) = c("ALL","UP")
for( i in 1:ncol(tab)){
  chi2 = chisq.test(tab[,i])
  print(colnames(tab)[i])
  print(chi2$p.value)
}

for( i in 1:length(pALL)){
  chi2 = chisq.test(UP[i], p = pALL[i] )
  print(chi2$p.value)
}

# H0 : La proportion de gènes IES+ et - est la même que parmis tous les gènes
# Plus p-value est petite, plus la théorie et l'observation diffèrent

# if(pv < 1e-200) {
#   signif="****"
# } else {
#   if(pv < 1e-100) {
#     signif="***"
#   } else {   
#     
#     if(pv < 1e-20) {
#       signif="**"
#     } else {
#       if(pv < 1e-10) {
#         signif="*"
#       }