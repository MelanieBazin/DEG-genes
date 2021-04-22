data_directories = c("CTIP", "CTIP_CTRL" , "ICL7", "KU80c" , "ND7", "PGM", "XRCC4", "XRCC4_CTRL" )

CTIP = c("ID", "T0", "T5.5", "T12.5", "T25", "Veg")
CTIP_CTRL = c("ID", "T0", "T5", "T10", "T20", "T30", "Veg")
ICL7 = c("ID", "T0", "T5", "T10", "T20", "T35", "T50", "Veg")
KU80c = c("ID", "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
ND7 = c("ID", "T0", "T5", "T10", "T20", "T30", "T40", "Veg")
PGM = c("ID", "T2", "T5", "T10", "T20", "T30", "T40", "Veg")
XRCC4 = c("ID", "T2", "T7", "T22", "T32","Veg")
XRCC4_CTRL = c("ID", "T2", "T7", "T22", "T32","Veg")
name = list(CTIP, CTIP_CTRL , ICL7, KU80c , ND7, PGM, XRCC4, XRCC4_CTRL )

dir.create("1_EXPRESSION")

for (i in 1:length(data_directories)){
  list= list.files(data_directories[i])
  table = NULL
  for(j in list){
    tab = read.table(paste0(data_directories[i],"/",j))
    table = cbind(table, tab[,2])
  }
  table = cbind(as.character(tab[,1]), table)
  colnames(table) = name[[i]]
  write.table(table,paste0("1_EXPRESSION/",data_directories[i],".tab"),sep="\t", row.names=F,quote=F)
}

