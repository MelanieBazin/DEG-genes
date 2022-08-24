####
# Definition of analyse variable
####
options(stringsAsFactors = FALSE)
library(stringr)

# Open the table with gene annotation
annotation = read.table("./DATA/My_annotation.tab",header=T,sep="\t")

#### Definition of the time course that will be analysed together ####
tout = sub(".tab","",list.files("./DATA/EXPRESSION/", pattern = ".tab"))

rnai_list = list(
  HiSeqvsNextSeq = tout[which(is.element(tout,c("EZL1bis","ICL7bis","EZL1", "ICL7" )))],
  analyseDE = tout[which(!is.element(tout,c("EZL1bis","EZL1" )))]
)
rm(tout)

#### Pseudo-replicate groups definition ####
cluster = list(
  ICL7 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  ICL7bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
  EZL1 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  EZL1bis = c(rep("EARLY",1),rep("INTER",1),rep("LATE",2)),
  
  ND7_K = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  PGM = c(rep("VEG",1),rep("INTER",3),rep("LATE",3)),
  KU80c = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  
  ND7_L = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",3)),
  CTIP = c(rep("VEG",1),rep("EARLY",1),rep("INTER",2),rep("LATE",1)),
  ND7_X = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",3)),
  XRCC4 = c(rep("VEG",1),rep("EARLY",1),rep("INTER",1),rep("LATE",3))
)


#### Definition of the different color set ####
# "darkorange1"
cluster_color = c("VEG" = "mediumpurple4",
                  "EARLY" = "deepskyblue",
                  "INTER" = "chartreuse3",
                  "LATE" = "red")

method_color = c("HiSeq" = "chartreuse4",
                 "NextSeq" = "blue4",
                 "both" = "mediumturquoise")

profile_color = c("Early peak" = "purple3",
                  "Intermediate peak" = "red2",
                  "Late peak" = "chartreuse4",
                  "Early repression" = "dodgerblue3",
                  "Late induction" = "deeppink",
                  "Late repression" = "darkorange",
                  "none" = "snow3")


#### Definition autogamy timing list for all time course #####
tabs = list.files("./DATA/EXPRESSION", pattern = ".tab")
timing_list = as.list(tabs)
for (k in 1:length(tabs)){
  table = read.table(paste0("./DATA/EXPRESSION/",tabs[k]), header = T, row.names = 1)
  names(timing_list)[k] = gsub(".tab", "", tabs[k])
  timing_list[[k]] = colnames(table)
}
rm(table, tabs,k)
