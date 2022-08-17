options(stringsAsFactors = FALSE)

ies_color = c("red", "orange", "chartreuse3", "dodgerblue", "snow3")

timing = read.table("./DATA/Timing_IES.csv", header = T, sep = ";")

mtFL1 = read.table("./DATA/IES_mtFL1.txt", header = T, sep = "\t")
mtFL1 = merge(mtFL1, timing[,1:2], by = "ID")
mtFL1 = mtFL1[which(mtFL1$SIGNIFICANT == TRUE),]

propotion_timing = cbind(table(timing$GROUP_NAME), table(mtFL1$GROUP_NAME))
colnames(propotion_timing) = c("ALL", "mtFL1")

propotion_timing[,1]= propotion_timing[,1]/sum(propotion_timing[,1])*100 
propotion_timing[,2]= propotion_timing[,2]/sum(propotion_timing[,2])*100 

propotion_timing = propotion_timing[c(5,1:4),]

png("./Timing_IES_mtFL1.png")
par(mar = c(5, 5, 4, 10))
barplot(propotion_timing, 
        main = "IES excsion timing",
        ylab = "% of IES",
        col = ies_color, 
        names.arg = paste(colnames(propotion_timing), c(nrow(timing),nrow(mtFL1)), sep = "\n"),
        legend.text = T,
        args.legend = list(x = "topright",
                           inset = c(- 0.5, 0)))
dev.off()


