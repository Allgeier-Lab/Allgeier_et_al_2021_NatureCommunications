rm(list = ls())


g <- read.table("data/modified/Allgeier_FullStoicDataset_woPalmyra",stringsAsFactors=F, header=T) #these are fish only data 
g.tab <- data.frame( gen.sp = as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ),
                     count = tapply(g$FG,g$gen.sp,length),
                     FG = (tapply(g$FG,g$gen.sp,unique)),
                     family = tapply(g$family,g$gen.sp,unique),#[ as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ) ])),
                     region = tapply(g$region,g$gen.sp, unique) , stringsAsFactors = F)
ggs <- g.tab[order(g.tab$region),]

write.csv(data.frame(ggs), file = "data/output/MainTable.csv")
