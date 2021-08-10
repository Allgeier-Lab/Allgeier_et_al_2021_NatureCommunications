rm(list = ls())
# Load data ---------
R2 <- read.csv(file = "data/output/AllStoicTraitModels_R2.glmm" , header = T)[,-1]
pl <- read.csv(file = "data/output/AllStoicTraitModels_P.Lambda" , header = T)[,-1]
dic <- read.csv(file = "data/output/AllStoicTraitModels_DIC" , header = T)[,-1]
gs <- readRDS(file = "data/output/AllStoicTraitModels_ParamOutputs",refhook=NULL)


stoic.trait <- list( c('nexc.m', "bodyn" ),
                     c('pexc.m', "bodyp" ),
                     c("npexc", "bodyn","bodyp"),
                     c("bodyn", "bodyp", "bodyc"),
                     c("bodyp", "bodyn", "bodyc"),
                     c("bodyc", "bodyp"),
                     c("bodynp", "bodyc"),
                     c("bodycn","bodyp"),
                     c("bodycp","bodyn") )
stoic.trait2 <- list( c("Body N" ),
                     c("Body P" ),
                     c("Body N","Body P"),
                     c("Body P"),
                     c("Body N","Body C"),
                     c("Body P"),
                     c("Body C"),
                     c("Body P"),
                     c("Body N") )
stoic.name <- c("Exc N",
                "Exc P",
                "Exc N:P",
                "Body %N",
                "Body %P",
                "Body %C",
                "Body N:P",
                "Body C:N",
                "Body C:P")

# PLOT SET UP BS -----------
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = c(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))) )
colspicked = col_vector[c(3,5,6,7,8,9,11,12,13)]

y.names <- c("predator","omnivore","invertivore","herbivore",
             expression(paste(delta^{15},"N")), expression(paste(delta^{13},"C")), 
             "region","mass","Body %P","Body %N", "Body %C" )

y.locations <-c( 1,#pred
                 2,#omni
                 3,#invert
                 4,#herb/det
                 4.85,5.15,#n15
                 5.85,6.15,#c13
                 7,#mass
                 8,#region
                 9,#bodyp,
                 10, #bodyn
                 11)#bodyc) 
pches <- c(18)

y.ticks <- seq(1,length(y.names),1)
x.labs1 <- c("","","","","Prop. Variance Explained","","","","")
x.labs2 <- c("","","","","Standardized Effect Size","","","","")

# ####
png(filename="output plots/Figure 1.png", 
    units="in", 
    width=12.5, 
    height=7, 
    res=200)


par(mfrow = c(length(stoic.trait),2), 
    mgp = c(2.2,.8,0),
    family = "serif",
    cex.lab = 1.4,
    cex.axis = 1.0 ,
    oma = c(4,7,3,2))
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18), nrow = 5, ncol = 18, byrow=TRUE) )

#####################################
#Top Row "Variance Explained" plots##
#####################################

for(j in 1:length(stoic.trait)){
  
  ci.u <- c( pl[,c("L.CIu")][j], R2[,c("R2c.CIu")][j], R2[,c("R2m.CIu")][j] )
  ci.l <- c( pl[,c("L.CIl")][j], R2[,c("R2c.CIl")][j],  R2[,c("R2m.CIl")][j] )
  means <- c( pl[,c("Lambda")][j], R2[,c("R2c")][j],  R2[,c("R2m")][j] )
  ys <- c(0.15,0.5,0.75)
  
  par(mar = c(1,1,1.5,0))
  plot(seq(0,1,1),seq(0,1,1),
       ylim = c(0,1),#max(combo$means + combo$sds)*1.1), 
       xlim = c(0,1), 
       ylab = "" , 
       xlab = "",
       col = "white",
       cex.axis = 1.2,
       yaxt = "n",
       cex.main = 1.7,
       main = stoic.name[j])
  mtext(x.labs1[j],1,cex = 1.3, line = 2.6)
  
  if(j == 1) { axis(2,at = ys,tick = TRUE, labels = c("Phylogeny-only", "Full Model", "Ecology" ), las = 1, cex.axis = 1.3) }
  
  segments(ci.u,ys, ci.l, ys, col = colspicked[j],lwd = 1.5 )
  points(means, ys, pch = pches , col = "gray30", cex = 2)
  points(means, ys, pch = pches , col = colspicked[j], cex = 1.7)
  
  arrows(-10,0.3,150,0.3,lty = 1, col = "black", angle = 0)
  
  expression(paste(delta,"DIC",round( (dic[j,1] - dic[j,2] ),1)))
  dd = round( (dic[j,1] - dic[j,2] ))

  legend(.1, .9,
                 c( expression(paste(Delta,"DIC =")), abs(dd) ),
                 bty="n", ncol = 2, xjust =0 , cex = 1.2, yjust = .5, x.intersp = -.51)
  
}

#####################################
#Bottom Row fixed effects plots######
#####################################

for(n in 1:length(stoic.trait)){
  #n=1
  par(mar = c(3,1,3,0))
  combo <- gs[[n]][,2:length(gs[[n]][1,])]
  combo[is.na(combo)] <- 0
  combo$ord <- length(combo[,1]):1
  combo <- combo[order(combo$ord),]
  
  xmin <- -1
  xmax <- 1
  
  #first - deal with the pchs - the values that are sign are  filled and those not, are not
  cexes <- c(2.4, #pred
             2.4,#omni
             2.4,#invert
             2.4,#herb/det
             2,2,#n15
             2,2,#c13
             2,#region
             2.4, #mass
             2.4, #body something
             2.4, #body something
             2.4 )#bodyotherthing
  
  pches <- c(18, #pred
             18,#omni
             18,#invert
             18,#herb/det
             16,15,#n15
             16,15,#c13
             15,#region
             18, #mass
             18, #body something
             18, #body something
             18 )#bodyotherthing
  pches.e <- c(5, #pred
               5,#omni
               5,#invert
               5,#herb/det
               1,0,#n15
               1,0,#c13
               0,#region
               5, #mass
               5, #body something
               5, #body something
               5 )#bodyotherthing
  
  if ( length(stoic.trait[[n]]) == 2 ) {
    runner <- paste(stoic.trait[[n]][1],stoic.trait[[n]][2] , sep = "_" )
    } else if ( length(stoic.trait[[n]]) == 3 ) {
    runner <- paste(stoic.trait[[n]][1],stoic.trait[[n]][2],stoic.trait[[n]][3] , sep = "_" )
    }
  
  plot(seq(xmin,xmax,(xmax + abs(xmin))/ (length(y.names))),
       seq(1,length(y.names)+1,1),
       ylim = c(0.75,length(y.names)),
       xlim = c(xmin,xmax), 
       ylab = "" , 
       xlab = "",
       col = "white",
       yaxt = "n",
       cex.axis = 1.2)
    
  if(n == 1) { axis(2,at = y.ticks,tick = TRUE, labels = c(y.names), las = 1, cex.axis = 1.6) }
  cols.now <- ifelse(combo$median == 0, "white", colspicked[n])
  cols.now1 <- ifelse(combo$upper>0 & combo$lower<0, "white",cols.now)
  cols.nowg <- ifelse(combo$median == 0, "white", "gray30")
  segments(combo$upper,y.locations, combo$lower, y.locations, col = cols.now,lwd = 1.5 )
  points(combo$median, y.locations, pch = pches , col = cols.nowg, cex = cexes+.4)
  points(combo$median, y.locations, pch = pches , col = cols.now1, cex = cexes)
  arrows(0,0,0,150,lty = 3, col = "gray50", angle = 0)
 
  if(n == 1) { legend("topleft", c("Caribbean", "Polynesia"), pch = c(16,15),
                      ncol = 1,bty = "n", cex = 1.4) }
  mtext(x.labs2[n],1,cex = 1.3, line = 2.8)
}
mtext("A)",2,cex = 1.3, line = 2.8, las = 2, adj = 46, padj = -24.3)
mtext("B)",2,cex = 1.3, line = 2.8, las = 2, adj = 48.5, padj = -15)

dev.off()