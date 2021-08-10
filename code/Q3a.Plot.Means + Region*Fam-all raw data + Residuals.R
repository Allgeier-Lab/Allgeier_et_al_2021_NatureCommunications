rm(list = ls())
library(ape)
library(grDevices)
library(TeachingDemos)
library(RColorBrewer)
library(gplots)
greens <- brewer.pal(7,"Greens")[2:6]
purps <- brewer.pal(7,"Purples")[2:6]
grays <- brewer.pal(7,"Greys")[1:6]
cols <- brewer.pal(length(fam.order),"Set3")


stoic.name <- c("Exc N",
                "Exc P",
                "Exc N:P",
                "Body %N",
                "Body %P",
                "Body %C",
                "Body N:P",
                "Body C:N",
                "Body C:P")
stoic.trait.names <- c( ('nexc.m' ),
                        ('pexc.m' ),
                        ("npexc"),
                        ("bodyn"),
                        ("bodyp"),
                        ("bodyc"), 
                        ("bodynp"),
                        ("bodycn"),
                        ("bodycp") )

#modeled data
raws <- readRDS(file = "data/output/RawMeansPerParamPerFam.PerStoicTrait-notPhyloCorrected-FamNoFam",refhook=NULL)
resid.list <- readRDS(file = "data/output/ResidOutput_Region*Fam*Phylo-NoPhylo-fulldataset",refhook=NULL)
aov.list <- readRDS(file = "data/output/AOVSummaryOutput_Region*Fam*Phylo-NoPhylo-fulldataset",refhook=NULL)
mcmc.list <- readRDS(file = "data/output/MCMCSummaryOutput_Region*Fam*Phylo-NoPhylo-fulldataset",refhook=NULL)

#raw data
g <- read.table("data/modified/Allgeier_FullStoicDataset_woPalmyra",stringsAsFactors=F, header=T) 
g <- g[g$family != "Scorpaenidae",]
g.tab <- data.frame( gen.sp = as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ),
                     FG = (tapply(g$FG,g$gen.sp,unique)),
                     family = tapply(g$family,g$gen.sp,unique),
                     region = tapply(g$region,g$gen.sp, unique) , stringsAsFactors = F)
g <- merge(g.tab,g)

#This is a bit sloppy but I do this to just get the family names in order
famnums <- do.call(rbind, lapply(raws, function(x) {dim(x)[2]}))
max.num.fams <- max(famnums)/2
fam.o <- colnames(raws[[which(famnums == max(famnums))[1]]])
#this is a zinger of a line - creates the family order for the plots
fam.order <- do.call(rbind,strsplit(do.call(rbind,strsplit(fam.o,":"))[,2], "fam"))[c(TRUE,FALSE)]

png(filename="output plots/Figure 3.png",
    units="in",
    width=12.5,
    height=7,
    res=200)

par(mfrow = c(length(stoic.trait.names),1), 
    mgp = c(2.2,.8,0),
    family = "serif",
    cex.lab = 1.4,
    cex.axis = 1.0 ,
    oma = c(4,8,4,2))
layout(matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9, 
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                #10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18,18,
                19,19,20,20,21,21,22,22,23,23,24,24,25,25,26,26,27,27), nrow = 6, ncol = 18, byrow=TRUE) )
cex.azees <- 1.2

#######################
###A) PLOT 1-9 MEANS###
#######################

g <- read.table("data/modified/Allgeier_FullStoicDataset_woPalmyra",stringsAsFactors=F, header=T) 
g.tab <- data.frame( gen.sp = as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ),
                     FG = (tapply(g$FG,g$gen.sp,unique)),
                     family = tapply(g$family,g$gen.sp,unique),
                     region = tapply(g$region,g$gen.sp, unique) , stringsAsFactors = F)
g <- merge(g.tab,g)

xmaxs <- c(200,27,1800,14,9,62,27,8.5,120)

for(k in 1:length(stoic.trait.names)){
  
  par(mar = c(2,1,0,0))
  
  candp <- na.omit(g[,c("region",stoic.trait.names[k])])
  c.now <- candp[candp[,1] == "Caribbean",2]
  p.now <- candp[candp[,1] == "Moorea",2]
  
  Transparency<- '55'
  
  #Density curves for all three habitats
  d1<-density(c.now,adjust=1.5)
  d2<-density(p.now,adjust=1.5)
  
  usethis <- summary(c.now,p.now)[6]
  
  if( max(d1$y)>max(d2$y) ) { jesus <- d1 } else { jesus <-  d2 } 
  
  
  if(k == 1 ) {
    plot( jesus$x, jesus$y ,
          xlim = c(min(c.now,p.now), xmaxs[k] ),
          col = "white",
          ylab = "Density",
          ann=F, cex.lab = 1.7, cex.axis = cex.azees ,
          las = 1)
  } else {
    plot( jesus$x, jesus$y ,
          xlim = c(min(c.now,p.now), xmaxs[k] ),
          col = "white",
          yaxt = "n", ylab = "Density",
          ann=F, cex.lab = 1.7, cex.axis = cex.azees ,
          las = 1)
  }
  
  polygon(x=c(d1$x,0),y=c(d1$y,0),col=paste0(greens[5],Transparency),border = greens[5])
  polygon(x=c(d2$x,0),y=c(d2$y,0),col=paste0(purps[5],Transparency),border = purps[5])
  
  if(k == 1 | k == 5 | k == 7 | k == 9) { legend("topright", "**",text.col="black",ncol=1,bty="n",cex=2,xjust=1 )}
}

mtext(c("Density"), side = 2, col = "black", 
      outer = TRUE, adj = .98, padj = -1.3,line = c(1), cex= 1.2) 

#######################
###B) fam means     ###
#######################

thismuch <- c(0.87,1.1)
add2 <- c(-.1,.1)
Transparency<- '25'
for(stoic in 1:length(stoic.name)){

  par(mar = c(3,1,1,0))
  aov.sign <- aov.list[[stoic]]
  mcmc.sign <- mcmc.list[[stoic]]
  
  #determine which are sign MCMC
  colnames(mcmc.sign) <- c("mean","ci.l","ci.u","samp","p")
  mcmcs <- data.frame(mcmc.sign[c(TRUE,FALSE),],
                      mcmc.sign[c(FALSE,TRUE),])
  
  mcmcs$sig <- ifelse(mcmcs$mean > mcmcs$mean.1 & mcmcs$ci.l > mcmcs$ci.u.1 |
           mcmcs$mean < mcmcs$mean.1 & mcmcs$ci.u < mcmcs$ci.l.1, 
         1, 0 )
  mcmcs$fam <- unlist(strsplit(rownames(mcmcs),"fam"))[c(FALSE,TRUE)]
  
  #determine which are sign AOVMCMC
  colnames(aov.sign) <- c("mean","ci.l","ci.u","samp","p")
  rovs <- data.frame(aov.sign[c(TRUE,FALSE),],
                      aov.sign[c(FALSE,TRUE),])
  
  rovs$sig <- ifelse(rovs$mean > rovs$mean.1 & rovs$ci.l > rovs$ci.u.1 |
                        rovs$mean < rovs$mean.1 & rovs$ci.u < rovs$ci.l.1, 
                      1, 0 )
  rovs$fam <- unlist(strsplit(rownames(rovs),"fam"))[c(FALSE,TRUE)]
  
  #get the raw data out
  print(stoic.trait.names[stoic])
  
  #isolate the fams and get their raw means
  ind.ns <- na.omit(data.frame( gen.sp = g$gen.sp, 
                                param = g[ ,c(stoic.trait.names[stoic]) ],
                                region = g$region,
                                fam = g$family ) )
                                
  curr.fam <- names(which( ifelse( table(ind.ns$fam,ind.ns$region)[,1] > 1 & table(ind.ns$fam,ind.ns$region)[,2] > 2, TRUE, FALSE) == TRUE ))
  order.to.use <-  which( rev(fam.order) %in% rev(curr.fam) )
  
  #identify the families that have data for the trait for each region
  gg <- g[which(g$family %in% curr.fam),]  
  m <- na.omit( data.frame(t(tapply(gg[,c(stoic.trait.names[stoic])],list(gg$region, gg$family),mean,na.rm =T))) )
  sd <- na.omit( data.frame(t(tapply(gg[,c(stoic.trait.names[stoic])],list(gg$region, gg$family),sd,na.rm =T))) )
  
  #Plot one - currently plotting actual means - might go to estimates
  means.c <- m[,1]
  sds.c <- sd[,1]
  means.p <- m[,2]
  sds.p <- sd[,2]
  
  #maxes and mins
  minsies <- min( c( (means.c - sds.c), (means.p - sds.p) ) ) *thismuch[1]
  maxies <- max( c( (means.c + sds.c), (means.p + sds.p) ) ) *thismuch[2]
  
  plot(1:10,1:10, col = "white", pch = 16, 
         ylim = c(0,length(fam.order)),
         xlim = c( minsies, maxies ),
         las = 1, cex = 1.6, yaxt = 'n', xlab = "", 
         ylab = "", cex.lab = 1.4, cex.axis = cex.azees, xaxs = "i",yaxs="i") 
  
  if(stoic == 1) { 
    axis(2,at = 1:length(fam.order)-0.5, labels = rev(fam.order), las = 1, cex.axis = 1.4)
    }
  
  #create sign polygons
  poly.dat <- data.frame(poly1 = seq(0, max.num.fams-1, 1),
                         y = rep(minsies,max.num.fams),
                         poly2 = seq(1, max.num.fams, 1),
                         x = rep(maxies, max.num.fams, 1) )
  
  for(q in 1:length(order.to.use) ){
    if(rev(rovs$sig)[q] == 1) {
      polygon(x = c(minsies,minsies,maxies,maxies) , 
              y = c(order.to.use[q]-1,order.to.use[q],order.to.use[q],order.to.use[q]-1),
              col=paste0(grays[6],Transparency),
              border = paste0(grays[1],Transparency))
    }
    if(mcmcs$sig[q] == 1){
      text( maxies*.9, order.to.use[q]-.25 , "*", cex = 2.5, col = "black")
    }
  }
  
  print(mcmcs$sig)
  
  arrows(  c(means.c+sds.c), rev(order.to.use - 0.5 + add2[1]), c( ifelse( means.c-sds.c > 0, means.c-sds.c, 0 )),rev(order.to.use - 0.5 + add2[1]),
           angle = 1, length = 0, lwd = 2, col = greens[4] )
  
  points( means.c, rev(order.to.use - 0.5 + add2[1]),col =  "black", pch = 16, cex = 1.7)
  points( means.c, rev(order.to.use - 0.5 + add2[1]),col =  greens[4], pch = 16, cex = 1.6)
  
  arrows(  c(means.p+sds.p), rev(order.to.use - 0.5 + add2[2]), c( ifelse( means.p-sds.p > 0, means.p-sds.p, 0 ) ),rev(order.to.use - 0.5 + add2[2]),
           angle = 1, length = 0, lwd = 2, col = purps[4] )
  points( means.p, rev(order.to.use - 0.5 + add2[2]),col = "black", pch = 15, cex = 1.6)
  points( means.p, rev(order.to.use - 0.5 + add2[2]),col = purps[4], pch = 15, cex = 1.5)
  
  arrows(  c(minsies), c(1:length(fam.order)),  c(maxies), c(1:length(fam.order)),
          angle = 1, length = 0, lwd = 1, col = "gray70",lty = 3 )
  
}
 
#######################
###C) resids       ###
#######################

##These are the residuals
par(mar = c(2,1,0,0))
for(stoic in 1:length(stoic.name)){
  #Plot 2 - Residuals
  resy <- resid.list[[stoic]]

  #Density curves for all three habitats
  d1<-density(resy$resids.cond,adjust=1)
  d2<-density(resy$resids.aov,adjust=1)
  
  usethis <- summary(resy$resids.cond,resy$resids.aov)[6]
  
  if( max(d1$y)>max(d2$y) ) { jesus <- d1 } else { jesus <-  d2 } 
  
  if(stoic == 1 ) {
    plot( jesus$x, jesus$y ,
          xlim = c(min(resy$resids.cond,resy$resids.aov), usethis ),
          col = "white",
          ylab = "",
          ann=F, cex.lab = 1.2, cex.axis = 1.2,
          las = 1)
    
  } else {
    plot( jesus$x, jesus$y ,
          xlim = c(min(resy$resids.cond,resy$resids.aov), usethis ),
          col = "white",
          ylab = "",yaxt = "n",
          ann=F, cex.lab = 1.2, cex.axis = 1.2,
          las = 1)
  }
  
  polygon(x=c(d2$x,0),y=c(d2$y,0),col=paste0(grays[6],Transparency),border = "black")
  polygon(x=c(d1$x,0),y=c(d1$y,0),col=paste0(grays[1],Transparency),border = "gray30")
  
  mtext(stoic.name[stoic], side = 1, cex = 1.3, padj = 2.2)
  
}

mtext(c("Density"), side = 2, col = "black", 
      outer = TRUE, adj = .05, padj = -1.3,line = c(1), cex= 1.2) 

mtext( c("Caribbean = "), outer = TRUE, side = 1, col=c( "black" ), cex=1.5, adj = .3, padj = -37)
mtext(c("-"), side = 3, col = c(greens[4]), 
      outer = TRUE, adj = .395, padj = .29,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(greens[4]), 
      outer = TRUE, adj = .395, padj = .38,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(greens[4]), 
      outer = TRUE, adj = .395, padj = .47,line = c(1), cex= 6)

mtext( c("Polynesia = "), outer = TRUE, side = 1, col=c( "black" ), cex=1.5, adj = .675, padj = -37) 
mtext(c("-"), side = 3, col = c(purps[4]), 
      outer = TRUE, adj = .73, padj = .29,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(purps[4]), 
      outer = TRUE, adj = .73, padj = .38,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(purps[4]), 
      outer = TRUE, adj = .73, padj = .47,line = c(1), cex= 6)

mtext( c("A)"), outer = TRUE, side = 2, col= "black", cex=1.2, adj = .6, padj = -23, las = 2)
mtext( c("B)"), outer = TRUE, side = 2, col= "black", cex=1.2, adj = .6, padj = -14.5, las = 2) 
mtext( c("C)"), outer = TRUE, side = 2, col= "black", cex=1.2, adj = .6, padj = 14.5, las = 2) 

dev.off()


