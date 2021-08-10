
rm(list = ls())
library(ape)
library(grDevices)
library(TeachingDemos)
library(siar)
library(RColorBrewer)
greens <- brewer.pal(7,"Greens")[2:6]
purps <- brewer.pal(7,"Purples")[2:6]
grays <- brewer.pal(7,"Greys")[2:6]

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

################################
#bring in data and organize ####
################################
fulldata <- read.csv("data/output/All_BK.Lambda.AllSampleSizes_Boot=1000.csv",stringsAsFactors=F, header=T) 
fulldata$stoic.name <- rep(stoic.name,each = 100)
cdata <- read.csv("data/output/Caribbean_BK.Lambda.AllSampleSizes_Boot=1000.csv",stringsAsFactors=F, header=T) 
cdata$stoic.name <- rep(stoic.name,each = 100)
mdata <- read.csv("data/output/Moorea_BK.Lambda.AllSampleSizes_Boot=1000.csv",stringsAsFactors=F, header=T) 
mdata$stoic.name <- rep(stoic.name,each = 100)


fulldata$stoic.num <- rep(seq(1:9), each=100)
cdata$stoic.num <- rep(seq(1:9), each=100)
mdata$stoic.num <- rep(seq(1:9), each=100)

fulldata$lp.test <- ifelse(fulldata$lambda.p<0.05,1, 0)
cdata$lp.test <- ifelse(cdata$lambda.p<0.05,1, 0)
mdata$lp.test <- ifelse(mdata$lambda.p<0.05,1, 0)

fulldata$bk.test <- ifelse(fulldata$bk.p<0.05,1, 0)
cdata$bk.test <- ifelse(cdata$bk.p<0.05,1, 0)
mdata$bk.test <- ifelse(mdata$bk.p<0.05,1, 0)

#compile for BK
sig.bk <- rbind(tapply(fulldata$bk.test,fulldata$stoic.name,sum),
                tapply(cdata$bk.test,cdata$stoic.name,sum),
                tapply(mdata$bk.test,mdata$stoic.name,sum) )

rownames(sig.bk) <- c("full","cdata","mdata")

apply(sig.bk,1,sum)

alldata.bk = list(fulldata[which(fulldata$bk.test == 1), ] ,
                   cdata[ which( cdata$bk.test == 1), ] ,
                   mdata[ which( mdata$bk.test == 1), ] )

#this checks to make sure the data are lining up through all the filtering
sum(apply(sig.bk,1,sum)) == dim(do.call(rbind.data.frame,alldata.bk))[1]

#compile for Lambda
sig.lp <- rbind(tapply(fulldata$lp.test,fulldata$stoic.name,sum),
                tapply(cdata$lp.test,cdata$stoic.name,sum),
                tapply(mdata$lp.test,mdata$stoic.name,sum) )

rownames(sig.lp) <- c("full","cdata","mdata")

apply(sig.lp,1,sum)

alldata.lp = list(fulldata[which(fulldata$lp.test  == 1), ] ,
                  cdata[ which( cdata$lp.test  == 1), ] ,
                  mdata[ which( mdata$lp.test  == 1), ] )

#this checks to make sure the data are lining up through all the filtering
sum(apply(sig.lp,1,sum)) == dim(do.call(rbind.data.frame,alldata.lp))[1]
################################

#ANOVAs OF THE P. LAMBDAs####
#subsetted data 
aov.bk.mat <- matrix(ncol = 4, nrow = 9, dimnames = list(stoic.name,c("dataset","full-C","P-C","P-full" )) )
aov.cor.mat <- matrix(ncol = 4, nrow = 9, dimnames = list(stoic.name,c("dataset","full-C","P-C","P-full"  )) )

for(h in 1:9){
  #data = alldatalist[[i]]
  now.bk <- rbind( alldata.bk[[1]][alldata.bk[[1]]$stoic.num == h,],
                alldata.bk[[2]][alldata.bk[[2]]$stoic.num == h,],
                alldata.bk[[3]][alldata.bk[[3]]$stoic.num == h,] )
  now.bk$X <- as.factor (c(rep(c("full"), dim(alldata.bk[[1]][alldata.bk[[1]]$stoic.num == h,])[1] ), 
                       rep(c("C"), dim(alldata.bk[[2]][alldata.bk[[2]]$stoic.num == h,])[1] ), 
                       rep(c("P"), dim(alldata.bk[[3]][alldata.bk[[3]]$stoic.num == h,])[1] )) )
  
  now.lp <- rbind( alldata.lp[[1]][alldata.lp[[1]]$stoic.num == h,],
                   alldata.lp[[2]][alldata.lp[[2]]$stoic.num == h,],
                   alldata.lp[[3]][alldata.lp[[3]]$stoic.num == h,] )
  now.lp$X <- as.factor (c(rep(c("full"), dim(alldata.lp[[1]][alldata.lp[[1]]$stoic.num == h,])[1] ), 
                           rep(c("C"), dim(alldata.lp[[2]][alldata.lp[[2]]$stoic.num == h,])[1] ), 
                           rep(c("P"), dim(alldata.lp[[3]][alldata.lp[[3]]$stoic.num == h,])[1] )) )
  
  
  bks <- aov(bk~X, data = now.bk)
  bktk <- TukeyHSD(bks)
  aov.bk.mat[h,] <-  ifelse( c(summary(bks)[[1]][1,5], bktk[[1]][,4] ) <0.05, "**","NS")
  
  ls <- aov(lambda~X, data = now.lp)
  ltk <- TukeyHSD(ls)
  aov.cor.mat[h,] <-  ifelse( c(summary(ls)[[1]][1,5], ltk[[1]][,4] ) <0.05, "**","NS")
  
}


bks <- do.call(rbind.data.frame,alldata.bk)

bks$X <- as.factor (c(rep(c("full"), dim(alldata.bk[[1]])[1] ), 
                         rep(c("C"), dim(alldata.bk[[2]])[1] ), 
                         rep(c("P"), dim(alldata.bk[[3]])[1] ) ) ) 
#for BK
c.f <- list()
tingsmon <- c("full","C","P")
for(i in 1:3){
  
  fullk <- aov(bk~stoic.name, data = bks[bks$X == tingsmon[i],])
  summary(fullk)
  c.full <- data.frame(TukeyHSD(fullk)$stoic.name)
  c.full$sig <- ifelse(c.full[,4]<0.05,1,0)
  c.f[[i]] <- c.full[c.full$sig == 1,]
  
}

#for P.L.
c.f <- list()
tingsmon <- c("full","C","P")
for(i in 1:3){
  
  fullk <- aov(lambda~stoic.name, data = bks[bks$X == tingsmon[i],])
  summary(fullk)
  c.full <- data.frame(TukeyHSD(fullk)$stoic.name)
  c.full$sig <- ifelse(c.full[,4]<0.05,1,0)
  c.f[[i]] <- c.full[c.full$sig == 1,]
  
}

# END OF ANOVAS ####

#FINAL PLOTS
######################
###PLOT 1-9 MEANS#####
######################

g <- read.table("data/modified/Allgeier_FullStoicDataset_woPalmyra",stringsAsFactors=F, header=T) #these are fish only data 
g.tab <- data.frame( gen.sp = as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ),
                     FG = (tapply(g$FG,g$gen.sp,unique)),
                     family = tapply(g$family,g$gen.sp,unique),#[ as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ) ])),
                     region = tapply(g$region,g$gen.sp, unique) , stringsAsFactors = F)
g <- merge(g.tab,g)

#quartz(width = 12, height = 8)
png(filename="output plots/Figure 2.png", 
    units="in", 
    width=12, 
    height=8,
    #pointsize=12, 
    res=200)
par(mfrow = c(2,2), mar = c(1,5,0,2), oma = c(6,6,2,3), mgp = c(3,1,0),
    family = "serif")
layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,
                2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4),
                nrow = 6, ncol = 36, byrow=TRUE) )

xmaxs <- c(200,27,1800,14,9,62,27,8.5,120)


######################
###PLOT A - Pagels#####
######################

gs <- readRDS(file = "data/output/AllStoicTraitModels_ParamOutputs",refhook=NULL)

coloreslist <- list(grays, greens, purps)
alldatalist <- list(fulldata, cdata, mdata)
add <- c(0,0.3,0.6)


param.data <- list(alldata.bk, alldata.lp)
param.names <- c( "K-statistic",expression(paste("Pagel's ", lambda, sep = " ")) )
curr.param <- c("bk", "lambda")
sig <- list(sig.bk,sig.lp)
bk.or.pl <- c("bk","lambda")

mm <- 2
curr.data <- param.data[[mm]]
dataall <- do.call(rbind.data.frame,curr.data)
dmeans <- tapply( dataall[,c(curr.param[mm])], dataall$stoic.num, mean )
dsds <- tapply( dataall[,c(curr.param[mm])], dataall$stoic.num, sd)
dses <- dsds/sqrt(tapply(dataall[,c( curr.param[mm]) ], dataall$stoic.num, length))

par(mar = c(1,0,3,0))

thismuch <- c(0.92, 1.08)

plot(1:10,1:10, col = "white", pch = 16, ylim = c(min(dmeans-dses)*thismuch[1],max(dmeans+ dses)*thismuch[2]) ,
     xlim = c(1,10),las = 1, cex = 1.7, xaxt = 'n', xlab = "", ylab = param.names[mm], cex.lab = 2.2, cex.axis = 1.5,
     xaxs = "i")
axis(1,at = 1:9+0.5,labels = stoic.name, cex.axis = 1.9, las = 1)

mtext(param.names[mm], side = 2, col = "black", 
      outer = TRUE, adj = .6, padj = -1.5,line = c(1), cex= 1.4) 

save.the.means <- matrix(ncol = length(stoic.name), nrow = 3)
add2 <- c(0.25,0.5,0.75)
for(i in 1:3){
  
  data = curr.data[[i]]
  colores = coloreslist[[i]]
  
  means <- tapply(data[,c( curr.param[mm]) ], data$stoic.num, mean)
  sds <- tapply(data[,c( curr.param[mm]) ], data$stoic.num, sd)
  ses <- sds/sqrt(tapply(data[,c( curr.param[mm]) ], data$stoic.num, length))
  
  arrows( c(1:9+add2[i]), c(means+ses), c(1:9+add2[i]), c(means-ses),
          angle = 1, length = 0, lwd = 6, col = colores[5])
  arrows( c(1:9+add2[i]), c(means+ses), c(1:9+add2[i]), c(means-ses),
          angle = 1, length = 0, lwd = 4, col = colores[3])
  
  points(1:9+add2[i], means, col =  colores[5], pch = 16, cex = 2.1)
  points(1:9+add2[i], means, col = colores[3], pch = 16, cex = 1.8)
  
  save.the.means[i,] <- means 
  
}

arrows(c(1:8+1), c(rep(-10,8)), c(1:8+1), c(rep(100,8)),col = "black" )

mtext(c("Full Dataset = "), side = 3, col = c("black"), 
      outer = TRUE, adj = .2, padj = 2.7,line = c(1), cex= 1.5)
mtext(c("-"), side = 3, col = c(grays[4]), 
      outer = TRUE, adj = .32, padj = .81,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(grays[4]), 
      outer = TRUE, adj = .32, padj = .9,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(grays[4]), 
      outer = TRUE, adj = .32, padj = .99,line = c(1), cex= 6)

mtext(c("Caribbean = "), side = 3, col = c("black"), 
      outer = TRUE, adj = .5, padj = 2.7,line = c(1), cex= 1.5)
mtext(c("-"), side = 3, col = c(greens[4]), 
      outer = TRUE, adj = .58, padj = .81,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(greens[4]), 
      outer = TRUE, adj = .58, padj = .9,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(greens[4]), 
      outer = TRUE, adj = .58, padj = .99,line = c(1), cex= 6)

mtext(c("Polynesia = "), side = 3, col = c("black"), 
      outer = TRUE, adj = .77, padj = 2.7,line = c(1), cex= 1.5)
mtext(c("-"), side = 3, col = c(purps[4]), 
      outer = TRUE, adj = .83, padj = .81,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(purps[4]), 
      outer = TRUE, adj = .83, padj = .9,line = c(1), cex= 6)
mtext(c("-"), side = 3, col = c(purps[4]), 
      outer = TRUE, adj = .83, padj = .99,line = c(1), cex= 6)

########################################################  
############Plot B i & ii C13 and N15 ############################# 
######################################################## 

stand <- function(X) { (X-mean(X,na.rm=T))/(sd(X,na.rm=T)) } #standardize variables
thisun <- g[which( !is.na(g[,c("n15")] & g[,c("c13")]) ), ]
thisun$c13s <- stand(thisun$c13)
thisun$n15s <- stand(thisun$n15)

Transparency<- '55'

for(j in 1:2){ 
  
  if(j == 1) { par(mar = c(2,0,3,4) ) } else { par( mar = c(2,3,3,2) ) }
  
  if(j == 1){
    
    c.now <- thisun[thisun$region == "Caribbean",]$c13s
    p.now <- thisun[thisun$region == "Moorea",]$c13s
    
  } else {
    
    c.now <- thisun[thisun$region == "Caribbean",]$n15s#( na.omit(raw[raw$region == "Caribbean",c(stoic.trait.names[k])]) )
    p.now <- thisun[thisun$region == "Moorea",]$n15s#( na.omit(raw[raw$region == "Moorea",c(stoic.trait.names[k])]) )
    
  }
  
  d1<-density(c.now,adjust=1)
  d2<-density(p.now,adjust=1)
  
  if( max(d1$y)>max(d2$y) ) { jesus <- d1 } else { jesus <-  d2 } 
  
  if(j == 1){
    
    plot(jesus$x,jesus$y, col = "white",
         xlab=expression(delta^13~C),
         ylab="",cex.lab=1.8, 
         cex.axis=1.5,asp=0,las=1,xlim=c( min(c.now,p.now)*.9, max(c.now,p.now)*1.1 ) )
  } else {
    
    plot(jesus$x,jesus$y, col = "white",
         xlab=expression(delta^15~N),
         ylab="",cex.lab=1.8, 
         cex.axis=1.5,asp=0,las=1,xlim=c( min(c.now,p.now)*.9, max(c.now,p.now)*1.1 ) )
    
  }
  
  #fill polygons from density distributions with transparent colors
  polygon(x=c(d1$x,0),y=c(d1$y,0),col=paste0(greens[5],Transparency),border = greens[5])
  polygon(x=c(d2$x,0),y=c(d2$y,0),col=paste0(purps[5],Transparency),border = purps[5])
  
  if(j == 1) {  legend(-2,.55,legend=c(ifelse(ks.test(c.now, p.now)$p.value<0.05, "ks = **","t = NS")) ,text.col="black",ncol=1,bty="n",cex=2,xjust=1)
  } else {  legend(4.5,.66,legend=c(ifelse(ks.test(c.now, p.now)$p.value<0.05, "ks = **","t = NS")) ,text.col="black",ncol=1,bty="n",cex=2,xjust=1) }
  
  
}

mtext(expression(delta^13~C), side = 1, col = "black", 
      outer = TRUE, adj = .12, padj = .1,line = c(1), cex= 1.4) 

mtext(expression(delta^15~N), side = 1, col = "black", 
      outer = TRUE, adj = .51, padj = .1,line = c(1), cex= 1.4) 

mtext("Density", side = 2, col = "black", 
      outer = TRUE, adj = .12, padj = -2.5,line = c(1), cex= 1.4) 

mtext("Density", side = 2, col = "black", 
      outer = TRUE, adj = .12, padj = 23.5,line = c(1), cex= 1.4) 

#########################  
#Plot B iii SEAs############################# 
######################### 

par(mar = c(2,4.5,3,0))
boot <- 500

seacx <- list()
seapx <- list()

seacy <- list()
seapy <- list()

seas <- matrix(ncol = 4, nrow = boot, dimnames = list(NULL, c("C", "P","n15","c13") ) )

howmuch <- .5*min(table(thisun$region))

#for SEA plot
plot(NULL,NULL,col=c(1,2,3),type="p",
     xlab="",
     ylab="",
     cex.lab=2,
     cex.axis=1.5,asp=0,las=1,xlim=c(-1.2,1.2),ylim=c(-1.6,1.6))

mtext(expression(delta^13~C), side = 1, col = "black", 
      outer = TRUE, adj = .88, padj = .3,line = c(1), cex= 1.4) 

mtext(expression(delta^15~N), side = 2, col = "black", 
      outer = TRUE, adj = .14, padj = 35,line = c(1), cex= 1.4) 

for(its in 1:boot){
  
  print(its)
  
  cs <- thisun[thisun$region == "Caribbean",]
  c <- cs[sample(1:length(cs[,1]),howmuch,replace = FALSE),]
  
  ps <- thisun[thisun$region == "Moorea",]
  p <- ps[sample(1:length(cs[,1]),howmuch,replace = FALSE),]
  
  SEc <- standard.ellipse(c$c13s,c$n15s,steps=1)
  SEp <- standard.ellipse(p$c13s,p$n15s,steps=1)
  
  #for standard elipse
  seacx[[its]] <- SEc$xSEA
  seapx[[its]] <- SEp$xSEA
  seacy[[its]] <- SEc$ySEA
  seapy[[its]] <- SEp$ySEA
  
  #For SEA
  seas[its,] <- c(SEc$SEA, SEp$SEA,
                  ifelse(ks.test(c$c13s, p$c13s)$p.value<0.05, 1,0),
                  ifelse(ks.test(c$n15s, p$n15s)$p.value<0.05, 1,0) ) 
  lines(SEc$xSEAc,SEc$ySEAc,col=greens[1],lty=1,lwd=.5)  
  lines(SEp$xSEAc,SEp$ySEAc,col=purps[1],lty=1,lwd=.5) 
  
  #For TA
  #seas[its,] <- c(SEc$TA, SEp$TA)
  #lines(SEc$xcoords,SEc$ycoords,col=greens[1],lty=1,lwd=.5)  
  #lines(SEp$xcoords,SEp$ycoords,col=purps[1],lty=1,lwd=.5) 
  
  #  lines(SE$xCEA,SE$yCEA,col=colors[j],lty=2,lwd=3)  
  #points(median(SEc$xSEAc),median(SEc$ySEAc),col=greens[1],pch=16,cex = 1)
  #points(median(SEp$xSEAc),median(SEp$ySEAc),col=purps[1],pch=16,cex = 1)
  
}

sum(ifelse( (seas[,2] - seas[,1])>0, 0, 1 ) )
legend("topright",legend=ifelse(t.test( as.numeric(seas[,1]), as.numeric(seas[,2]) )$p.value<0.05, "**","NS"),text.col="black",ncol=1,bty="n",cex=2,xjust=1)

#add in the output from all KS test for the isotopes
mtext(paste( round(sum(as.numeric(seas[,4]))/500, 2)*100, "%" ), side = 1, col = "black", 
      outer = TRUE, adj = .59, padj = -10.5,line = c(1), cex= 1.4) 

mtext(paste( round(sum(as.numeric(seas[,3]))/500, 2)*100, "%" ), side = 1, col = "black", 
      outer = TRUE, adj = .025, padj = -10.5,line = c(1), cex= 1.4) 


#########################  
#Panel Labels######
mtext("A)", side = 2, col = "black", 
      outer = TRUE, adj = 1, padj = -21,line = c(3), cex= 1.4, las = 1) 
mtext("B)", side = 2, col = "black", 
      outer = TRUE, adj = 1, padj = 9.3,line = c(3), cex= 1.4, las = 1) 


mtext("i)", side = 2, col = "black", 
      outer = TRUE, adj = -1.8, padj = 11.6,line = c(3), cex= 1.2, las = 1) 
mtext("ii)", side = 2, col = "black", 
      outer = TRUE, adj = -23, padj = 11.6,line = c(3), cex= 1.2, las = 1) 
mtext("iii)", side = 2, col = "black", 
      outer = TRUE, adj = -33.5, padj = 11.6,line = c(3), cex= 1.2, las = 1) 

dev.off()

    ######################### 
    ######################### 

# Tests for ~75%, and 90% of the number of individuals in the Caribbean ####

seacx <- list()
seapx <- list()

seacy <- list()
seapy <- list()

seas <- matrix(ncol = 4, nrow = boot, dimnames = list(NULL, c("C", "P","n15","c13") ) )

howmuch <- c( .75*min(table(thisun$region)), .9*min(table(thisun$region)) )

quartz(width = 7, height = 3)
par(mfrow = c(1,2))

for (gg in 1:2){
  
  #for SEA plot
  plot(NULL,NULL,col=c(1,2,3),type="p",
       #xlab="",
       #ylab="",
       cex.lab=2,
       cex.axis=1.5,asp=0,las=1,xlim=c(-1.2,1.2),ylim=c(-1.6,1.6))
  
  mtext(expression(delta^13~C), side = 1, col = "black", 
        outer = TRUE, adj = .88, padj = .3,line = c(1), cex= 1.4) 
  
  mtext(expression(delta^15~N), side = 2, col = "black", 
        outer = TRUE, adj = .14, padj = 35,line = c(1), cex= 1.4) 
  
  for(its in 1:boot){
    #its =1 
    print(its)
    
    cs <- thisun[thisun$region == "Caribbean",]
    c <- cs[sample(1:length(cs[,1]),howmuch[gg],replace = FALSE),]
    
    ps <- thisun[thisun$region == "Moorea",]
    p <- ps[sample(1:length(cs[,1]),howmuch[gg],replace = FALSE),]
    
    SEc <- standard.ellipse(c$c13s,c$n15s,steps=1)
    SEp <- standard.ellipse(p$c13s,p$n15s,steps=1)
    
    #for standard elipse
    seacx[[its]] <- SEc$xSEA
    seapx[[its]] <- SEp$xSEA
    seacy[[its]] <- SEc$ySEA
    seapy[[its]] <- SEp$ySEA
    
    #For SEA
    seas[its,] <- c(SEc$SEA, SEp$SEA,
                    ifelse(ks.test(c$c13s, p$c13s)$p.value<0.05, 1,0),
                    ifelse(ks.test(c$n15s, p$n15s)$p.value<0.05, 1,0) ) 
    
    lines(SEc$xSEAc,SEc$ySEAc,col=greens[1],lty=1,lwd=.5)  
    lines(SEp$xSEAc,SEp$ySEAc,col=purps[1],lty=1,lwd=.5) 
    
  }
  
  sum(ifelse( (seas[,2] - seas[,1])>0, 0, 1 ) )
  legend("center",legend=ifelse(t.test( as.numeric(seas[,1]), as.numeric(seas[,2]) )$p.value<0.05, "**","NS"),text.col="black",ncol=1,bty="n",cex=2,xjust=1)
  
  #add in the output from all KS test for the isotopes
  legend("bottomleft", paste("c13 = ", round(sum(as.numeric(seas[,4]))/500, 2)*100, "%" ), cex= 1.4, bty = "n") 
  legend("topleft", paste("n15 = ", round(sum(as.numeric(seas[,4]))/500, 2)*100, "%" ), cex= 1.4, bty = "n") 
  
  
}

