#This tests how similar the results are from the non-phlyo corrected and phylo correct analysis

rm(list = ls())
#library(data.table)
#library(gplots)
#library(picante)
#library(dplyr)
library(fishtree)
library(operators)
library("ape")
library("geiger")
library(phytools)
library(MCMCglmm)
library(MuMIn)

#Second Test
#On average are the low P Caribbean fish more closely related to the low P Poly fish?
###############################################################
# Load functions####
###############################################################
stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) } #standardize variables
center <- function(X) {(X-mean(X,na.rm=T))}

#fucntion to quantify strength of random effects (phylogeny)
p.lambda <- function(mod){
  (r.effects <- ( mod$VCV[,1] )/ ( mod$VCV[,1] + mod$VCV[,2] )) # could use rowSums(m1$VCV))
  ( r.prop.mode <<- posterior.mode( r.effects ) )
  ( r.prop <<- median( r.effects ) )
  ( r.prop.mean <<- mean( r.effects ) )
  ( r.propCI <<- HPDinterval(r.effects, 0.95) )
  print(lambda <<- matrix(dimnames = list(NULL, c("Lambda","L.CIl","L.CIu")),
                          c(r.prop[[1]],c(r.propCI)), ncol = 3, nrow = 1,byrow = T) ) 
  
}

#mmf is a MCMCGlmm model object
#mmF = m1
#prints [1]R2m [2]R2m.cil [3]R2m.ciu [4]R2c [5]R2c.cil [6]R2c.ciu
#gnerates variance explained + CI for three things:
#1) FE vs all varienace + error = R2m
#2) FE + RE vs all variance + error = R2c
#3) RE vs all variance + error = R2RE
R.squaredGLMM <- function(mmF){ 
  
  if ( length(mmF$VCV[1,]) == 2 ) {
    runner <- mmF$VCV[,1]+mmF$VCV[,2]
  } else if ( length(mmF$VCV[1,]) == 3 ) {
    runner <- mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3] 
  } else if ( length(mmF$VCV[1,]) == 4 ) {
    runner <- mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3]+mmF$VCV[,4] }
  
  if ( length(mmF$VCV[1,]) == 2 ) {
    numerator <- mmF$VCV[,1]
  } else if ( length(mmF$VCV[1,]) == 3 ) {
    numerator <- mmF$VCV[,1]+mmF$VCV[,2]
  } else if ( length(mmF$VCV[1,]) == 4 ) {
    numerator <- mmF$VCV[,1]+mmF$VCV[,2]+mmF$VCV[,3] }
  
  mVarF <- var(as.vector(apply(mmF$Sol,2,mean) %*% t(mmF$X)))
  # MCMCglmm - marginal
  #print(mVarF/(mVarF+sum(apply(mmF$VCV,2,mean))))
  
  # with crebile intervals
  vmVarF<-numeric(dim(mmF$Sol)[1])
  for(i in 1:dim(mmF$Sol)[1]){
    Var<-var(as.vector(mmF$Sol[i,] %*% t(mmF$X)))
    vmVarF[i]<-Var}
  R2m1<-vmVarF/(vmVarF+runner)
  #mean(R2m)
  ( R2m.mode <<- round(posterior.mode(R2m1)[[1]],3) )
  ( R2m <<- round(median(R2m1),3) )
  ( R2m.CIl <<- round(HPDinterval(R2m1)[1],3) )
  ( R2m.CIu <<- round(HPDinterval(R2m1)[2],3) )
  # R2GLMM(c) - conditional R2GLMM for full model
  # MCMCglmm - conditional
  #(mVarF+sum(apply(mmF$VCV,2,mean)[-3]))/(mVarF+sum(apply(mmF$VCV,2,mean)))
  # Varieance explained by RE + FE - alternative with crebile intervals
  R2c1<-(vmVarF+numerator)/(vmVarF+runner)
  #mean(R2c1)
  ( R2c.mode <<- round(posterior.mode(R2c1)[[1]],3) )
  ( R2c <<- round(median(R2c1),3) )
  ( R2c.CIl <<- round(HPDinterval(R2c1)[1],3) )
  ( R2c.CIu <<- round(HPDinterval(R2c1)[2],3) )
  R2RE1 <- (numerator)/(vmVarF+runner)
  ( R2RE.mode <<- round(posterior.mode(R2RE1)[[1]],3) )
  ( R2RE <<- round(median(R2RE1),3) )
  ( R2RE.CIl <<- round(HPDinterval(R2RE1)[1],3) )
  ( R2RE.CIu <<- round(HPDinterval(R2RE1)[2],3) )
  print(allz <<- matrix(dimnames = list(NULL, c("R2m","R2m.CIl","R2m.CIu","R2c","R2c.CIl","R2c.CIu","R2RE","R2RE.CIl","R2RE.CIu")),
                        c(R2m,R2m.CIl,R2m.CIu,R2c,R2c.CIl,R2c.CIu,R2RE,R2RE.CIl,R2RE.CIu), ncol = 9, nrow = 1,byrow = T) )
}


###############################################################
# Load data and set up g.tab####
###############################################################
#Read the tree 
sp.tree.full <- read.newick("data/raw/actinopt_12k_treePL.tre") 

#g <- read.table("Allgeier_FullStoicDataset",stringsAsFactors=F, header=T)
g <- read.table("data/modified/Allgeier_FullStoicDataset_woPalmyra",stringsAsFactors=F, header=T) #%>% #these are fish only data 
g <- g[g$family != "Scorpaenidae",]#remove lionfish because they are invasive
g.tab <- data.frame( gen.sp = as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ),
                     FG = (tapply(g$FG,g$gen.sp,unique)),
                     family = tapply(g$family,g$gen.sp,unique),#[ as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ) ])),
                     region = tapply(g$region,g$gen.sp, unique) , stringsAsFactors = F)
g <- merge(g.tab,g)

#I need to standardize or center the n15 values and c13 values
g$keyID <- paste(g$gen.sp,g$nexc.m,g$mass,g$bodycn,g$pexc.m,g$npexc,g$bodyp,g$bodyc,sep = "_")
length(unique(g$keyID)) == dim(g)[1]

s.c <- g[g$region == "Caribbean",][,c("keyID","n15","c13")]
s.m <- g[g$region == "Moorea",][,c("keyID","n15","c13")]

s.c$c.n15 <- stand(s.c$n15)
s.c$c.c13 <- stand(s.c$c13)
s.m$c.n15 <- stand(s.m$n15)
s.m$c.c13 <- stand(s.m$c13)

g <- merge(g,data.frame(rbind(s.c,s.m)) )

stoic.name <- c("Exc N",
                "Exc P",
                "Exc NP",
                "Body %N",
                "Body %P",
                "Body %C",
                "Body NP",
                "Body CN",
                "Body CP")
stoic.trait.names <- c( ('nexc.m' ),
                        ('pexc.m' ),
                        ("npexc"),
                        ("bodyn"),
                        ("bodyp"),
                        ("bodyc"), 
                        ("bodynp"),
                        ("bodycn"),
                        ("bodycp") )

standardize <- c("no")
centerize <- c("no")
MCMCglmm.updateable<- updateable(MCMCglmm)

# ####
namethemthangs <- c("AICc.null", "AICc.n15", "AICc.c13","AICc.nofam",
                    "AICc.n15nofam", "AICc.c13nofam","AICc.fam", "AICc.famregion",#"AICc.seabnofam",
                    
                    "DICc.null", "DICc.n15", "DICc.c13", "DICc.nofam",
                    "DICc.n15nofam", "DICc.c13nofam", "DICc.fam", "DICc.fam*region",#"DICc.seabnofam", 
                    
                    "lambda.null", "lambda.n15", "lambda.c13", "lambda.nofam",
                    "lambda.nofamn15", "lambda.nofamc13", "lambda.fam", "lambda.fam*region", #"lambda.nofamseab",
                    
                    "R2m.null", "R2m.n15", "R2m.c13", "R2m.nofam",
                    "R2m.n15nofam", "R2m.c13nofam", "R2m.fam", "R2m.fam*region",
                    
                    "R2c.null", "R2c.n15", "R2c.c13",  "R2c.nofam",
                    "R2c.n15nofam", "R2c.c13nofam", "R2c.fam", "R2c.fam*region")
                    
fam.reg.names.use <-c("Acanthuridae_Caribbean", "Acanthuridae_Moorea","Chaetodontidae_Caribbean","Chaetodontidae_Moorea","Holocentridae_Caribbean", "Holocentridae_Moorea",
                  "Labridae_Caribbean", "Labridae_Moorea" , "Lutjanidae_Caribbean","Lutjanidae_Moorea", "Mullidae_Caribbean", "Mullidae_Moorea",         
                  "Pomacentridae_Caribbean","Pomacentridae_Moorea", "Scaridae_Caribbean","Scaridae_Moorea", "Serranidae_Caribbean","Serranidae_Moorea",       
                  "Tetraodontidae_Caribbean", "Tetraodontidae_Moorea" )
g$fam.reg.names <- paste(g$family,g$region,sep = "_")#create name to pull these above out

param.list.stoic <- list()
iterations <- 1
mcmc.its <- 80000
nburn <- 10000
resid.list <- list()
aov.list <- list()
mcmc.list <- list()

basicstats.mat <- matrix(nrow = length(stoic.name), ncol = 6, dimnames = list(NULL,c("Poly.ind", "Carib.ind","Poly.sp", "Carib.sp", "Poly.fam", "Carib.fam")))

quartz()
par(mfrow = c(3,3))

for(n in 1:length(stoic.name)) {
 
  print(stoic.trait.names[n])

  ind.ns <- na.omit(data.frame( gen.sp = g$gen.sp, 
                                 param = g[ ,c(stoic.trait.names[n]) ],
                                 region = g$region,
                                 fam = g$family,
                                 c13 = g$c.c13,
                                 n15 = g$c.n15 )  )
  
  fam.keep <- names(which( ifelse( table(ind.ns$fam,ind.ns$region)[,1] > 1 & table(ind.ns$fam,ind.ns$region)[,2] > 1, TRUE, FALSE) == TRUE ))
  
  #isolate the families that are in both regions
  yeyo <- droplevels(ind.ns[ind.ns$fam %in% fam.keep,])
  basicstats.mat[n,] <- c( dim(yeyo[yeyo$region == "Moorea",])[1],
                           dim(yeyo[yeyo$region == "Caribbean",])[1],
                           length(unique(yeyo[yeyo$region == "Moorea",]$gen.sp)),
                           length(unique(yeyo[yeyo$region == "Caribbean",]$gen.sp)),
                           length(unique(yeyo[yeyo$region == "Moorea",]$fam)),
                           length(unique(yeyo[yeyo$region == "Caribbean",]$fam)) )
  
  #matrix to put the full model output in
  stacks <- matrix(ncol = length(namethemthangs), nrow = iterations, dimnames = list(NULL, namethemthangs) )
  #list to store the params per iterations
  param.list.per.its <- list()
  
  # ####
  its <- 1
    print( paste("iteration = ",its, " param = ", stoic.name[n]) )
    
    #getting full phylo data now
    sp.list <- as.character(yeyo$gen.sp)
    #Drop any tips from the tree that are not in the species list - note that name.check will still fail for this, even though it has dropped all tips not in the data.  
    sp.tree<-drop.tip(sp.tree.full, setdiff(sp.tree.full$tip.label, sp.list))
    
    #see if there is something missing 
    if( length(setdiff(sp.list, sp.tree$tip.label)) > 0 ) { " DOLLY PARDON !"}
    
    sp.tree.nu <- force.ultrametric(sp.tree, method="nnls")
    
    Ainv<-inverseA(sp.tree.nu, nodes = "TIPS")$Ainv#this converts it to the variance/covariance matrix inverse matrix of shared phyloegnetic history
    
    # ####
   
    prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002)) #THIS SHOULD BE WIDE OPEN PRIORS
    #prior1.2 <- list(G = list(G1 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000) ),
    #               R = list(V=1, fix=1))
    
    null<-MCMCglmm.updateable(log(param+1) ~ -1 + region:fam, 
                            random=~gen.sp, ginverse=list(gen.sp=Ainv), 
                            data=yeyo, prior = prior1.1, verbose=FALSE, singular.ok = FALSE,
                            nitt=mcmc.its,thin=50,burnin=nburn,
                            pr = TRUE )#need to generate conditional preds
    
    summary(null)
    
    preds.marg <- predict.MCMCglmm(null)
    preds.cond <- predict.MCMCglmm(null, marginal = NULL)
    resids.marg <- log(yeyo$param+1) - preds.marg  
    resids.cond <- log(yeyo$param+1) - preds.cond
    
    #generate outputs for AOV
    
    
    #prior1.2 <- list(G = list(G1 = list(V = 1, nu = 0.002)) ) #THIS SHOULD BE WIDE OPEN PRIORS
    null.aov <- MCMCglmm.updateable(log(param+1) ~ -1 + region:fam, 
                        #random=~gen.sp, ginverse=list(gen.sp=Ainv), 
                        data=yeyo,  verbose=FALSE, singular.ok = FALSE,
                        nitt=mcmc.its,thin=50,burnin=nburn,
                        pr = TRUE )
    summary(null.aov)
    
    preds.aov <- predict.MCMCglmm(null.aov)
    resids.aov <- log(yeyo$param+1) - preds.aov  
    
    #null.aov <- aov(log(param+1) ~ -1 + region:fam, data=yeyo)
    #resids.aov <- resid(null.aov)
    
    print(summary(null))
    print(summary(null.aov))
    
    resid.list[[n]] <- data.frame(resids.marg, resids.cond, resids.aov)
    aov.list[[n]] <- summary(null.aov)[[5]]
    mcmc.list[[n]] <- summary(null)[[5]]
    plot(resids.marg, resids.aov)
}
    
saveRDS(resid.list, file = "data/output/ResidOutput_Region*Fam*Phylo-NoPhylo-fulldataset", 
        ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL) #save is as a list  
saveRDS(aov.list, file = "data/output/AOVSummaryOutput_Region*Fam*Phylo-NoPhylo-fulldataset", 
        ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL) #save is as a list  
saveRDS(mcmc.list, file = "data/output/MCMCSummaryOutput_Region*Fam*Phylo-NoPhylo-fulldataset", 
        ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL) #save is as a list  
write.csv(data.frame(basicstats.mat),file="data/output/numb.fams")
