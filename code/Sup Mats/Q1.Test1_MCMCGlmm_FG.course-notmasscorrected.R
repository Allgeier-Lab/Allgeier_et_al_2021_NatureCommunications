#QUESTION 1: Not mass-corrected excretion
rm(list = ls())
library(dplyr)
library(fishtree)
library(operators)
library("ape")
library("geiger")
library(phytools)
library(MCMCglmm)
library(MuMIn)
library(data.table)

#######################
# functions needed ####
#######################

stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) } #standardize variables
#fucntion to quantify strength of random effects (phylogeny)
p.lambda <- function(mod){
  (r.effects <- ( mod$VCV[,1] )/ ( mod$VCV[,1] + mod$VCV[,2] )) 
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


# ####

g <- read.table("data/modified/Allgeier_FullStoicDataset_woPalmyra",stringsAsFactors=F, header=T)#these are fish only data 
stoic.trait <- list( c('nexc', "bodyn" ),
                     c('pexc', "bodyp" ),
                    c("npexc", "bodyn","bodyp"),
                    c("bodyn", "bodyp"), #used n and p not c because the correlation between C and P is 0.6 and the correlation between N and P is greater than N and C
                    c("bodyp", "bodyn", "bodyc"),
                    c("bodyc", "bodyp"), #based on the old stoic paper
                    c("bodynp", "bodyc"),
                    c("bodycn","bodyp"),
                    c("bodycp","bodyn") )

#output list - for all params and two matrices
allstoic.params <-list()
R2.glmm.mat <- matrix(ncol = 9,nrow = length(stoic.trait), dimnames = list(NULL,c("R2m","R2m.CIl","R2m.CIu","R2c","R2c.CIl","R2c.CIu","R2RE","R2RE.CIl","R2RE.CIu")))
Lambda.mat <- matrix(ncol = 3,nrow = length(stoic.trait), dimnames = list(NULL,c("Lambda","L.CIl","L.CIu")))
DIC.mat <- matrix(ncol = 2,nrow = length(stoic.trait), dimnames = list(NULL,c("Full","Intercept")))
#where the plot outputs go
#setwd("/Users/jeallg/My Work Biatch/Fish Datasets/Caribbean+Pacific Fish Stoic/Evolution of fish pee/R output/StoicAnalysis/traceplots")
nitts <- 80000
nburn <- 10000

for(i in 1:length(stoic.trait)){
   
   #load in the correct stoic data - and pull NAs
   traitstouse <- c("mass","n15","c13",stoic.trait[[i]])
   m <- g[-unique(which(is.na(g[,c(traitstouse)]), arr.ind=TRUE)[,1]),]#identify and remove rows with NA
   
   #standardize all data
   now <- data.frame(region = as.factor(m$region),
                     FG = as.factor(m$FG.course),
                     fam = as.factor(m$family), #can't write family because of the mcmcmglmm calls
                     genus = as.factor(m$genus),
                     gen.sp = as.character(m$gen.sp),
                     m[,traitstouse],
                     apply(m[,traitstouse], 2, function(x) stand(as.numeric(as.character(x))) ) )
   colnames(now)[ (dim(now)[2] - length(traitstouse)+1) : dim(now)[2]  ] <- c(paste("s",c(traitstouse),sep = ""))
   #no fucking clue why this won't convert this to character.
   now$gen.sp <- as.character(now$gen.sp)
   
   ############################################################
   ##########     BRING IN THE PHYLOGENY       ################
   ############################################################
   #Make the phylogeny and force it to be ultrametric
   sp.list <- as.character(unique(now$gen.sp))
   #Read the tree that was downloaded
   sp.tree.full <- read.newick("data/raw/actinopt_12k_treePL.tre") 
   #Drop any tips from the tree that are not in the species list - note that name.check will still fail for this, even though it has dropped all tips not in the data.  
   #I think this is because there are repeated values for each species
   sp.tree<-drop.tip(sp.tree.full, setdiff(sp.tree.full$tip.label, sp.list))
   
   #If you want to know what's missing - right now, you get 183 species from your dataset
   if( length(setdiff(sp.list, sp.tree$tip.label)) > 0 ) { " DOLLY PARDON !"}
   
   #At this pont, sp.tree is essentially ultrametric, but there are some minor differences in branch length, which is fine, you can force it to be ultrametric
   sp.tree.nu <- force.ultrametric(sp.tree, method="nnls")
   obj<-name.check(sp.tree.nu,sp.list) 
   
   #sp.tree <- force.ultrametric(sp.tree.nu, method=c("nnls")) #double check to see what the fuck nnls and extend mean - send to Brian
   sum(sp.tree.nu$tip.label %in% sp.list ) #check if all are there 
   sp.tree$Nnode# equals n-1 number of species (here 164)
   is.ultrametric(sp.tree.nu)
   tips<-sp.tree.nu$tip.label
   #using the now data 
   Ainv<-inverseA(sp.tree.nu, nodes = "TIPS")$Ainv#this converts it to the variance/covariance matrix inverse matrix of shared phyloegnetic history
   
   #this is so I can use dredge
   MCMCglmm.updateable<- updateable(MCMCglmm)
   
   #for the prior term you need a prior (G1, G2, etc.) for any random effect and a proir for the variance structure (R)
   s.trait <- paste("s",stoic.trait[[i]][1],sep = "")
   
   if( length( paste("+ s",stoic.trait[[i]][-1], sep = "") ) == 2 ) {
    stoic.preds <- paste( "+ s",stoic.trait[[i]][2], " + s",stoic.trait[[i]][3], sep = "") 
   } else { stoic.preds <-  paste("+ s",stoic.trait[[i]][-1], sep = "") }
   
   #keep mass in the equation for excretion   
   mods <- formula(paste(s.trait, "~"," -1 +  smass + FG + region + region:sn15 + region:sc13 ",stoic.preds)) 
   
   
   prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1, nu = 0.002)) #THIS SHOULD BE WIDE OPEN PRIORS
   m1<-MCMCglmm.updateable(mods, 
                           random=~gen.sp, ginverse=list(gen.sp=Ainv), 
                           data=now, prior = prior1.1, verbose=TRUE,
                           nitt=nitts,thin=50,burnin=nburn) #prior = prior.biv,
   
   
   #fixed <- paste(stoic.trait[[i]][1],"FE","pdf", sep = ".")
   #pdf( file= paste("data/output/plots/",fixed,sep = "" ) )
   plot(m1$Sol)
   
   #dev.off()
   #random <- paste(stoic.trait[[i]][1],"RE","pdf", sep = ".")
   #pdf( file= paste("data/output/plots/",random,sep = "" ))
   plot(m1$VCV)
   #dev.off()
   
   param.order <-  c("sbodyc",
                     "sbodyn",
                     "sbodyp",
                     "smass",
                     "regionMoorea",  
                     "regionMoorea:sc13", 
                     "regionCaribbean:sc13",
                     "regionMoorea:sn15", 
                     "regionCaribbean:sn15",
                     "FGherbivore.detritivore",
                     "FGinvertivore",
                     "FGomnivore", 
                     "FGpredator")                                                                      
   
   poopoo <- data.frame(params = names(posterior.mode(m1$Sol)),
                        mode = posterior.mode(m1$Sol),
                        HPDinterval(m1$Sol),
                        mean = apply(m1$Sol,2,mean),
                        median = apply(m1$Sol,2,median))
   pp <- poopoo[param.order,]
   
   prior2 <- list(G = list(G1 = list(V = 1, nu = 2, alpha.mu = 0, alpha.V = 10000) ),
                  R = list(V=1, fix=1))
   
   mods2 <- formula(paste(s.trait, "~","1") )
   m2<-MCMCglmm.updateable(mods2, 
                           random=~gen.sp, ginverse=list(gen.sp=Ainv), 
                           data=now, prior = prior1.1, verbose=TRUE,
                           nitt=nitts,thin=50,burnin=nburn) 
   
   #matrix of all the conditional R2 and the CIs for the full model with the FEs
   R2.glmm.mat[i,] <- c ( R.squaredGLMM(m1) )
   #matrix of the conditional R2 of an intercept only model (which is functionally lambda) with CIs 
   Lambda.mat[i,] <-  c( p.lambda(m2) )
   DIC.mat[i,] <- c(DIC(m1),DIC(m2))
   allstoic.params[[i]] <- poopoo[param.order,]
   names(allstoic.params)[[i]] <- stoic.trait[[i]][1]
   
 }


write.csv(data.frame(R2.glmm.mat), file = "data/output/Sup Mats/AllStoicTraitModels_R2.glmm-notmasscorrected")
write.csv(Lambda.mat,file = "data/output/Sup Mats/AllStoicTraitModels_P.Lambda-notmasscorrected")
write.csv(DIC.mat,file = "data/output/Sup Mats/AllStoicTraitModels_DIC-notmasscorrected")
saveRDS(allstoic.params, file = "data/output/Sup Mats/AllStoicTraitModels_ParamOutputs-notmasscorrected", 
        ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL) #save is as a list   


