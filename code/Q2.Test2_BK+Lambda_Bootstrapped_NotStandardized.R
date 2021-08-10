###phylo TEST YALL!##
#this code generates Mantel tests - which I am not using, Phylo corr test - which I am not using and pagels lambda and bloomberg k - both of which I do use 
rm(list = ls())
library(data.table)
library(gplots)
library(picante)
library(phytools)
library(caper)

###############################################################
##########################Functions############################
###############################################################
stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) } #standardize variables
center <- function(X) {(X-mean(X,na.rm=T))}
#make the matrices
make.mat <- function(data,columnofinterest){
  sp.ord <- unique(data$to.sp)
  new.mat <- matrix(ncol = length(sp.ord), nrow = length(sp.ord), dimnames = list(sp.ord, sp.ord))
  for(i in 1:length(sp.ord)){
    new.mat[sp.ord[i],which(colnames(new.mat) %in% data[data$to.sp == sp.ord[i],"from.sp"])] <- data[data$to.sp == sp.ord[i],columnofinterest]
  }
  diag(new.mat) <- 0
  return(new.mat)
}
###############################################################

###############################################################
# Load data and set up g.tab####
###############################################################

g <- read.table("data/modified/Allgeier_FullStoicDataset_woPalmyra",stringsAsFactors=F, header=T) #these are fish only data 
g.tab <- data.frame( gen.sp = as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ),
                     FG = (tapply(g$FG,g$gen.sp,unique)),
                     family = tapply(g$family,g$gen.sp,unique),
                     region = tapply(g$region,g$gen.sp, unique) , stringsAsFactors = F)
g <- merge(g.tab,g)
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
###############################################################
#because the below metrics are calculated on means per species I iterate through the data varying the number of individuals required per species for the analysis
#that is what 'sample size' is all about
#sample.size <- 5
boot <- 1000

alldata <- list(g,
                g[g$region == "Moorea",],
                g[g$region == "Caribbean",])
names(alldata) <- c("All", "Moorea", "Caribbean")

for(region in 1:length(alldata)){

  gg <- alldata[[region]]
  
  All.sample.sizes <- list()
  
    for(n in 1:length(stoic.name)) {
    
    final.output <- matrix(ncol = 6, nrow = boot, 
                         dimnames = list(NULL,  c("lambda","lambda.p","bk","bk.p","sample.size","n.species")))
  
      for(size in 1:boot){
        
      print(size)
      print(stoic.trait.names[n])
      
      #get the means per species
      ggg <- gg[sample(1:length(gg[,1]),replace = TRUE),]
      
      if(standardize == "yes") { ggg[ ,c(stoic.trait.names[n]) ] <- stand( ggg[ ,c(stoic.trait.names[n]) ] ) }
      if(centerize == "yes") { ggg[ ,c(stoic.trait.names[n]) ] <- center( ggg[ ,c(stoic.trait.names[n]) ] ) }
      
        tot.list <- list()
        
        #isolate species and trait
        now <- data.frame( gen.sp = names(tapply(ggg[,c(stoic.trait.names[n])],list(ggg$gen.sp),mean,na.rm =T)),
                           param = matrix(tapply(ggg[,c(stoic.trait.names[n])],ggg$gen.sp,mean,na.rm =T)),
                           se = matrix(tapply(ggg[,c(stoic.trait.names[n])],ggg$gen.sp,sd,na.rm =T)/
                                         sqrt(sapply(tapply(ggg[,c(stoic.trait.names[n])],ggg$gen.sp,unique,na.rm =T),length))) )
        now <- now[!is.nan(now$param),] 
        
        #for species with no SD I am giving the mean SD
        now$se <- ifelse(is.na(now$se), mean(na.omit(now$se)), now$se)
        
        ############################################################
        ##########     BRING IN THE PHYLOGENY       ################
        #Make the phylogeny and force it to be ultrametric
        sp.list <- as.character(unique(now$gen.sp))
        #Read the tree that was downloaded
        sp.tree.full <- read.newick("data/raw/actinopt_12k_treePL.tre") 
        #Drop any tips from the tree that are not in the species list - note that name.check will still fail for this, even though it has dropped all tips not in the data.  
        sp.tree<-drop.tip(sp.tree.full, setdiff(sp.tree.full$tip.label, sp.list))
        
        #If you want to know what's missing - right now, you get 183 species from your dataset
        if( length(setdiff(sp.list, sp.tree$tip.label)) > 0 ) { " PUCK SKIP !"}
        
        #At this pont, sp.tree is essentially ultrametric, but there are some minor differences in branch length, which is fine, you can force it to be ultrametric
        sp.tree.nu <- force.ultrametric(sp.tree, method="nnls")
        ############################################################
        tot.list <- list()
        
        to.use <- now[,2]
        to.use.se <- now[,3]
        names(to.use) <- now$gen.sp
        names(to.use.se) <- now$gen.sp
        #without sampling error - don't need because I am bootstrapping
        bk <- phylosig(sp.tree.nu,to.use,  test = TRUE, method = "K", nsim = 999)
        l <-  phylosig(sp.tree.nu,to.use,  test = TRUE, method = "lambda", nsim = 999)
        
        #with sampling error
        #bk <- phylosig(sp.tree.nu,to.use, se = to.use.se, test = TRUE, method = "K", nsim = 999)
        #l <-  phylosig(sp.tree.nu,to.use, se = to.use.se, test = TRUE, method = "lambda", nsim = 999)
        
        final.output[size,]<- c(l$lambda,
                             l$P,
                             bk$K,
                             bk$P,
                             size,
                             dim(now)[1])
        
    }
    
    All.sample.sizes[[n]] <- final.output
      
    }
    
  subfile <- paste(names(alldata)[region],"_BK.Lambda.AllSampleSizes","_Boot=",boot,".csv",sep = "")
  
  write.csv(do.call(rbind.data.frame,All.sample.sizes), file = paste("data/output/",subfile, sep = "") ) 
  
}

    




