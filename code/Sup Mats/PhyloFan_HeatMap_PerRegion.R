rm(list= ls())
# Load data ####
g <- read.table("data/modified/Allgeier_FullStoicDataset_woPalmyra",stringsAsFactors=F, header=T)#these are fish only data 

g.tab <- data.frame( gen.sp = as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ),
                     family = tapply(g$family,g$gen.sp,unique),#[ as.character(c(names(tapply(g$FG,g$gen.sp,unique))) ) ])),
                     region = tapply(g$region,g$gen.sp, unique) , stringsAsFactors = F)
g <- merge(g1,g.tab)

library(phytools)
library(geiger)

stand <- function(X) { (X-mean(X,na.rm=T))/(2*sd(X,na.rm=T)) } #standardize variables

############################################################
##########     BRING IN THE PHYLOGENY       ################
#Make the phylogeny and force it to be ultrametric
sp.list <- as.character(unique(g$gen.sp))
#Read the tree that was downloaded
sp.tree.full <- read.newick("data/raw/actinopt_12k_treePL.tre") 
#Drop any tips from the tree that are not in the species list - 
#note that name.check will still fail for this, even though it has dropped all tips not in the data.  
sp.tree<-drop.tip(sp.tree.full, setdiff(sp.tree.full$tip.label, sp.list))

#If you want to know what's missing - right now, you get 183 species from your dataset
if( length(setdiff(sp.list, sp.tree$tip.label)) > 0 ) { " FUCK SHIT !"}

sp.tree.nu <- force.ultrametric(sp.tree, method="nnls")

trait.to.use <- c("nexc.m" , "pexc.m", "npexc", "bodyc", "bodyn", "bodyp", "bodycn", "bodycp", "bodynp")
trait.to.use.proppa <- c("Exc N" , "EXC P", "Exc NP", "Body C", "Body N", "Body P", "Body CN", "Body CP", "Body NP")


for(i in 1:2) {
  
  now <- g[g$region == unique(g$region)[i],]
  
  pdf( paste( "output plots/Sup Mats/", "FanPhylogeny_AllTraits_",unique(g$region)[i],".pdf", sep = "" ))
  
  for(ii in 1:length(trait.to.use)){
    
    trait <- tapply(now[,trait.to.use[ii]], now$gen.sp,mean, na.rm = T, drop = T)[!is.nan(tapply(now[,trait.to.use[ii]], now$gen.sp,mean, na.rm = T, drop = T))]
    
    #this makes sure the tree matches how I culled the species pool
    obj<-name.check(sp.tree,trait)
    treetouse <- drop.tip(sp.tree,  obj$tree_not_data )
    
    #this plots the phylogeny and the traits
    #likelihood reconstruction of the ancestral states and interpolates in between them
    contMap(treetouse, trait, fsize = .5, lwd = 1.5, outline = F, type = "fan") #direction = "upwards")
    legend("topright", paste(unique(g$region)[i],trait.to.use.proppa[ii], sep = " "), bty = "n", cex = 1.5)
    
  }
 
  dev.off()
  
}



