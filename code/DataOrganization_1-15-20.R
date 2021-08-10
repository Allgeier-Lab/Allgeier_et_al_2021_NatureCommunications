rm(list = ls())
library(operators)
library(dplyr)

'%!in%' <- function(x,y)!('%in%'(x,y))

#Moorea data
m1 <- read.csv("data/raw/Moorea_Compiled_2016_2017_Update-12-2019.csv",header=T,sep=",",stringsAsFactors = F) 
m1 <- m1[m1$island == "Moorea",]
m1$spID <- paste(m1$genus,m1$species,sep = "_")
m1$region <-  "Moorea"

mfg <- read.csv("data/raw/Moorea_FunctionaGroups.csv",header=T,sep=",") #%>%
mfg$spID <- paste(mfg$genus,mfg$species,sep = "_")

m <- droplevels(merge(m1,mfg[,c("spID","FG")]))

#Fix the names to match the Caribbean dataset
names(m)[c(8,9,10,19,20,21,25,26,27,28,29)] <- c("length","mass","dmass","nexc_min", "pexc_min", "urea_min","bodyn","bodyc","bodyp", "n15","c13")

#get excretion per hour to match Caribbean data
m$nexc <- as.numeric(m$nexc_min)*60
m$pexc <- as.numeric(as.character(m$pexc_min))*60
#get out the nonsensical data (determined based on Allgeier et al. 2015 PNAS)
m$pexc[m$pexc<.002] <- NA
m$nexc[m$nexc<.02] <- NA

#bring in Caribbean data
c <- droplevels(read.csv("data/raw/Caribbean_Compiled_2009-2013.csv",header=T,sep=",") )
c <- droplevels(c[c$v.i == "v",]) #get rid of the invertebrates
c <- c[c$isotope_only != 1,] #get rid of the samples that I didn't do excretion on
c$region = "Caribbean"

#get out all values that are questionable based on Allgeier et al. 2015 PNAS
c$pexc[c$pexc<.002] <- NA
c$nexc[c$nexc<.02] <- NA

head(m)
head(c)

#isolate the columns needed for analysis
ofinterest <- c("region", "family","genus","species",
                "mass","nexc","pexc","bodyn", "bodyp", 
                "bodyc", "n15", "c13", "FG")#add in functional group don't have length for some reason
b1 <- rbind(m[,ofinterest],c[,ofinterest])

#add in the super course FG
fg.c <- read.csv("data/raw/FG.course.csv",header=T,sep=",") 

b <- droplevels(merge(b1,fg.c))
table(b$region,b$FG)#check the FGs are lined up yo
table(b$region,b$FG.course)

#organizing the specific values
#these are transformed mass-specific rates
b$pexc.m <-   (b$pexc)/(b$mass) 
hist(b$pexc.m)
b$nexc.m <-  ( b$nexc/b$mass ) 
hist(b$nexc.m)
b$npexc<- (b$nexc/b$pexc * 31/14) 
hist(b$npexc)

#body stoic
b$bodycn <- b$bodyc/b$bodyn * 14/12 
hist(b$bodycn)
b$bodycp <- b$bodyc/b$bodyp * 31/12
hist(b$bodycp)
b$bodynp <- b$bodyn/b$bodyp * 31/14
hist(b$bodynp)

#Here I do some major name cleaning 
#there are issues in name linkages across all the data so here I clean all that up in a very unimpressive way
b$gen.sp <-  paste(b$genus,b$species,sep = "_")
b$family <- ifelse(b$family == "Blennidae", "Blenniidae", b$family)

b$gen.sp <- ifelse(b$gen.sp == "Haemulon_plumieri", "Haemulon_plumierii", 
                     ifelse ( b$gen.sp == "Eucinostmus_melanopterus", "Eucinostomus_melanopterus",
                              ifelse ( b$gen.sp == "Tylosurus_crocodilus", "Tylosurus_crocodilus_crocodilus", 
                                       ifelse( b$gen.sp == "Chrysiptera_leucopoma", "Chrysiptera_brownriggii",
                                               ifelse( b$gen.sp == "Chrysiptera_leucopoma", "Chrysiptera_brownriggii",
                                                       ifelse( b$gen.sp == "Myripristis_sp", "Myripristis_kuntee", 
                                                               ifelse( b$gen.sp == "Scarus_iserti", "Scarus_iseri",
                                                                       ifelse( b$gen.sp == "Holocentrus_coruscum", "Sargocentron_coruscum",
                                                                               ifelse( b$gen.sp == "Psudeopeneus_maculatus", "Pseudupeneus_maculatus",
                                                                                       ifelse( b$gen.sp == "paraclinus_nigripinnis", "Paraclinus_nigripinnis",
                                                                                               ifelse( b$gen.sp == "labrisomus_nuchipinnis", "Labrisomus_nuchipinnis",
                                                                                                       ifelse( b$gen.sp == "paraclinus_acuminatus", "Pareques_acuminatus", 
                                                                                                               ifelse(b$gen.sp == "Nicholsina_usta", "Nicholsina_usta_usta",
                                                       b$gen.sp ) ) ) ) ) ) ) ) ) ) ) ) )

###########################################################################
#these fish are missing from the phylogenetic tree all together and were not replaced for various reasons
###########################################################################

#Sparisoma_sp. - this sucks to lose these but I am not sure what else to do, I just can justify giving these individuals a sp name
#Scarus_coeruleus - the don't have and nothing seems close
#Malacoctenus_macropus Malacoctenus macropus
#"Carcharhinus_melanopterus"
#"Paraclinus_nigripinnis"    
#"Syngnathus floridae"

#these need to be looked into more
#Stegastes_punctatus
#Epinephelus_hexagonatus
#Exallias_brevis

#these are missing from the tree but I am replaing them with a close species or just correcting the name

#Cosmocampus_brachycephalus	Cosmocampus brachycephalus - replace with  Cosmocampus albirostris
#Calamus_bajonado replace with Calamus nodosus
#Syngnathus_springeri replacing with Syngnathus floridae
#Nicholsina_usta replace with Nicholsina_usta_usta
b$gen.sp <- ifelse(b$gen.sp == "Cosmocampus_brachycephalus", "Cosmocampus_albirostris", #this one is fine
                   ifelse(b$gen.sp == "Calamus_bajonado", "Calamus_nodosus", #this one is fine
                          ifelse(b$gen.sp == "Syngnathus_springeri", "Syngnathus floridae", #this one is fine
                                 b$gen.sp ) ) )

#sp names that I need to remove for various reasons
b <- b[b$gen.sp %!in% c("Sparisoma_sp.", "Scarus_coeruleus", "Malacoctenus_macropus", 
                        "Stegastes_punctatus", "sp_sp", "Epinephelus_hexagonatus",
                        "Exallias_brevis","Carcharhinus_melanopterus","Paraclinus_nigripinnis","Syngnathus floridae"), ]

###############

#these are species-specific ecological 'functional' traits from Moullet et al. 2014 PNAS
w <- read.csv("data/raw/traits_gaspar.csv",stringsAsFactors=F, header=T)  
w$gen.sp <- paste(w$Genus,w$Species,sep = "_")

#there are species in our dataset that are not in gaspar
#[1] "Pomachromis_fuscidorsalis"       "Pristiapogon_exostigma"          "Pristiapogon_kallopterus"        "Kyphosus_sectatrix"              "Nicholsina_usta_usta"            "Trinectes_inscriptus"           
#[7] "Rhinesomus_triqueter"            "Anguilla_rostrata"               "Tylosurus_crocodilus_crocodilus" "Atherinomorus_stipes"           

#match the above fish that are in our dataset and not in gaspar with similar ecological fishes 
w$gen.sp <- ifelse(w$gen.sp == "Pomachromis_richardsoni", "Pomachromis_fuscidorsalis", 
                   ifelse ( w$gen.sp ==  "Ostorhinchus_luteus", "Pristiapogon_exostigma", 
                            ifelse ( w$gen.sp == "Ostorhinchus_leslie", "Pristiapogon_kallopterus", #
                                     ifelse( w$gen.sp == "Tylosurus_crocodilus crocodilus", "Tylosurus_crocodilus_crocodilus", 
                                             ifelse( w$gen.sp == "Kyphosus_sectator", "Kyphosus_sectatrix",
                                                     ifelse( w$gen.sp == "Nicholsina_usta usta", "Nicholsina_usta_usta", 
                                                             ifelse( w$gen.sp == "Soleichthys_heterorhinos", "Trinectes_inscriptus", 
                                                                     ifelse( w$gen.sp == "Lactophrys_trigonus", "Rhinesomus_triqueter", 
                                                                             ifelse( w$gen.sp == "Gymnothorax_hubbsi", "Anguilla_rostrata", 
                                                                                     ifelse( w$gen.sp == "Atherion_elymus", "Atherinomorus_stipes", 
                                                                                                                    w$gen.sp ) ) ) ) ) ) ) ) ) ) 


suckkamothaclukka <- w[,c("Class","Order","Home.Range","Activity","Schooling","Level.water","gen.sp")]
g <- merge(suckkamothaclukka, b)
#test if the data transferred correctly
if(dim(g)[1] == dim(b)[1]) {"HELL YES"} else {"DAAAAG"}

#checking for any duplicate rows - there are 4 in the moorea dataset for some reason
g$keyID <- paste(g$gen.sp,g$nexc.m,g$mass,g$bodycn,g$pexc.m,g$npexc,g$bodyp,g$bodyc,sep = "_")
length(unique(g$keyID)) == dim(g)[1]
g <- g[-which(duplicated(g$keyID)),] #remove them


write.table(g,file = "data/modified/Allgeier_FullStoicDataset_woPalmyra")



