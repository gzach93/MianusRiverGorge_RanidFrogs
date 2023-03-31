#Packages
library(vegan)
library(ggplot2)
library(viridis)
library(multcompView)
library(plyr)
library(tidyr)
library(dplyr)
library(reshape2)
library(ggpubr)
library(cowplot)
library(ggnewscale)
library(pairwiseAdonis)
library(lsmeans)
library(MASS)
library(indicspecies)
library(EcolUtils)


#Load in Data
#ASV Data
asv <- read.csv('Data/feature-table-with-taxonomy.csv')

#Taxonomy data and create columns for Domain down
tax <- data.frame(asv[,79])
tax$ID <- asv$X.OTU.ID

#Remove Taxonomy from asv file to run dissimlarity matrices
asv <- asv[,-79] 

#Sample Information
frog <- read.csv('Data/Frog.csv')

#Alpha Diversity Measurements
alpha <- read.csv('Data/AlphaDiversityMeasurements.csv')


#Format Data with Samples as Rows
col.name <- asv$X.OTU.ID
asv <- asv[,-1]
asv <- data.frame(t(asv))
colnames(asv) <- col.name

#Change Columns to Numeric
for (i in 1:1079) {
  asv[,i] <- as.numeric(asv[,i])
}

#Create ID Column in ASV Data Frame
asv$ID <- row.names(asv)

#Order the two data sets by the sample ID
frog <- frog[order(frog$ind),]
asv <- asv[order(asv$ID),]
alpha <- alpha[order(alpha$X),]

#Combine sample information with alpha diversity measurements for later
alpha$site <- frog$location
alpha$species <- frog$species
alpha$infection <- frog$infect

#Take out environmental swabs -- 
#just to have amphibian samples only for PERMANOVA and betadisper
asv$species <- frog$species
asv$site <- frog$location

asv.amph <- asv[which(asv$species == 'WF' | asv$species == 'GF' | asv$species == 'PF' |
                        asv$species == 'GF' | asv$species == 'Tad'),]

asv.amph$species <- as.character(asv.amph$species)
asv.amph$species <- as.factor(asv.amph$species)

# ALPHA DIVERSITY by site -- removed environmental samples 
alpha.amph <- alpha[which(alpha$species == 'WF' | alpha$species == 'GF' | alpha$species == 'PF' | alpha$species == 'GF' | alpha$species == 'Tad') ,]

alpha.amph$species <- as.character(alpha.amph$species)
alpha.amph$species <- as.factor(alpha.amph$species)

#NMDS
asv.nmds.bray <- metaMDS(asv[,-c(1080:1082)], distance = 'bray',
                         k = 2, try = 50)
asv.nmds.jac <- metaMDS(asv[,-c(1080:1082)], distance = 'jaccard', 
                        k = 2, try = 50, binary = T)

#Merge NMDS points with Frog data -- Bray
pt.nmds.bray <- cbind(unname(asv.nmds.bray$points[,1]),
                      unname(asv.nmds.bray$points[,2]))

colnames(pt.nmds.bray) <- c('bray_x.nmds', 'bray_y.nmds')
frog <- cbind(frog, pt.nmds.bray)

#Merge NMDS points with Frog data -- Jaccard
pt.nmds.jac <- cbind(unname(asv.nmds.jac$points[,1]),
                     unname(asv.nmds.jac$points[,2]))

colnames(pt.nmds.jac) <- c('jac_x.nmds', 'jac_y.nmds')
frog <- cbind(frog, pt.nmds.jac)


#### PERMANOVA
#Bray-Curtis by species
perm_species_bray <- adonis2(asv.amph[,-c(1080:1082)] ~ asv.amph$species, method = 'bray', permutations = 999)
perm_species_bray

bray.perm.pair <- pairwise.adonis(asv.amph[,-c(1080:1082)], asv.amph$species)
bray.perm.pair

#Brary by site -- removed tadpoles and vernal pool
bray.perm.site <- adonis2(asv.amph[-which(asv.amph$site == 'Vernal Pool'),-c(1080:1082)]~asv.amph[-which(asv.amph$site == 'Vernal Pool'), 'site'], method = 'bray', permutations = 999)
bray.perm.site

bray.perm.pair.site <- pairwise.adonis(asv.amph[-which(asv.amph$site == 'Vernal Pool'),-c(1080:1082)], asv.amph[-which(asv.amph$site == 'Vernal Pool'), 'site'])
bray.perm.pair.site

#Jaccard by species
perm_species_jaccard <- adonis2(asv.amph[,-c(1080:1082)] ~ asv.amph$species, method = 'jaccard', permutations = 999, binary = T)
perm_species_jaccard

jaccard.perm.pair.species <- pairwise.adonis(asv.amph[,-c(1080:1082)], asv.amph[, 'site'])
jaccard.perm.pair.site

#Jaccard by site -- removed tadpoles and vernal pool
jaccard.perm.site <- adonis2(asv.amph[-which(asv.amph$site == 'Vernal Pool'),-c(1080:1082)]~asv.amph[-which(asv.amph$site == 'Vernal Pool'), 'site'], method = 'jaccard', permutations = 999, binary = T)
jaccard.perm.site

jaccard.perm.pair.site <- pairwise.adonis(asv.amph[-which(asv.amph$site == 'Vernal Pool'),-c(1080:1082)], asv.amph[-which(asv.amph$site == 'Vernal Pool'), 'site'])
jaccard.perm.pair.site


#BETA DISPER MULTIVARIATE HOMOGENEITY OF DISPERSION
#By Species
asv.amph.bray <-vegdist(asv.amph[,-c(1080:1082)], method = 'bray', binary = F)
disp.species.bray <- betadisper(asv.amph.bray, asv.amph$species, type = "centroid")
disp.species.bray


asv.amph.jac <-vegdist(asv.amph[,-c(1080:1082)], method = 'jaccard', binary = T)
disp.species.jac <- betadisper(asv.amph.jac,asv.amph$species, type = "centroid")
disp.species.jac

#By site
asv.site.bray <-vegdist(asv.amph[-which(asv.amph$site == 'Vernal Pool'),-c(1080:1082)], method = 'bray', binary = F)
disp.site.bray <- betadisper(asv.site.bray,  asv.amph[-which(asv.amph$site == 'Vernal Pool'), 'site'], type = "centroid")
disp.site.bray


asv.site.jac <-vegdist(asv.amph[-which(asv.amph$site == 'Vernal Pool'),-c(1080:1082)], method = 'jaccard', binary = T)
disp.site.jac <- betadisper(asv.site.jac,  asv.amph[-which(asv.amph$site == 'Vernal Pool'), 'site'], type = "centroid")
disp.site.jac

#Check for differences in dispersion amoung groups
#Species/Lifestage
permutest(disp.species.bray, permutations = 999, pairwise = FALSE) 
permutest(disp.species.bray, permutations = 999, pairwise = TRUE) 


permutest(disp.species.jac, permutations = 999, pairwise = FALSE) 
permutest(disp.species.jac, permutations = 999, pairwise = TRUE) 

#Site - vernal pool
permutest(disp.site.bray, permutations = 999, pairwise = FALSE) 
permutest(disp.site.jac, permutations = 999, pairwise = FALSE) 

#Infection GLM
alpha.pf <- alpha.amph[which(alpha.amph$species == 'PF'),]

inf.mod.shannon <- glm(infection~shannon, alpha.pf, family = 'binomial')
summary(inf.mod.shannon)

inf.mod.faith <- glm(infection~faith, alpha.pf, family = 'binomial')
summary(inf.mod.faith)

inf.mod.rich <- glm(infection~observed_OTU, alpha.pf, family = 'binomial')
summary(inf.mod.rich)


frog.amph <- na.omit(frog)

#Dominate Taxa
new.asv <- asv
new.asv <- new.asv[order(new.asv$ID),]
frog <- frog[order(frog$ind),]
new.asv$species <- frog$species

new.asv <- new.asv[,-1080]

asv.ag <- data.frame(t(aggregate(new.asv[,-1080], 
                                 by = list(new.asv$species), FUN = mean)))
colnames(asv.ag) <- asv.ag[1,]
asv.ag <- asv.ag[-c(1),]
asv.ag$asv <- rownames(asv.ag)

#edit tax data frame
tax <- separate(data = tax, col = asv...79.,
                sep = ";", into = c('Kingdom', 
                                    'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
                fill = "right")


#Edit Tax names
tax$Kingdom <- substr(tax$Kingdom, start = 6, stop = nchar(tax$Kingdom))
tax$Phylum <- substr(tax$Phylum, start = 7, stop = nchar(tax$Phylum))
tax$Class <- substr(tax$Class, start = 7, stop = nchar(tax$Class))
tax$Order <- substr(tax$Order, start = 7, stop = nchar(tax$Order))
tax$Family <- substr(tax$Family, start = 7, stop = nchar(tax$Family))
tax$Genus <- substr(tax$Genus, start = 7, stop = nchar(tax$Genus))
tax$Species <- substr(tax$Species, start = 7, stop = nchar(tax$Species))

#asv.ag <- asv.ag[-1079,]
asv.ag$family <- NA

for(i in 1:nrow(asv.ag)){
  row <- which(asv.ag$asv[i] == tax$ID)
  asv.ag$family[i] <- as.character(tax$Family[row])
}

asv.ag[,-c(10,11)] <- sapply(asv.ag[,-c(10,11)], as.numeric)

asv.ag.fam <- aggregate(asv.ag[,-c(10,11)], by = list(asv.ag$family), FUN = sum)

percent.fam <- data.frame(matrix(ncol= 10, nrow = nrow(asv.ag.fam)))

for(i in 2:ncol(asv.ag.fam)){
  percent.fam[,i] <- asv.ag.fam[,i]/sum(asv.ag.fam[,i])
}

percent.fam$X1 <- asv.ag.fam$Group.1
colnames(percent.fam) <- colnames(asv.ag.fam)


percent.fam[which(percent.fam$UP < .05 & 
                    percent.fam$VP < .05 &
                    percent.fam$Bog < .05 &
                    percent.fam$Tun < .05 &
                    percent.fam$Ann < .05 &
                    percent.fam$GF < .05 &
                    percent.fam$Tad < .05 &
                    percent.fam$WF < .05 &
                    percent.fam$PF < .05), 'Group.1'] <- 'Other'

melt.percent.fam <- melt(percent.fam)

melt.percent.fam$variable <- factor(melt.percent.fam$variable,levels = c("Ann", "Bog", "Tun", "UP", "VP", 
                                                                         "GF", "PF", "WF", "Tad"))


#Finding Dominate family and what Percentage they make up
#Tadpole
which(percent.fam$Tad > .1)
percent.fam[c(17,21),'Group.1']
percent.fam$Tad[c(17,21)]

#GF
which(percent.fam$GF > .1)
percent.fam[c(21),'Group.1']
percent.fam$GF[c(21)]

#WF
which(percent.fam$WF > .1)
percent.fam[c(21,124),'Group.1']
percent.fam$WF[c(21,124)]

#PF
which(percent.fam$PF > .1)
percent.fam[c(21,40, 101),'Group.1']
percent.fam$PF[c(21,40,101)]

#Adding Information for Plotting
frog$type <- NA
frog$type <- 'Amphibian'
frog[which(frog$species == 'VP' | frog$species == 'UP' | frog$species == 'Bog' |
             frog$species == 'Tun' | frog$species == 'Ann'),'type'] <- 'Environment'


#frog <- frog[-which(frog$type == 'Environment'),]

frog$species <- as.factor(frog$species)
frog1 <-  frog[which(frog$type == 'Environment'),]
frog2 <-  frog[-which(frog$type == 'Environment'),]



#Alpha Diversity GLMs
rich.sp <- glm.nb(observed_OTU~species, alpha.amph)
summary(rich.sp)
emmeans(rich.sp, specs = pairwise~species)

shan.sp <- glm(shannon~species, alpha.amph, family = Gamma)
summary(shan.sp)
emmeans(shan.sp, specs = pairwise~species)


rich.site <- glm.nb(observed_OTU~site, alpha.amph)
summary(rich.site)
emmeans(rich.site, specs = pairwise~site)

shan.site <- glm(shannon~site, alpha.amph, family = Gamma)
summary(shan.site)
emmeans(shan.site, specs = pairwise~site)

#Site Residual Models
rich.site.res <- lm(rich.sp$residuals~alpha.amph$site, alpha.amph)
summary(rich.site.res)

shan.site.res <- lm(shan.sp$residuals~alpha.amph$site)
summary(shan.site.res)

#How are Species Distrupted Across Site 1 -- difference in subsites?
table(alpha.amph$site, alpha.amph$species)


#indicator species anaylsis -- what in tadpole samples is difference than adults amphibians

ind.sp <- multipatt(asv.amph[,-c(1080:1082)], asv.amph$species, control = how(nperm=999))
summary(ind.sp)

xx <- ind.sp$sign
uni.tad <- xx[which(xx$s.Tad == 1 & xx$p.value < 0.05 & xx$s.GF == 0  & xx$s.PF == 0 & xx$s.WF == 0),]

uni.tad$Family <- NA
uni.tad$Phylum <- NA

for(i in 1:nrow(uni.tad)){
  uni.tad$Family[i] <- tax[which(tax$ID == rownames(uni.tad)[i]), 'Family']
}


for(i in 1:nrow(uni.tad)){
  uni.tad$Phylum[i] <- tax[which(tax$ID == rownames(uni.tad)[i]), 'Phylum']
}


uni.tad$Family <- as.factor(uni.tad$Family)
uni.tad$Phylum <- as.factor(uni.tad$Phylum)


count(uni.tad,Family)
count(uni.tad,Phylum)















