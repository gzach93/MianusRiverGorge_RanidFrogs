#Libraries
library(phyloseq)
library(emmeans)
library(tidyverse)

#Load in Data
#ASV Data
asv <- read.csv('Data/feature-table-with-taxonomy.csv')

#Taxonomy data and create columns for Domain down
tax <- data.frame(asv[,79])
tax$ID <- asv$X.OTU.ID

#Remove Taxonomy from asv file to run dissimlarity matrices
asv <- asv[,-79] 


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

otu_table(as.matrix(asv), taxa_are_rows = FALSE)


test <- read.csv('~/Documents/2020_ZG18A_ASV_Qiime2/ZG18A_Final_TABLE_24May/feature-table.csv')
taxa <- read.csv('~/Documents/2020_ZG18A_ASV_Qiime2/ZG18A_Final_TABLE_24May/feature-table.csv')

test <- phyloseq(otu_table(test[,-1], taxa_are_rows = TRUE))



##############################################################
################## Infection Model ###########################
##############################################################

#Load in Data
frog.data <- read.csv('~/Desktop/FullFrog_Data.csv')
colnames(frog.data) <- c('year', 'site', 'subsite', 'ID', 'species',
                         'mass', 'svl', 'infection', 'total.length', 'head.length',
                         'tail.length', 'notes')


#Infection Model - Species Adults
sp.inf <- glm(infection ~ species,
              frog.data[-which(frog.data$species == 'L. tadpole'),], 
              family = binomial);summary(sp.inf)
emmeans(sp.inf, specs = pairwise~species)

table(frog.data[-which(frog.data$species == 'L. tadpole'),'species'], 
      frog.data[-which(frog.data$species == 'L. tadpole'),'infection'])

#Infection Model - sites Adults
site.inf <- glm(infection ~ site,
              frog.data[-which(frog.data$species == 'L. tadpole'),], 
              family = binomial);summary(site.inf)
emmeans(site.inf, specs = pairwise~site)

table(frog.data[-which(frog.data$species == 'L. tadpole'),'site'], 
      frog.data[-which(frog.data$species == 'L. tadpole'),'infection'])

#Infection Model - sites Adults
site.inf.tad <- glm(infection ~ site,
                frog.data[which(frog.data$species == 'L. tadpole'),], 
                family = binomial);summary(site.inf.tad)
emmeans(site.inf.tad, specs = pairwise~site)

table(frog.data[which(frog.data$species == 'L. tadpole'),'site'], 
      frog.data[which(frog.data$species == 'L. tadpole'),'infection'])




#SITE PERMANOVA -- REMOVE TADPOLES
#Residuals vs site

