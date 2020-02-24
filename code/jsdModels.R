library(tidyverse)
library(magrittr)
library(phyloseq)
library(mvabund)
source("code/Rfunctions.R")

#set seed for reproducability
set.seed(32453)

# load phyloseq data
(phy <- loadPhyloseq())

##########################
## Heterogeneous effort ##
##########################

# load bias correction factors - estimated from mock communities
bias <- read.csv("output/tabs/bias.csv")

# There is unequal effort for both species (sequencing bias, rows) and samples (sampling depth, columns). 
# To account for both we will create a matrix, multipling the two and then take the log (because the log link 
# function in the negative binomial glm we will use below). Using this as an offset will be equal to modeling 
# the proportional abundance after correcting for sequencing bias while retaining the error structure of the 
# count data. For sequencing depth we will use the sample sums after correcting for the sequencing bias for 
# each species so that the effective proportional abundances reflect the bias corrections.
effort <- (sample_sums(unbias(phy, bias)) %*% t(bias$Bhat[match(taxa_names(phy),bias$Taxon)])) %>% log 

##########################
### Fit genotype model ###
##########################

# Make mvabund object 
mvDat <- otu_table(phy) %>% data.frame  %>% mvabund

# Fit joint-species model for genotype effect
mv.full <- manyglm(mvDat ~ Genotype*Treatment, 
                   offset=effort, 
                   family="negative.binomial",
                   data=data.frame(sample_data(phy)))

# Check model assumptions
#plot(mv.full)
#meanvar.plot(mvDat~sample_data(phy)$Treatment)
#meanvar.plot(mvDat~sample_data(phy)$Genotype)

# Test with anova.manyglm 
# Using unstructured correlation matrix and wald tests.  
# Including univariate test with adjustment for multiple testing.  
mv.anova <- anova(mv.full, nBoot=4999, p.uni="adjusted", cor.type="R", test="wald")

# Save results 
mv.anova$table %>% write.csv("output/tabs/mv.genotype.csv")
saveRDS(mv.anova, "output/rds/mv.genotype.rds")

########################
### Fit region model ###
########################

#' # Fit joint-species model for genotype effect
mv.region <- manyglm(mvDat ~ Region*Treatment, 
                     offset=effort, 
                     family="negative.binomial",
                     data=data.frame(sample_data(phy)))

#' ## Check model assumptions
plot(mv.region)
meanvar.plot(mvDat~sample_data(phy)$Region)

#' ## Test with anova.manyglm 
#+ cache=T, results='asis'
mv.region.anova <- anova(mv.region, nBoot=4999, p.uni="adjusted", cor.type="R", test="wald")

mv.region.anova$table %>% write.csv("output/tabs/mv.region.csv")
saveRDS(mv.region.anova, "output/rds/mv.region.rds")

