# Produce ordinations of community data using dbRDA.

library(tidyverse)
library(magrittr)
library(phyloseq)
library(vegan)
library(ggvegan)
source("code/Rfunctions.R")
source("code/colors.R")

# load phyloseq data
(phy <- loadPhyloseq())

# Because the ordinations are primarily for visualization and there is a fair amount of noise in the full data set,
# we will first merge the replicate samples for each Genotype / Treatment combinatiton (after applying the sequencing 
# bias correction). 
merged.phy <- phy %>% unbias(.,read.csv("output/tabs/bias.csv")) %>% 
  merge_samples(.,paste(sample_data(.)$Genotype,sample_data(.)$Treatment))
# Repair sample data
sample_data(merged.phy) %<>% .[,c('Genotype','Region','Treatment')]
sample_data(merged.phy)$Genotype <- sample_names(merged.phy) %>% strsplit(" ") %>% sapply(., `[`, 1)
sample_data(merged.phy)$Region <- sample_data(merged.phy)$Genotype %>% strsplit("-") %>% sapply(., `[`, 1)
sample_data(merged.phy)$Treatment <- sample_names(merged.phy) %>% strsplit(" ") %>% sapply(., `[`, 2)

# Calculate comunity dissimilarity using Jensen-Shannon distance (sqrt of Jensen-Shannon divergence).
# We will also first apply a Wisconsin double standardization.
otu_table(merged.phy) %<>% data.frame %>% wisconsin %>% otu_table(taxa_are_rows = F)
merged.dist <- merged.phy %>% phyloseq::distance("jsd") %>% sqrt

########################
### Treatment effect ###
########################

# fit dbRDA constrained by treatment, conditioned on genotype
(dbrda.treatment <- capscale(merged.dist~Condition(Genotype)+Treatment, data=data.frame(sample_data(merged.phy))))
# get variance explained by the ordination axes
dbrda.treatment.r2 <- (eigenvals(dbrda.treatment)/sum(eigenvals(dbrda.treatment))) %>% multiply_by(100) %>% round(1)
# plot ordination
(ordination.treatment <- fortify(dbrda.treatment,display="wa") %>%
    filter(Score=="sites") %>% 
    bind_cols(data.frame(sample_data(merged.phy))) %>%
    group_by(Treatment) %>%
    summarize(n=n(),x=mean(CAP1),y=mean(CAP2),sd1=sd(CAP1),sd2=sd(CAP2)) %>%
    ggplot(aes(x=x,y=y,fill=Treatment)) +
    geom_errorbar(aes(ymin=y-sd2,ymax=y+sd2),width=0)+
    geom_errorbarh(aes(xmin=x-sd1,xmax=x+sd1),height=0)+
    geom_point(size=5,shape=21)+
    annotate("text", x = -Inf, y = Inf, label = "(a)",hjust=-0.4, vjust=1.75, size=5)+
    scale_fill_manual(values=pal.treatment)+
    coord_fixed(eigenvals(dbrda.treatment)[2]/eigenvals(dbrda.treatment)[1])+
    labs(x=paste0("dbRDA1 [",dbrda.treatment.r2[1],"%]"),
         y=paste0("dbRDA2 [",dbrda.treatment.r2[2],"%]"))+
    ggthemes::theme_few()+
    theme(legend.justification = "left"))

#######################
### Genotype effect ###
#######################

# fit dbRDA constrained by genotype, conditioned on treatment
(dbrda.genotype <- capscale(merged.dist~Condition(Treatment)+Genotype, data=data.frame(sample_data(merged.phy))))
# get variance explained by the ordination axes
dbrda.genotype.r2 <- (eigenvals(dbrda.genotype)/sum(eigenvals(dbrda.genotype))) %>% multiply_by(100) %>% round(1)
#' plot ordination
(ordination.genotype <- fortify(dbrda.genotype,display="wa") %>%
    filter(Score=="sites") %>% 
    bind_cols(data.frame(sample_data(merged.phy))) %>%
    group_by(Genotype) %>%
    summarize(n=n(),x=mean(CAP1),y=mean(CAP2),sd1=sd(CAP1),sd2=sd(CAP2)) %>%
    left_join(sample_data(phy) %>%
                data.frame %>%
                select(Genotype,Region) %>%
                unique) %>%
    ggplot(aes(x=x,y=y,fill=Region)) +
    geom_errorbar(aes(ymin=y-sd2,ymax=y+sd2),width=0)+
    geom_errorbarh(aes(xmin=x-sd1,xmax=x+sd1),height=0)+
    geom_point(size=5,shape=21)+
    annotate("text", x = -Inf, y = Inf, label = "(b)",hjust=-0.4, vjust=1.75, size=5)+
    scale_fill_manual(values=pal.region)+
    coord_fixed(eigenvals(dbrda.genotype)[2]/eigenvals(dbrda.genotype)[1])+
    labs(x=paste0("dbRDA1 [",dbrda.genotype.r2[1],"%]"),
         y=paste0("dbRDA2 [",dbrda.genotype.r2[2],"%]"))+
    ggthemes::theme_few()+
    theme(legend.justification = "left"))

# Make combined plot for pub
gridExtra::grid.arrange(ordination.treatment,ordination.genotype,nrow=2,heights=c(1.28,1)) %>%
  ggsave("output/figs/Fig.1.pdf",.,width=18,height=14,units="cm")

#####################
### Joint effects ###
#####################

# fit dbRDA with genotype * treatment interaction. 
# The ordination will explain all variance because there is only one point for each genotype / treatment combination. 
# For visualization purposes this should be fine.
(dbrda.interaction <- capscale(merged.dist~Treatment*Genotype, data=data.frame(sample_data(merged.phy))))
#' get variance explained by the ordination axes
dbrda.interaction.r2 <- (eigenvals(dbrda.interaction)/sum(eigenvals(dbrda.interaction))) %>% multiply_by(100) %>% round(1)
#' plot ordination
(ordination.interaction <- fortify(dbrda.interaction,display="wa") %>%
    filter(Score=="sites") %>% 
    bind_cols(data.frame(sample_data(merged.phy))) %>%
    ggplot(aes(x=CAP1,y=CAP2)) +
    geom_point(size=4,aes(fill=Treatment,shape=Region))+
    scale_fill_manual(values=pal.treatment)+
    scale_shape_manual(values=c(21,23))+
    coord_fixed(eigenvals(dbrda.interaction)[2]/eigenvals(dbrda.interaction)[1])+
    labs(x=paste0("dbRDA1 [",dbrda.interaction.r2[1],"%]"),
         y=paste0("dbRDA2 [",dbrda.interaction.r2[2],"%]"))+
    guides(fill=guide_legend(order=1,override.aes=list(shape=21)))+
    ggthemes::theme_few()+
    theme(legend.justification = "left"))
ggsave("output/figs/Fig.S6.pdf",width=20,height=12.5,units="cm")



