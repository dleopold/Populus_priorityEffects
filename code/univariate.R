# Make figues showing variation in individual species across treatents as well as univarite statistical results underlying mvabund models.

library(tidyverse)
library(magrittr)
library(mvabund)
library(phyloseq)
library(ggbeeswarm)
library(foreach)
library(MASS)
library(ggtext)

source("code/Rfunctions.R")
source("code/colors.R")

# load phyloseq data
(phy <- loadPhyloseq())

##################################
### Plot proportional abundace ###
##################################

phy %>% phy_logit %>%
  otu_table %>% data.frame %>%
  bind_cols(sample_data(phy) %>%
              data.frame %>%
              dplyr::select("Treatment","Genotype","Region")) %>%
  pivot_longer(1:8,"Taxon",values_to="RelAbund") %>%
  mutate(Treatment=substr(Treatment,1,3)) %>%
  ggplot(aes(x=Treatment,y=RelAbund,fill=Region))+
  scale_fill_manual(values=pal.region,name="Host\ngenotype")+
  geom_violin(scale="width",draw_quantiles=0.5)+
  ggbeeswarm::geom_quasirandom(aes(group=Region),shape=19,alpha=0.3,size=0.5,dodge.width=1,width=0.18,na.rm=T, show.legend = F)+
  facet_wrap(~Taxon,scales="free",nrow=2)+
  labs(y="logit( proportional abundance )",
       x="Initial colonist treatment")+
  ggthemes::theme_few()+
  theme(legend.title=element_text(size=16),
        legend.text=element_text(size=14),
        legend.position = "right",
        strip.text = element_text(size=14,face="italic"),
        axis.text.x = element_text(size=12,face="italic"),
        axis.title.y = element_text(size=16, margin=margin(r = 5)),
        axis.title.x = element_text(size=16, margin=margin(t = 10)))
ggsave("output/figs/Fig.S5.pdf",width=28,height=14,units="cm")

###################################
### Univariate tests as heatmap ###
###################################

# In addition to indicating "significance" from the joint models, we will fit individual negative
# binomial glms and extract pseudo R2 values as a indicator of (relative) effect sizes.

# create effort variable
bias <- read.csv("output/tabs/bias.csv")
effort <- (sample_sums(unbias(phy, bias)) %*% t(bias$Bhat[match(taxa_names(phy),bias$Taxon)])) %>% log 
colnames(effort) <- taxa_names(phy)
effort.long <- effort %>% data.frame %>%
  mutate(SampID=sample_names(phy)) %>%
  pivot_longer(1:8,"Taxon",values_to = "effort")

# Fit individual nb.glm models and extract pseudo R2s
uni.r2 <- otu_table(phy) %>% data.frame %>%
  rownames_to_column("SampID") %>%
  mutate(Genotype=sample_data(phy)$Genotype,
         Treatment=sample_data(phy)$Treatment,
         Region=sample_data(phy)$Region) %>%
  pivot_longer(2:9,"Taxon",values_to="abund") %>%
  left_join(effort.long) %>%
  split(.$Taxon) %>%
  map_df( function(x){
    full.genotype <- glm.nb(abund~Treatment*Genotype+offset(effort),data=x)
    full.region <- glm.nb(abund~Treatment*Genotype+offset(effort),data=x)
    list(Treatment=rsq::rsq.partial(full.genotype,update(full.genotype, .~. -Treatment/Genotype),type='n',adj = T)$partial.rsq,
         Genotype=rsq::rsq.partial(full.genotype,update(full.genotype, .~. -Genotype/Treatment),type='n',adj = T)$partial.rsq,
         'Genotype:Treatment'=rsq::rsq.partial(full.genotype,update(full.genotype, .~. -Genotype:Treatment),type='n',adj = T)$partial.rsq,
         Region=rsq::rsq.partial(full.region,update(full.region, .~. -Region/Treatment),type='n',adj = T)$partial.rsq) %>% return}, .id="Taxon") %>%
  pivot_longer(., 2:5,"predictor",values_to = "partialR2") %>%
  mutate(partialR2=ifelse(partialR2<0,0,partialR2)) 

# Get univariate p-values from mvabund model 
mv.genotype <- readRDS("output/rds/mv.genotype.rds")
mv.region <- readRDS("output/rds/mv.region.rds")
uni.pvals <- mv.genotype$uni.p %>% data.frame %>% .[-1,] %>%
  rbind(mv.region$uni.p %>% data.frame %>% .[2,]) %>%
  rownames_to_column("predictor") %>%
  pivot_longer(2:9,"Taxon",values_to = "pval") %>%
  mutate(stars=gtools::stars.pval(pval) %>% gsub(".","+",.,fixed=T))

# Merge univariate data for plotting
uni.dat <- full_join(uni.r2,uni.pvals) 
uni.dat$predictor %<>% factor(levels=c("Region","Genotype","Treatment","Genotype:Treatment"))
uni.dat$Taxon %<>% factor(levels=phy %>% phy_logit %>% taxa_sums() %>% sort(decreasing = F) %>% names)

# make plot
vec_fontface <- ifelse(levels(uni.dat$Taxon) %in% c(unique(sample_data(phy)$Treatment)),"bold.italic","italic")
ggplot(uni.dat,aes(x=predictor,y=Taxon,fill=partialR2))+
  geom_tile(color="black",size=0.2)+
  scale_fill_gradientn(expression(plain("Partial ")*italic("r")^2),
                       colours = viridisLite::viridis(256, begin=0.45, end=0.95,alpha=0.98, direction=-1,option = "A"))+
  geom_text(aes(label=stars),size=6)+
  scale_x_discrete(labels = c("Region","Genotype","Treatment","Genotype x\nTreatment"))+
  #scale_y_discrete(limits = rev(levels(partialR2$taxa)))+
  ggthemes::theme_few()+
  theme(panel.grid.minor = element_line(color="black",size=1),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size=12,face=rev(vec_fontface)),
        axis.text.x =  element_text(size=12,angle=0),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))
ggsave("output/figs/Fig.2.pdf",width=18,height=20,units="cm")
