# Effects of genotype x arrival order interactions on rust severity
# Analyses conducte seperately for resistent and susceptible genotypes (likely representing varaition in major gene resistance)

library(tidyverse)
library(magrittr)
library(betareg)
library(emmeans)
library(lmtest)
library(ggthemes)
library(patchwork)
library(ggtext)

source("code/Rfunctions.R")
source("code/colors.R")

# load rust data and subset to inoculated plants
rust <- loadRust() %>% filter(Treatment!="Control") 

# subset susceptible and resistant genotypes
susceptible <- read.csv("output/tabs/susceptibility.csv", stringsAsFactors = F) 
susceptible$Genotype %<>% factor(.,levels=.[order(susceptible$pctLesion)])

rust.s <-  rust %>%
  filter(Genotype %in% filter(susceptible,Cluster=="susceptible")$Genotype) %>%
  mutate(Genotype=factor(Genotype, levels=levels(susceptible$Genotype)),
         weights=nleaf*length(nleaf)/sum(nleaf)) #proportionality weight for betagreg model

rust.r <-  rust %>%
  filter(Genotype %in% filter(susceptible,Cluster=="resistant")$Genotype) %>%
  mutate(Genotype=factor(Genotype, levels=levels(susceptible$Genotype)),
         weights=nleaf*length(nleaf)/sum(nleaf)) #proportionality weight for betagreg model

#############################
### Betaregression models ###
#############################

# Susceptible genotypes
sink("output/tabs/betareg.susceptible.txt"); print("Susceptible genotypes")
lesion.s <- betareg(pctLesion ~ Genotype*Treatment, weights=weights, data=rust.s)
lesion.s.noInt <- betareg(pctLesion ~ Genotype+Treatment, weights=weights, data=rust.s)
lrtest(lesion.s, lesion.s.noInt) #test interation
lrtest(lesion.s.noInt, .~. -Treatment) #test treatment 
lrtest(lesion.s.noInt, .~. -Genotype) #test genotype

uridinia.s <- betareg(pctRust ~ Genotype*Treatment, weights=weights, data=rust.s)
uridinia.s.noInt <- betareg(pctRust ~ Genotype+Treatment, weights=weights, data=rust.s)
lrtest(uridinia.s, uridinia.s.noInt) #test interation
lrtest(uridinia.s.noInt, .~. -Treatment) #test treatment 
lrtest(uridinia.s.noInt, .~. -Genotype) #test genotype
sink()

# Resistant genotypes - precision parameter modeled as function of genotype to account for heteroskedasticity
sink("output/tabs/betareg.resistant.txt"); print("Resistant genotypes")
lesion.r <- betareg(pctLesion ~ Genotype*Treatment|Genotype, weights=weights, data=rust.r)
lesion.r.noInt <- betareg(pctLesion ~ Genotype+Treatment|Genotype, weights=weights, data=rust.r)
lrtest(lesion.r, lesion.r.noInt) #test interation
lrtest(lesion.r.noInt, .~. -Treatment) #test treatment 
lrtest(lesion.r.noInt, .~. -Genotype) #test genotype

uridinia.r <- betareg(pctRust ~ Genotype*Treatment|Genotype, weights=weights, data=rust.r)
uridinia.r.noInt <- betareg(pctRust ~ Genotype+Treatment|Genotype, weights=weights, data=rust.r)
lrtest(uridinia.r, uridinia.r.noInt) #test interation
lrtest(uridinia.r.noInt, .~. -Treatment) #test treatment 
lrtest(uridinia.r.noInt, .~. -Genotype) #test genotype
sink()

####################
### Make figures ###
####################

dat.fig <- emmeans(lesion.s, ~Treatment|Genotype) %>% data.frame %>%
  mutate(group="Susceptible",
         Genotype=factor(Genotype,levels=levels(susceptible$Genotype)),
         upper.CL=asymp.UCL,
         lower.CL=asymp.LCL) %>%
  bind_rows(
    emmeans(lesion.r, ~Treatment|Genotype) %>% data.frame %>%
      mutate(group="Resistant",
             Genotype=factor(Genotype,levels=levels(susceptible$Genotype)),
             upper.CL=asymp.UCL,
             lower.CL=asymp.LCL)
  )

ggplot(dat.fig,aes(x=Treatment,y=emmean,fill=Treatment)) +
  geom_hline(data=filter(susceptible,Cluster=="susceptible"),
             aes(yintercept=pctLesion),linetype="dotted")+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0)+
  geom_point(shape=21,size=4)+
  geom_point(data=rust.s,
             aes(y=pctLesion,alpha=nleaf),shape=21,size=2,
             show.legend = F)+
  scale_alpha(range=c(0.4,1))+
  labs(y="Disease severity (% rust lesion)")+
  scale_fill_manual(values=pal.treatment)+
  scale_y_continuous(labels = function(x) x*100)+
  facet_wrap(~group+Genotype,nrow=2)+
  theme_few()+
  theme(strip.background = element_blank(), 
        axis.text.x = element_blank(),
        strip.text = element_text(size=11),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.text = element_text(face="italic"))




# Susceptible genotypes
(fig.a <- emmeans(lesion.s, ~Treatment|Genotype) %>% data.frame %>%
    mutate(Genotype=factor(Genotype,levels=levels(susceptible$Genotype)),
           upper.CL=asymp.UCL,
           lower.CL=asymp.LCL) %>%
    ggplot(aes(x=Treatment,y=emmean,fill=Treatment)) +
    geom_hline(data=filter(susceptible,Cluster=="susceptible"),
               aes(yintercept=pctLesion),linetype="dotted")+
    geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0)+
    geom_point(shape=21,size=3)+
    geom_point(data=rust.s,
               aes(y=pctLesion,alpha=nleaf),shape=21,size=1.5,
               show.legend = F)+
    scale_alpha(range=c(0.4,1))+
    labs(y="Susceptible genotypes<br>leaf rust (% lesion)")+
    scale_fill_manual("Initial colonist:",values=pal.treatment)+
    scale_y_continuous(labels = function(x) x*100)+
    facet_wrap(~Genotype,nrow=1)+
    theme_few()+
    theme(strip.background = element_blank(), 
          axis.text.x = element_blank(),
          strip.text = element_text(size=11),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_markdown(),
          legend.text = element_text(size=12,face="italic")))

# Resistant genotypes
(fig.b <- emmeans(lesion.r.noInt, ~Treatment|Genotype) %>% data.frame %>%
  mutate(Genotype=factor(Genotype,levels=levels(susceptible$Genotype)),
         upper.CL=asymp.UCL,
         lower.CL=asymp.LCL) %>%
  ggplot(aes(x=Treatment,y=emmean,fill=Treatment)) +
  geom_hline(data=filter(susceptible,Cluster=="resistant"),
             aes(yintercept=pctLesion),linetype="dotted")+
  geom_errorbar(aes(ymin=lower.CL,ymax=upper.CL),width=0)+
  geom_point(shape=21,size=3)+
  geom_point(data=rust.r,
             aes(y=pctLesion,alpha=nleaf),shape=21,size=1.5,
             show.legend = F)+
  scale_alpha(range=c(0.4,1))+
  labs(y="Resistant genotypes<br>leaf rust (% lesion)")+
  scale_y_continuous(limits=c(0,0.4),labels = function(x) x*100)+
  scale_fill_manual("Initial colonist:",values=pal.treatment)+
  facet_wrap(~Genotype,nrow=1)+
  theme_few()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=11),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(size=12),
        legend.text = element_text(size=12,face="italic")))

# Combine into multipanel figure
fig.a + fig.b +
  plot_layout(ncol=1,guides = 'collect') +
  plot_annotation(tag_levels = 'a',tag_prefix = "(",tag_suffix = ")") &
  theme(plot.tag = element_text(size=12,face="bold"),
        legend.position = 'bottom')
ggsave("output/figs/Fig.3.pdf",width=25,height=18,units="cm")
ggsave("MS/figs/Fig.3.jpg",width=25,height=18,units="cm") 
