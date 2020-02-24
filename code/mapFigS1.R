# make map showing origins of P. trichocarpa genotypes

library(ggmap)
library(ggrepel)
library(cowplot)
library(magick)

points <- read.csv("data/GenotypeOrigins.csv",header=T,as.is=T)
ptri <- png::readPNG("data/500px-Populus_trichocarpa_range_map.svg.png")

map <- make_bbox(points$Longitude,points$Latitud,f=0.35) %>%
  get_map(zoom=7, maptype = "watercolor",source="osm",force=F) %>%
  ggmap()+
  geom_point(data=points,aes(x=Longitude,y=Latitude,fill=Region,shape=Region),size=3)+
  scale_shape_manual(values=c(23,21))+
  scale_fill_manual(values = c("#fb832d","#016392"))+
  geom_label_repel(data=points,aes(x=Longitude,y=Latitude,label=Genotype),
                   point.padding = 0.5,box.padding = 0.6)+
  labs(x="Longitude",y="Latitude")

ggdraw()+
  draw_plot(map)+
  draw_image("data/500px-Populus_trichocarpa_range_map.svg.png",scale=1,width=0.32,x=0.52,y=0.24)
ggsave("output/figs/Figs.S1.pdf",width=20,height=16,units="cm")

  
