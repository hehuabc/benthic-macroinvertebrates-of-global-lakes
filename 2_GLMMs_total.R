# 1. prepare work-----------------------
# set your own work place
setwd("C:/data") # set your own

# loading packages
library(dplyr)
library(plyr)
library(tidyr)
library(Rmisc)
library(ggplot2)
library(ggpubr)
library(data.table)
library(glmmTMB)
library(ggeffects)

# 2. upload data -------
data <- read.csv("diversity.csv",header = T) 

# separate warm/cool lakes
data$climate <- factor(data$region,levels = c("Tropical","Subtropical",
                                              "Mediterranean","Temperate"),
                       labels = c("Warm","Warm","Warm","Cool"))


# 3. fit Sn, SPIE and interaction ---------------------------------------------
# log10-transformed data
data2 <-  mutate(data,s = log10(S),
                Area = log10(lakearea),
                Depth = log10(waterdepth),
                TP = (log10(tp)),
                TN = (log10(tn)),
                s_n = log10(S_n),
                s_pie = log10(S_PIE),
                density = log10(N)) %>% setDT()


# fit_sn
fit_sn <- glmmTMB(s_n ~ Area + (Area|source),
                     family = gaussian(link = "identity"),
                     data2)
summary(fit_sn)

# get fitted lines from predict 
data_sn <- as.data.frame(ggpredict(fit_sn,terms = "Area"))


# fit_spie
fit_spie <- glmmTMB(s_pie ~ Area + (Area|source),
                  family = gaussian(link = "identity"),
                  data2)
summary(fit_spie)

data_spie <- as.data.frame(ggpredict(fit_spie,terms = "Area"))


# fit area * climate interaction
# sn
fit_sn_int <- glmmTMB(s_n ~ Area * climate + (Area|source),
                    family = gaussian(link = "identity"),
                    data2)
summary(fit_sn_int)

data_sn_int <- as.data.frame(ggpredict(fit_sn_int,terms = c("Area","climate")))

# Spie
fit_spie_int <- glmmTMB(s_pie ~ Area * climate + (Area|source),
                      family = gaussian(link = "identity"),
                      data2)
summary(fit_spie_int)

data_spie_int <- as.data.frame(ggpredict(fit_spie_int,
                                         terms = c("Area","climate")))


# 4. check Moran I for spatial autocorrelation ----
library(SoDA)
library(spdep)
plot.xy <- geoXY(data2$long, data2$lat,unit=1000)# unit is 1km
k1<-knn2nb(knearneigh(plot.xy,k=1))
dsts <- unlist(nbdists(k1, plot.xy))
quant95 <- as.numeric(quantile(dsts,0.95))
nbsr1<-dnearneigh(as.matrix(plot.xy),0,quant95) # calculate the neighbors distance (based on lags)
# sn
moran_sn<-moran.test(resid(fit_sn),
                     listw=nb2listw(nbsr1,zero.policy=TRUE),
                     zero.policy=TRUE)
moran_sn

# spie
moran_spie<-moran.test(resid(fit_spie),
                       listw=nb2listw(nbsr1,zero.policy=TRUE),
                       zero.policy=TRUE)
moran_spie

# spie-int
moran_int_spie <-moran.test(resid(fit_spie_int),
                          listw=nb2listw(nbsr1,zero.policy=TRUE),
                          zero.policy=TRUE)
moran_int_spie

# sn-int
moran_int_sn <-moran.test(resid(fit_sn_int),
                            listw=nb2listw(nbsr1,zero.policy=TRUE),
                            zero.policy=TRUE)
moran_int_sn


# 5. fit productivity and predation strength  ----
# TP 
## here we used data soure as a random intercept,as the AIC value is lower 
### compared with when used source as random slopes, 
### Moreover, there are convergence problems if we include source as random slopes
fit_tp_int <- glmmTMB(TP ~ Area*climate + (1|source),
                  family = gaussian(link = "identity"),
                  data2)
summary(fit_tp_int)

data_tp_int <- as.data.frame(ggpredict(fit_tp_int,terms = c("Area","climate")))


# tn
fit_tn_int <- glmmTMB(TN ~ Area * climate  + (1|source),
                  family = gaussian(link = "identity"),
                  data2)
summary(fit_tn_int)

data_tn_int <- as.data.frame(ggpredict(fit_tn_int,terms = c("Area","climate")))

# density
# here we used both TP and TN as random intercepts, as zoobenthos density was
# determined by both predation (top-down) and lake productivity (bottom-up)
fit_n_int <- glmmTMB(density ~ Area * climate + (1|TP) + (1|TN),
                 family = gaussian(link = "identity"),
                 data2)

summary(fit_n_int)

data_n_int <- as.data.frame(ggpredict(fit_n_int,terms = c("Area","climate")))

# 6. plot Fig. 2-3  ---------------
## define a theme for plots
mytheme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        axis.title = element_text(size = 8),
        panel.border=element_rect(size = 0.5),
        legend.title = element_blank(),
        legend.text = element_text(size = 8),
        strip.text = element_text(size = 8),
        aspect.ratio = 2/2)

## choose studies included lakes > 6 (n = 11)
study <- filter(data2,source == "bob"|
                  source == "caiyongjiu"|
                  source == "Chrysoula"|
                  source == "erik"|
                  source == "heikki"|
                  source == "Qianming_Dou"|
                  source == "wanghaijun"|
                  source == "zhangyou"|
                  source == "David"|
                  source == "Joseline"|
                  source == "Javier")
## plot Sn
plot_sn <-  ggplot(data= data2,aes(Area,s_n)) +
  geom_point(aes(color = climate),size = 1.5,alpha = 0.2) + 
  geom_smooth(data = study,aes(group = source,color = climate),
              method = "lm",
              se = F,
              size = 0.3)  +
  geom_line(data = data_sn, aes(x, predicted),linetype = "dashed",size = 1) +
  geom_ribbon(data = data_sn, aes(x, predicted,ymin=conf.low, ymax=conf.high), 
              alpha=0.2) +  
  mytheme +
  labs(x = "", 
       y =  expression(italic(S)[n]~(Rarefied~richness))) +
  theme(legend.position = c(0.8,0.2),
        legend.background = element_blank()) +
  scale_color_manual(values = c("darkorange1",
                                "mediumpurple")) +
  annotate("text",x = 1, y = 2, label = "Slope = 0.042 (-0.236 - 0.227)")
  
plot_sn

## plot Spie
plot_spie <-  ggplot(data= data2,aes(Area,s_pie)) +
  geom_point(aes(color = climate),size = 1.5,alpha = 0.2) + 
  geom_smooth(data = study,aes(group = source,color = climate),
              method = "lm",
              se = F,
              size = 0.3)  +
  geom_line(data = data_spie, aes(x, predicted),linetype = "dashed",size = 1) +
  geom_ribbon(data = data_spie, aes(x, predicted,ymin=conf.low, ymax=conf.high), 
              alpha=0.2) +  
  mytheme +
  labs(x = "", 
       y =  expression(italic(S)[PIE]~(evenness))) +
  theme(legend.position = "none",
        legend.background = element_blank()) +
  scale_color_manual(values = c("darkorange1",
                                "mediumpurple"))+
  annotate("text",x = 1, y = 1.5, label = "Slope = 0.026 (-0.084 - 0.191)")

plot_spie

## plot sar_interaction
plot_sar_int <-  ggplot(data = data_sn_int, aes(x, predicted))+
  geom_line(aes(color = group),linetype = "solid",size = 1) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high,fill = group), 
              alpha=0.2) +  
  mytheme +
  labs(x = expression(Area~(km^{2})), 
       y =   expression(italic(S)[n]~(Rarefied~richness))) +
  theme(legend.position = "none",
        legend.background = element_blank()) +
  scale_color_manual(values = c("darkorange1","mediumpurple"))+
  scale_fill_manual(values = c("darkorange1","mediumpurple")) +
  annotate("text",x = 1,y = 1.5,
           label = expression('Interaction':~italic(P)== '0.022'~(italic(Z)=="-2.287")))
plot_sar_int

## plot spie_interaction
plot_spie_int <-  ggplot(data = data_spie_int, aes(x, predicted))+
  geom_line(aes(color = group),size = 1,
            linetype = "dashed") +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high,fill = group), 
              alpha=0.2) +  
  mytheme +
  labs(x = expression(Area~(km^{2})), 
       y =   expression(italic(S)[PIE]~(evenness))) +
  theme(legend.position = "none",
        legend.background = element_blank()) +
  scale_color_manual(values = c("darkorange1","mediumpurple"))+
  scale_fill_manual(values = c("darkorange1","mediumpurple")) +
  annotate("text",x = 1,y = 0.92,
           label = expression('Interaction':~italic(P)== '0.099'~(italic(Z)=="-1.646")))
plot_spie_int

## Figure 2 arrange
tiff("Figure2.tiff",res = 900,compression = "lzw",
     width = 16,height = 16,units = "cm")
ggarrange(plot_sn,
          plot_spie,
          plot_sar_int,
          plot_spie_int,
          ncol = 2,nrow = 2,
          align = "hv",
          labels = "auto")
dev.off()


## plot tn,tp,density
plot_n <-  ggplot(data= data2, 
                  aes(Area,density))+ 
  geom_point(aes(color = climate),
             size = 2,alpha = 0.2)+
  geom_line(data = data_n_int, aes(x, predicted,color = group),size = 1) +
  geom_ribbon(data = data_n_int,aes(x, predicted,ymin=conf.low, ymax=conf.high,
                  fill = group), show.legend = F,
              alpha=0.2) +  
  mytheme +
  labs(x = expression(Area~(km^{2})), 
       y =  expression(Density~(ind~m^{-2}))) +
  theme(legend.position = c(0.5,0.1),
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  scale_color_manual(values = c("darkorange1",
                                "mediumpurple")) +
  scale_fill_manual(values = c("darkorange1",
                               "mediumpurple")) +
  annotate("text",x = 1,y = 5.2,
           label = expression('Interaction':~italic(P)== '0.019'~(italic(Z)=="-2.33")))
plot_n


plot_area_tp <-  ggplot(data= data2, aes(Area,TP))+ 
  geom_point(aes(color = climate),
             size = 1.5,alpha = 0.2) +
  geom_line(data = data_tp_int, aes(x, predicted,color = group),size = 1,
            linetype= "dashed") +
  geom_ribbon(data = data_tp_int,aes(x, predicted,ymin=conf.low, ymax=conf.high,
                                fill = group), show.legend = F,
              alpha=0.2) +
  mytheme +
  labs(x = expression(Area~(km^{2})), 
       y =  expression(TP~(mg~L^{-1}))) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.direction = "vertical") +
  scale_color_manual(values = c("darkorange1",
                                "mediumpurple")) +
  scale_fill_manual(values = c("darkorange1",
                                "mediumpurple")) +
  annotate("text",x = 1,y = 0.8,
           label = expression('Interaction':~italic(P)== '0.121'~(italic(Z)=="1.549")))
plot_area_tp

plot_area_tn <-  ggplot(data= data2, aes(Area,TN))+ 
  geom_point(aes(color = climate),
             size = 1.5,alpha = 0.2) +
  geom_line(data = data_tn_int, aes(x, predicted,color = group),size = 1,
            linetype= "dashed") +
  geom_ribbon(data = data_tn_int,aes(x, predicted,ymin=conf.low, ymax=conf.high,
                                 fill = group), show.legend = F,
              alpha=0.2) +
  mytheme +
  labs(x = expression(Area~(km^{2})), 
       y =  expression(TN~(mg~L^{-1}))) +
  theme(legend.position = "none",
        strip.text = element_blank(),
        legend.background = element_blank(),
        legend.direction = "vertical") +
  scale_color_manual(values = c("darkorange1",
                                "mediumpurple"))+
  scale_fill_manual(values = c("darkorange1",
                               "mediumpurple")) +
  annotate("text",x = 1,y = 1.3,
           label = expression('Interaction':~italic(P)== '0.837'~(italic(Z)=="0.205")))
plot_area_tn

# Figure 3 arrange
tiff("Figure3.tiff",res = 1200,compression = "lzw",width = 30,height = 10,
     units = "cm")
ggarrange(plot_area_tp, 
           plot_area_tn,
           plot_n, 
           align = "hv",
           ncol = 3,labels = "auto")

dev.off()



# 7. plot world map_Fig. 1 ----
require(maps)
library(maptools)
library(ggplot2)
data <- read.csv("diversity.csv",header = T)

data$region <- factor(data$region,levels = c("Tropical","Subtropical",
                                             "Mediterranean","Temperate"))

coord <- select(data,lat,long,region)

mapworld <- borders(database = "world", 
                    colour="gray70", 
                    fill = "gray70",
                    size = 0) #world map
mp_world <- ggplot() + mapworld  

# Figure 1-site map
tiff("map.tiff",
     compression = "lzw",
     res = 900,
     width = 19,
     height = 12,
     units = "cm")
mp_world + geom_point(data = coord, aes(x = long, 
                                        y = lat,
                                        color = region),
                      size = 2,alpha = 0.3,shape = 16) +
  labs(x="Longtitude (¡ãE)",y="Latitude (¡ãN)") + 
  scale_x_continuous(breaks = c(-180,-120,-60,0,60,120,180)) +
  scale_y_continuous(breaks = c(-90,-60,-30,0,30,60,90)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "white"),
        legend.position = c(0.15,0.3),
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10)) +
  scale_color_manual(values = c("red","darkorange1",
                                "steelblue","mediumpurple"))
dev.off()
