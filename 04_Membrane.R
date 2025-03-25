library(RColorBrewer)
library(ggplot2)
library(ggalt)

d_CEF_1 <-  read.csv("CEF_potential.csv")
head(d_CEF_1)

ggplot(d_CEF_1, aes(x=Time, y=Mean, color=Surfactant))+
  theme_bw()+
  scale_y_continuous(limits = c(0.5, 1.3))+
  geom_point(size=4, alpha=0.8)+
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width=0)+
  geom_xspline(spline_shape = 1.0)+
  theme(strip.text.x = element_text(size=16))+
  labs(x=expression("Time (min)"),
       y=expression("Fold change of potential"),
       title = expression("CEF"),
       color=NULL)+
  scale_color_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#C24E44", "SDS"="#2A5196", "TX-100"="#7367BE"),
                     breaks = c("Ctrl", "DTAC", "SDS", "TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size=20, face = "bold", vjust = -0.5))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(size=16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF_potential.png", width=600/90, height=390/90, dpi=600, unit="in")



d_TET_1 <- read.csv("TET_potential.csv")
head(d_TET_1)

ggplot(d_TET_1, aes(x=Time, y=Mean, color=Surfactant))+
  theme_bw()+
  scale_y_continuous(limits = c(0.5, 1.3))+
  geom_point(size=4, alpha=0.8)+
  geom_errorbar(aes(ymin = Mean-SD, ymax = Mean+SD), width=0)+
  geom_xspline(spline_shape = 1.0)+
  theme(strip.text.x = element_text(size=16))+
  labs(x=expression("Time (min)"),
       y=expression("Fold change of potential"),
       title = expression("TET"),
       color=NULL)+
  scale_color_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#C24E44", "SDS"="#2A5196", "TX-100"="#7367BE"),
                     breaks = c("Ctrl", "DTAC", "SDS", "TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size=20, face = "bold", vjust = -0.5))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(size=16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET_potential.png", width=600/90, height=390/90, dpi=600, unit="in")



library(sf)
library(tidyverse)
library(ggpubr)
library(scales)

d_CEF_2 <- read.csv("CEF_membrane.csv")
head(d_CEF_2)

residuals <- (aov(Membrane.permeability ~ Surfactant, data=d_CEF_2))$residuals
shapiro.test(residuals)
bartlett.test(Membrane.permeability ~ Surfactant, data=d_CEF_2)
aov.1 <- aov(Membrane.permeability ~ Surfactant, data=d_CEF_2)
summary(aov.1)
TukeyHSD(aov.1)

ggplot(d_CEF_2, aes(Surfactant, Membrane.permeability, fill= Surfactant))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1.5))+
  geom_jitter(aes(Surfactant, Membrane.permeability), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.95)+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.625, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Fold change of permeability"),
       title=expression("CEF"))+
  scale_fill_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#C24E44", "SDS"="#2A5196", "TX-100"="#7367BE"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF_permeability.png", width=565/90, height=360/90, dpi=600, unit="in")



residuals <- (aov(Membrane.hydrophobicity ~ Surfactant, data=d_CEF_2))$residuals
shapiro.test(residuals)
bartlett.test(Membrane.hydrophobicity ~ Surfactant, data=d_CEF_2)
aov.2 <- aov(Membrane.hydrophobicity ~ Surfactant, data=d_CEF_2)
summary(aov.2)
TukeyHSD(aov.2)

ggplot(d_CEF_2, aes(Surfactant, Membrane.hydrophobicity, fill= Surfactant))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1.25))+
  geom_jitter(aes(Surfactant, Membrane.hydrophobicity), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.95)+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.625, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Fold change of hydrophobicity"),
       title=expression("CEF"))+
  scale_fill_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#C24E44", "SDS"="#2A5196", "TX-100"="#7367BE"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF_hydrophobicity.png", width=565/90, height=360/90, dpi=600, unit="in")



d_TET_2 <- read.csv("TET_membrane.csv")
head(d_TET_2)

residuals <- (aov(Membrane.permeability ~ Surfactant, data=d_TET_2))$residuals
shapiro.test(residuals)
bartlett.test(Membrane.permeability ~ Surfactant, data=d_TET_2)
aov.3 <- aov(Membrane.permeability ~ Surfactant, data=d_TET_2)
summary(aov.3)
TukeyHSD(aov.3)

ggplot(d_TET_2, aes(Surfactant, Membrane.permeability, fill= Surfactant))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 3), labels = number_format(accuracy = 0.1))+
  geom_jitter(aes(Surfactant, Membrane.permeability), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.95)+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.625, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Fold change of permeability"),
       title=expression("TET"))+
  scale_fill_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#C24E44", "SDS"="#2A5196", "TX-100"="#7367BE"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET_permeability.png", width=565/90, height=360/90, dpi=600, unit="in")



residuals <- (aov(Membrane.hydrophobicity ~ Surfactant, data=d_TET_2))$residuals
shapiro.test(residuals)
bartlett.test(Membrane.hydrophobicity ~ Surfactant, data=d_TET_2)
aov.4 <- aov(Membrane.hydrophobicity ~ Surfactant, data=d_TET_2)
summary(aov.4)
TukeyHSD(aov.4)

ggplot(d_TET_2, aes(Surfactant, Membrane.hydrophobicity, fill= Surfactant))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1.25))+
  geom_jitter(aes(Surfactant, Membrane.hydrophobicity), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.95)+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.625, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Fold change of hydrophobicity"),
       title=expression("TET"))+
  scale_fill_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#C24E44", "SDS"="#2A5196", "TX-100"="#7367BE"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET_hydrophobicity.png", width=565/90, height=360/90, dpi=600, unit="in")