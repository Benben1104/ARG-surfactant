d3 <- read.csv("TET_membrane.csv")
head(d3)

residuals <- (aov(Membrane.permeability ~ Surfactant, data=d3))$residuals
shapiro.test(residuals)  #正态分布检验#
bartlett.test(Membrane.permeability ~ Surfactant, data=d3)  #方差齐性检验#
aov.1 <- aov(Membrane.permeability ~ Surfactant, data=d3)
summary(aov.1)
TukeyHSD(aov.1)

library(sf)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(RColorBrewer)
library(ggpubr)
library(scales)

ggplot(d3, aes(Surfactant, Membrane.permeability, fill= Surfactant))+
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
  scale_fill_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#203f99", "SDS"="#C31419", "TX-100"="#824CA2"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET_permeability.png", width=565/90, height=360/90, dpi=600, unit="in")



residuals <- (aov(Membrane.hydrophobicity ~ Surfactant, data=d3))$residuals
shapiro.test(residuals)  #正态分布检验#
bartlett.test(Membrane.hydrophobicity ~ Surfactant, data=d3)  #方差齐性检验#
aov.2 <- aov(Membrane.hydrophobicity ~ Surfactant, data=d3)
summary(aov.2)
TukeyHSD(aov.2)


ggplot(d3, aes(Surfactant, Membrane.hydrophobicity, fill= Surfactant))+
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
  scale_fill_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#203f99", "SDS"="#C31419", "TX-100"="#824CA2"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET_hydrophobicity.png", width=565/90, height=360/90, dpi=600, unit="in")



d4 <- read.csv("CEF_membrane.csv")
head(d4)

residuals <- (aov(Membrane.permeability ~ Surfactant, data=d4))$residuals
shapiro.test(residuals)  #正态分布检验#
bartlett.test(Membrane.permeability ~ Surfactant, data=d4)  #方差齐性检验#
aov.3 <- aov(Membrane.permeability ~ Surfactant, data=d4)
summary(aov.3)
TukeyHSD(aov.3)

ggplot(d4, aes(Surfactant, Membrane.permeability, fill= Surfactant))+
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
  scale_fill_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#203f99", "SDS"="#C31419", "TX-100"="#824CA2"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF_permeability.png", width=565/90, height=360/90, dpi=600, unit="in")




residuals <- (aov(Membrane.hydrophobicity ~ Surfactant, data=d4))$residuals
shapiro.test(residuals)  #正态分布检验#
bartlett.test(Membrane.hydrophobicity ~ Surfactant, data=d4)  #方差齐性检验#
aov.4 <- aov(Membrane.hydrophobicity ~ Surfactant, data=d4)
summary(aov.4)
TukeyHSD(aov.4)

ggplot(d4, aes(Surfactant, Membrane.hydrophobicity, fill= Surfactant))+
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
  scale_fill_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#203f99", "SDS"="#C31419", "TX-100"="#824CA2"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(legend.position = "none")+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF_hydrophobicity.png", width=565/90, height=360/90, dpi=600, unit="in")

