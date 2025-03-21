d_SDS <- read.csv("CEF+SDS.csv")
head(d_SDS)
d_SDS$Group <- factor(d_SDS$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Intracellular ~ Group, data=d_SDS))$residuals
shapiro.test(residuals)
bartlett.test(Intracellular ~ Group, data=d_SDS)
aov.1 <- aov(Intracellular ~ Group, data=d_SDS)
summary(aov.1)
TukeyHSD(aov.1)

library(multcompView)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

ggplot(d_SDS, aes(Group, Intracellular, color=Group))+
  geom_point(stat="summary", fun="mean", alpha=1.0, size=9)+
  geom_errorbar(stat="summary",
                fun.min=function(x) mean(x)-sd(x),
                fun.max=function(x) mean(x)+sd(x),
                width=0, size=3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(5*10^-11, 5*10^-10))+
  labs(x=NULL,
       y=expression("Intracellular CEF (ng"~cell^"-1"*")"),
       title=expression("CEF+SDS"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.5))+
  scale_fill_brewer(palette = "OrRd")+
  scale_color_brewer(palette = "OrRd")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
  
ggsave("CEF+SDS_in.png", width=690/90, height=420/90, dpi=600, unit="in")

  

d_DTAC <- read.csv("CEF+DTAC.csv")
head(d_DTAC)
d_DTAC$Group <- factor(d_DTAC$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Intracellular ~ Group, data=d_DTAC))$residuals
shapiro.test(residuals)
bartlett.test(Intracellular ~ Group, data=d_DTAC)
aov.2 <- aov(Intracellular ~ Group, data=d_DTAC)
summary(aov.2)
TukeyHSD(aov.2)

ggplot(d_DTAC, aes(Group, Intracellular, color=Group))+
  geom_point(stat="summary", fun="mean", alpha=1.0, size=9)+
  geom_errorbar(stat="summary",
                fun.min=function(x) mean(x)-sd(x),
                fun.max=function(x) mean(x)+sd(x),
                width=0, size=3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(1.5*10^-10, 6*10^-10))+
  labs(x=NULL,
       y=expression("Intracellular CEF (ng"~cell^"-1"*")"),
       title=expression("CEF+DTAC"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.5))+
  scale_fill_brewer(palette = "Blues")+
  scale_color_brewer(palette = "Blues")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF+DTAC_in.png", width=690/90, height=420/90, dpi=600, unit="in")



d_TX <- read.csv("CEF+TX-100.csv")
head(d_TX)
d_TX$Group <- factor(d_TX$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Intracellular ~ Group, data=d_TX))$residuals
shapiro.test(residuals)
bartlett.test(Intracellular ~ Group, data=d_TX)
aov.3 <- aov(Intracellular ~ Group, data=d_TX)
summary(aov.3)
TukeyHSD(aov.3)

ggplot(d_TX, aes(Group, Intracellular, color=Group))+
  geom_point(stat="summary", fun="mean", alpha=1.0, size=9)+
  geom_errorbar(stat="summary",
                fun.min=function(x) mean(x)-sd(x),
                fun.max=function(x) mean(x)+sd(x),
                width=0, size=3.25)+
  theme_bw()+
  scale_y_continuous(limits = c(5*10^-11, 5*10^-10))+
  labs(x=NULL,
       y=expression("Intracellular CEF (ng"~cell^"-1"*")"),
       title=expression("CEF+TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.5))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF+TX-100_in.png", width=690/90, height=420/90, dpi=600, unit="in")

