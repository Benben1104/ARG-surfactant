d_SDS <- read.csv("E. coli_TET+SDS.csv")
head(d_SDS)
d_SDS$Group <- factor(d_SDS$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_SDS))$residuals
shapiro.test(residuals)  #正态分布检验#
bartlett.test(Relative.value ~ Group, data=d_SDS)  #方差齐性检验#
aov.1 <- aov(Relative.value ~ Group, data=d_SDS)
summary(aov.1)
TukeyHSD(aov.1)

library(multcompView)
library(dplyr)
df_p1 <- (TukeyHSD(aov.1)$Group)[,4]
let1 <- multcompLetters(df_p1, compare="<", threshold=0.05, Letters=letters, reversed = F)
let1

library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(scales)
ggplot(d_SDS, aes(Group, Relative.value, fill=Group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 5), labels = number_format(accuracy = 0.1))+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.65, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
      y=expression("Relative ROS level"),
      title=expression("TET+SDS"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.25))+
  scale_fill_brewer(palette = "OrRd")+
  scale_color_brewer(palette = "OrRd")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("ROS_TET+SDS.png", width=650/90, height=420/90, dpi=600, unit="in")

  


 

d_DTAC <- read.csv("E. coli_TET+DTAC.csv")
head(d_DTAC)
d_DTAC$Group <- factor(d_DTAC$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_DTAC))$residuals
shapiro.test(residuals)  #正态分布检验#
bartlett.test(Relative.value ~ Group, data=d_DTAC)  #方差齐性检验#
aov.2 <- aov(Relative.value ~ Group, data=d_DTAC)
summary(aov.2)
TukeyHSD(aov.2)

df_p2 <- (TukeyHSD(aov.2)$Group)[,4]
let2 <- multcompLetters(df_p2, compare="<", threshold=0.05, Letters=letters, reversed = F)
let2

ggplot(d_DTAC, aes(Group, Relative.value, fill=Group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 2.5), labels = number_format(accuracy = 0.1))+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.65, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Relative ROS level"),
       title=expression("TET+DTAC"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.25))+
  scale_fill_brewer(palette = "Blues")+
  scale_color_brewer(palette = "Blues")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("ROS_TET+DTAC.png", width=650/90, height=420/90, dpi=600, unit="in")





d_TX <- read.csv("E. coli_TET+TX-100.csv")
head(d_TX)
d_TX$Group <- factor(d_TX$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_TX))$residuals
shapiro.test(residuals)  #正态分布检验#
bartlett.test(Relative.value ~ Group, data=d_TX)  #方差齐性检验#
aov.3 <- aov(Relative.value ~ Group, data=d_TX)
summary(aov.3)
TukeyHSD(aov.3)

df_p3 <- (TukeyHSD(aov.3)$Group)[,4]
let3 <- multcompLetters(df_p3, compare="<", threshold=0.05, Letters=letters, reversed = F)
let3

ggplot(d_TX, aes(Group, Relative.value, fill=Group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 8), labels = number_format(accuracy = 0.1))+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.65, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Relative ROS level"),
       title=expression("TET+TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.25))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("ROS_TET+TX-100.png", width=650/90, height=420/90, dpi=600, unit="in")

