library(multcompView)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(scales)

d_CEF_DTAC <- read.csv("CEF+DTAC.csv")
head(d_CEF_DTAC)
d_CEF_DTAC$Group <- factor(d_CEF_DTAC$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_CEF_DTAC))$residuals
shapiro.test(residuals)
bartlett.test(Relative.value ~ Group, data=d_CEF_DTAC)
aov.1 <- aov(Relative.value ~ Group, data=d_CEF_DTAC)
summary(aov.1)
TukeyHSD(aov.1)

ggplot(d_CEF_DTAC, aes(Group, Relative.value, fill=Group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1.5))+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.65, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Relative ROS level"),
       title=expression("CEF+DTAC"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.25))+
  scale_fill_brewer(palette = "Blues")+
  scale_color_brewer(palette = "Blues")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF+DTAC_ROS.png", width=650/90, height=420/90, dpi=600, unit="in")



d_CEF_SDS <- read.csv("CEF+SDS.csv")
head(d_CEF_SDS)
d_CEF_SDS$Group <- factor(d_CEF_SDS$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_CEF_SDS))$residuals
shapiro.test(residuals)
bartlett.test(Relative.value ~ Group, data=d_CEF_SDS)
aov.2 <- aov(Relative.value ~ Group, data=d_CEF_SDS)
summary(aov.2)
TukeyHSD(aov.2)

ggplot(d_CEF_SDS, aes(Group, Relative.value, fill=Group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 1.5))+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.65, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Relative ROS level"),
       title=expression("CEF+SDS"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.25))+
  scale_fill_brewer(palette = "OrRd")+
  scale_color_brewer(palette = "OrRd")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF+SDS_ROS.png", width=650/90, height=420/90, dpi=600, unit="in")



d_CEF_TX <- read.csv("CEF+TX-100.csv")
head(d_CEF_TX)
d_CEF_TX$Group <- factor(d_CEF_TX$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_CEF_TX))$residuals
shapiro.test(residuals)
bartlett.test(Relative.value ~ Group, data=d_CEF_TX)
aov.3 <- aov(Relative.value ~ Group, data=d_CEF_TX)
summary(aov.3)
TukeyHSD(aov.3)

ggplot(d_CEF_TX, aes(Group, Relative.value, fill=Group))+
  theme_bw()+
  scale_y_continuous(limits = c(0, 2.5))+
  stat_summary(geom="col", fun="mean",
               position="dodge", width = 0.65, alpha=0.7)+
  stat_summary(geom="errorbar",
               fun.min=function(x) mean(x)-sd(x),
               fun.max=function(x) mean(x)+sd(x),
               width=0.15)+
  labs(x=NULL,
       y=expression("Relative ROS level"),
       title=expression("CEF+TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.25))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF+TX-100_ROS.png", width=650/90, height=420/90, dpi=600, unit="in")



d_TET_DTAC <- read.csv("TET+DTAC.csv")
head(d_TET_DTAC)
d_TET_DTAC$Group <- factor(d_TET_DTAC$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_TET_DTAC))$residuals
shapiro.test(residuals)
bartlett.test(Relative.value ~ Group, data=d_TET_DTAC)
aov.4 <- aov(Relative.value ~ Group, data=d_TET_DTAC)
summary(aov.4)
TukeyHSD(aov.4)

ggplot(d_TET_DTAC, aes(Group, Relative.value, fill=Group))+
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

ggsave("TET+DTAC_ROS.png", width=650/90, height=420/90, dpi=600, unit="in")



d_TET_SDS <- read.csv("TET+SDS.csv")
head(d_TET_SDS)
d_TET_SDS$Group <- factor(d_TET_SDS$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_TET_SDS))$residuals
shapiro.test(residuals)
bartlett.test(Relative.value ~ Group, data=d_TET_SDS)
aov.5 <- aov(Relative.value ~ Group, data=d_TET_SDS)
summary(aov.5)
TukeyHSD(aov.5)

ggplot(d_TET_SDS, aes(Group, Relative.value, fill=Group))+
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

ggsave("TET+SDS_ROS.png", width=650/90, height=420/90, dpi=600, unit="in")



d_TET_TX <- read.csv("TET+TX-100.csv")
head(d_TET_TX)
d_TET_TX$Group <- factor(d_TET_TX$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Relative.value ~ Group, data=d_TET_TX))$residuals
shapiro.test(residuals)
bartlett.test(Relative.value ~ Group, data=d_TET_TX)
aov.6 <- aov(Relative.value ~ Group, data=d_TET_TX)
summary(aov.6)
TukeyHSD(aov.6)

ggplot(d_TET_TX, aes(Group, Relative.value, fill=Group))+
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

ggsave("TET+TX-100_ROS.png", width=650/90, height=420/90, dpi=600, unit="in")