d_SDS <- read.csv("CEF+SDS.csv")
head(d_SDS)
d_SDS$Group <- factor(d_SDS$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_SDS))$residuals
shapiro.test(residuals)
bartlett.test(Cell.density ~ Group, data=d_SDS)
aov.1 <- aov(Cell.density ~ Group, data=d_SDS)
summary(aov.1)
TukeyHSD(aov.1)

library(multcompView)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(scales)

ggplot(d_SDS, aes(Group, Cell.density), color=Group)+
  geom_boxplot(aes(fill=Group),size=1.35, width=0.75, outlier.fill="grey2",outlier.color="grey2",alpha=0.6)+
  theme_bw()+
  scale_y_log10()+
  scale_y_continuous(limits = c(1200000000, 2000000000), labels = scales::scientific)+
  geom_jitter(aes(Group, Cell.density), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.75)+
  labs(x=NULL,
       y=expression("Cell density (cells"~mL^"-1"*")"),
       title=expression("CEF+SDS"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.5))+
  scale_fill_brewer(palette = "OrRd")+
  scale_color_brewer(palette = "OrRd")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))
  
ggsave("CEF+SDS_density.png", width=690/90, height=420/90, dpi=600, unit="in")



d_DTAC <- read.csv("CEF+DTAC.csv")
head(d_DTAC)
d_DTAC$Group <- factor(d_DTAC$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_DTAC))$residuals
shapiro.test(residuals)  #正态分布检验#
bartlett.test(Cell.density ~ Group, data=d_DTAC)  #方差齐性检验#
aov.2 <- aov(Cell.density ~ Group, data=d_DTAC)
summary(aov.2)
TukeyHSD(aov.2)

ggplot(d_DTAC, aes(Group, Cell.density), color=Group)+
  geom_boxplot(aes(fill=Group),size=1.35, width=0.75, outlier.fill="grey2",outlier.color="grey2",alpha=0.6)+
  theme_bw()+
  scale_y_log10()+
  scale_y_continuous(limits = c(1100000000, 1600000000), labels = scales::scientific)+
  geom_jitter(aes(Group, Cell.density), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.75)+
  labs(x=NULL,
       y=expression("Cell density (cells"~mL^"-1"*")"),
       title=expression("CEF+DTAC"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.5))+
  scale_fill_brewer(palette = "Blues")+
  scale_color_brewer(palette = "Blues")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF+DTAC_density.png", width=690/90, height=420/90, dpi=600, unit="in")



d_TX <- read.csv("CEF+TX-100.csv")
head(d_TX)
d_TX$Group <- factor(d_TX$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_TX))$residuals
shapiro.test(residuals)
bartlett.test(Cell.density ~ Group, data=d_TX)
aov.3 <- aov(Cell.density ~ Group, data=d_TX)
summary(aov.3)
TukeyHSD(aov.3)

ggplot(d_TX, aes(Group, Cell.density), color=Group)+
  geom_boxplot(aes(fill=Group),size=1.35, width=0.75, outlier.fill="grey2",outlier.color="grey2",alpha=0.6)+
  theme_bw()+
  scale_y_log10()+
  scale_y_continuous(limits = c(1300000000, 2090000000), labels = scales::scientific)+
  geom_jitter(aes(Group, Cell.density), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.75)+
  labs(x=NULL,
       y=expression("Cell density (cells"~mL^"-1"*")"),
       title=expression("CEF+TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2.5))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))


ggsave("CEF+TX-100_density.png", width=690/90, height=420/90, dpi=600, unit="in")
