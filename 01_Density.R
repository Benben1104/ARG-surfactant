library(multcompView)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(scales)

d_CEF_DTAC <- read.csv("CEF+DTAC.csv")
head(d_CEF_DTAC)
d_CEF_DTAC$Group <- factor(d_CEF_DTAC$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_CEF_DTAC))$residuals
shapiro.test(residuals)
bartlett.test(Cell.density ~ Group, data=d_CEF_DTAC)
aov.1 <- aov(Cell.density ~ Group, data=d_CEF_DTAC)
summary(aov.1)
TukeyHSD(aov.1)

ggplot(d_CEF_DTAC, aes(Group, Cell.density), color=Group)+
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
  theme(axis.title.y = element_text(size=30, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "OrRd")+
  scale_color_brewer(palette = "OrRd")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF+DTAC_density.png", width=690/90, height=420/90, dpi=600, unit="in")



d_CEF_SDS <- read.csv("CEF+SDS.csv")
head(d_CEF_SDS)
d_CEF_SDS$Group <- factor(d_CEF_SDS$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_CEF_SDS))$residuals
shapiro.test(residuals)
bartlett.test(Cell.density ~ Group, data=d_CEF_SDS)
aov.2 <- aov(Cell.density ~ Group, data=d_CEF_SDS)
summary(aov.2)
TukeyHSD(aov.2)

ggplot(d_CEF_SDS, aes(Group, Cell.density), color=Group)+
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
  theme(axis.title.y = element_text(size=30, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Blues")+
  scale_color_brewer(palette = "Blues")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF+SDS_density.png", width=690/90, height=420/90, dpi=600, unit="in")



d_CEF_TX <- read.csv("CEF+TX-100.csv")
head(d_CEF_TX)
d_CEF_TX$Group <- factor(d_CEF_TX$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_CEF_TX))$residuals
shapiro.test(residuals)
bartlett.test(Cell.density ~ Group, data=d_CEF_TX)
aov.3 <- aov(Cell.density ~ Group, data=d_CEF_TX)
summary(aov.3)
TukeyHSD(aov.3)

ggplot(d_CEF_TX, aes(Group, Cell.density), color=Group)+
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
  theme(axis.title.y = element_text(size=30, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))


ggsave("CEF+TX-100_density.png", width=690/90, height=420/90, dpi=600, unit="in")



d_TET_DTAC <- read.csv("TET+DTAC.csv")
head(d_TET_DTAC)
d_TET_DTAC$Group <- factor(d_TET_DTAC$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_TET_DTAC))$residuals
shapiro.test(residuals)
bartlett.test(Cell.density ~ Group, data=d_TET_DTAC)
aov.4 <- aov(Cell.density ~ Group, data=d_TET_DTAC)
summary(aov.4)
TukeyHSD(aov.4)

ggplot(d_TET_DTAC, aes(Group, Cell.density), color=Group)+
  geom_boxplot(aes(fill=Group),size=1.35, width=0.75, outlier.fill="grey2",outlier.color="grey2",alpha=0.6)+
  theme_bw()+
  scale_y_log10()+
  scale_y_continuous(limits = c(20000000, 40000000))+
  geom_jitter(aes(Group, Cell.density), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.75)+
  labs(x=NULL,
       y=expression("Cell density (cells"~mL^"-1"*")"),
       title=expression("TET+DTAC"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=30, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "OrRd")+
  scale_color_brewer(palette = "OrRd")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET+DTAC_density.png", width=690/90, height=420/90, dpi=600, unit="in")



d_TET_SDS <- read.csv("TET+SDS.csv")
head(d_TET_SDS)
d_TET_SDS$Group <- factor(d_TET_SDS$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_TET_SDS))$residuals
shapiro.test(residuals)
bartlett.test(Cell.density ~ Group, data=d_TET_SDS)
aov.5 <- aov(Cell.density ~ Group, data=d_TET_SDS)
summary(aov.5)
TukeyHSD(aov.5)

ggplot(d_TET_SDS, aes(Group, Cell.density), color=Group)+
  geom_boxplot(aes(fill=Group),size=1.35, width=0.75, outlier.fill="grey2",outlier.color="grey2",alpha=0.6)+
  theme_bw()+
  scale_y_log10()+
  scale_y_continuous(limits = c(20000000, 35000000))+
  geom_jitter(aes(Group, Cell.density), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.75)+
  labs(x=NULL,
       y=expression("Cell density (cells"~mL^"-1"*")"),
       title=expression("TET+SDS"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=30, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Blues")+
  scale_color_brewer(palette = "Blues")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET+SDS_density.png", width=690/90, height=420/90, dpi=600, unit="in")



d_TET_TX <- read.csv("TET+TX-100.csv")
head(d_TET_TX)
d_TET_TX$Group <- factor(d_TET_TX$Group, level=c("Ctrl", "0.02", "0.04", "0.2", "1", "5", "10"))

residuals <- (aov(Cell.density ~ Group, data=d_TET_TX))$residuals
shapiro.test(residuals)
bartlett.test(Cell.density ~ Group, data=d_TET_TX)
aov.6 <- aov(Cell.density ~ Group, data=d_TET_TX)
summary(aov.6)
TukeyHSD(aov.6)

ggplot(d_TET_TX, aes(Group, Cell.density), color=Group)+
  geom_boxplot(aes(fill=Group),size=1.35, width=0.75, outlier.fill="grey2",outlier.color="grey2",alpha=0.6)+
  theme_bw()+
  scale_y_log10()+
  scale_y_continuous(limits = c(15000000, 35000000))+
  geom_jitter(aes(Group, Cell.density), width =0.22, 
              size = 3.75, shape = 21, stroke=0.01, fill="black", alpha=0.75)+
  labs(x=NULL,
       y=expression("Cell density (cells"~mL^"-1"*")"),
       title=expression("TET+TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0.4, angle = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.y = element_text(size=30, face = "bold", vjust = 2))+
  scale_fill_brewer(palette = "Purples")+
  scale_color_brewer(palette = "Purples")+
  theme(legend.position = "none")+
  theme(title = element_text(size = 18, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET+TX-100_density.png", width=690/90, height=420/90, dpi=600, unit="in")
