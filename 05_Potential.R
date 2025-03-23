library(RColorBrewer)
library(ggplot2)
library(ggalt)

d1 <-  read.csv("CEF_potential.csv")
head(d1)

ggplot(d1, aes(x=Time, y=Mean, color=Surfactant))+
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
  scale_color_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#203f99", "SDS"="#C31419", "TX-100"="#824CA2"),
                     breaks = c("Ctrl", "DTAC", "SDS", "TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size=20, face = "bold", vjust = -0.5))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(size=16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("CEF_potential.png", width=600/90, height=390/90, dpi=600, unit="in")



d2 <- read.csv("TET_potential.csv")
head(d2)

ggplot(d2, aes(x=Time, y=Mean, color=Surfactant))+
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
  scale_color_manual(values = c("Ctrl"="#B3D9CE", "DTAC"="#203f99", "SDS"="#C31419", "TX-100"="#824CA2"),
                     breaks = c("Ctrl", "DTAC", "SDS", "TX-100"))+
  theme(axis.text.x=element_text(size = 18, face = "bold", vjust = 0))+
  theme(axis.text.y = element_text(size = 18, face = "bold"))+
  theme(axis.title.x = element_text(size=20, face = "bold", vjust = -0.5))+
  theme(axis.title.y = element_text(size=20, face = "bold", vjust = 2))+
  theme(title = element_text(size = 16, face = "bold"))+
  theme(legend.text = element_text(size=16, face = "bold"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"))

ggsave("TET_potential.png", width=600/90, height=390/90, dpi=600, unit="in")