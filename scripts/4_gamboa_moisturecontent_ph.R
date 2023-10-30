library(ggplot2)
library(dplyr)
library(reshape2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)

data <-read.csv("panama_moisture_content_ph.csv")
head(data)

data$Site <-as.factor(data$Site)
data$Treatment <- factor(data$Treatment, 
                             levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)"))


#remake as boxplots probably

p_ph <- ggplot(data, aes(x=Treatment, ph)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = data, aes(color=Site), width = 0.1, size = 3) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 10, hjust = .7)) +
  stat_compare_means(method = "anova", label.y = 8.2, label.x = 0.8) +
  labs(title="Refuse pile pH", x="Sample Type", y="pH")
p_ph


p_mc <- ggplot(data, aes(x=Treatment, y=percent_water)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = data, aes(color=Site), width = 0.1, size = 3) +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(angle = 10, hjust = .7)) +
  stat_compare_means(method = "anova", label.y = 550, label.x = 0.8) +
  labs(title="Refuse pile moisture content", x="Sample Type", y="Moisture Content (%)")

p_mc

p <- ggarrange(p_ph, p_mc, common.legend = TRUE, legend = "bottom", labels = "AUTO",
               font.label = list(size=20))
p

ggsave("ph_mc.tiff", device = "tiff", dpi = 700)
ggsave("ph_mc.png", device = "png", dpi = 700)
ggsave("ph_mc.pdf", device = "pdf", dpi = 700)


#stats
#https://www.scribbr.com/statistics/anova-in-r/
data_stats <- data %>%
  group_by(Treatment) %>%
  summarise(count=n(),
            mean=mean(percent_water),
            sd=sd(percent_water))
head(data_stats)



data_stats2 <- data %>%
  group_by(Type) %>%
  summarise(count=n(),
            meanph = mean(ph),
            sdph = sd(ph),
            meanwc=mean(percent_water),
            sdwc=sd(percent_water))
head(data_stats2)



#save data into sep vectors
exp <- subset(data, Type == "experimental")
nat <- subset(data, Type == "natural")


ph_exp <- exp$ph
ph_nat <- nat$ph
wc_exp <- exp$percent_water
wc_nat <- nat$percent_water

ttest_ph <- t.test(ph_exp, ph_nat, var.equal = TRUE)
ttest_wc <- t.test(wc_exp, wc_nat, var.equal = TRUE)


ttest_ph
ttest_wc

one.way <- aov(data=data, percent_water ~ Treatment)
summary(one.way)

two.way <- aov(data=data, percent_water ~ Treatment + Type)
summary(two.way)
