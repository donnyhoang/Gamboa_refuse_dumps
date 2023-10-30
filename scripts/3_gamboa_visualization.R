library(vegan)
library(ggplot2)
library(stringr)
#library(plyr) #was only needed for calcing div. messes with dplyr later
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(usedist)
library(hillR)
library(cowplot)



#references
#https://rpubs.com/an-bui/vegan-cheat-sheet
#https://sites.google.com/site/mb3gustame/home
#https://www.flutterbys.com.au/stats/tut/tut13.2.html

metadata <- read.csv("panama_phyloseq_metadata_genus.csv", header=TRUE)
otu <- read.csv("panama_phyloseq_otu_genus.csv", header=TRUE, check.names = FALSE)

colnames(metadata)[1] = "Sample"

metadata$Site<- as.character(metadata$Site)
metadata$Replicate<- as.character(metadata$Replicate)
metadata$Subreplicate<- as.character(metadata$Subreplicate)

metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "top","Top (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "bottom","Bottom (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil(?!_)","Soil (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil_X","Soil only"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil_net","Soil & Net"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "X_net","Net only"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "X_X","No additions"))

metadata$Treatment <- factor(metadata$Treatment, 
                         levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)","Soil (environmental)"))


#attach otu to metadata
data <- cbind(metadata,otu)
#drop columns
data <-data[-c(10)]

data$Treatment <- factor(data$Treatment, 
                         levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)","Soil (environmental)"))

###ALPHA DIVERSITY###


## Go with Hill numbers ###


comm <- data[,c(1,10:1056)]
row.names(comm) <- comm[,1]
comm <- comm[,-c(1)]

#h0 is richness
#h1 is exponent of shannon entropy -> "number of common ASVs"
#h2 is inverse simpson

h0 <- as.data.frame(hill_taxa(comm, q = 0, MARGIN = 1, base = exp(1)))
h0$Sample <- row.names(h0)
h1 <- as.data.frame(hill_taxa(comm, q = 1, MARGIN = 1, base = exp(1)))
h1$Sample <- row.names(h1)
h2 <- as.data.frame(hill_taxa(comm, q = 2, MARGIN = 1, base = exp(1)))
h2$Sample <- row.names(h2)


hill_total <- merge(h0, h1, by = "Sample")
hill_total <- merge(hill_total, h2, by = "Sample")
colnames(hill_total) <- c("Sample", "Hill0","Hill1", "Hill2")

hill_total <- merge(hill_total, metadata, by = "Sample")

#write.csv(hill_total, "panama_hill_diversity_metadata.csv", row.names = FALSE)

hill_summary <- hill_total %>%
  group_by(Treatment) %>%
  summarise(meanHill1 = mean(Hill1),
            sdHill1 = sd(Hill1),
            minHill1 = min(Hill1),
            maxHill1 = max(Hill1),
            meanHill2 = mean(Hill2),
            sdHill2 = sd(Hill2),
            meanHill0 = mean(Hill0),
            sdHill0 = sd(Hill0),
            minHill0 = min(Hill0),
            maxHill0 = max(Hill0),)
stats_data0 <- function(y) {
  return(data.frame(
    y=350,
    label = paste('n=', length(y), '\n')
  ))
}

h0 <- ggplot(hill_total, aes(x = Treatment, y = Hill0, fill = Treatment)) +
  geom_boxplot() +
  stat_summary(fun.data = stats_data0,
               geom = "text") +
  theme_classic(base_size=15) +
  theme(#aspect.ratio=1,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(method="anova", label.x = 1, label.y = 370) +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")) + 
  labs(y="Number of ASVs") + 
  ggtitle("Hill 0 Diversity \n Richness")
h0

stats_data1 <- function(y) {
  return(data.frame(
    y=130,
    label = paste('n=', length(y), '\n')
  ))
}

h1 <- ggplot(hill_total, aes(x = Treatment, y = Hill1, fill = Treatment)) +
  geom_boxplot() +
  stat_summary(fun.data = stats_data1,
               geom = "text") +
  theme_classic(base_size=15) +
  theme(#aspect.ratio=1,
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  stat_compare_means(method="anova", label.x = 1, label.y = 140) +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")) + 
  labs(y="Number of Common ASVs") + 
  ggtitle("Hill 1 Diversity \n Exponential of Shannon's Entropy")
h1

stats_data2 <- function(y) {
  return(data.frame(
    y=85,
    label = paste('n=', length(y), '\n')
  ))
}


h2 <- ggplot(hill_total, aes(x = Treatment, y = Hill2, fill = Treatment)) +
  geom_boxplot() +
  stat_summary(fun.data = stats_data2,
               geom = "text") +
  theme_classic(base_size=15) +
  theme(#aspect.ratio=1,
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
  stat_compare_means(method="anova", label.x = 1, label.y = 95) +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")) + 
  labs(y="Number of Very Common ASVs") + 
  ggtitle("Hill 2 Diversity \n Inverse-Simpson Index")
h2


p1 <- h0 + theme(legend.position = "none")
p2 <- h1 + theme(legend.position = "none")
p3 <- h2 + theme(legend.position = "none")
p4 <- get_legend(h0)

philltemp <- plot_grid(p1, p2, p3, p4, rel_widths = c(1,1,1,1), labels = c("A","B","C"), label_size = 20)
philltemp


ggsave("hill_diversity.tiff", device=tiff, dpi = 700)
ggsave("hill_diversity.png", device=png, dpi = 700)
ggsave("hill_diversity.pdf", device=pdf, dpi = 700)


######### Ordination


#https://rstudio-pubs-static.s3.amazonaws.com/694016_e2d53d65858d4a1985616fa3855d237f.html
#http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization
#https://chrischizinski.github.io/rstats/adonis/


metadata <- read.csv("panama_phyloseq_metadata_genus.csv", header=TRUE)
otu <- read.csv("panama_phyloseq_otu_genus.csv", header=TRUE, row.names=1, check.names = FALSE)

colnames(metadata)[1] = "Sample"

metadata$Site<- as.character(metadata$Site)
metadata$Replicate<- as.character(metadata$Replicate)
metadata$Subreplicate<- as.character(metadata$Subreplicate)

metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "top","Top (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "bottom","Bottom (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil(?!_)","Soil (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil_X","Soil only"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil_net","Soil & Net"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "X_net","Net only"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "X_X","No additions"))

metadata$Treatment <- factor(metadata$Treatment, 
                             levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)","Soil (environmental)"))





str(otu)
attach(metadata)



#ordinate 
set.seed(2)
ordination.model <- metaMDS(otu, distance='bray', k=2)

#extract coordinates so I can plot with ggplot2
nmds_coords <- ordination.model$points

#attach to metadata so I only have to refer to one object
data <- cbind(metadata, nmds_coords)

#write.csv(data, "metadata_vegan_NMDS_coords.csv", row.names = FALSE)

#stress 0.1351673 ... Procrustes: rmse 0.008222721  max resid 0.08123369 

############## ordination stats ##############

#https://chrischizinski.github.io/rstats/adonis/

dist <- vegdist(otu) #create distance matrix
attach(metadata) #attach metadata


#anosim #close to 1 means high dissimilarity
anoTreatment <- anosim(dist, Treatment) #similarity analysis
anoSite <- anosim(dist, Site)
anoSoil <-anosim(dist, Soil)
anoNet <- anosim(dist, Net)
anoType <- anosim(dist,Type)

summary(anoTreatment) # see summary
summary(anoSite)
summary(anoSoil)
summary(anoNet)
summary(anoType)



#adonis #how much of variation explained by particular variable
adTreatment <- adonis2(dist ~ Treatment, data=metadata, permutations=99, method="bray")
adSite <- adonis2(dist ~ Site, data=metadata,  permutations=99, method="bray")
adSoil <- adonis2(dist ~ Soil, data=metadata,  permutations=99, method="bray")
adNet <- adonis2(dist ~ Net, data=metadata,  permutations=99, method="bray")
adType <- adonis2(dist ~ Type, data=metadata,  permutations=99, method="bray")

adTreatment
adSite
adSoil
adNet
adType



Variable <- c("Treatment", "Site", "Soil", "Net","Experimental vs. Environmental")
ANOSIM <- c(0.6104, 0.2012, 0.1735, 0.01153, 0.7635)
ano_p <- c(0.001,0.001,0.001,0.344,0.001)
PERMANOVA <- c(0.4321,0.12363,0.07747,0.07144,0.18856)
permanova_p <- c(0.01,0.01,0.01,0.01,0.01)

ord_stat_table <- data.frame(Variable, ANOSIM, ano_p, PERMANOVA, permanova_p)
colnames(ord_stat_table) <- c("Variable", "ANOSIM", "ANOSIM \n p-value", "PERMANOVA", "PERMANOVA \n p-value")

p_stat_table <- ggtexttable(ord_stat_table, rows = NULL,
                            theme = ttheme("light"))
p_stat_table



# "#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d"
# "#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf"
# "#8dd3c7","#ffffb3","#bebada","#80b1d3","#fdb462","#b3de69","#fccde5"

data <- read.csv("metadata_vegan_NMDS_coords.csv", header = TRUE)
data$Site<- as.character(data$Site)
data$Replicate<- as.character(data$Replicate)
data$Subreplicate<- as.character(data$Subreplicate)
data$Treatment <- factor(data$Treatment, 
                         levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)","Soil (environmental)"))



plot <- ggplot() +
  geom_point(data=data, aes(x=MDS1, y=MDS2, color=Treatment, shape=Site),
             size=3) +
  annotate(geom="text", x = -1.25, y = 1.25, label = "Stress = 0.135", size = 6) +
  theme_classic(base_size = 20) +
  theme(aspect.ratio=1) +
  xlab("NMDS1") +
  ylab("NMDS2") +
  ggtitle("Biplot of Bray-Curtis Distance Matrix") +
  scale_color_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")) +
  scale_shape_manual(values=c(16, 17, 15))
plot


ordp1 <- plot
ordp2 <- p_stat_table 


ordplot <- plot_grid(ordp1, ordp2, ncol = 1, rel_heights = c(3,1), align = "h", labels = "AUTO", label_size = 20)
ordplot


ggsave("ordination.tiff", device = "tiff", dpi = 700)
ggsave("ordination.png", device = "png", dpi = 700)
ggsave("ordination.pdf", device = "pdf", dpi = 700)





######Abundance Bar graph######

metadata <- read.csv("panama_phyloseq_metadata_class.csv", header=TRUE)
otu <- read.csv("panama_phyloseq_otu_class.csv", header=TRUE, check.names = FALSE)



#order samples with separate file that I worked with in excel
#order <-read.csv("disturbance_total_order_top15.csv", header=TRUE)
#metadata$Sample_Name <- factor(metadata$Sample_Name, levels=unique(order$Sample_Name))

metadata$Site<- as.character(metadata$Site)
metadata$Replicate<- as.character(metadata$Replicate)
metadata$Subreplicate<- as.character(metadata$Subreplicate)


#metadata <- metadata %>%  mutate(Treatment = str_replace(Treatment, "top","Top"))
#metadata <- metadata %>% mutate(Treatment = str_replace(Treatment, "bottom","Bottom"))

metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "top","Top (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "bottom","Bottom (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil(?!_)","Soil (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil_X","Soil only"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil_net","Soil & Net"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "X_net","Net only"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "X_X","No additions"))

metadata$Treatment <- factor(metadata$Treatment, 
                             levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)","Soil (environmental)"))
colnames(metadata)[1] = "Sample"


#attach otu to metadata
data <- cbind(metadata,otu)
#drop columns
data <-data[-c(10)]


data <- melt(data,id=c("Sample","Site","Type","Initial","Treatment","Soil","Net",
                       "Replicate","Subreplicate"))
head(data)


### let's get top...11 list. Since Set3 is up to 12 and I need 1 for "other" #####


top <- data %>%
  group_by(variable) %>%
  summarise(sum = sum(value))
top <- top %>% arrange(desc(sum))
top$variable
top_list <- c("Alphaproteobacteria","Gammaproteobacteria","Bacteroidia",
              "Actinobacteria","Clostridia","Planctomycetes",
              "Polyangia","Verrucomicrobiae","Bacilli",
              "Vicinamibacteria","Acidobacteriae")
#for silva v128 
#top_list <- c("Alphaproteobacteria","Gammaproteobacteria","Actinobacteria","Clostridia","Betaproteobacteria","Deltaproteobacteria","Sphingobacteriia","Planctomycetacia","Flavobacteriia","Bacilli","Subgroup_6")

#get top microbes
data_filt <- filter(data, variable %in% top_list)
#get other microbes to condense into "other"
other_filt <- filter(data, !variable %in% top_list)
other_filt <- other_filt %>%
  group_by(Sample, Site, Type, Initial, Treatment, Soil, Net, Replicate, Subreplicate) %>%
  summarise(value=sum(value))
other_filt$variable <- "Other"

data_table <- rbind(data_filt, other_filt)
data_table$value <- (data_table$value)*100


#calc avg abundane of "other"
other_filt <- filter(data, !variable %in% top_list)
other_filt_summary <- other_filt %>%
  group_by(variable) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            mode=mode(value),
            median=median(value))



######## summaries for easy reference


data_table_expsum_temp <- subset(data_table, Type == "experimental")
data_table_expsum <- data_table_expsum_temp %>%
  group_by(Treatment, variable) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            min=min(value),
            max=max(value))
data_table_expsum2 <- data_table_expsum_temp %>%
  group_by(variable) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            min=min(value),
            max=max(value))
data_table_envsum_temp <- subset(data_table, Type == "environmental")
data_table_envsum <- data_table_envsum_temp %>%
  group_by(Treatment,variable) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            min=min(value),
            max=max(value))
data_table_envsum_temp2 <- subset(data_table_envsum_temp, Treatment == "Soil (environmental)")
data_table_envsum2 <- data_table_envsum_temp2 %>%
  group_by(Site,variable) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            min=min(value),
            max=max(value))



length(unique(data_table$variable))
getPalette = colorRampPalette(brewer.pal(12, "Set3"))

##actually.....let's write out Set3 manually so I can arrange the colors.
##make red and yellow lower abundance groups bc those colors hurt my eyes lol
##make "other" grey

pal <- c("#8DD3C7","#FFFFB3","#BEBADA","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F","#FB8072","#D9D9D9")

data1 <- subset(data_table, Site == "1")
data2 <- subset(data_table, Site == "2")
data3 <- subset(data_table, Site == "3")



p1 <-ggplot(data1, aes(y=value,x=Sample,fill=variable, color=variable)) + 
  geom_bar(stat="identity") + 
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous("Relative Abundance (%)", sec.axis = sec_axis( ~ . * 1 , name = "Site 1", labels = NULL, breaks = NULL)) + 
  labs(fill='Class') +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  ggtitle("Relative Abundance of Bacterial Classes") +
  guides(color = "none", fill = guide_legend(ncol=10)) +
  facet_grid( ~ Treatment, scales = "free", switch = "y")

p1

p2 <-ggplot(data2, aes(y=value,x=Sample,fill=variable, color=variable)) + 
  geom_bar(stat="identity") + 
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous("Relative Abundance (%)", sec.axis = sec_axis( ~ . * 1 , name = "Site 2", labels = NULL, breaks = NULL)) +  
  labs(fill='Class') +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  guides(color = "none", fill = guide_legend(ncol=10)) +
  facet_grid( ~ Treatment, scales = "free", switch = "y")
p3 <-ggplot(data3, aes(y=value,x=Sample,fill=variable, color=variable)) + 
  geom_bar(stat="identity") + 
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        #strip.background = element_blank(),
        #strip.text.x = element_blank(),
        legend.position = "bottom") +
  scale_y_continuous("Relative Abundance (%)", sec.axis = sec_axis( ~ . * 1 , name = "Site 3", labels = NULL, breaks = NULL)) + 
  xlab("Samples") +
  labs(fill='Class') +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  guides(color = "none", fill = guide_legend(ncol=10)) +
  facet_grid( ~ Treatment, scales = "free", switch = "y")



ggarrange(p1, p2, p3,
          ncol = 1,
          common.legend = TRUE,
          legend = "bottom",
          labels = "AUTO",
          font.label=(list(size = 20)))



ggsave("abundace_bar.tiff", device = "tiff", dpi = 700)
ggsave("abundace_bar.png", device = "png", dpi = 700)
ggsave("abundace_bar.pdf", device = "pdf", dpi = 700)

###### BOXPLOT of ENVIRONMENTAL ####################

data_table_nat_dump <- subset(data_table, Type == "environmental")
#data_table_nat_dump <- filter(data_table_nat_dump, !Treatment == "Soil (environmental)")


data_table_nat_dump$Treatment <- factor(data_table_nat_dump$Treatment, levels=c("Top (environmental)","Bottom (environmental)", "Soil (environmental)"))


data_table_nat_dump_summary <- data_table_nat_dump %>%
  group_by(Treatment,variable) %>%
  summarise(mean=mean(value),
            sd=sd(value),
            count = n())


pal <- c("#8DD3C7","#FFFFB3","#BEBADA","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F","#FB8072","#D9D9D9")
label <- c("Top", "Bottom", "Soil")

p <-ggboxplot(data_table_nat_dump, x="Treatment",y="value", fill="variable", lwd=1) +
  theme_classic(base_size=15)+
  theme(#aspect.ratio=1,
        legend.position = "none",
        strip.background = element_blank()) +
  scale_fill_manual(values=pal) + 
  scale_x_discrete(labels = label) +
  labs(y="Relative Abundance of Bacterial Classes (%)", x="Refuse Dump Layer") +
  stat_compare_means(method = "anova") +
  facet_wrap( . ~ variable , scales = "free", ncol = 3, nrow = 4) +
  ggtitle("Abundance of Bacterial Classes in Environmental Samples")

p

ggsave("nat_soil_box.tiff", device = "tiff", dpi = 700)
ggsave("nat_soil_box.png", device = "png", dpi = 700)
ggsave("nat_soil_box.pdf", device = "pdf", dpi = 700)



######## BOXPLOT of EXPERIMENTAL #######################

data_table_exp <- subset(data_table, Type == "experimental")


data_table_exp$Treatment <- factor(data_table_exp$Treatment, levels=c("Soil only","Net only","Soil & Net","No additions"))


pal <- c("#8DD3C7","#FFFFB3","#BEBADA","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#BC80BD","#CCEBC5","#FFED6F","#FB8072","#D9D9D9")

p <-ggboxplot(data_table_exp, x="Treatment",y="value", fill="variable", lwd=1) +
  theme_classic(base_size=15)+
  theme(#aspect.ratio=1,
    legend.position = "none",
    strip.background = element_blank(),
    axis.text.x = element_text(angle=25, hjust = 1)) +
  scale_fill_manual(values=pal) + 
  labs(y="Relative Abundance of Bacterial Classes (%)", x="Treatment") +
  stat_compare_means(method = "t.test", ref.group = "Soil only", label= "p.signif", vjust=0.2) +
  facet_wrap( . ~ variable , scales = "free", ncol = 3)+
  ggtitle("Abundance of Bacterial Classes in Experimental Treatments")

p

ggsave("exp_box.tiff", device = "tiff", dpi = 700)
ggsave("exp_box.png", device = "png", dpi = 700)
ggsave("exp_box.pdf", device = "pdf", dpi = 700)


########## focus on gamma
###need to pull in counts file where different taxonomic levels are preserved

metadata <- read.csv("panama_phyloseq_metadata_all.csv", header=TRUE)
otu <- read.csv("panama_phyloseq_otu_all.csv", header=TRUE, check.names = FALSE)



metadata$Site<- as.character(metadata$Site)
metadata$Replicate<- as.character(metadata$Replicate)
metadata$Subreplicate<- as.character(metadata$Subreplicate)


#metadata <- metadata %>%  mutate(Treatment = str_replace(Treatment, "top","Top"))
#metadata <- metadata %>% mutate(Treatment = str_replace(Treatment, "bottom","Bottom"))

metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "top","Top (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "bottom","Bottom (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil(?!_)","Soil (environmental)"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil_X","Soil only"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "soil_net","Soil & Net"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "X_net","Net only"))
metadata <- metadata %>%
  mutate(Treatment = str_replace(Treatment, "X_X","No additions"))

metadata$Treatment <- factor(metadata$Treatment, 
                             levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)","Soil (environmental)"))
colnames(metadata)[1] = "Sample"


#attach otu to metadata
data <- cbind(metadata,otu)
#drop columns
data <-data[-c(10)]


data <- melt(data,id=c("Sample","Site","Type","Initial","Treatment","Soil","Net",
                       "Replicate","Subreplicate"))
head(data)

data[c("Kingdom","Phylum","Class","Order","Family","Genus")] <- str_split_fixed(data$variable, "_", 6)
data <-data[-c(10)]
head(data)

data$percent <- (data$value)*100



data$Treatment <- factor(data$Treatment, 
                             levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)","Soil (environmental)"))



#specific Orders



data_gamma <- subset(data, Class == "Gammaproteobacteria" & Type == "environmental" & Soil == "FALSE")
length(unique(data_gamma$Order))




##Gamma

top <- data_gamma %>%
  group_by(Order) %>%
  summarise(sum = sum(percent))
top <- top %>% arrange(desc(sum))
top$Order



top_list <- c("Pseudomonadales", "Enterobacterales","Burkholderiales","Xanthomonadales")

#get top microbes
data_filt <- filter(data_gamma, Order %in% top_list)
data_filt <- select(data_filt, -c("Kingdom","Phylum","Class","Family","Genus"))
#get other microbes to condense into "other"
other_filt <- filter(data_gamma, !Order %in% top_list)
other_filt <- other_filt %>%
  group_by(Sample, Site, Type, Initial, Treatment, Soil, Net, Replicate, Subreplicate,value, percent) %>%
  summarise(percent=sum(percent))
other_filt$Order<- "Other"

data_gamma <- rbind(data_filt, other_filt)
data_gamma$Order <- factor(data_gamma$Order, levels=c("Pseudomonadales", "Enterobacterales","Burkholderiales",
                                                      "Xanthomonadales", "Other"))


pal <- c("#e78ac3","#8da0cb","#a6d854","#e5c494","#D9D9D9")
p_gam<-ggplot(data_gamma, aes(y=percent,x=Sample,fill=Order, color=Order)) + 
  geom_bar(stat="identity") + 
  theme_classic(base_size=15) + 
  theme(axis.text.x=element_blank(),
        strip.background = element_blank(),
        legend.position = "bottom") +
  #scale_y_continuous("Relative Abundance", sec.axis = sec_axis( ~ . * 1 , name = "Site 1", labels = NULL, breaks = NULL)) + 
  labs(fill='Order') +
  ylab("Relative Abundance of Gammaproteobacteria Orders (%)") +
  xlab("Samples") +
  scale_fill_manual(values = pal) +
  scale_color_manual(values = pal) +
  #ggtitle("Relative Abundance Gammapr") +
  guides(color = "none", fill = guide_legend(ncol=10)) +
  facet_grid( ~ Treatment, scales = "free", switch = "y")

p_gam


data_gamma_summary <- data_gamma %>%
  group_by(Sample, Treatment, Order) %>%
  summarise(value=sum(percent))

labels <- c("Top", "Bottom")

p_gam_box <-ggboxplot(data_gamma_summary, x="Treatment",y="value", fill="Order", lwd=1) +
  theme_classic(base_size=15)+
  theme(#aspect.ratio=1,
    legend.position = "none",
    #axis.text.x = element_text(angle=45, hjust = 1),
    strip.background = element_blank()
    ) +
  scale_fill_manual(values=c("#e78ac3","#8da0cb","#a6d854","#e5c494","#D9D9D9")) + 
  scale_x_discrete(labels = labels) +
  labs(y="Relative Abundance (%)", x="Refuse Dump Layer") +
  stat_compare_means(method = "t.test", hjust = -0.2, label = "p.signif") +
  facet_wrap( . ~ Order , scales = "free", ncol = 5) +
  ggtitle("Relative Abundance of Gammaproteobacteria Orders")

p_gam_box




ggsave("gamma_box.tiff", device = "tiff", dpi = 700)
ggsave("gamma_box.png", device = "png", dpi = 700)
ggsave("gamma_box.pdf", device = "pdf", dpi = 700)



################################################################################
## old alpha div, ignore

#menhinick

menhinick <- function(x) {
  sum(x>0)/sqrt(sum(x))
}


richness <- ddply(data,~Sample,function(x){
  data.frame(Richness=menhinick(x[,10:1056]>0))
})
#rich_list<-list(metadata,richness)


shan <- ddply(data,~Sample,function(x){
  data.frame(Shannon=diversity(x[,10:1056],index="shannon"))
})
#shan_list<-list(metadata,shan)



invsimpson <- ddply(data,~Sample,function(x){
  data.frame(InverseSimpson=diversity(x[,10:1056],index="invsimpson"))
})
#in_list<-list(metadata,invsimpson)

#evenness
pielou <-ddply(data,~Sample,function(x){
  data.frame(Pielou=exp(diversity(x[,10:1056], index="shannon"))/sum(x[,10:1056]>0))
})
#pie_list<-list(metadata,pielou)

metadata <- merge(metadata, richness, by="Sample")
metadata <- merge(metadata, shan, by="Sample")
metadata <- merge(metadata, invsimpson, by="Sample")
metadata <- merge(metadata, pielou, by="Sample")

#write.csv(metadata, "panama_phyloseq_metadata_diversity_indices.csv", row.names = FALSE)




#### Diversity plots ############


metadata <- read.csv("panama_phyloseq_metadata_diversity_indices.csv", header=TRUE)
metadata$Treatment <- factor(metadata$Treatment, 
                             levels=c("Soil only","Net only","Soil & Net","No additions","Top (environmental)","Bottom (environmental)","Soil (environmental)"))


# "#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d"
# "#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf"
# "#8dd3c7","#ffffb3","#bebada","#80b1d3","#fdb462","#b3de69","#fccde5"

pshan <-ggboxplot(metadata, x="Treatment",y="Shannon", fill="Treatment", lwd=1) +
  theme_classic(base_size=20)+
  theme(aspect.ratio=1,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(method="anova",label.x = 1.5, label.y = 5) +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")) + 
  labs(y="Shannon Diversity Index") + 
  ggtitle("Shannon Diversity")

pshan


pinv <-ggboxplot(metadata, x="Treatment",y="InverseSimpson", fill="Treatment", lwd=1) +
  theme_classic(base_size=20)+
  theme(aspect.ratio=1,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(method="anova",label.x = 1.5, label.y = 85) +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")) + 
  labs(y="Inverse Simpson Index") + 
  ggtitle("Inverse Simpson Diversity")

pinv


ppie <-ggboxplot(metadata, x="Treatment",y="Pielou", fill="Treatment", lwd=1) +
  theme_classic(base_size=20)+
  theme(aspect.ratio=1,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(method="anova",label.x = 1.5, label.y = 0.5) +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")) + 
  labs(y="Pielou's Evenness Index") + 
  ggtitle("Evenness")

ppie


prich <-ggboxplot(metadata, x="Treatment",y="Richness", fill="Treatment", lwd=1) +
  theme_classic(base_size=20)+
  theme(aspect.ratio=1,
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  stat_compare_means(method="anova",label.x = 1.5, label.y = 19) +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d")) + 
  labs(y="Menhinick's Index") + 
  ggtitle("Richness")

prich

ggarrange(pshan, pinv, ppie, prich, 
          common.legend = TRUE,
          legend = "bottom")

#ggsave("diversity.png")
#ggsave("diversity.pdf")

shan_summary <- metadata %>%
  group_by(Treatment) %>%
  summarise(mean=mean(Shannon),
            sd=sd(Shannon),
            count=n())


invsim_summary <- metadata %>%
  group_by(Treatment) %>%
  summarise(mean=mean(InverseSimpson),
            sd=sd(InverseSimpson),
            count=n())


even_summary <- metadata %>%
  group_by(Treatment) %>%
  summarise(mean=mean(Pielou),
            sd=sd(Pielou),
            count=n())


rich_summary <- metadata %>%
  group_by(Treatment) %>%
  summarise(mean=mean(Richness),
            sd=sd(Richness),
            count=n())




#centroid analysis, ignore



dist_2cent_treatment <- dist_to_centroids(dist, Treatment)
colnames(dist_2cent_treatment)[1] <- "Sample"
dist_2cent_treatment <- merge(dist_2cent_treatment, metadata, by = "Sample")
dist_2cent_treatment <- dist_2cent_treatment %>%
  filter(CentroidGroup == Treatment)
dist_2cent_treatment$Variable <- "Treatment"

dist_2cent_site <- dist_to_centroids(dist, Site)
colnames(dist_2cent_site)[1] <- "Sample"
dist_2cent_site <- merge(dist_2cent_site, metadata, by = "Sample")
dist_2cent_site <- dist_2cent_site %>%
  filter(CentroidGroup == Site)
dist_2cent_site$Variable <- "Site"



dist2cent <- rbind(dist_2cent_treatment, dist_2cent_site)
dist2cent$CentroidGroup <- as.factor(dist2cent$CentroidGroup)


p_cent <- ggplot(data=dist2cent, aes(x=CentroidGroup, y = CentroidDistance)) +
  geom_jitter(data=dist2cent,width = 0.2, size = 1.5) +
  geom_boxplot(outlier.shape = NA, lwd = 1, alpha = 0) +
  theme_classic(base_size=18) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(y= "Distance from Centroid") + 
  stat_compare_means(method = "anova") +
  facet_grid(~Variable, scales = "free") +
  ggtitle("Dispersion of Samples by Variable")
p_cent
