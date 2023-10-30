library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(vegan)
library(QsRutils)
library(dplyr)
#references
#http://www.hiercourse.com/docs/microbial/02_alphaDiversity.html
#https://rpubs.com/an-bui/vegan-cheat-sheet

#load data
tax <- read.csv("panama_dada2_taxonomy.csv", header=TRUE, row.names=1)
counts <- read.csv("panama_dada2_counts.csv", header=TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("panama_metadata.csv", header=TRUE, row.names= 1)


#change things to character (not numeric)
metadata$Site<- as.character(metadata$Site)
metadata$Replicate<- as.character(metadata$Replicate)
metadata$Subreplicate<- as.character(metadata$Subreplicate)


#create phyloseq ob
count.ps <- otu_table(as.matrix(counts),taxa_are_rows = TRUE) #define our count table as our otu table
tax.ps <- tax_table(as.matrix(tax)) #define our taxonomy table as a taxonomy table for phyloseq to use
meta.ps <- sample_data(metadata) #define our metadata as sample information
ps <- phyloseq(count.ps, tax.ps, meta.ps) #put them all together into a phyloseq object



ps <- subset_samples(ps, sample_sums(ps) >=10000)




##take top 20 abundant, transform to proportion
#top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50] 
ps <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
#ps <- prune_taxa(top20, ps)

#write otu_table were asv is replaced with Class 
#https://github.com/joey711/phyloseq/issues/1559
ps <- tax_glom(ps, "Class")
otu_class <- otu_table(ps)
otu_class <- t(otu_class)
df_class <- as.data.frame(otu_class)
colnames(df_class) <- as.data.frame(tax_table(ps))$Class
df_class[is.na(df_class)] <- 0

newmeta <- vegansam(ps)
#write data
write.csv(newmeta, "panama_phyloseq_metadata_class_silva128.csv", row.names=TRUE)
write.csv(df_class, "panama_phyloseq_otu_class_silva128.csv", row.names=TRUE)



#write otu_table were asv is replaced with order
#https://github.com/joey711/phyloseq/issues/1559
ps <- tax_glom(ps, "Order")
otu_ord <- otu_table(ps)
otu_ord <- t(otu_ord)
df_ord <- as.data.frame(otu_ord)
colnames(df_ord) <- as.data.frame(tax_table(ps))$Order
df_ord[is.na(df_ord)] <- 0
newmeta <- vegansam(ps)
#write data
write.csv(newmeta, "panama_phyloseq_metadata_order.csv", row.names=TRUE)
write.csv(df_fam, "panama_phyloseq_otu_order.csv", row.names=TRUE)

#write otu_table were asv is replaced with genus
#https://github.com/joey711/phyloseq/issues/1559
ps <- tax_glom(ps, "Genus")
otu_genus <- otu_table(ps)
otu_genus <- t(otu_genus)
df_genus <- as.data.frame(otu_genus)
colnames(df_genus) <- as.data.frame(tax_table(ps))$Genus
df_genus[is.na(df_genus)] <- 0
newmeta <- vegansam(ps)
#write data
write.csv(newmeta, "panama_phyloseq_metadata_genus.csv", row.names=TRUE)
write.csv(df_genus, "panama_phyloseq_otu_genus.csv", row.names=TRUE)

#write otu_table were asv is replaced with Phylum 
#https://github.com/joey711/phyloseq/issues/1559
ps <- tax_glom(ps, "Phylum")
otu_phylum <- otu_table(ps)
otu_phylum <- t(otu_phylum)
df_phylum <- as.data.frame(otu_phylum)
colnames(df_phylum) <- as.data.frame(tax_table(ps))$Phylum
df_phylum[is.na(df_phylum)] <- 0
newmeta <- vegansam(ps)
#write data
write.csv(newmeta, "panama_phyloseq_metadata_phylum.csv", row.names=TRUE)
write.csv(df_phylum, "panama_phyloseq_otu_phylum.csv", row.names=TRUE)

#write otu_table were asv is maintains phylogeny levels
otu_out <- otu_table(ps)
otu_out <- t(otu_out)
df_out <- as.data.frame(otu_out)

tax$Tax <- paste(tax$Kingdom, tax$Phylum, tax$Class, tax$Order, tax$Family, tax$Genus, sep="_")
colnames(df_out) <- as.data.frame(t(tax$Tax))

df_out[is.na(df_out)] <- 0

newmeta <- vegansam(ps)
#write data
write.csv(newmeta, "panama_phyloseq_metadata_all.csv", row.names=TRUE)
write.csv(df_out, "panama_phyloseq_otu_all.csv", row.names=TRUE)










