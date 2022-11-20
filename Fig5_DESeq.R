#This script first removes ASVs with <1% prevalence
#Then combines taxa at the genus/family level 
#Then uses DESeq2 to call differential abundance at the genus and family levels 
#modified from preprocess_vst.R

library("phyloseq")
library("ggplot2")
library("plyr")
library("DESeq2")
library("vegan")
library("viridis")

##load in phyloseq data 
load("ASV_physeq_silva.rda")
load("sample_info.rda")
load("count_tab.rda")

##############################################################################################################################################
#Preprocessing for DESEq
# Compute prevalence (how many samples each taxa appears in) of each feature, store as data.frame
prevdf = apply(X = otu_table(ASV_physeq),
               MARGIN = ifelse(taxa_are_rows(ASV_physeq), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ASV_physeq),
                    tax_table(ASV_physeq))

#a plot of ASV taxa abundances 
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ASV_physeq, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ASV_physeq),color=Phylum)) + 
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) + 
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") + 
  facet_wrap(~Phylum) + theme_bw() + guides(color = FALSE, size = FALSE)

# Define prevalence threshold as 1% of total samples
prevalenceThreshold = 0.01 * nsamples(ASV_physeq)
prevalenceThreshold

# Remove taxa with <1% prevalence
# This removed 1662 taxa; we have a lot of low count taxa! 
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ASV_physeq_2 = prune_taxa(keepTaxa, ASV_physeq)
##########################################################################################################
#significance testing at the genera level 
#merge ASV_physeq_2 at the genus level 
ASV_physeq_3 <- tax_glom(ASV_physeq_2, taxrank="Genus")

#convert to deseq object
#This will model using non-continuous Bin_Years, but will control for patient
#Note to self: this is right. If you don't believe it, check the bioconductor vingette 'analyzing RNA-seq data with DESeq2', 'Multi-factor designs' section
ASV_deseq_genus <- phyloseq_to_deseq2(ASV_physeq_3, ~Patient + Bin_EarlyvLate)

#workaround to deal with 0s
cts <- counts(ASV_deseq_genus)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(sum(log(row[row != 0]))/length(row)))
ASV_deseq_genus <- estimateSizeFactors(ASV_deseq_genus, geoMeans=geoMeans)

#deseq standard analysis
ASV_deseq_genus <- DESeq(ASV_deseq_genus)

##########################################################################################################

# pulling out our results table, we specify the object, the p-value we are going to use to filter our results 
#and what contrast we want to consider by first naming the column, then the two groups we care about
deseq_res_EvL <- results(ASV_deseq_genus, alpha=0.01, contrast=c("Bin_EarlyvLate","Late","Early"))

# we can get a glimpse at what this table currently holds with the summary command
summary(deseq_res_EvL)

# stitch together with the ASV's taxonomic annotations for a quick look at both together
deseq_res_EvL_tax <- cbind(as(deseq_res_EvL, "data.frame"), as(tax_table(ASV_physeq_2)[row.names(deseq_res_EvL), ], "matrix"))

# and sort that table by the baseMean column
sum_EvL <- deseq_res_EvL_tax[order(deseq_res_EvL_tax$baseMean, decreasing=T), ]

write.csv(sum_EvL, "sum_EvL_genus_all.csv")

# subset this table to only include these that pass our specified significance level
sum_EvL_sig <- sum_EvL[which(sum_EvL$padj < 0.01), ]
write.csv(sum_EvL_sig, "sum_EvL_genus.csv")

##########################################################################################################
#I did some manual work here to identify whether significantly altered taxa were relevant to CF
#Under CFChange column
#classified taxa as more CF like with age ("CF"), less CF-like with age ("Normal"), or not CF associated ("None")

##########################################################################################################
#Making graphs
#DESeq2 visualization methods from phyloseq tutorial
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

#read in manually annotated csv
EvL_plot <- read.csv("sum_EvL_genus_update.csv", row.names = 1)

# Genus order EvL
x = tapply(EvL_plot$log2FoldChange, EvL_plot$Genus, function(x) max(x))
x = sort(x, TRUE)
EvL_plot$Genus = factor(as.character(EvL_plot$Genus), levels=names(x))

gen_EvL <- ggplot(EvL_plot, aes(x=log2FoldChange, y=Genus, color=CFChange)) + 
  geom_point(size=5) + 
  scale_color_viridis(discrete=TRUE, labels = c("More CF-like with age", "Not CF-associated",
                                                "Less CF-like with age")) +
  labs(x="Log2 Fold Change", y="Genus",color="CF-Associated Changes")+
  geom_vline(xintercept = 0) +
  theme(axis.text=element_text(size=8), 
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 8, face = "italic"),
        axis.title=element_text(size=14,face="bold"), 
        legend.text=element_text(size=12), 
        legend.title=element_text(size=14, face="bold"))+
  geom_label(label="Decreased with Age", x=-15, y=2, label.size = 0.4, color = "black", size =3) +
  geom_label(label="Increased with Age", x= 15, y=2, label.size = 0.4, color = "black", size =3)

#using tiff() and dev.off
tiff('Fig5.tiff', units="in", width=9, height=10, res=300)
gen_EvL
dev.off()

##########################################################################################################
#Relative abundance and prevalence plots for supplement Fig S3
#Blautia, Bacteroides, and Akkermansia 

##load in phyloseq data 
load("ASV_physeq_silva.rda")
load("sample_info.rda")

#Transform to relative abundance 
ASV_physeq = transform_sample_counts(ASV_physeq, function(x) (x / sum(x))*100 )

phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Genus"))

# making a vector of family names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="Genus"))[,6])
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

dim(phyla_counts_tab)

#make a different file for manipulating
major_taxa_for_plot <- as.data.frame(phyla_counts_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,
                                  "Sample"=row.names(sample_info_tab), 
                                  "Age"=sample_info_tab$Bin_Months, 
                                  "Patient"=sample_info_tab$Patient, 
                                  "Year" =sample_info_tab$Bin_Years,
                                  "BinEvL" = sample_info_tab$Bin_EarlyvLate,
                                  stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)
########################################################################################
#subset to keep taxa of interest
Ecoli <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Escherichia.Shigella", ]

#ribbon plot
Ecoli_1 <- ggplot(Ecoli, aes(x=Age_Days, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), size =.6, show.legend = FALSE) + theme_bw() +
  geom_point(size=1) +
  ggtitle("Escherichia/Shigella") +
  labs(x = "Age (Days)", y = "Relative Abundance (%)") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))
Ecoli_1

#relative abundance plot 
EcoliEvL <- Ecoli[Ecoli$BinEvL != "Middle", ]

#boxplot
Ecoli_2 <- ggplot(EcoliEvL, aes(x=Age_Days, y=Proportion)) + 
  geom_boxplot(aes(x=BinEvL, y=Proportion), size=.3) + theme_bw() +
  ggtitle("") +
  scale_y_continuous(trans='log10') +
  labs(x = "Age", y = "Relative Abundance (%)") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))+
  scale_x_discrete(labels=c("<6 months", "2+ years"))
Ecoli_2

library(ggpubr)
tiff('FigS3A.tiff', units="in", width=6.5, height=4.2, res=300)
ggarrange(Ecoli_1,Ecoli_2)
dev.off()

########################################################################################
#subset to keep taxa of interest
Blautia <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Blautia", ]
head(Blautia)

#ribbon plot
Blautia_1 <- ggplot(Blautia, aes(x=Age_Days, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), size =.6, show.legend = FALSE) + theme_bw() +
  geom_point(size=1) +
  ggtitle("Blautia") +
  labs(x = "Age (Days)", y = "Relative Abundance (%)") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))
Blautia_1

#relative abundance plot 
BlautiaEvL <- Blautia[Blautia$BinEvL != "Middle", ]

#boxplot
Blautia_2 <- ggplot(BlautiaEvL, aes(x=Age_Days, y=Proportion)) + 
  geom_boxplot(aes(x=BinEvL, y=Proportion), size=.3) + theme_bw() +
  ggtitle("") +
  scale_y_continuous(trans='log10') +
  labs(x = "Age", y = "Relative Abundance (%)") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))+
  scale_x_discrete(labels=c("<6 months", "2+ years"))
Blautia_2

tiff('FigS3C.tiff', units="in", width=6.5, height=4.2, res=300)
ggarrange(Blautia_1,Blautia_2)
dev.off()

########################################################################################
#subset to keep taxa of interest
Bacteroides <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Bacteroides", ]
head(Bacteroides)

#ribbon plot
Bacteroides_1 <- ggplot(Bacteroides, aes(x=Age_Days, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), size =.6, show.legend = FALSE) + theme_bw() +
  geom_point(size=1) +
  ggtitle("Bacteroides") +
  labs(x = "Age (Days)", y = "Relative Abundance (%)") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))
Bacteroides_1

#relative abundance plot 
BacteroidesEvL <- Bacteroides[Bacteroides$BinEvL != "Middle", ]

#boxplot
Bacteroides_2 <- ggplot(BacteroidesEvL, aes(x=Age_Days, y=Proportion)) + 
  geom_boxplot(aes(x=BinEvL, y=Proportion), size=.3) + theme_bw() +
  ggtitle("") +
  scale_y_continuous(trans='log10') +
  labs(x = "Age", y = "Relative Abundance (%)") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))+
  scale_x_discrete(labels=c("<6 months", "2+ years"))
Bacteroides_2

tiff('FigS3B.tiff', units="in", width=6.5, height=4.2, res=300)
ggarrange(Bacteroides_1,Bacteroides_2)
dev.off()


