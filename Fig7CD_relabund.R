#Crohn's Dysbiosis Index (CDI)
##For Fig 6C-D, need to make a new, smaller phyloseq object
##load in phyloseq data
#Add several libraries that will be used - primary ones for analysis are phyloseq and DESeq2
#others are primarily for visualization
library("phyloseq")
library("ggplot2")
library(dplyr)
library(tidyr)

##load in count, tax, and sample info data##
count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")


tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))



sample_info_tab <- read.csv("MDindex_late_short.csv", header=T, row.names=1)

#transform info metadata to character; helps with downstream analysis
sample_info_tab$Patient <- as.character(sample_info_tab$Patient)

##remove unannotated samples from count_tab and sample_info_tab
count_tab <- count_tab[,colnames(count_tab) %in% rownames(sample_info_tab)]
sample_info_tab <- sample_info_tab[(rownames(sample_info_tab) %in% colnames(count_tab)),]

##make a phyloseq object##
#transform data (?) into phyloseq object
#for some reason, if you don't do this, phyloseq will not recognize "Month" header
sample_info_tab_phy <- sample_data(sample_info_tab)
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)

#make phyloseq object
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)
ASV_physeq

###############################################################################################################
#This does some basic data cleanup
#Looks at counts per ASV, as well as counts per sample, and filters by # of counts
#The way it's currently written, only removes 1 sample

#any taxa with zero counts?
any(taxa_sums(ASV_physeq) == 0)
sum(taxa_sums(ASV_physeq) == 0)

#removes taxa with zero counts 
ASV_physeq = prune_taxa(taxa_sums(ASV_physeq) > 0, ASV_physeq)

#remove taxa that aren't annotated at the phylum level
ASV_physeq <- subset_taxa(ASV_physeq, !is.na(Phylum))

#how many reads per sample, and what is the distribution
readsumsdf = data.frame(nreads = sort(taxa_sums(ASV_physeq), TRUE), sorted = 1:ntaxa(ASV_physeq), 
                        type = "ASVs")

readsumsdf = rbind(readsumsdf, data.frame(nreads = sort(sample_sums(ASV_physeq), 
                                                        TRUE), sorted = 1:nsamples(ASV_physeq), type = "Samples"))
#plot the number of reads per sample 
title = "Read counts per ASV or sample"
p = ggplot(readsumsdf, aes(x = sorted, y = nreads)) + geom_bar(stat = "identity")
p + ggtitle(title) + ylab("Number of Reads") + xlab("") +
  scale_y_log10() + facet_wrap(~type, 1, scales = "free") + theme_bw()

#mean counts per sample
mean(colSums(otu_table(ASV_physeq)))
###############################################################################################################
##remove unannotated samples from count_tab and sample_info_tab
count_tab[,colnames(count_tab) %in% rownames(sample_info_tab)]

count_tab <- count_tab[,colnames(count_tab) %in% rownames(sample_info_tab)]
sample_info_tab <- sample_info_tab[(rownames(sample_info_tab) %in% colnames(count_tab)),]
###############################################################################################################
#convert everthing to to % relative abundance
ASV_physeq = transform_sample_counts(ASV_physeq, function(x) (x / sum(x))*100 )
###############################################################################################################
#calculate relative abundance of each phylum 
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Phylum"))

# making a vector of family names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="Phylum"))[,2])
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
                                  "Early" =sample_info_tab$EarlyCDHi,
                                  stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)

#average relative abundance by person
x <- major_taxa_for_plot.g2 %>% 
  group_by(MTaxa, Early, Patient) %>%
  summarize(mean_relAb = (mean(Proportion, na.rm = TRUE)))

library(viridis)

ggplot(x, aes(x=Patient, y=mean_relAb, fill=MTaxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  labs(x="Early High Dysbiosis Index", y="Relative Abundance (%)") +
  scale_fill_viridis(discrete=TRUE, name = "Phylum")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=12), legend.title=element_text(size=16, face="bold"))

##relative abundance
x <- x %>% 
  group_by(MTaxa, Early) %>%
  summarize(mean_relAb = (mean(mean_relAb, na.rm = TRUE)))

x <- x[x$MTaxa=="Proteobacteria" |
         x$MTaxa=="Firmicutes"|
         x$MTaxa=="Bacteroidota"|
         x$MTaxa=="Actinobacteriota"|
         x$MTaxa=="Verrucomicrobiota",]

write.csv(x, "Rel_abund_phy_CDI.csv")

phy_abund <- ggplot(x, aes(x=Early, y=mean_relAb, fill=MTaxa)) +
  geom_bar(width=0.6, stat="identity") +
  theme_bw() +
  labs(x="Early High CDI", y="Relative Abundance (%)") +
  scale_fill_viridis(discrete=TRUE, name = "Phylum",
                     labels = c("Actinobacteria",
                                "Bacteroidetes",
                                "Firmicutes",
                                "Proteobacteria",
                                "Verrucomicrobiota"))+
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=12, face="bold"))
phy_abund

