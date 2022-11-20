#Add libraries that will be used
library("phyloseq")
packageVersion("phyloseq") #1.32.0
library("ggplot2") 
packageVersion("ggplot2") #3.3.6

##load in count, tax, and sample info data##
count_tab <- read.table("ASVs_counts.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")


tax_tab <- as.matrix(read.table("ASVs_taxonomy.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))



sample_info_tab <- read.csv("metadata3yrs.csv", header=T, row.names=1)

#transform info metadata to character; helps with downstream analysis
sample_info_tab$Patient <- as.character(sample_info_tab$Patient)
sample_info_tab$Bin_Years <- as.character(sample_info_tab$Bin_Years)
sample_info_tab$Bin_Change <- as.character(sample_info_tab$Bin_Change)
#sample_info_tab$Year <- as.character(sample_info_tab$Year)

##remove unannotated samples from count_tab and sample_info_tab
count_tab <- count_tab[,colnames(count_tab) %in% rownames(sample_info_tab)]
sample_info_tab <- sample_info_tab[(rownames(sample_info_tab) %in% colnames(count_tab)),]

#merge Sinfo and sample info tab
Sinfo <- read.csv("Sinfo.csv", header=T)
sample_info_tab$Sample <- row.names(sample_info_tab)
sample_info_tab <- merge(sample_info_tab, Sinfo)
rownames(sample_info_tab) <- sample_info_tab$Sample

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
#It only removes 1 sample

#any taxa with zero counts? Yes, 21
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

#prune the dataset at 10^4
#This removes 1 sample only
ASV_physeq = prune_samples(sample_sums(ASV_physeq) > 10000, ASV_physeq)
ASV_physeq

#mean counts per sample
mean(colSums(otu_table(ASV_physeq)))
###############################################################################################################
##remove unannotated samples from count_tab and sample_info_tab
count_tab[,colnames(count_tab) %in% rownames(sample_info_tab)]

count_tab <- count_tab[,colnames(count_tab) %in% rownames(sample_info_tab)]
sample_info_tab <- sample_info_tab[(rownames(sample_info_tab) %in% colnames(count_tab)),]
###############################################################################################################
#save to rda file
save(ASV_physeq, file = "ASV_physeq_silva.rda")
save(sample_info_tab, file = "sample_info.rda")
save(count_tab, file = "count_tab.rda")
save(tax_tab, file = "tax_tab.rda")

ASV_physeq
all_metadata <- as.data.frame(sample_data(ASV_physeq))
write.csv (all_metadata, "all_metadata.csv")
