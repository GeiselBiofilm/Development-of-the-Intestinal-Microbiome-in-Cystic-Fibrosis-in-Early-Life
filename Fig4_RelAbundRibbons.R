#Add several libraries that will be used - primary ones for analysis are phyloseq and DESeq2
#others are primarily for visualization
library("phyloseq")
library("vegan")
library("DESeq2")
library("ggplot2")
library("dendextend")
library("tidyr")
library("viridis")
library("reshape")
library("gplots")
library("ggpubr")
library("dplyr")

load("ASV_physeq_silva.rda")
load("sample_info.rda")

###############################################################################################################
#convert everything to to % relative abundance
ASV_physeq = transform_sample_counts(ASV_physeq, function(x) (x / sum(x))*100 )

###############################################################################################################
#calculate relative abundance of each Phylum 
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Phylum"))

# making a vector of Family names to set as row names
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
                                  stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)

###############################################################################################################
#average relative abundance at the phylum level for all time points
#find mean of each Taxa from each person 
avg_by_patient <- major_taxa_for_plot.g2%>%                                        
  group_by(MTaxa, Patient) %>%                        
  summarise_at(vars(Proportion),
               list(mean_relab = mean))
head(avg_by_patient)

#average by taxa
avg_total <- avg_by_patient%>%                                        
  group_by(MTaxa) %>%                        
  summarise_at(vars(mean_relab),
               list(mean_relab = mean))
avg_total

#average relative abundance at the phylum level for 6months or younger 
##0-6mo relative abundance
major_taxa_6mo <- subset(major_taxa_for_plot.g2, Age_Days < 183)

#find mean of each Taxa from each person 
avg_mo6_by_patient <- major_taxa_6mo %>%                                        
  group_by(MTaxa, Patient) %>%                        
  summarise_at(vars(Proportion),
               list(mean_relab = mean))
head(avg_mo6_by_patient)

#average by taxa
avg_mo6 <- avg_mo6_by_patient%>%                                        
  group_by(MTaxa) %>%                        
  summarise_at(vars(mean_relab),
               list(mean_relab_6mo = mean))
avg_mo6

#average relative abundance at the phylum level for >183 days to 500 days 
major_taxa_middle1 <- subset(major_taxa_for_plot.g2, Age_Days >= 183 & Age_Days <= 500)

#find mean of each Taxa from each person 
avg_middle_by_patient1 <- major_taxa_middle1 %>%                                        
  group_by(MTaxa, Patient) %>%                        
  summarise_at(vars(Proportion),
               list(mean_relab = mean))
head(avg_middle_by_patient1)

#average by taxa
avg_middle1 <- avg_middle_by_patient1%>%                                        
  group_by(MTaxa) %>%                        
  summarise_at(vars(mean_relab),
               list(mean_relab_middle1 = mean))
avg_middle1

#average relative abundance at the phylum level for 500 days to 730 days 
major_taxa_middle2 <- subset(major_taxa_for_plot.g2, Age_Days > 500 & Age_Days <= 730)

#find mean of each Taxa from each person 
avg_middle_by_patient2 <- major_taxa_middle2 %>%                                        
  group_by(MTaxa, Patient) %>%                        
  summarise_at(vars(Proportion),
               list(mean_relab = mean))
head(avg_middle_by_patient2)

#average by taxa
avg_middle2 <- avg_middle_by_patient2%>%                                        
  group_by(MTaxa) %>%                        
  summarise_at(vars(mean_relab),
               list(mean_relab_middle2 = mean))
avg_middle2


##2+ year relative abundance
major_taxa_2yr <- subset(major_taxa_for_plot.g2, Age_Days > 730)

#find mean of each Taxa from each person 
avg_2yr_by_patient <- major_taxa_2yr %>%                                        
  group_by(MTaxa, Patient) %>%                        
  summarise_at(vars(Proportion),
               list(mean_relab = mean))
head(avg_2yr_by_patient)

#average by taxa
avg_2yr <- avg_2yr_by_patient%>%                                        
  group_by(MTaxa) %>%                        
  summarise_at(vars(mean_relab),
               list(mean_relab_2yr = mean))
avg_2yr

df_list <- list(avg_total, avg_mo6, avg_middle1, avg_middle2, avg_2yr)
relabund_table <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list, accumulate=FALSE)

head(relabund_table)

write.csv(relabund_table, "Avg_Rel_Abundances.csv")

###############################################################################################################
#Phylum level line graph 
#subset on phyla of interest
new_gplot <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa =="Proteobacteria" |
                                      major_taxa_for_plot.g2$MTaxa =="Firmicutes" |
                                      major_taxa_for_plot.g2$MTaxa =="Bacteroidota" |
                                      major_taxa_for_plot.g2$MTaxa =="Actinobacteriota" |
                                      major_taxa_for_plot.g2$MTaxa =="Verrucomicrobiota",]

#ribbon plot
Fig3A <- ggplot(new_gplot, aes(x=Age_Days, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), size =2) + theme_bw() +
  labs(x = "Age (Days)", y = "Relative Abundance (%)") +
  scale_color_viridis(discrete=TRUE, name = "Phylum", labels = c("Actinobacteria",
                                                                 "Bacteroidetes",
                                                                 "Firmicutes",
                                                                 "Proteobacteria",
                                                                 "Verrucomicrobiota")) +
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=14,face="bold"), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=14, face="bold"),
        legend.position = c(0.5, 0.55),
        plot.margin=unit(c(0.4,1,0.4,0.4),"cm"),
        legend.key.size = unit(0.7, 'cm'))
Fig3A

###############################################################################################################
#subset Proteobacteria
Phy_prot = subset_taxa(ASV_physeq, Phylum=="Proteobacteria")
count_tab <- as.data.frame(otu_table(Phy_prot))

phyla_counts_tab <- otu_table(tax_glom(Phy_prot, taxrank="Family"))

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(Phy_prot, taxrank="Family"))[,5]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

#turn NAs into their own group and add this row to our Phylum count table:
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(phyla_and_unidentified_counts_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)

#check the dimensions of this table at this point
dim(phyla_and_unidentified_counts_tab)

# here, we'll only keep rows (taxa) that make up greater than 10% in any
temp_filt_major_taxa_proportions_tab <- as.data.frame(phyla_and_unidentified_counts_tab[apply(phyla_and_unidentified_counts_tab, 1, max) > 10, ])

# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
rownames(temp_filt_major_taxa_proportions_tab)

# though each of the filtered taxa made up less than 2% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(phyla_and_unidentified_counts_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(filt_major_taxa_proportions_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
head(major_taxa_for_plot.g2)

####################################################################################################################
#make a new sheet with just the taxa of interest
Prot_of_int <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Enterobacteriaceae", ]
dim(Prot_of_int)
###############################################################################################################
xProt1 <- major_taxa_for_plot.g2 %>% 
  group_by(MTaxa, Patient) %>% 
  summarize(mean_relAb = (mean(Proportion, na.rm = TRUE)))

xProt <- xProt1 %>% 
  group_by(MTaxa) %>% 
  summarize(mean_relAb = (mean(mean_relAb, na.rm = TRUE)))
xProt

sum(xProt$mean_relAb)

xProt$PhylaPct <- ((xProt$mean_relAb)/(sum(xProt$mean_relAb)))*100
xProt

FigS6D <- ggplot(xProt, aes(x="", y=PhylaPct, fill=MTaxa)) +
  geom_bar(stat="identity", width=1, color = "white") +
  coord_polar("y", start=0) + theme_void() + scale_fill_viridis(discrete=TRUE, name = "Family")+ 
  ggtitle("Proteobacteria") + theme(#legend.position = "bottom", 
                                    text = element_text(size = 16))

FigS6D

###############################################################################################################
#subset Actinobacteria
Act_prot = subset_taxa(ASV_physeq, Phylum=="Actinobacteriota")
count_tab <- as.data.frame(otu_table(Act_prot))

phyla_counts_tab <- otu_table(tax_glom(Act_prot, taxrank="Family"))

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(Act_prot, taxrank="Family"))[,5]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

#turn NAs into their own group and add this row to our Phylum count table:
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(phyla_and_unidentified_counts_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)

#check the dimensions of this table at this point
dim(phyla_and_unidentified_counts_tab)

# here, we'll only keep rows (taxa) that make up greater than 5% in any
temp_filt_major_taxa_proportions_tab <- as.data.frame(phyla_and_unidentified_counts_tab[apply(phyla_and_unidentified_counts_tab, 1, max) > 10, ])

# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
rownames(temp_filt_major_taxa_proportions_tab)

# though each of the filtered taxa made up less than 2% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(phyla_and_unidentified_counts_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(filt_major_taxa_proportions_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
head(major_taxa_for_plot.g2)

####################################################################################################################
#make a new sheet with just the taxa of interest
Act_of_int <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Bifidobacteriaceae", ]
dim(Act_of_int)
###############################################################################################################
xAct1 <- major_taxa_for_plot.g2 %>% 
  group_by(MTaxa, Patient) %>% 
  summarize(mean_relAb = (mean(Proportion, na.rm = TRUE)))

xAct <- xAct1 %>% 
  group_by(MTaxa) %>% 
  summarize(mean_relAb = (mean(mean_relAb, na.rm = TRUE)))
sum(xAct$mean_relAb)

xAct$PhylaPct <- ((xAct$mean_relAb)/(sum(xAct$mean_relAb)))*100
xAct

FigS6C <- ggplot(xAct, aes(x="", y=PhylaPct, fill=MTaxa)) +
  geom_bar(stat="identity", width=1, color = "white") +
  coord_polar("y", start=0) + theme_void() + scale_fill_viridis(discrete=TRUE, name = "Family") + 
  ggtitle("Actinobacteria") + theme(#legend.position = "bottom",
                                    text = element_text(size = 16))
FigS6C
###############################################################################################################
#subset Firmicutes
Firm_prot = subset_taxa(ASV_physeq, Phylum=="Firmicutes")
count_tab <- as.data.frame(otu_table(Firm_prot))

phyla_counts_tab <- otu_table(tax_glom(Firm_prot, taxrank="Family"))

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(Firm_prot, taxrank="Family"))[,5]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

#turn NAs into their own group and add this row to our Phylum count table:
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(phyla_and_unidentified_counts_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)

#check the dimensions of this table at this point
dim(phyla_and_unidentified_counts_tab)

# here, we'll only keep rows (taxa) that make up greater than 5% in any
temp_filt_major_taxa_proportions_tab <- as.data.frame(phyla_and_unidentified_counts_tab[apply(phyla_and_unidentified_counts_tab, 1, max) > 10, ])

# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
rownames(temp_filt_major_taxa_proportions_tab)

# though each of the filtered taxa made up less than 2% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(phyla_and_unidentified_counts_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(filt_major_taxa_proportions_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
head(major_taxa_for_plot.g2)

####################################################################################################################
#make a new sheet with just the taxa of interest
Firm_of_int <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Lachnospiraceae", ]
dim(Firm_of_int)
###############################################################################################################
xFirm1 <- major_taxa_for_plot.g2 %>% 
  group_by(MTaxa, Patient) %>% 
  summarize(mean_relAb = (mean(Proportion, na.rm = TRUE)))

xFirm <- xFirm1 %>% 
  group_by(MTaxa) %>% 
  summarize(mean_relAb = (mean(mean_relAb, na.rm = TRUE)))
sum(xFirm$mean_relAb)

xFirm$PhylaPct <- ((xFirm$mean_relAb)/(sum(xFirm$mean_relAb)))*100
xFirm

FigS6E <- ggplot(xFirm, aes(x="", y=PhylaPct, fill=MTaxa)) +
  geom_bar(stat="identity", width=1, color = "white") +
  coord_polar("y", start=0) + theme_void() + scale_fill_viridis(discrete=TRUE, name = "Family")+ 
  ggtitle("Firmicutes") + theme(text = element_text(size = 16))

FigS6E
###############################################################################################################
#subset Bacteroidetes
Bact_prot = subset_taxa(ASV_physeq, Phylum=="Bacteroidota")
count_tab <- as.data.frame(otu_table(Bact_prot))

phyla_counts_tab <- otu_table(tax_glom(Bact_prot, taxrank="Family"))

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(Bact_prot, taxrank="Family"))[,5]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

#turn NAs into their own group and add this row to our Phylum count table:
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(phyla_and_unidentified_counts_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)

#check the dimensions of this table at this point
dim(phyla_and_unidentified_counts_tab)

# here, we'll only keep rows (taxa) that make up greater than 5% in any
temp_filt_major_taxa_proportions_tab <- as.data.frame(phyla_and_unidentified_counts_tab[apply(phyla_and_unidentified_counts_tab, 1, max) > 10, ])

# checking how many we have that were above this threshold
dim(temp_filt_major_taxa_proportions_tab) 
rownames(temp_filt_major_taxa_proportions_tab)

# though each of the filtered taxa made up less than 2% alone, together they
# may add up and should still be included in the overall summary
# so we're going to add a row called "Other" that keeps track of how much we
# filtered out (which will also keep our totals at 100%)
filtered_proportions <- colSums(phyla_and_unidentified_counts_tab) - colSums(temp_filt_major_taxa_proportions_tab)
filt_major_taxa_proportions_tab <- rbind(temp_filt_major_taxa_proportions_tab, "Other"=filtered_proportions)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(filt_major_taxa_proportions_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
head(major_taxa_for_plot.g2)

####################################################################################################################
#make a new sheet with just the taxa of interest
Bact_of_int <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Bacteroidaceae", ]
dim(Bact_of_int)
###############################################################################################################
xBact1 <- major_taxa_for_plot.g2 %>% 
  group_by(MTaxa, Patient) %>% 
  summarize(mean_relAb = (mean(Proportion, na.rm = TRUE)))

xBact <- xBact1 %>% 
  group_by(MTaxa) %>% 
  summarize(mean_relAb = (mean(mean_relAb, na.rm = TRUE)))
sum(xBact$mean_relAb)

xBact$PhylaPct <- ((xBact$mean_relAb)/(sum(xBact$mean_relAb)))*100
xBact

FigS6A <- ggplot(xBact, aes(x="", y=PhylaPct, fill=MTaxa)) +
  geom_bar(stat="identity", width=1, color = "white") +
  coord_polar("y", start=0) + theme_void() + scale_fill_viridis(discrete=TRUE, name = "Family")+ 
  ggtitle("Bacteroidetes") + theme(#legend.position = "bottom",
                                   text = element_text(size = 16))

FigS6A
###############################################################################################################
#subset Verrucomicrobiota
Verr_prot = subset_taxa(ASV_physeq, Phylum=="Verrucomicrobiota")
count_tab <- as.data.frame(otu_table(Verr_prot))

phyla_counts_tab <- otu_table(tax_glom(Verr_prot, taxrank="Family"))

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(Verr_prot, taxrank="Family"))[,5]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

#turn NAs into their own group and add this row to our Phylum count table:
unclassified_tax_counts <- colSums(count_tab) - colSums(phyla_counts_tab)
phyla_and_unidentified_counts_tab <- rbind(phyla_counts_tab, "Unclassified"=unclassified_tax_counts)

#make a copy of our table that's safe for manipulating
major_taxa_for_plot <- as.data.frame(phyla_and_unidentified_counts_tab)

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

# take a look at the new table
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_tab$Age_Days,"Sample"=row.names(sample_info_tab), "Age"=sample_info_tab$Bin_Months, "Patient"=sample_info_tab$Patient, "Year" =sample_info_tab$Bin_Years, stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
head(major_taxa_for_plot.g2)

####################################################################################################################
#make a new sheet with just the taxa of interest
Verr_of_int <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Akkermansiaceae", ]
dim(Verr_of_int)
###############################################################################################################
xVerr1 <- major_taxa_for_plot.g2 %>% 
  group_by(MTaxa, Patient) %>% 
  summarize(mean_relAb = (mean(Proportion, na.rm = TRUE)))

xVerr <- xVerr1 %>% 
  group_by(MTaxa) %>% 
  summarize(mean_relAb = (mean(mean_relAb, na.rm = TRUE)))
sum(xVerr$mean_relAb)

xVerr$PhylaPct <- ((xVerr$mean_relAb)/(sum(xVerr$mean_relAb)))*100
xVerr

FigS6B <- ggplot(xVerr, aes(x="", y=PhylaPct, fill=MTaxa)) +
  geom_bar(stat="identity", width=1, color = "white") +
  coord_polar("y", start=0) + theme_void() + scale_fill_viridis(discrete=TRUE, name = "Family")+ 
  ggtitle("Verrucomicrobiota") + theme(#legend.position = "bottom", 
                                       text = element_text(size = 16))
FigS6B
###############################################################################################################
dim(Prot_of_int)
dim(Verr_of_int)
dim(Bact_of_int)
dim(Act_of_int)
dim(Firm_of_int)

#this is how you made the faceted line plots subset on specific phyla 
z <- rbind(Act_of_int, Bact_of_int, Firm_of_int, Prot_of_int, Verr_of_int)

#reorder levels; this makes the plot in the right order 
z$MTaxa <- factor(z$MTaxa, levels = c("Bifidobacteriaceae", "Bacteroidaceae", "Lachnospiraceae", 
                                  "Enterobacteriaceae", "Akkermansiaceae"))

#ribbon plot
Fig3B <- ggplot(z, aes(x=Age_Days, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), size=2) + theme_bw() +
  labs(x = "Age (Days)", y = "Relative Abundance (%)") +
  scale_color_viridis(discrete=TRUE, name = "Phylum", 
                      labels = c("Actinobacteria",
                                 "Bacteroidetes",
                                 "Firmicutes",
                                 "Proteobacteria",
                                 "Verrucomicrobiota")) +
  facet_wrap(.~MTaxa) + theme(axis.text = element_text(size = 10),
                              axis.title=element_text(size=14, face="bold"),
                              legend.text=element_text(size=10), 
                              legend.title=element_text(size=14, face="bold"),
                              legend.position = c(.85, 0.25),
                              plot.margin=unit(c(0.2,.6,0.2,0.2),"cm"),
                              panel.spacing = unit(.6, "lines"),
                              legend.key.size = unit(0.7, 'cm'))

Fig3B  


tiff("Figure3.tiff", units="in", width=6, height=9, res=600)
ggarrange(Fig3A,Fig3B, labels = c("A.","B."), ncol=1)
dev.off()

