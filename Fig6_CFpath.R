#load packages 
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

#load phyloseq object and metadata
load("ASV_physeq_silva.rda")
load("sample_info.rda")

#########################################################################################################
##CF Pathogen Relative abundance##
ASV_physeq = transform_sample_counts(ASV_physeq, function(x) (x / sum(x))*100 )

#pull out only taxa corresponding to CF pathogens 
CF_pathogens = subset_taxa(ASV_physeq, Genus=="Pseudomonas" |
                             Genus=="Staphylococcus" | 
                             Genus=="Veillonella" |
                             Genus == "Prevotella" |
                             Genus == "Prevotella_7" |
                             Genus == "Prevotella_9" |
                             Genus == "Staphylococcus"|
                             Genus == "Neisseria"|
                             Genus == "Haemophilus"|
                             Genus == "Gemella" |
                             Genus == "Achromobacter" |
                             Genus == "Stenotrophomonas" |
                             Genus == "Burkholderia" |
                             Genus == "Streptococcus")

#merge ASVs by Genus
phyla_counts_tab <- otu_table(tax_glom(CF_pathogens, taxrank="Genus"))

# making a vector of Family names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(CF_pathogens, taxrank="Genus"))[,6])
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
                                  "Bin" = sample_info_tab$Bin_EarlyvLate,
                                  stringsAsFactors=F)

# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)

#Add a prevalence column
major_taxa_for_plot.g2$Prevalence[major_taxa_for_plot.g2$Proportion > 0] <- 1
major_taxa_for_plot.g2$Prevalence[major_taxa_for_plot.g2$Proportion == 0] <- 0
############################################################################################################
#Table S3A
#Calculate overall proportions & prevalence
#average by patient first 
All_prop_bp <- major_taxa_for_plot.g2 %>%
  group_by(MTaxa, Patient) %>%
  summarise(
    Proportion = mean(Proportion),
    Prevalence = mean(Prevalence)
  )

#average by Taxa
All_prop <- All_prop_bp %>%
  group_by(MTaxa) %>%
  summarise(
    sd = sd(Proportion, na.rm = TRUE),
    Proportion = mean(Proportion),
    Prevalence = mean(Prevalence*100)
  )
All_prop

write.csv(All_prop, "CFpath_all_summary.csv")
############################################################################################################
#Table S3B, Figure 6A-B
#Calculate proportion & prevalence, grouped by 3-month age bins
#average by patient first
Age_prop_bp <- major_taxa_for_plot.g2 %>%
  group_by(Age, Patient, MTaxa) %>%
  summarise(
    Proportion = mean(Proportion),
    Prevalence = mean(Prevalence)
  )

#average by taxa
Age_prop <- Age_prop_bp %>%
  group_by(Age, MTaxa) %>%
  summarise(
    sd = sd(Proportion, na.rm = TRUE),
    Proportion = mean(Proportion),
    Prevalence = mean(Prevalence*100)
  )
Age_prop

write.csv(Age_prop, "CFpath_age_summary.csv")
####################################################################################################################
#subset to keep taxa of interest
Veillonella <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Veillonella", ]
head(Veillonella)

#ribbon plot
Veillonella_1 <- ggplot(Veillonella, aes(x=Age_Days, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), size =.6, show.legend = FALSE) + theme_bw() +
  geom_point(size=1) +
  ggtitle("Veillonella") +
  labs(x = "Age (Days)", y = "Relative Abundance (%)") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold.italic", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))
Veillonella_1

#relative abundance plot 
VeillonellaEvL <- Veillonella[Veillonella$Bin != "Middle", ]

#makes a summary sheet that include the mean and SD from the sheet of
#Proportion, grouped by Bin
Veillonella_sum_bp <- VeillonellaEvL %>%
  group_by(Bin, Patient) %>%
  summarise(
    Proportion = mean(Proportion)
  )

Veillonella_sum <- Veillonella_sum_bp %>%
  group_by(Bin) %>%
  summarise(
    sd = sd(Proportion, na.rm = TRUE),
    Proportion = mean(Proportion)
  )

Veillonella_sum 

#makes a plot with bars, standard deviation, and points 
Veillonella_2 <- ggplot(VeillonellaEvL, aes(Bin, Proportion)) +
  geom_col(data = Veillonella_sum, fill = NA, color = "black") +  #adds the bar graph
  geom_jitter( position = position_jitter(0.1), color = "black") + #makes the dots more visible
  geom_errorbar( aes(ymin = Proportion, ymax = Proportion+sd), #adds erros bars
                 data = Veillonella_sum, width = 0.2) +
  theme_bw() + #makes it less ugly 
  ggtitle("") + #adds a title (empty)
  labs(x = "Age", y = "Relative Abundance (%)") +
  annotate(geom="text",x=1.5, y=37.5, label="p = 1.09e-6") +
  theme(axis.text=element_text(size=8), #phont sizes & face
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))+
  scale_x_discrete(labels=c("<6 months", "2+ years")) #labels bars
Veillonella_2

tiff('Fig6C.tiff', units="in", width=6.5, height=4.2, res=300)
ggarrange(Veillonella_1,Veillonella_2)
dev.off()

####################################################################################################################
#subset to keep taxa of interest
Streptococcus <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa == "Streptococcus", ]
head(Streptococcus)

#ribbon plot
Streptococcus_1 <- ggplot(Streptococcus, aes(x=Age_Days, y=Proportion)) + 
  geom_smooth(aes(color = MTaxa), size =.6, show.legend = FALSE) + theme_bw() +
  geom_point(size=1) +
  ggtitle("Streptococcus") +
  labs(x = "Age (Days)", y = "Relative Abundance (%)") +
  theme(axis.text=element_text(size=8), 
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold.italic", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))
Streptococcus_1

#relative abundance plot 
StreptococcusEvL <- Streptococcus[Streptococcus$Bin != "Middle", ]

#makes a summary sheet that include the mean and SD from the sheet of
#Proportion, grouped by Bin
Streptococcus_sum_bp <- StreptococcusEvL %>%
  group_by(Bin, Patient) %>%
  summarise(
    Proportion = mean(Proportion)
  )

Streptococcus_sum <- Streptococcus_sum_bp %>%
  group_by(Bin) %>%
  summarise(
    sd = sd(Proportion, na.rm = TRUE),
    Proportion = mean(Proportion)
  )

Streptococcus_sum 

#makes a plot with bars, standard deviation, and points 
Streptococcus_2 <- ggplot(StreptococcusEvL, aes(Bin, Proportion)) +
  geom_col(data = Streptococcus_sum, fill = NA, color = "black") +  #adds the bar graph
  geom_jitter( position = position_jitter(0.1), color = "black") + #makes the dots more visible
  geom_errorbar( aes(ymin = Proportion, ymax = Proportion+sd), #adds erros bars
                 data = Streptococcus_sum, width = 0.2) +
  theme_bw() + #makes it less ugly 
  ggtitle("") + #adds a title (empty)
  labs(x = "Age", y = "Relative Abundance (%)") + #labels axes
  annotate(geom="text",x=1.5, y=37.5, label="ns") +
  theme(axis.text=element_text(size=8), #phont sizes & face
        axis.title=element_text(size=12,face="bold"), 
        plot.title = element_text(face = "bold", size = 14),
        plot.margin = margin(.4,.4,.4,.4, "cm"))+
  scale_x_discrete(labels=c("<6 months", "2+ years")) #labels bars
Streptococcus_2

tiff('Fig6D.tiff', units="in", width=6.5, height=4.2, res=300)
ggarrange(Streptococcus_1, Streptococcus_2)
dev.off()
