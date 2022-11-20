library("dplyr")
library("reshape")
library("dendextend")
library("nlme")
library("tidyr")
library("ggplot2")
library("viridis")
library("ggpubr")
library("phyloseq")

#####################################################################################################################################################################
#load important files
load("ASV_physeq_silva.rda")
load("sample_info.rda")

#####################################################################################################################################################################
#Recording alpha-diversity values in the metadata
divIndex <- estimate_richness(ASV_physeq, measures = c("Simpson", "Shannon", "Fisher", "InvSimpson"))

sample_info_tab$Sample <- row.names(sample_info_tab)
divIndex$Sample <- row.names(divIndex)
sample_info_new <- merge(sample_info_tab, divIndex)
head(sample_info_new)

#####################################################################################################################################################################
#pull out only taxa corresponding to Bacteroidetes & Firmicutes
Bt_Fi = subset_taxa(ASV_physeq, Phylum=="Bacteroidota" |
                             Phylum=="Firmicutes")


#merge ASVs by Phylum 
phyla_counts_tab <- otu_table(tax_glom(Bt_Fi, taxrank="Phylum"))

# making a vector of family names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(Bt_Fi, taxrank="Phylum"))[,2])
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)
dim(phyla_counts_tab)

##Making the Ratio Table##
#transpose the count table and change to a data frame
phyla_counts <- as.data.frame(t(phyla_counts_tab))
#calculate Bacteroidetes/Firmicutes ratio & re-transpose the table 
phyla_counts$BtoFratio <- (phyla_counts$Firmicutes)/(phyla_counts$Bacteroidota)
major_taxa_for_plot <- as.data.frame(t(phyla_counts))

# add a column of the taxa names so that it is within the table
# transform the table into narrow, or long, format
major_taxa_for_plot$MTaxa <- row.names(major_taxa_for_plot)
major_taxa_for_plot.g <- gather(major_taxa_for_plot, Sample, Proportion, -MTaxa)

#major_taxa_for_plot.g <- major_taxa_for_plot.g[major_taxa_for_plot.g$Sample == "BtoFratio",]
head(major_taxa_for_plot.g)

# now we want a table with metadata for each sample
sample_info_for_merge<-data.frame("Age_Days"=sample_info_new$Age_Days,
                                  "Sample"=sample_info_new$Sample, 
                                  "Age"=sample_info_new$Bin_Months, 
                                  "Patient"=sample_info_new$Patient, 
                                  "Year" =sample_info_new$Bin_Years,
                                  "Shannon" =sample_info_new$Shannon,
                                  stringsAsFactors=F)


# and here we are merging this table with the plotting table we just made
major_taxa_for_plot.g2 <- merge(major_taxa_for_plot.g, sample_info_for_merge)
dim(major_taxa_for_plot.g2)
major_taxa_for_plot.g2 <- major_taxa_for_plot.g2[major_taxa_for_plot.g2$MTaxa =="BtoFratio",]

####################################################################################################################
#Some data prep for doing stats
sample_info <- subset(major_taxa_for_plot.g2, !(Proportion == "Inf"))
sample_info_0 <- sample_info[sample_info$Year == "0",]
sample_info_1 <- sample_info[sample_info$Year == "1",]
sample_info_2 <- sample_info[sample_info$Year == "2",]
sample_info_3 <- sample_info[sample_info$Year == "3",]

#all points, Age is not significantly correlated with B/F ratio 
SDIlme <- lme(Proportion ~ Age_Days, data=sample_info, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #Not sig, 0.65460935 

#Year 0-1
SDIlme <- lme(Proportion ~ Age_Days, data=sample_info_0, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #significant, 0.03970828

#Year 1-2
SDIlme <- lme(Proportion ~ Age_Days, data=sample_info_1, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #not significant, 0.3781653

#Year 2
SDIlme <- lme(Proportion ~ Age_Days, data=sample_info_2, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #not significant, 0.3001484 

#Year 3
SDIlme <- lme(Proportion ~ Age_Days, data=sample_info_3, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #not significant, 0.1686674 

####################################################################################################################

#Fig 2A
Fig2A <- ggplot(sample_info, aes(x = Age_Days, y = Proportion, color = Year)) +
  geom_point(size=2.5) +
  geom_smooth(method='lm',formula=y~x) +
  theme_bw() +
  labs(x="Age (Days)", y="Firmicutes/Bacteroidetes Ratio",color="Age (Years)") +
  scale_color_viridis(discrete=TRUE, labels = c("0-1", ">1-2", ">2-3",">3-4")) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=16), legend.title=element_text(size=16, face="bold"),
        plot.margin = margin(0.3, 0.3, 0.3, 1.2, "cm"))+
  scale_y_continuous(trans='log10')
Fig2A

dat_text <- data.frame(label = c("p = 0.040", "p = 0.378", "p = 0.300", "p = 0.169"), 
                       Year = c("0", "1", "2", "3"), x = c(182.5,547.5,912.5,1277.5), 
                       y = c(0.25, 0.25, 0.25, 0.25))


Fig2A <- Fig2A + geom_text(show.legend = FALSE, 
                           data = dat_text, 
                           mapping = aes(x = x, y = y, label = label))
Fig2A

####################################################################################################################
#all points, Shannon is marginally significantly correlated with B/F ratio 
SDIlme <- lme(Proportion ~ Shannon, data=sample_info, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #Marginally sig, 0.069052124

SDIlme <- lme(Proportion ~ Shannon, data=sample_info_0, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #not significant, 0.3908545

SDIlme <- lme(Proportion ~ Shannon, data=sample_info_1, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #marginally significant, 0.05599482

SDIlme <- lme(Proportion ~ Shannon, data=sample_info_2, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #not significant, 0.1672219 

SDIlme <- lme(Proportion ~ Shannon, data=sample_info_3, random=~1|Patient)
summary(SDIlme)$tTable[,"p-value"] #not significant, 0.6145129 


####################################################################################################################
#graph Fig 2B
Fig2B <- ggplot(sample_info, aes(x = Proportion, y = Shannon, color = Year)) +
  geom_point(size=2.5, show.legend = FALSE) +
  theme_bw() + geom_smooth(method='lm',formula=y~x, show.legend = FALSE) + facet_wrap(~Year,ncol = 4) +
  labs(x="Shannon Diversity", y="Firmicutes/Bacteroidetes Ratio") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16,face="bold"),
        plot.margin = margin(0.3, 1, 0.3, 1.2, "cm"))+
  scale_x_continuous(trans='log10') + 
  scale_color_viridis(discrete=TRUE)
Fig2B


dat_text <- data.frame(label = c("p = 0.391", "p = 0.056", "p = 0.167", "p = 0.615"), 
                       Year = c("0", "1", "2", "3"), x = c(100,100,100,100), 
                       y = c(0.5, 0.5, 0.5, 0.5))

Fig2B <- Fig2B + geom_text(show.legend = FALSE,
                           data = dat_text, mapping = aes(x = x, y = y, label = label))

Fig2B

#######################################################################################################
tiff("Fig2_FBRatio.tiff", units="in", width=7, height=8, res=300)
ggarrange(Fig2A, Fig2B, ncol=1, labels = c("A.","B."),
          font.label = list(size = 20))
dev.off()
#######################################################################################################


