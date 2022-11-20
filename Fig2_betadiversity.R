library("phyloseq")
library("ggplot2")
library("plyr")
library("DESeq2")
library("vegan")
library("viridis")
library("ggpubr")

packageVersion("phyloseq") #1.32.0
packageVersion("ggplot2") #3.3.6
packageVersion("plyr") #1.8.7
packageVersion("DESeq2") #1.28.1
packageVersion("vegan") #2.5.7
packageVersion("viridis") #0.6.2
packageVersion("ggpubr") #0.4.0

##load in phyloseq data 
load("ASV_physeq_silva.rda")
load("sample_info.rda")
load("count_tab.rda")

#set seed for reproducability 
set.seed(11271989)
########Pre-processing######################################################################################################################################
#Preprocessing for both ordination & DESEq
#Most of this is visualizations, only pre-processing is removing ASVs with <1% prevalence

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
# This removed 2190 taxa; we have a lot of low count taxa! 
keepTaxa = rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ASV_physeq = prune_taxa(keepTaxa, ASV_physeq)

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
  facet_wrap(~Phylum) + theme(legend.position="none") +
  theme_bw() + guides(color = FALSE, size = FALSE)

########Figure2A all######################################################################################################################################
#tests the contribution and significance of each variable after
#accounting for the contribution of all other variables 

#Remove unknowns
ASV_physeq_sm <- subset_samples(ASV_physeq, !(Gender == "Unknown"))
ASV_physeq_sm <- subset_samples(ASV_physeq_sm, !(preterm == "Unknown"))
ASV_physeq_sm <- subset_samples(ASV_physeq_sm, !(Pancreatic_sufficiency == "Unknown"))
ASV_physeq_sm <- subset_samples(ASV_physeq_sm, !(delivery == "Unknown"))

#bray distance
dist = phyloseq::distance(ASV_physeq_sm, method="bray")

#test Age age
adon.age <- adonis(dist ~ sample_data(ASV_physeq_sm)$Gender +
                     sample_data(ASV_physeq_sm)$Pancreatic_sufficiency +
                     sample_data(ASV_physeq_sm)$breastfed_ever +
                     sample_data(ASV_physeq_sm)$RecentAbx +
                     sample_data(ASV_physeq_sm)$preterm +
                     sample_data(ASV_physeq_sm)$delivery +
                     sample_data(ASV_physeq_sm)$RecentSaSwab +
                     sample_data(ASV_physeq_sm)$RecentPaSwab +
                     sample_data(ASV_physeq_sm)$Age_Days,
                   strata = sample_data(ASV_physeq_sm)$Patient)
age <- as.data.frame(adon.age$aov.tab)[9,]

#test Gender
adon.gender <- adonis(dist ~ sample_data(ASV_physeq_sm)$Age_Days +
                        sample_data(ASV_physeq_sm)$Pancreatic_sufficiency +
                        sample_data(ASV_physeq_sm)$breastfed_ever +
                        sample_data(ASV_physeq_sm)$RecentAbx +
                        sample_data(ASV_physeq_sm)$preterm +
                        sample_data(ASV_physeq_sm)$delivery +
                        sample_data(ASV_physeq_sm)$RecentSaSwab +
                        sample_data(ASV_physeq_sm)$RecentPaSwab +
                        sample_data(ASV_physeq_sm)$Gender,
                      strata = sample_data(ASV_physeq_sm)$Patient)
gender <- as.data.frame(adon.gender$aov.tab)[9,]

#test PS
adon.PI <- adonis(dist ~ sample_data(ASV_physeq_sm)$Age_Days +
                    sample_data(ASV_physeq_sm)$Gender +
                    sample_data(ASV_physeq_sm)$breastfed_ever +
                    sample_data(ASV_physeq_sm)$RecentAbx +
                    sample_data(ASV_physeq_sm)$preterm +
                    sample_data(ASV_physeq_sm)$delivery +
                    sample_data(ASV_physeq_sm)$RecentSaSwab +
                    sample_data(ASV_physeq_sm)$RecentPaSwab +
                    sample_data(ASV_physeq_sm)$Pancreatic_sufficiency,
                  strata = sample_data(ASV_physeq_sm)$Patient)
PI <- as.data.frame(adon.PI$aov.tab)[9,]

#test breast feeding
adon.BF <- adonis(dist ~ sample_data(ASV_physeq_sm)$Age_Days +
                    sample_data(ASV_physeq_sm)$Gender +
                    sample_data(ASV_physeq_sm)$Pancreatic_sufficiency +
                    sample_data(ASV_physeq_sm)$RecentAbx +
                    sample_data(ASV_physeq_sm)$preterm +
                    sample_data(ASV_physeq_sm)$delivery +
                    sample_data(ASV_physeq_sm)$RecentSaSwab +
                    sample_data(ASV_physeq_sm)$RecentPaSwab +
                    sample_data(ASV_physeq_sm)$breastfed_ever,
                  strata = sample_data(ASV_physeq_sm)$Patient)
BF <- as.data.frame(adon.BF$aov.tab)[9,]

#test total abx prior to sample collection 
adon.abx <- adonis(dist ~ sample_data(ASV_physeq_sm)$Age_Days +
                     sample_data(ASV_physeq_sm)$Gender +
                     sample_data(ASV_physeq_sm)$Pancreatic_sufficiency +
                     sample_data(ASV_physeq_sm)$breastfed_ever +
                     sample_data(ASV_physeq_sm)$preterm +
                     sample_data(ASV_physeq_sm)$delivery +
                     sample_data(ASV_physeq_sm)$RecentSaSwab +
                     sample_data(ASV_physeq_sm)$RecentPaSwab +
                     sample_data(ASV_physeq_sm)$RecentAbx,
                   strata = sample_data(ASV_physeq_sm)$Patient)
abx <- as.data.frame(adon.abx$aov.tab)[9,]

#test preterm
adon.preterm <- adonis(dist ~ sample_data(ASV_physeq_sm)$Age_Days +
                         sample_data(ASV_physeq_sm)$Gender +
                         sample_data(ASV_physeq_sm)$Pancreatic_sufficiency +
                         sample_data(ASV_physeq_sm)$breastfed_ever +
                         sample_data(ASV_physeq_sm)$RecentAbx +
                         sample_data(ASV_physeq_sm)$delivery +
                         sample_data(ASV_physeq_sm)$RecentSaSwab +
                         sample_data(ASV_physeq_sm)$RecentPaSwab +
                         sample_data(ASV_physeq_sm)$preterm,
                       strata = sample_data(ASV_physeq_sm)$Patient)
preterm <- as.data.frame(adon.preterm$aov.tab)[9,]

#test delivery
adon.delivery <- adonis(dist ~ sample_data(ASV_physeq_sm)$Age_Days +
                          sample_data(ASV_physeq_sm)$Gender +
                          sample_data(ASV_physeq_sm)$Pancreatic_sufficiency +
                          sample_data(ASV_physeq_sm)$breastfed_ever +
                          sample_data(ASV_physeq_sm)$RecentAbx +
                          sample_data(ASV_physeq_sm)$preterm +
                          sample_data(ASV_physeq_sm)$RecentSaSwab +
                          sample_data(ASV_physeq_sm)$RecentPaSwab +
                          sample_data(ASV_physeq_sm)$delivery,
                        strata = sample_data(ASV_physeq_sm)$Patient)
delivery <- as.data.frame(adon.delivery$aov.tab)[9,]

#test Sa
adon.Sa <- adonis(dist ~ sample_data(ASV_physeq_sm)$Age_Days +
                    sample_data(ASV_physeq_sm)$Gender +
                    sample_data(ASV_physeq_sm)$Pancreatic_sufficiency +
                    sample_data(ASV_physeq_sm)$breastfed_ever +
                    sample_data(ASV_physeq_sm)$RecentAbx +
                    sample_data(ASV_physeq_sm)$preterm +
                    sample_data(ASV_physeq_sm)$delivery +
                    sample_data(ASV_physeq_sm)$RecentPaSwab +
                    sample_data(ASV_physeq_sm)$RecentSaSwab,
                  strata = sample_data(ASV_physeq_sm)$Patient)
Sa <- as.data.frame(adon.Sa$aov.tab)[9,]

#test Recent Pa
adon.Pa <- adonis(dist ~ sample_data(ASV_physeq_sm)$Age_Days +
                    sample_data(ASV_physeq_sm)$Gender +
                    sample_data(ASV_physeq_sm)$Pancreatic_sufficiency +
                    sample_data(ASV_physeq_sm)$breastfed_ever +
                    sample_data(ASV_physeq_sm)$RecentAbx +
                    sample_data(ASV_physeq_sm)$preterm +
                    sample_data(ASV_physeq_sm)$delivery +
                    sample_data(ASV_physeq_sm)$RecentSaSwab +
                    sample_data(ASV_physeq_sm)$RecentPaSwab)
Pa <- as.data.frame(adon.Pa$aov.tab)[9,]

################################################################################
#makes a graph of adonis (permanova) data 
#combine statistical data into one dataframe 
adon.all <- rbind(Sa,Pa, delivery, preterm, abx, BF, gender, PI, age)

#give each condition better names 
adon.all$condition <- c("OP S. aureus", "OP P. aeruginosa", "Delivery Mode", 
                        "Gestation", 
                        "Recent Antibiotics", 
                        "Breastfeeding", "Sex", 
                        "PS", "Age")

adon.all$Significance <- adon.all$`Pr(>F)`
adon.all$Significance[adon.all$Significance < 0.05] <- "Significant"
adon.all$Significance[adon.all$Significance < 0.1] <- "Marginal Significance, p <0.1"
adon.all$Significance[adon.all$Significance < 1] <- "Not Significant"


Fig2A <- ggplot(adon.all, aes(x = reorder(condition, R2), y = R2, fill = Significance)) + 
  geom_bar(stat = "identity", color = "black") +
  coord_flip() +
  theme_bw() +
  scale_fill_viridis(discrete=TRUE) +
  labs(y = "R-squared", x ="") +
  scale_fill_viridis(discrete = TRUE) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=16), legend.title=element_text(size=16, face="bold"),
        plot.margin = margin(.3, .3, .3, .8, "cm"))
Fig2A

tiff("Fig2A_Resub.tiff", units="in", width=8.5, height=7, res=600)
Fig2A
dev.off()

########Figure2B age######################################################################################################################################
#Figure 1B & Figure 1 export
#Note that Fig1A object must be present to export full figure; this was generated in Fig1A_S1A_alphadiversity.R 

#ordination 
vst_pcoa <- ordinate(ASV_physeq_sm, method="MDS", distance="bray")

#plotting
Fig2B <- plot_ordination(ASV_physeq_sm, vst_pcoa, color = "Bin_Years") + 
  geom_point(size=2.5) + 
  stat_ellipse(type = "t") +
  theme_bw() + 
  scale_color_viridis(discrete=TRUE, labels = c("0-1", ">1-2", ">2-3",">3-4")) +
  labs(col= "Age (Years)") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=16), legend.title=element_text(size=16, face="bold"),
        plot.margin = margin(.3, .3, .3, .8, "cm"))
Fig2B

tiff("Fig2B_Resub.tiff", units="in", width=8.5, height=7, res=600)
Fig2B
dev.off()
########Age Statistics######################################################################################################################################
##PERMANOVA##
#Statistics by age bins
#dummy variable 
sample_data(ASV_physeq_sm)$Bin_Years_0 <- ifelse(sample_data(ASV_physeq_sm)$Bin_Years == "0",1,0)
sample_data(ASV_physeq_sm)$Bin_Years_1 <- ifelse(sample_data(ASV_physeq_sm)$Bin_Years == "1",1,0)
sample_data(ASV_physeq_sm)$Bin_Years_2 <- ifelse(sample_data(ASV_physeq_sm)$Bin_Years == "2",1,0)
sample_data(ASV_physeq_sm)$Bin_Years_3 <- ifelse(sample_data(ASV_physeq_sm)$Bin_Years == "3",1,0)

#year 0 vs others

#year 0 is significantly different from yr1 when adjusting for others (p=0.001)
adonis(dist ~ sample_data(ASV_physeq_sm)$Bin_Years_2 +
         sample_data(ASV_physeq_sm)$Bin_Years_3 +
         sample_data(ASV_physeq_sm)$Bin_Years_1 +
         sample_data(ASV_physeq_sm)$Bin_Years_0,
       strata = sample_data(ASV_physeq_sm)$Patient)

#year 0 is significantly different from yr2 when adjusting for others (p=0.001)
adonis(dist ~ sample_data(ASV_physeq_sm)$Bin_Years_1 +
         sample_data(ASV_physeq_sm)$Bin_Years_3 +
         sample_data(ASV_physeq_sm)$Bin_Years_2 +
         sample_data(ASV_physeq_sm)$Bin_Years_0,
       strata = sample_data(ASV_physeq_sm)$Patient)

#yr 0 is different from 3 when adjusting for others (p = 0.001)
adonis(dist ~ sample_data(ASV_physeq_sm)$Bin_Years_1 +
         sample_data(ASV_physeq_sm)$Bin_Years_2 +
         sample_data(ASV_physeq_sm)$Bin_Years_3 +
         sample_data(ASV_physeq_sm)$Bin_Years_0,
       strata = sample_data(ASV_physeq_sm)$Patient)

#year 1 vs others

#year 1 is not significantly different from yr2 when adjusting for others (p=0.255)
adonis(dist ~ sample_data(ASV_physeq_sm)$Bin_Years_0 +
         sample_data(ASV_physeq_sm)$Bin_Years_3 +
         sample_data(ASV_physeq_sm)$Bin_Years_2 +
         sample_data(ASV_physeq_sm)$Bin_Years_1,
       strata = sample_data(ASV_physeq_sm)$Patient)

#yr1 is different from 3 when adjusting for others (p=0.008)
adonis(dist ~ sample_data(ASV_physeq_sm)$Bin_Years_0 +
         sample_data(ASV_physeq_sm)$Bin_Years_2 +
         sample_data(ASV_physeq_sm)$Bin_Years_3 +
         sample_data(ASV_physeq_sm)$Bin_Years_1,
       strata = sample_data(ASV_physeq_sm)$Patient)

#year 2 vs others

#yr2 is not different from 3 when adjusting for others (p = 0.184)
adonis(dist ~ sample_data(ASV_physeq_sm)$Bin_Years_0 +
         sample_data(ASV_physeq_sm)$Bin_Years_1 +
         sample_data(ASV_physeq_sm)$Bin_Years_3 +
         sample_data(ASV_physeq_sm)$Bin_Years_2,
       strata = sample_data(ASV_physeq_sm)$Patient)

########Heatmap#########################################################################
#Remove unknowns

#convert everything to to % relative abundance
ASV_physeq_sm_relab = transform_sample_counts(ASV_physeq_sm, function(x) (x / sum(x))*100 )

ASV_physeq_gen <- tax_glom(ASV_physeq_sm_relab, taxrank="Genus")
ASV_physeq_gen <- prune_taxa(names(sort(taxa_sums(ASV_physeq_gen),TRUE)[1:10]), ASV_physeq_gen)

x <- as.data.frame(otu_table(ASV_physeq_gen))
y <- as.data.frame(sample_data(ASV_physeq_gen))
z <- as.data.frame(tax_table(ASV_physeq_gen))

identical(row.names(x), row.names(z))
row.names(x) <- z$Genus

# Define colors for each levels of qualitative variables
# Define gradient color for continuous variable (mpg)
col = list(Sex = c("Male" = "black", "Female" = "gray"), 
           PS = c("Insufficient" = "black", "Sufficient" = "gray"),
           Delivery = c("Csection" = "black", "Vaginal" = "gray"), 
           BF = c("Yes" = "black", "No" = "gray"),
           Age =c("0" = "#440154", "1" = "#31688e", "2" = "#35b779", "3" = "#fde725"))

library(ComplexHeatmap)
# Create the heatmap annotation
ha <- HeatmapAnnotation(Sex = y$Gender,
                        PS = y$Pancreatic_sufficiency,
                        Delivery = y$delivery,
                        BF = y$breastfed_ever,
                        Age = y$Bin_Years, col = col)



# Combine the heatmap and the annotation
Fig2C <- Heatmap(as.matrix(x), name = "% Rel Abund", bottom_annotation = ha, 
        show_column_names = FALSE)
Fig2C

tiff("Fig2C_Resub.tiff", units="in", width=12, height=6, res=600)
Fig2C
dev.off()
