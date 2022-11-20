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
library("viridis")
library("ggeffects")
library("ggplot2")
library("ggpubr")

load("ASV_physeq_silva.rda")
load("count_tab.rda")
load("sample_info.rda")

###############################################################################################################
##plot grouped relative abundances##

#merge taxa at the genus level
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Genus"))

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="Genus"))[,6]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

#keep only the names we want
#these are taxa up in Crohn's 
names = c("Veillonella", "Escherichia-Shigella", "Fusobacterium", "Haemophilus")
phyla_counts_up <- subset(phyla_counts_tab, rownames(phyla_counts_tab) %in% names)
dim(phyla_counts_up)

#these are taxa down in Crohn's
#The following are not included because they will be captured at the family level"
  #Ruminococcaceae -> Faecalibacterium, Oscillospira, Ruminococcus
  #Lachnospiraceae -> Dorea, Coprococcus
names = c("Parabacteroides", "Bilophila","Dialister", "Sutterella", "Bacteroides")

phyla_counts_down <- subset(phyla_counts_tab, rownames(phyla_counts_tab) %in% names)
rownames(phyla_counts_down)
dim(phyla_counts_down)

#do this again at the family level
phyla_counts_tab <- otu_table(tax_glom(ASV_physeq, taxrank="Family"))

# making a vector of phyla names to set as row names
phyla_tax_vec <- as.vector(tax_table(tax_glom(ASV_physeq, taxrank="Family"))[,5]) 
rownames(phyla_counts_tab) <- as.vector(phyla_tax_vec)

names = c("Ruminococcaceae", "Erysipelotrichaceae", "Rikenellaceae",
          "Lachnospiraceae")

phyla_counts_down_family <- subset(phyla_counts_tab, rownames(phyla_counts_tab) %in% names)
rownames(phyla_counts_down_family)
dim(phyla_counts_down_family)
all_counts_down <- rbind(phyla_counts_down,phyla_counts_down_family)
dim(all_counts_down)

MDindex <- log(colSums(phyla_counts_up)/colSums(all_counts_down), 10)
MDindex <- as.data.frame(MDindex)

sample_info_tab <- as.data.frame(sample_data(ASV_physeq))

identical(rownames(sample_info_tab), rownames(MDindex))
sample_info_tab <- sample_info_tab[order(rownames(sample_info_tab)),]
identical(rownames(sample_info_tab),rownames(MDindex))

#Adding sample info to MDindex
#Not doing this the other way around because impossible to convert sample_info_tab to a dataframe (probably an easier way to do this...)
MDindex$Patient <- sample_info_tab$Patient
MDindex$Age_Days <- sample_info_tab$Age_Days
MDindex$Age_Years <- sample_info_tab$Age_Years
MDindex$Bin_Years <- sample_info_tab$Bin_Years
MDindex$Bin_Months <- sample_info_tab$Bin_Months
MDindex$Bin_Months <- sample_info_tab$Bin_Months
MDindex$Sample <- sample_info_tab$Sample
MDindex$Gender <- sample_info_tab$Gender
MDindex$Genotype <- sample_info_tab$Genotype
MDindex$Pancreatic_sufficiency <- sample_info_tab$Pancreatic_sufficiency
MDindex$Preterm <- sample_info_tab$preterm
MDindex$Delivery <- sample_info_tab$delivery
MDindex$Breastfeeding <- sample_info_tab$breastfed_ever
MDindex$TotalAbxPrior <- sample_info_tab$TotalAbxPrior
MDindex$RecentAbx <- sample_info_tab$RecentAbx
MDindex$RecentPa <- sample_info_tab$RecentPaSwab
MDindex$RecentSa <- sample_info_tab$RecentSaSwab


#Calculate diversity and add to MDindex
divIndex <- estimate_richness(ASV_physeq, split = TRUE, measures = c("Shannon", "Simpson"))
identical(rownames(divIndex), rownames(MDindex))

MDindex$Shannon <- divIndex$Shannon

#write csv for making Table S1
write.csv(MDindex, "Table_S1.csv")

sample_info_new <- MDindex

#check new table
head(sample_info_new)
class(sample_info_new)

#individual patient graphs with line at 1 = severe crohn's
#Fig S5A
FigS5A <- ggplot(sample_info_new, aes(x=Age_Days, y=MDindex)) + geom_point() + 
  ylim(-4,4) + xlim(0,1500)+ theme_bw() + 
  geom_hline(yintercept=1, linetype="dashed", color = "red")+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=16,face="bold"),
        panel.spacing = unit(.5,"lines"))+
  labs(y="Crohn's Dysbiosis Index") +
  facet_wrap(~Patient)
FigS5A

tiff("FigS4A_CDIbyPatient.tiff", units="in", width=10, height=6, res=300)
FigS5A
dev.off()

#################################################################################################
#nlme for analyzing correlation between SDI and age in days
library("nlme")
sample_info <- subset(sample_info_new, !(MDindex == "-Inf")) 
sample_info <- subset(sample_info, !(MDindex == "Inf")) 

#Linear model stats for Fig 6A
SDIlme <- lme(Shannon ~ MDindex, data=sample_info, random=~1|Patient)
summary(SDIlme)
summary(SDIlme)$tTable[,"p-value"] #2.388022e-22

#Linear model stats for Fig 6B
SDIlme <- lme(Age_Years ~ MDindex, data=sample_info, random=~1|Patient)
summary(SDIlme)
summary(SDIlme)$tTable[,"p-value"] #5.904115e-15

write.csv(sample_info_new, "MDindex.csv")

#########################################################################################
#Figures with correlations 
#Shannon vs Crohn's
#Fig 6A
Fig7A <- ggplot(sample_info_new, aes(x=Shannon, y=MDindex)) +
  geom_point(aes(colour=Bin_Years)) + 
  theme_bw() + 
  scale_color_viridis(discrete = TRUE, name = "Age (Years)",
                      labels = c("0-1",">1-2",">2-3","3+")) +
  labs(x = "Shannon Diversity Index", y = "Crohn's Dysbiosis Index") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=12, face="bold"))+
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  annotate(geom="text",x=1.5, y=-3, label="p = 2.39e-22")
Fig7A

#Age vs. MD index for all patients with line at 1=severe crohn's
#Figure 6B
Fig7B <- ggplot(sample_info_new, aes(x=Age_Days, y=MDindex)) + 
  geom_point(aes(color = Bin_Years)) + 
  ylim(-4,4) + xlim(0,1500)+ 
  theme_bw() + 
  geom_hline(yintercept=1, linetype="dashed", color = "red") +
  labs(x = "Age (Days)", y = "Crohn's Dysbiosis Index") +
  theme(axis.text=element_text(size=12), 
        axis.title=element_text(size=12,face="bold"), 
        legend.text=element_text(size=10), 
        legend.title=element_text(size=12, face="bold")) +
  scale_color_viridis(discrete=TRUE, name = "Age (Years)", labels = c("0-1",">1-2",">2-3","3+"))+
  annotate(geom="text",x=300, y=-3, label="p = 5.90e-15")
Fig7B

#########################################################################################
#Use exported file MDIndex.csv to manually chose samples
#For Early High CDI ('Yes') vs No Early High CDI ('No')
#1. Patient had samples in both <2year and 2+ year range
#2. Binned sample by whether patient had any early samples w/CD index >0.5
#3. Chose the latest available sample for each patient so there is only 1 sample per patient 
#########################################################################################
set.seed(11271989)
#Read new file back in 
compMD <- read.csv("MDindex_late_short_v2.csv")

#Significance testing by t-tests
t.test(Shannon ~ EarlyCDHi, data = compMD) #Not sig, p=0.3343
t.test(MDindex ~ EarlyCDHi, data = compMD) #Not sig, p=0.1251
t.test(Age_Years ~ EarlyCDHi, data = compMD) #Not sig, p=0.1825

t.test(Age_Years ~ Gender, data = compMD) #Sig p = 0.01112
t.test(Age_Years ~ breastfed_ever, data = compMD) #NS p = 0.8843
t.test(Age_Years ~ delivery, data = compMD) #NS p = 0.07631
t.test(Age_Years ~ preterm, data = compMD) #NS p = 0.6651
t.test(Age_Years ~ Pancreatic_sufficiency, data = compMD) #NS p = 0.7843

chisq.test(compMD$Gender, compMD$EarlyCDHi, 
           simulate.p.value = TRUE) #Sig, p=0.02599, but throws a warning due to small sample size #v2, p=1.397e-6

chisq.test(compMD$breastfed_ever, compMD$EarlyCDHi, 
           simulate.p.value = TRUE) #p = 0.4128

chisq.test(compMD$delivery, compMD$EarlyCDHi, 
           simulate.p.value = TRUE) #p = 0.4058

chisq.test(compMD$preterm, compMD$EarlyCDHi, 
           simulate.p.value = TRUE) #p = 0.4128

chisq.test(compMD$Pancreatic_sufficiency, compMD$EarlyCDHi, 
           simulate.p.value = TRUE) #ns p = 0.2224


#Figures
FigS4B <- ggplot(compMD, aes(x=EarlyCDHi, y=Shannon)) + geom_boxplot() +
  labs(x="Early Dysbiosis", y="Shannon Diversity Index") + theme_bw() + ylim(0, 4) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) +
  annotate(geom="text",x=1.5, y=1.25, label="p = 0.3343")
FigS4B



FigS4C <- ggplot(compMD, aes(x=EarlyCDHi, fill=Gender)) + geom_bar(position = "fill") +
  labs(x="Early Dysbiosis", y="Proportion") + 
  theme_bw() + 
  scale_fill_viridis(discrete=TRUE, name = "Sex") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14)) 
FigS4C

#breastfeeding
FigS4D <- ggplot(compMD, aes(x=EarlyCDHi, fill=breastfed_ever)) + geom_bar(position = "fill") +
  labs(x="Early Dysbiosis", y="Proportion") + 
  theme_bw() + 
  scale_fill_viridis(discrete=TRUE, name = "Breastfed") +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14))
FigS4D

#delivery 
FigS4E <- ggplot(compMD, aes(x=EarlyCDHi, fill=delivery)) + geom_bar(position = "fill") +
  labs(x="Early Dysbiosis", y="Proportion") + 
  theme_bw() + 
  scale_fill_viridis(discrete=TRUE, name = "Delivery", 
                     labels = c("C-section", "Vaginal")) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
FigS4E

#preterms had a lower rate of early dysbiosis 
FigS4F <- ggplot(compMD, aes(x=EarlyCDHi, fill=preterm)) + geom_bar(position = "fill") +
  labs(x="Early Dysbiosis", y="Proportion") + theme_bw() + 
  scale_fill_viridis(discrete=TRUE, name = "Gestation", 
                     labels = c("Full Term", "Preterm")) +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=14))
FigS4F

#not significant, there are only 2 pancreatic sufficient 
FigS4G <- ggplot(compMD, aes(x=EarlyCDHi, fill=Pancreatic_sufficiency)) + 
  geom_bar(position = "fill") +
  labs(x="Early Dysbiosis", y="Proportion") + 
  theme_bw() + 
  scale_fill_viridis(discrete=TRUE, name = "PS") +
  theme(axis.text=element_text(size=14), 
        axis.title=element_text(size=14))
FigS4G

#FigS5 B-G
library(cowplot)
tiff("FigS4BG_CatCharacter.tiff", units="in", width=10.5, height=7, res=300)
plot_grid(FigS4B, FigS4C,FigS4D, FigS4E, FigS4F, FigS4G, align = "v", axis = "lr")
dev.off()
