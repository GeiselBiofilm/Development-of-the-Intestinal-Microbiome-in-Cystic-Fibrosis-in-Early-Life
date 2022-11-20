##### Libraries and data loading ################################################################################################################################################################
#Add several libraries that will be used
#primary ones for analysis are phyloseq and DESeq2
#others are for visualization
library("phyloseq") 
library("ggplot2")
library("ggpubr")
library("viridis")
library("nlme")
library("ggeffects")
library("cowplot")

packageVersion("phyloseq") #1.32.0
packageVersion("ggplot2") #3.3.6
packageVersion("ggpubr") #0.4.0
packageVersion("viridis") #0.6.2
packageVersion("nlme") #3.1.151
packageVersion("ggeffects") #1.1.3
packageVersion("cowplot") #1.1.1

##load in phyloseq data 
load("ASV_physeq_silva.rda")
load("sample_info.rda")


##### Recording alpha-diversity values in the metadata ###############################################################################################################################################################
#Recording alpha-diversity values in the metadata
#This is necessary for all graphs
divIndex <- estimate_richness(ASV_physeq, measures = c("Shannon"))

#sample_info_tab$Sample <- row.names(sample_info_tab)
divIndex$Sample <- row.names(divIndex)
sample_info_new <- merge(sample_info_tab, divIndex)
head(sample_info_new)

##### Stats for single model ######################################################################################################################
#Statistics for demographics/exposures; single model approach 

#remove unknowns 
sample_info <- as.data.frame(subset(sample_info_new, preterm != "Unknown" &
                                      Gender != "Unknown" & 
                                      Genotype != "Unknown" & 
                                      delivery != "Unknown" & 
                                      Pancreatic_sufficiency != "Unknown" &
                                      breastfed_ever != "Unknown"))
sample_info$RecentSaSwab[sample_info$RecentSaSwab == "Unknown"] <- NA
sample_info$RecentPaSwab[sample_info$RecentPaSwab == "Unknown"] <- NA

#testing which antibiotic measurement to use
#by seeing whether either/both interact with Age in days

#recent Abx does not significantly interact with age, p = 0.9521
SDIlme <- lme(Shannon ~ RecentAbx + Age_Days + RecentAbx * Age_Days,
              na.action = na.omit, data=sample_info, random=~1|Patient)
summary(SDIlme)

#total abx significantly interacts, p = 0.0107
SDIlme <- lme(Shannon ~ TotalAbxPrior + Age_Days + TotalAbxPrior * Age_Days,
              na.action = na.omit, data=sample_info, random=~1|Patient,
              control = lmeControl(msMaxIter = 1000, msMaxEval = 1000,
                                   sing.tol = 1e-20))
summary(SDIlme)

#model with everything 
SDIlme <- lme(Shannon ~ RecentPaSwab + 
                RecentSaSwab + 
                RecentAbx + 
                preterm + 
                Gender + 
                delivery + 
                Pancreatic_sufficiency + 
                breastfed_ever +
                Age_Days,
              na.action = na.omit, data=sample_info, random=~1|Patient,
              control = lmeControl(msMaxIter = 1000, msMaxEval = 1000,
                                   sing.tol = 1e-20))
summary(SDIlme)

#sequentially removed: PS, breastfeeding, Pa swab, delivery, Gender, Sa swab
SDIlme <- lme(Shannon ~ RecentAbx + 
                preterm +
                Age_Days,
              na.action = na.omit, data=sample_info, random=~1|Patient,
              control = lmeControl(msMaxIter = 1000, msMaxEval = 1000,
                                   sing.tol = 1e-20))

summary(SDIlme)
summary(SDIlme)$tTable[,"p-value"] #Age p=3.633e-19 #preterm p=3.087365e-2, RecentAbx p=5.0642e-2

######Figure 1 plots #########################################################################################################
#plot for Figure 1A
Fig1A <- ggplot(sample_info_new, aes(x = Age_Days, y = Shannon, color = Bin_Years)) +
  geom_point(size=2.5) +
  geom_smooth(method='lm',formula=y~x) +
  theme_bw() +
  labs(x="Age (Days)", y="Shannon Diversity Index",color="Age (Years)") +
  scale_color_viridis(discrete=TRUE, labels = c("0-1", ">1-2", ">2-3",">3-4")) +
  annotate(geom="text",x=1200, y=1.25, label="p = 3.63e-19") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=16), legend.title=element_text(size=16, face="bold"),
        plot.margin = margin(.3, .3, .3, .8, "cm"))
Fig1A

###gestation###
#linear model
Fig1B <- ggplot(sample_info, aes(x = Bin_Years, y = Shannon, color = preterm)) +
  geom_boxplot() +
  theme_bw() +
  labs(x="Age (Years)", y="Shannon Diversity Index",color="Gestation") +
  scale_color_viridis(discrete=TRUE, label = c("Full Term", "Premature")) + 
  annotate(geom="text",x=1, y=0.5, label="p = 0.0752") +
  annotate(geom="text",x=2, y=0.5, label="p = 0.6099") +
  annotate(geom="text",x=3, y=0.5, label="p = 0.0550") +
  annotate(geom="text",x=4, y=0.5, label="p = 0.3999") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=16), legend.title=element_text(size=16, face="bold"),
        plot.margin = margin(.3, .3, .3, .8, "cm"))
Fig1B

###gestation###
#linear model
Fig1C <- ggplot(sample_info, aes(x = Bin_Years, y = Shannon, color = RecentAbx)) +
  geom_boxplot() +
  theme_bw() +
  labs(x="Age (Years)", y="Shannon Diversity Index",color="Recent Antibiotics") +
  scale_color_viridis(discrete=TRUE) + 
  annotate(geom="text",x=1, y=0.5, label="p = 0.0006") +
  annotate(geom="text",x=2, y=0.5, label="p = 0.4908") +
  annotate(geom="text",x=3, y=0.5, label="p = 0.0503") +
  annotate(geom="text",x=4, y=0.5, label="p = 0.4793") +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=16), legend.title=element_text(size=16, face="bold"),
        plot.margin = margin(.3, .3, .3, .8, "cm"))
Fig1C

Fig1 <- ggarrange(Fig1A, Fig1B, Fig1C, ncol = 1, 
                  labels = c("A.","B.", "C."), align = "v")
Fig1

tiff("Fig1_Resub.tiff", units="in", width=7.5, height=10, res=600)
Fig1
dev.off()

##### Linear model for age in bins, Fig S1 #########################################################################################################
#Breaks data into single year age bins 
yr0 <- sample_info[sample_info$Bin_Years == "0",]
yr1 <- sample_info[sample_info$Bin_Years == "1",]
yr2 <- sample_info[sample_info$Bin_Years == "2",]
yr3 <- sample_info[sample_info$Bin_Years == "3",]

SDIlme0 <- lme(Shannon ~ RecentAbx + preterm + Age_Days,na.action = na.omit, data=yr0, random=~1|Patient)
SDIlme1 <- lme(Shannon ~ RecentAbx + preterm + Age_Days,na.action = na.omit, data=yr1, random=~1|Patient)
SDIlme2 <- lme(Shannon ~ RecentAbx + preterm + Age_Days,na.action = na.omit, data=yr2, random=~1|Patient)
SDIlme3 <- lme(Shannon ~ RecentAbx + preterm + Age_Days,na.action = na.omit, data=yr3, random=~1|Patient)

summary(SDIlme0) #Abx p=0.0006, preterm p=0.0752, age = 0.0091
summary(SDIlme1) #Abx p=0.4908, preterm p=0.6099, age = 0.0731
summary(SDIlme2) #Abx p=0.0503, preterm p=0.0550, age = 0.6414
summary(SDIlme3) #Abx p=0.4793, preterm p=0.3999, age = 0.1555

#Breaks data into two year age bins 
yr0_1 <- sample_info[sample_info$Bin_Years == "0" | 
                       sample_info$Bin_Years == "1", ]
yr2_3 <- sample_info[sample_info$Bin_Years == "2" |
                       sample_info$Bin_Years == "3",]

SDIlme0_1 <- lme(Shannon ~ RecentAbx + preterm + Age_Days,na.action = na.omit, data=yr0_1, random=~1|Patient)
SDIlme2_3 <- lme(Shannon ~ RecentAbx + preterm + Age_Days,na.action = na.omit, data=yr2_3, random=~1|Patient)

summary(SDIlme0_1)$tTable #Abx p=0.1374, preterm p=0.2088, age p=9.49e-19
summary(SDIlme2_3) #Abx p=0.0232, preterm p=0.0192, age p=0.7075
#########################################################################################################
FigS1A <- ggplot(yr0_1, aes(x = Age_Days, y = Shannon)) +
  geom_point(size=2.5) +
  geom_smooth(method='lm',formula=y~x) +
  theme_bw() +
  labs(x="Age (Days)", y="Shannon Diversity Index",color="Age (Years)") +
  scale_color_viridis(discrete=TRUE) +
  annotate(geom="text",x=365, y=0.4, label="p = 9.49e-13") +
  ylim(0,4) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=16), legend.title=element_text(size=16, face="bold"),
        plot.margin = margin(.3, .3, .3, .8, "cm"))
FigS1A

FigS1B <- ggplot(yr2_3, aes(x = Age_Days, y = Shannon)) +
  geom_point(size=2.5) +
  geom_smooth(method='lm',formula=y~x) +
  theme_bw() +
  labs(x="Age (Days)", y="Shannon Diversity Index",color="Age (Years)") +
  scale_color_viridis(discrete=TRUE) +
  annotate(geom="text",x=1095, y=0.4, label="p = 0.707") +
  ylim(0,4) +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold"), 
        legend.text=element_text(size=16), legend.title=element_text(size=16, face="bold"),
        plot.margin = margin(.3, .3, .3, .8, "cm"))
FigS1B

#FigS1 for export
FigS1 <- ggarrange(FigS1A, FigS1B, ncol = 1, labels = c("A.","B."))
FigS1

tiff("FigS1.tiff", units="in", width=7.5, height=7, res=600)
FigS1
dev.off()

