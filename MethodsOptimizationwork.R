rm(list = ls())
library(patchwork);library(rstatix);library(ggpubr)
source("~/Startup_script-Master2024.R")
setwd("~/Dropbox/Documents/Postdoc/OperationClapback/")

#save.image("MethodsOptimizationComplete.Rdata")
load("MethodsOptimizationComplete.Rdata")

Colpal_4 = c("lightgrey", "#8c6510", "#0D2B52", "#FF7399")


# Assembly stuff ----------------------------------------------------------


db = MakeSeqDB("combined_Virome_megahit_contigs.fna") %>% # 
  select(-sequences) %>% 
  #Testassembly is a testfile i made from the fasta assembly from /fs/project/PAS1117/guillermo/RNA_viromics_optimization/megahit_assemblies/C2_RNase_megahit_contigs.fna
  separate(., headers, into = c("Sample","B","C","D","E","Multi","Len"), sep = "_", remove = F) %>% 
  select(-E, -Len) %>% 
  unite(Misc, B,C,D) %>% 
  mutate(Multi = str_replace_all(Multi, "multi=",""),
         Sample = str_replace_all(Sample, ">","")) 
pull(seq_lengths) %>% summary()
# From here I have a dataframe that I can use map to execute acros all of the assembly files, each of which will be boiled down to summary statistics of the assembly lengths. I will write a short article for the results.
# For growth, I want to make this an executable script.



# Dataframe from megahit output. Information about the lengths of individual contigs and number of contigs for each treatment/sample
X = read.table("combined_Virome_megahit_contigs_headers.txt", header = F) %>% 
  #mutate(V1 = str_squish(V1)) %>% 
  mutate(V1 = str_replace_all(V1, "TARA_", "TARA"),
         V1 = str_replace_all(V1, "CONC_", "CONC"),
         V1 = str_replace_all(V1, "_RNase", "RNase"),
         #V1 = str_replace_all(V1, "\\+", "plus"),
         V1 = str_replace_all(V1, "_SRF", "SRF")) %>% 
  #mutate(V2 = str_count(V1, "_"))
  separate(V1, sep = "_", into = c("Sample_ID", "B", "C", "Flag", "Multi", "Len")) %>% 
  unite(SeqID,B,C) %>% 
  mutate(Len = as.numeric(str_replace_all(Len, "\\w+\\=", "")),
         Multi = as.numeric(str_replace_all(Multi, "\\w+\\=", "")),
         Flag = str_replace_all(Flag, "\\w+\\=", ""),
         Sample_ID = str_replace_all(Sample_ID, ">", "")) %>% 
  mutate(LibPrep = case_when(str_detect(Sample_ID, regex("RNA")) ~"Ovation",
                             str_detect(Sample_ID, regex("^TARA\\d{1,2}")) ~"TARA_Ovation",
                             str_detect(Sample_ID, regex("SISPA1[012]")) ~"SISPArandom",
                             str_detect(Sample_ID, regex("SISPARPC")) ~"Ctls_SISPArandom",
                             str_detect(Sample_ID, regex("SISPA[789][BE]")) ~"SISPApolyT",
                             #str_detect(Sample_ID, regex("polyT")) ~"SISPApolyT",
                             str_detect(Sample_ID, regex("polyTC[-+]$")) ~"Ctls_SISPApolyT",
                             str_detect(Sample_ID, regex("^C\\d{1}")) ~"Ctls_Ovation",
                             TRUE ~Sample_ID)) %>% 
  mutate(SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+")),
         Spikein = ifelse(SampDigit %in% c(10,11,12), "Spiked", "UnSpiked"),
         ExtMethod = case_when(str_detect(Sample_ID, ".+[D,E]$") ~"PS", 
                               str_detect(Sample_ID, ".+[A,B]$") ~"QVR",
                               str_detect(Sample_ID, "TARA") ~"QVR",
                               str_detect(Sample_ID, "^C") ~"Ctls_QVR",
                               str_detect(Sample_ID, "C[+-]$") ~"Ctls_QVR",
                               TRUE ~Sample_ID))
### I have identified the control samples. 
### Pairwise t.tests to assess significantly and exclude from other comparative plots (Undone)


#write_csv(X, file = "extable.csv")
#Make a boxplot median, quartiles, and outliers corresponding to each Sample_ID
X %>% 
  filter(!str_detect(LibPrep, "Ctls_.+")) %>% 
  ggplot(., aes(y = Sample_ID, x = Len, color = Spikein)) +
  geom_boxplot(outlier.size = .3, lwd = .3) +
  stat_summary(fun = "mean", shape = 23, size = .1) +
  labs(x = "Length",
       y = "Sample") +
  scale_x_continuous(breaks = c(seq(0,2500,500), seq(2500,18500,2000))) +
  geom_vline(xintercept = 500, linetype = "dashed", color = "red", alpha = .8) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1))

#write.csv(Len_summary, file = "Lensummary.csv")
X %>% pull(Len) %>% summary()


# Next task, Use Nest and Map and summary to get summary statistics for each
Len_summary = X %>%
  #filter(Sample_ID == "C2RNase") %>% pull(Len) %>% mean()
  group_by(Sample_ID) %>% 
  summarise(N_contigs = n(),
            MinLen = min(Len),
            MaxLen = max(Len),
            MeanLen = mean(Len),
            MedLen = median(Len),
            q1Len = quantile(Len, probs = 0.25),
            q3Len = quantile(Len, probs = 0.75),
            SdLen = sd(Len)) %>% 
  left_join(., X %>% group_by(Sample_ID, LibPrep, SampDigit, Spikein, ExtMethod) %>% count()) %>% 
  select(-SampDigit) %>% 
  mutate(Control = ifelse(str_detect(ExtMethod, "(Ctls_.+)|(denatSISPA)"), TRUE, FALSE))
write.csv(Len_summary, file = "Length_summary.csv")

Len_summary %>% data.frame()

### Make a names key
Metadata_key = Len_summary %>% 
  select(Sample_ID, Spikein, LibPrep, ExtMethod) %>% 
  data.frame() %>% 
  rbind(vOTUs_raw %>% group_by(Sample_ID, LibPrep, Spikein, ExtMethod) %>% summarise(total = mean(Len)) %>% filter(str_detect(Sample_ID, "(CONC)|(RNA10A)")) %>% mutate(Sample_ID = str_replace_all(Sample_ID, "-","")) %>%  select(-total)) 


### Make a pyramid barplot of the different extraction methods reflected across different lib preps
for_pyramid = Len_summary %>% 
  mutate(SampDigit = as.character(str_match_all(Sample_ID, "\\d+"))) %>% 
  select(Sample_ID, LibPrep, SampDigit, ExtMethod, N_contigs, Spikein, MeanLen) %>% pivot_longer(cols = c(N_contigs, MeanLen), names_to = "metric", values_to = "values") %>% 
  group_by(LibPrep, SampDigit, ExtMethod, metric, Spikein) %>% 
  summarise(MeanValue = mean(values)) %>% 
  ungroup() %>% 
  filter(!str_detect(ExtMethod, "Ctls_")) %>% 
  unite(ID1, LibPrep, SampDigit, remove = F) %>% 
  unite(ID2, ExtMethod, SampDigit, remove = F) %>% 
  filter(!str_detect(ID1, regex("(TARA)|(character)"))) %>% 
  arrange(SampDigit) %>% 
  mutate(MeanValue = ifelse(ExtMethod == "PS", MeanValue *-1,MeanValue))

for_pyramid %>% 
  filter(metric == "N_contigs") %>% 
  #mutate(ID2 = factor(ID1, levels = intermediatetable %>% distinct(ID1) %>% arrange())) %>% 
  ggplot(., aes(y = SampDigit, x = MeanValue, fill = LibPrep)) +
  geom_bar(stat = "identity", position = "dodge", show.legend = F) + 
  geom_vline(xintercept = 0, size = .1) +
  labs(y = "Sample",
       x = "Number of Contigs") +
  scale_fill_manual(name = "",
                    values = Colpal_4) +
  scale_x_continuous(breaks = seq(-30000,30000,10000)) +
  coord_cartesian(xlim = c(-30000,30000)) +
  theme_minimal() +
  facet_wrap(~Spikein) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  for_pyramid %>% 
  filter(metric == "MeanLen") %>% 
  #mutate(ID2 = factor(ID1, levels = intermediatetable %>% distinct(ID1) %>% arrange())) %>% 
  ggplot(., aes(y = SampDigit, x = MeanValue, fill = LibPrep)) +
  geom_bar(stat = "identity", position = "dodge") + 
  geom_vline(xintercept = 0, size = .1) +
  labs(y = "",
       x = "Average Contig Length") +
  scale_fill_manual(name = "",
                    values = Colpal_4) +
  scale_x_continuous(breaks = seq(-800,800,100)) +
  theme_minimal() +
  facet_wrap(~Spikein) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(last_plot(), filename = "Pyramid_prepext_AverageLength.pdf", height = 8, width = 8)

### overall does ps produce more and longer contigs?
for_pyramid %>% 
  ggplot(., aes(y = SampDigit, x = MeanValue, fill = ExtMethod)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Extraction",
                    values = c("grey", "steel blue")) +
  labs(y="Sample",
       x="Number of Contigs")
theme_minimal()
facet_wrap(~metric, scales = "free") 
###Jitter plot combining all of the sources that illustrates the number of contigs for the different library prep method
Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  ggplot(., aes(x = "Samples", y = N_contigs)) +
  geom_boxplot(width = .2, outlier.shape = 4, outlier.colour = "red", outlier.size = .4, color = "grey")+
  geom_jitter(width = .1,  aes(color = ExtMethod, shape = LibPrep)) +
  labs(x = "",
       y = "Number of Contigs") +
  scale_color_manual(values = c("grey","steel blue")) +
  scale_shape_manual(breaks = c("Ovation", "TARA_Ovation", "SISPApolyT", "SISPArandom"),
                     values = c(16,1,17,2)) +
  theme_minimal() +
  #ggsave(last_plot(), filename = "ASSemblyPrepExtJitter.pdf", height = 6, width = 6)
  
  # Jitterplot comparing the number of contigs produced from each sample. Colored by the extraction method (pink or blue) where shapes represent Library prep methods.
  # QVR produced the most contigs and PolyT priming outperformed random priming.
  
  ## Jitterplot of average length
  Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  ggplot(., aes(x = "Samples", y = MeanLen)) +
  geom_boxplot(width = .2, outlier.shape = 4, outlier.colour = "red", outlier.size = .4, color = "grey")+
  geom_jitter(width = .1,  aes(color = ExtMethod, shape = LibPrep), show.legend = F) +
  labs(x = "",
       y = "Mean Contig Length") +
  scale_color_manual(values = c("grey","steel blue")) +
  scale_shape_manual(breaks = c("Ovation", "TARA_Ovation", "SISPApolyT", "SISPArandom"),
                     values = c(16,1,17,2)) +
  theme_minimal()
ggsave(last_plot(), filename = "ASSemblyPrepExtJitter.pdf", height = 6, width = 6)
#### Again with facetted spikein because the results touted in the proposal show how well Ovation does with collecting spiked viruses
###Jitter plot combining all of the sources that illustrates the number of contigs for the different library prep method
Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  ggplot(., aes(x = "Samples", y = N_contigs)) +
  geom_boxplot(width = .2, outlier.shape = 4, outlier.colour = "red", outlier.size = .4, color = "grey")+
  geom_jitter(width = .1,  aes(color = ExtMethod, shape = LibPrep)) +
  labs(x = "",
       y = "Number of Contigs") +
  scale_color_manual(values = c("grey","steel blue")) +
  scale_shape_manual(breaks = c("Ovation", "TARA_Ovation", "SISPApolyT", "SISPArandom"),
                     values = c(16,1,17,2)) +
  theme_minimal() +
  facet_wrap(~Spikein) +
  #ggsave(last_plot(), filename = "ASSemblyPrepExtJitter.pdf", height = 6, width = 6)
  
  # Jitterplot comparing the number of contigs produced from each sample. Colored by the extraction method (pink or blue) where shapes represent Library prep methods.
  # QVR produced the most contigs and PolyT priming outperformed random priming.
  
  ## Jitterplot of average length
  Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  ggplot(., aes(x = "Samples", y = MeanLen)) +
  geom_boxplot(width = .2, outlier.shape = 4, outlier.colour = "red", outlier.size = .4, color = "grey")+
  geom_jitter(width = .1,  aes(color = ExtMethod, shape = LibPrep), show.legend = F) +
  labs(x = "",
       y = "Mean Contig Length") +
  scale_color_manual(values = c("grey","steel blue")) +
  scale_shape_manual(breaks = c("Ovation", "TARA_Ovation", "SISPApolyT", "SISPArandom"),
                     values = c(16,1,17,2)) +
  theme_minimal() + 
  facet_wrap(~Spikein)
ggsave(last_plot(), filename = "ASSemblyPrepExtJitter_spikefacet.pdf", height = 8, width = 8)


#Side by side boxplot
lookup_table = c("MaxLen" = "Max Length", "MeanLen" = "Average Length", "N_contigs" = "N.contigs", "MinLen" = "Min Length")
Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  pivot_longer(cols = c(N_contigs, MinLen, MaxLen, MeanLen), names_to = "variables", values_to = "values") %>% 
  mutate(variables = case_when(variables == "MaxLen" ~"Max Length", 
                               variables == "MeanLen" ~"Average Length",
                               variables == "N_contigs" ~"N.contigs",
                               variables == "MinLen" ~"Min Length")) %>% 
  ggplot(., aes(x = values, y = LibPrep)) +
  geom_boxplot(outlier.size = .8, outlier.shape = 3, outlier.color = "red", outlier.alpha = .8) +
  stat_summary(fun = mean, shape = 23, size = .2) +
  geom_jitter(width = .2) +
  labs(x = "Contig Length",
       y = "Library Prep") +
  theme_minimal() +
  facet_wrap(~variables, scales = "free", labeller = labeller(group = lookup_table))
ggsave(last_plot(), filename = "ASSemblyGroupSummary_boxplot.pdf", height = 8, width = 7)

#Plot Number of contigs
NcontigExt.plot = Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  mutate(ID = factor(Sample_ID, levels = Len_summary %>% arrange(N_contigs) %>%  pull(Sample_ID))) %>% 
  ggplot(., aes(x = N_contigs, y = ID, fill = ExtMethod, color = LibPrep)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(name = "Extraction",
                    values = c("white","black")) +
  labs(y = "Sample ID",
       x = "Number of Contigs") +
  theme_minimal()

NcontigLib.plot = Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  mutate(ID = factor(Sample_ID, levels = Len_summary %>% arrange(N_contigs) %>%  pull(Sample_ID))) %>% 
  ggplot(., aes(x = N_contigs, y = ID, fill = LibPrep, color = LibPrep)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(name = "Lib Prep",
                    breaks = c( "SISPApolyT", "SISPArandom", "Ovation", "TARA_Ovation"),
                    values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
  labs(y = "",
       x = "Number of Contigs") +
  theme_minimal() +
  theme(axis.text.y = element_blank())

NcontigExt.plot + NcontigLib.plot
ggsave(last_plot(), filename = "N_contigs_bars.pdf", height = 8, width = 7)
# QVR produced the largest number of contigs

#Plot average length
MeanLenExt.plot = Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  mutate(ID = factor(Sample_ID, levels = Len_summary %>% arrange(MeanLen) %>%  pull(Sample_ID))) %>% 
  ggplot(., aes(x = MeanLen, y = ID, fill = ExtMethod, color = LibPrep)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(name = "Extraction",
                    values = c("white","black")) +
  labs(y = "Sample ID",
       x = "Mean Contig Length") +
  theme_minimal()

MeanLenLib.plot = Len_summary %>% 
  filter(!Control,
         LibPrep != "denatSISPA") %>% 
  mutate(ID = factor(Sample_ID, levels = Len_summary %>% arrange(MeanLen) %>%  pull(Sample_ID))) %>% 
  ggplot(., aes(x = MeanLen, y = ID, fill = LibPrep, color = LibPrep)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(name = "Lib Prep",
                    breaks = c( "SISPApolyT", "SISPArandom", "Ovation", "TARA_Ovation"),
                    values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
  labs(y = "",
       x = "Mean Contig Length") +
  theme_minimal() +
  theme(axis.text.y = element_blank())

MeanLenExt.plot+MeanLenLib.plot
ggsave(last_plot(), filename = "AvgLength_bars.pdf", height = 8, width = 7)

NcontigExt.plot+NcontigLib.plot + MeanLenExt.plot+MeanLenLib.plot + patchwork::plot_layout(ncol = 2)


### Last step is to compare controls with samples using a t.test/ ANOVA or wilcox/ kruskal
#Conduct an Anova on the controls
TukResults = Len_Forplot %>% 
  aov(.$N_contigs ~ .$Colorme, data = .) %>% 
  TukeyHSD(ordered = TRUE)
Tuk_Df = as.data.frame(TukResults$`.$Colorme`)

write.csv(Tuk_Df, file = "TukeyHSD.csv")

TukMeanResults = Len_Forplot %>% 
  aov(.$MeanLen ~ .$Colorme, data = .) %>% 
  TukeyHSD(ordered = TRUE)
TukMean_Df = as.data.frame(TukMeanResults$`.$Colorme`)

write.csv(TukMean_Df, file = "TukeyHSD_MeanLen.csv")
library(corrplot)
cormat = Len_Forplot %>% column_to_rownames("Sample_ID") %>% select(-Colorme, -SdLen, -starts_with("q")) %>% cor()
corrplot(cormat)

cormat_ps = cor.mtest(Len_Forplot %>% column_to_rownames("Sample_ID") %>% select(-Colorme, -SdLen, -starts_with("q")))
#plot only correlations below alpha (p-value) below 0.01 of as circles containing the correlation value
corrplot(cormat,
         type = "lower",
         p.mat = cormat_ps$p,
         sig.level = c(0.0001, 0.001, 0.01),
         insig = "label_sig",
         addCoef.col = "black",
         #number.cex = 0.7,
         pch.cex = 0.7,
         method = "circle",
         order = "hclust",
         #tl.pos = "r", #text label position
         tl.srt = 0, #text label angle 
         cl.pos = "b" #color label position
)

#End of Assembly Work.
save.image("CBAssemblyWork.Rdata")
load("CBAssemblyWork.Rdata")


# VOTU stuff --------------------------------------------------------------


##### vOTU work
#Set up environment
rm(list = ls()) #clean
options(scipen = 999) #adjust number display
source("~/Startup_script-Master2024.R") #Load in all my heavy use tools
library(broom);library(reshape2);library(vegan);library(patchwork) #Load additional packages
setwd("~/Dropbox/Documents/Postdoc/OperationClapback/") #Set my working directory
#Decide color palette used for this analysis
display.jisho.pal("american")
Colpal_4 = c("#FFFFFF", "#8c6510", "#0D2B52", "#FF7399")
Colpal_5 = c("#0D2B52", "#8c6510", "#FF7399", "grey")
Colpal_6 = c("#0D2B52", "#8c6510", "#FF7399", "#FFFFFF")

## Produce the table of vOTUs
vOTUs_raw = read.table("vOTU_headers.txt", header = FALSE) %>% 
  mutate(V1 = str_replace_all(V1, "(\\d)_(\\d)", "\\1-\\2")) %>% #replace the underscore between k141 and the number with a dash.
  separate(., V1, into = c("Sample_ID", "vOTU_ID", "Flag", "Multi", "Len"), sep = "_") %>% #separate out the header
  mutate(Multi = as.numeric(str_replace_all(Multi, "\\w+\\=", "")), #clean the pieces
         Len = as.numeric(str_replace_all(Len, "\\w+\\=", "")),
         Flag = str_replace_all(Flag, "\\w+\\=", ""),
         Sample_ID = str_replace_all(Sample_ID, ">", ""),
         vOTU_ID = str_replace_all(vOTU_ID, "k141", "vOTU")) %>% 
  mutate(LibPrep = case_when(str_detect(Sample_ID, regex("RNA")) ~"Ovation",
                             str_detect(Sample_ID, regex("^TARA-\\d{1,3}")) ~"TARA_Ovation",
                             str_detect(Sample_ID, regex("SISPA1[012]")) ~"SISPArandom",
                             str_detect(Sample_ID, regex("SISPARPC")) ~"Ctls_SISPArandom",
                             str_detect(Sample_ID, regex("SISPA[789][BE]")) ~"SISPApolyT",
                             #str_detect(Sample_ID, regex("polyT")) ~"SISPApolyT",
                             str_detect(Sample_ID, regex("polyTC[-+]$")) ~"Ctls_SISPApolyT",
                             str_detect(Sample_ID, regex("^C\\d{1}")) ~"Ctls_Ovation",
                             TRUE ~Sample_ID)) %>% 
  mutate(SampDigit = as.numeric(str_match_all(Sample_ID, "\\d+")),
         Spikein = ifelse(SampDigit %in% c(10,11,12), "Spiked", "UnSpiked"),
         ExtMethod = case_when(str_detect(Sample_ID, ".+[D,E]$") ~"PS", 
                               str_detect(Sample_ID, ".+[A,B]$") ~"QVR",
                               str_detect(Sample_ID, "TARA") ~"QVR",
                               str_detect(Sample_ID, "^C") ~"Ctls_QVR",
                               str_detect(Sample_ID, "C[+-]$") ~"Ctls_QVR",
                               TRUE ~Sample_ID)) %>% 
  select(-Flag)

vOTUs_raw %>% distinct(LibPrep)
write.csv(vOTUs_raw, file = "raw_vOTU_analysis_Table.csv") #this is the rawest data table that will be used for analysis

Ordered = vOTUs_raw %>% 
  count(Sample_ID) %>% #Count the number of sequences corresponding to each sample ID. A measure of species richness.
  left_join(., vOTUs_raw %>% group_by(Sample_ID, LibPrep, Spikein, ExtMethod) %>% count() %>% select(-n), by = "Sample_ID") %>% 
  arrange(desc(Spikein),desc(n)) %>% pull(Sample_ID)

## Species richness
Rich.Ext = vOTUs_raw %>% distinct(Sample_ID)
count(Sample_ID) %>% #Count the number of sequences corresponding to each sample ID. A measure of species richness.
  left_join(., vOTUs_raw %>% group_by(Sample_ID, LibPrep, Spikein, ExtMethod) %>% count() %>% select(-n), by = "Sample_ID") %>% 
  mutate(Ordered = factor(Sample_ID, levels = rev(Ordered))) %>% 
  filter(!str_detect(ExtMethod, "(Ctl)|(99C)")) %>% 
  ggplot(., aes(y = Ordered, x = n, fill = ExtMethod)) + #Plot
  geom_bar(stat = "identity", color = "black") +
  geom_hline(yintercept = "SISPA7B", linetype = "dashed") +
  scale_fill_manual(name = "Extraction",
                    values = c("black","white")) +
  labs(y = "Sample",
       x = "Number of vOTUs") +
  theme_minimal()

Rich.Lib = vOTUs_raw %>% 
  count(Sample_ID) %>% #Count the number of sequences corresponding to each sample ID. A measure of species richness.
  left_join(., vOTUs_raw %>% group_by(Sample_ID, LibPrep, Spikein, ExtMethod) %>% count() %>% select(-n), by = "Sample_ID") %>% 
  mutate(Ordered = factor(Sample_ID, levels = rev(Ordered))) %>% 
  filter(!str_detect(ExtMethod, "(Ctl)|(99C)")) %>% 
  ggplot(., aes(y = Ordered, x = n, fill = LibPrep)) + #Plot
  geom_bar(stat = "identity", color = "black") +
  geom_hline(yintercept = "SISPA7B", linetype = "dashed") +
  scale_fill_manual(name = "Lib Prep",
                    breaks = c( "SISPApolyT", "SISPArandom", "Ovation", "TARA_Ovation"),
                    values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
  labs(y = "Sample",
       x = "Number of vOTUs") +
  theme_minimal()


## Summary statistics of vOTU sequence lengths
Lengths = vOTUs_raw %>% 
  count(Sample_ID)
vOTUs_raw %>% filter(LibPrep == "SISPArandom")
#1) Table
Length_summary_table = vOTUs_raw %>% 
  nest(Data = -Sample_ID) %>% 
  mutate(Summary_Stats = map(Data, ~summary(.$Len) %>% tidy())) %>% 
  select(-Data) %>% 
  unnest(Summary_Stats) %>% 
  left_join(., Lengths) %>% 
  left_join(., vOTUs_raw %>% group_by(Sample_ID, LibPrep, Spikein, ExtMethod) %>% count() %>% select(-n), by = "Sample_ID")

write.csv(Length_summary_table, file = "vOTU_Length_Summary_table.csv")
### barplot of mean length comparing the extmethods and other methods


OTUlen.Ext = Length_summary_table %>% 
  filter(!str_detect(ExtMethod, "(Ctl)|(99C)")) %>%
  #filter(!Spikein == "UnSpiked") %>% 
  mutate(Sample_ID = factor(Sample_ID, levels = c(vec2,vec))) %>% #vectors were custom made to stack unspiked atop of spiked
  ggplot(., aes(y = Sample_ID, x = mean, fill = ExtMethod)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(yintercept = "RNA7E", linetype = "dashed") +
  scale_fill_manual(name = "Extraction",
                    values = c("black","white")) +
  labs(x = "Mean vOTU Length",
       y = "Sample ID") +
  theme_minimal()

OTUlen.Lib = Length_summary_table %>% 
  filter(!str_detect(ExtMethod, "(Ctl)|(99C)")) %>%
  #filter(!Spikein == "UnSpiked") %>% 
  mutate(Sample_ID = factor(Sample_ID, levels = c(vec2,vec))) %>% #vectors were custom made to stack unspiked atop of spiked
  ggplot(., aes(y = Sample_ID, x = mean, fill = LibPrep)) +
  geom_bar(stat = "identity", color = "black") +
  geom_hline(yintercept = "RNA7E", linetype = "dashed") +
  scale_fill_manual(name = "Lib Prep",
                    breaks = c( "SISPApolyT", "SISPArandom", "Ovation", "TARA_Ovation"),
                    values = c("#a6cee3","#1f78b4","#b2df8a","#33a02c")) +
  labs(x = "Mean vOTU Length",
       y = "Sample ID") +
  theme_minimal()

Rich.Ext + Rich.Lib + OTUlen.Ext + OTUlen.Lib + patchwork::plot_layout(ncol = 2)

save.image("vOTUwork.Rdata")


# Community diveristy stuff -----------------------------------------------


######## Alpha diversity 
#Create vOTU table
vOTU_table = read.csv("Clapback_OTUtable.csv") %>% #select(contains("SISPApolyTC.")) %>% colSums()
  mutate(Contig = str_replace_all(Contig, "(\\d)_(\\d)", "\\1-\\2")) %>% #replace the underscore between k141 and the number with a dash.
  separate(., Contig, into = c("Sample_ID", "vOTU_ID", "Flag", "Multi", "Len"), sep = "_") %>% #separate out the header
  mutate(Multi = as.numeric(str_replace_all(Multi, "\\w+\\=", "")), #clean the pieces
         Len = as.numeric(str_replace_all(Len, "\\w+\\=", "")),
         Flag = str_replace_all(Flag, "\\w+\\=", ""),
         Sample_ID = str_replace_all(Sample_ID, ">", ""),
         vOTU_ID = str_replace_all(vOTU_ID, "k141", "vOTU")) %>% 
  unite(vOTUID, Sample_ID, vOTU_ID, sep = "_") %>% 
  select(-Flag, -Multi, -Len) %>% 
  column_to_rownames("vOTUID")


library(vegan)
div_InvSimp = vOTU_table %>% 
  diversity(., index = "invsimpson",2) %>% 
  data.frame("InvSimp" = .) %>% 
  rownames_to_column("Sample_ID")


#Calculate Sannon
div_Shannon = vOTU_table %>% 
  diversity(., index = "shannon",2) %>% 
  data.frame("Shannon" = .) %>% 
  rownames_to_column("Sample_ID")

#Calculate Species Richness
div_richness = vOTU_table %>% rownames_to_column("ID") %>%  
  pivot_longer(cols = -ID, values_to = "values", names_to = "variables") %>% 
  arrange(variables) %>% mutate(bin = ifelse(values > 0,1,0)) %>% 
  group_by(variables) %>% 
  summarise("Richness" = sum(bin)) %>% 
  rename("Sample_ID" = variables)

#Calculate Sequence Abundance
vOTU_abundances = vOTU_table %>% 
  colSums(.) %>% 
  data.frame("Abundance" = .) %>% 
  rownames_to_column("Sample_ID")

## Create a diveristy tables
Diversity_Table = div_richness %>% 
  left_join(., div_InvSimp) %>% 
  left_join(., div_Shannon) %>% 
  left_join(., vOTU_abundances) %>% 
  data.frame() %>% 
  mutate(Sample_ID = case_when(str_detect(Sample_ID, "SISPApolyTC\\.\\.") ~"SISPApolyTC+",
                               str_detect(Sample_ID, "SISPApolyTC\\.$") ~"SISPApolyTC-",
                               str_detect(Sample_ID, "SISPARPC\\.\\.") ~"SISPARPC+",
                               str_detect(Sample_ID, "SISPARPC\\.$") ~"SISPARPC-",
                               str_detect(Sample_ID, "SISPA_99C") ~"denatSISPA",
                               str_detect(Sample_ID, "RNAC\\.") ~"RNAC+",
                               str_detect(Sample_ID, "SISPARPC\\.") ~"SISPARPC-",
                               TRUE ~Sample_ID),
         Sample_ID = str_replace_all(Sample_ID, "_","")) %>% 
  left_join(., Metadata_key)

write.csv(Diversity_Table, file = "diversitytable.csv")

#Samples with less than 10 OTUs removed from further analysis
Remove = Diversity_Table %>% 
  filter(Richness < 10)
Remove %>% pull(Sample_ID)

#
Diversity_Table_Keep = Diversity_Table %>% 
  filter(Richness > 10)

#Long version with factors and grouping
Diversity_Table.long = Diversity_Table_Keep %>% 
  pivot_longer(cols = c(-Sample_ID, -Spikein, -LibPrep, -ExtMethod) , names_to = "variable", values_to = "values") %>% 
  mutate(variable = factor(variable, levels = c("Shannon", "InvSimp", "Richness", "Abundance")))

## Wilcox tables- Shannon
Shannon_pairPs = Diversity_Table.long %>% 
  filter(ExtMethod == "QVR"| ExtMethod == "PS",
         variable == "Shannon") %>% 
  split(.$ExtMethod) %>% 
  map(., #each dataframe
      ~pairwise_wilcox_test(values ~ LibPrep, data = .))

## Compare library preps after QVR
Shannon_QVR = Diversity_Table.long %>% 
  filter(ExtMethod == "QVR",
         variable == "Shannon") %>% 
  #Spikein == "Spiked") %>% 
  ggplot(., aes(x = LibPrep, y = values)) +
  geom_boxplot(aes(fill = LibPrep), show.legend = F) +
  geom_jitter()+
  #stat_pvalue_manual(Shannon_pairPs$QVR %>% 
  #                     add_xy_position()) +
  scale_fill_manual(values = jisho_picker("la_lakers")) +
  labs(y = "Shannon",
       x = "")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## compare library preps after PS
Shannon_PS = Diversity_Table.long %>% 
  filter(ExtMethod == "PS",
         variable == "Shannon") %>% 
  #Spikein == "Spiked") %>% 
  ggplot(., aes(x = LibPrep, y = values)) +
  geom_boxplot(aes(fill = LibPrep), show.legend = F) +
  geom_jitter()+
  #stat_pvalue_manual(Shannon_pairPs$PS %>% 
  #                     add_xy_position()) +
  scale_fill_manual(values = c("#fdb927", "#000000")) +
  labs(y = "Shannon",
       x = "")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


#### Simpson test
Simpson_pairPs = Diversity_Table.long %>% 
  filter(ExtMethod == "QVR"| ExtMethod == "PS",
         variable == "InvSimp") %>% 
  split(.$ExtMethod) %>% 
  map(., #each dataframe
      ~pairwise_wilcox_test(values ~ LibPrep, data = .))

InvSimp_QVR = Diversity_Table.long %>% 
  filter(ExtMethod == "QVR",
         variable == "InvSimp") %>% 
  #Spikein == "Spiked") %>% 
  ggplot(., aes(x = LibPrep, y = values)) +
  geom_boxplot(aes(fill = LibPrep)) +
  geom_jitter()+
  #stat_pvalue_manual(Simpson_pairPs$QVR %>% 
  #                     add_xy_position()) +
  scale_fill_manual(values = jisho_picker("la_lakers")) +
  labs(y = "Simpson",
       x = "")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## compare library preps after PS
InvSimp_PS = Diversity_Table.long %>% 
  filter(ExtMethod == "PS",
         variable == "InvSimp") %>% 
  #Spikein == "Spiked") %>% 
  ggplot(., aes(x = LibPrep, y = values)) +
  geom_boxplot(aes(fill = LibPrep)) +
  geom_jitter()+
  #stat_pvalue_manual(Simpson_pairPs$PS %>% 
  #                    add_xy_position()) +
  scale_fill_manual(values = c("#fdb927", "#000000")) +
  labs(y = "Simpson",
       x = "")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


#### Abundances
Abundance_pairPs = Diversity_Table.long %>% 
  filter(ExtMethod == "QVR"| ExtMethod == "PS",
         variable == "Abundance") %>% 
  split(.$ExtMethod) %>% 
  map(., #each dataframe
      ~pairwise_wilcox_test(values ~ LibPrep, data = .))

Abundance_QVR = Diversity_Table.long %>% 
  filter(ExtMethod == "QVR",
         variable == "Abundance") %>% 
  #Spikein == "Spiked") %>% 
  ggplot(., aes(x = LibPrep, y = values)) +
  geom_boxplot(aes(fill = LibPrep), show.legend = F) +
  geom_jitter()+
  #stat_pvalue_manual(Abundance_pairPs$QVR %>% 
  #                     add_xy_position()) +
  scale_fill_manual(values = jisho_picker("la_lakers")) +
  labs(y = "Abundance",
       x = "")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## compare library preps after PS
Abundance_PS = Diversity_Table.long %>% 
  filter(ExtMethod == "PS",
         variable == "Abundance") %>% 
  #Spikein == "Spiked") %>% 
  ggplot(., aes(x = LibPrep, y = values)) +
  geom_boxplot(aes(fill = LibPrep), show.legend = F) +
  geom_jitter()+
  #stat_pvalue_manual(Abundance_pairPs$PS %>% 
  #                     add_xy_position()) +
  scale_fill_manual(values = c("#fdb927", "#000000")) +
  labs(y = "Abundance",
       x = "")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

#### Sp richness
Richness_pairPs = Diversity_Table.long %>% 
  filter(ExtMethod == "QVR"| ExtMethod == "PS",
         variable == "Richness") %>% 
  split(.$ExtMethod) %>% 
  map(., #each dataframe
      ~pairwise_wilcox_test(values ~ LibPrep, data = .))

Richness_QVR = Diversity_Table.long %>% 
  filter(ExtMethod == "QVR",
         variable == "Richness") %>% 
  #Spikein == "Spiked") %>% 
  ggplot(., aes(x = LibPrep, y = values)) +
  geom_boxplot(aes(fill = LibPrep), show.legend = F) +
  geom_jitter()+
  #stat_pvalue_manual(Richness_pairPs$QVR %>% 
  #                     add_xy_position()) +
  scale_fill_manual(values = jisho_picker("la_lakers")) +
  labs(y = "Richness",
       x = "")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

## compare library preps after PS
Richness_PS = Diversity_Table.long %>% 
  filter(ExtMethod == "PS",
         variable == "Richness") %>% 
  #Spikein == "Spiked") %>% 
  ggplot(., aes(x = LibPrep, y = values)) +
  geom_boxplot(aes(fill = LibPrep), show.legend = F) +
  geom_jitter()+
  #stat_pvalue_manual(Richness_pairPs$PS %>% 
  #                     add_xy_position()) +
  scale_fill_manual(values = c("#fdb927", "#000000")) +
  labs(y = "Richness",
       x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

DivPlot_QVR = Shannon_QVR + InvSimp_QVR + Richness_QVR + Abundance_QVR + patchwork::plot_layout(nrow = 2)
DivPlot_PS = Shannon_PS + InvSimp_PS + Richness_PS + Abundance_PS + patchwork::plot_layout(nrow = 2)
ggsave(DivPlot_QVR, filename = "DivPlot_QVR.pdf", height = 8, width = 6)
ggsave(DivPlot_PS, filename = "DivPlot_PS.pdf", height = 8, width = 6)


getwd()
##Plot Overall diversity
Overall_Diversity = Diversity_Table.long %>% #distinct(LibPrep)
  ggplot(., aes(x = variable, y = values, color = LibPrep)) +
  geom_boxplot(width = .3, outlier.colour = "red", outlier.shape = 3, outlier.size = .4, outlier.alpha = .7) +
  geom_jitter(width = .1) +
  scale_color_manual(name = "",
                     values = c("#0D2B52", "#8c6510", "#FF7399", "grey")) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  theme(axis.text.x = element_blank())
ggsave(Overall_Diversity, filename = "OverallDiversity_combined.pdf", height = 8, width = 7)

##Boxplot of Group_diversity
Group_diversity = Diversity_Table.long %>%
  ggplot(., aes(y = Sample_ID, x = values)) +
  geom_boxplot(outlier.shape = 3, outlier.size = .5, outlier.color = "red", outlier.alpha = .7) +
  stat_summary(fun = mean,  size = .3, shape = 23) +
  geom_jitter() +
  labs(y = "",
       x = "Index/Values") +
  theme_minimal() +
  facet_wrap(~variable, scales = "free")
ggsave(Group_diversity, filename = "GroupDiversity.pdf", height = 8, width = 7)

##Shannon V Simpson
Sample_diversity = Diversity_Table %>% 
  pivot_longer(cols = -Sample_ID, names_to = "variable", values_to = "values") %>%
  mutate(Ordered = factor(Sample_ID, levels = rev(Diversity_Table %>% arrange(desc(Shannon)) %>% pull(Sample_ID))),
         Colorme = case_when(str_detect(Sample_ID, regex("RNA")) ~"RNAs",
                             str_detect(Sample_ID, regex("^TARA")) ~"TARA",
                             str_detect(Sample_ID, regex("SISPA")) ~"SISPAs",
                             str_detect(Sample_ID, regex("^C\\d{1}")) ~"Controls",
                             TRUE ~Sample_ID),
         variable = factor(variable, levels = c("Shannon", "InvSimp", "Richness", "Abundance"))) %>% 
  mutate(Colorme = factor(Colorme, levels = c("SISPAs", "RNAs", "TARA", "Controls"))) %>% 
  ggplot(., aes(x = values, y = Ordered, fill = Colorme)) +
  geom_bar(stat = "identity", color = "black") +
  labs(y = "Sample",
       x = "Index/Values") +
  scale_fill_manual(name = "",
                    values = c("#0D2B52", "#8c6510", "#FF7399", "#FFFFFF")) +
  facet_wrap(~variable, scales = "free") +
  theme_minimal()
ggsave(Sample_diversity, filename = "SampleDiversity.pdf", height = 8, width = 7)
# End pre processing and preliminary exploratory analyses


# Experiment Two reporting ------------------------------------------------


### Start distilled reporting scripts
#Create a function to calculate standard error for error bars.
se <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}
         
Diversity_Table %>% 
  filter(str_detect(Sample_ID, "(10[AB])|(11[AB])|(12[AB])"), LibPrep == "Ovation") %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, "RNA","RNA_"),
         Sample_ID = str_replace_all(Sample_ID, "(\\d{2})","\\1_")) %>% 
  separate(Sample_ID, into = c("c","a", "b"), sep = "_") %>% 
  select(-c, -Spikein, -LibPrep, -ExtMethod) %>% 
  pivot_longer(cols = c(-a,-b)) %>% 
  mutate(treatmt = case_when(str_detect(a, "10") ~"CU",
                             str_detect(a, "11") ~"CU+CsCl",
                             str_detect(a, "12") ~"Pellet+CsCl"),
         treatmt = factor(treatmt, levels = rev(c("CU", "CU+CsCl", "Pellet+CsCl")))) %>% 
  group_by(treatmt,name) %>% 
  summarise(avgvalue = mean(value),
            serror = se(value)) %>% 
  ggplot(., aes(x = avgvalue, y = treatmt, fill = treatmt)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(width = .1, size = .3, aes(xmin=avgvalue-serror, xmax=avgvalue+serror)) +
  scale_fill_manual(name = "",
                    breaks = c("CU", "CU+CsCl", "Pellet+CsCl"),
                    values = c("grey", "steelblue", "grey50")) +
  labs(x="Average Value",
       y="Treatment")+
  facet_wrap(~name, scales = "free") +
  theme_minimal()

# Barplot illustrating community diversity metrics for spiked experiment. 
X %>% distinct(Sample_ID)
## assembly stuff note from Lab Journal 13Nov24, missing RNA10A
#

Contig_table = X %>% data.frame() %>% 
  filter(str_detect(Sample_ID, "(10[AB])|(11[AB])|(12[AB])"), LibPrep == "Ovation") %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, "RNA","RNA_"),
         Sample_ID = str_replace_all(Sample_ID, "(\\d{2})","\\1_")) %>% 
  separate(Sample_ID, into = c("c","a", "b"), sep = "_") %>% 
  select(-c, -Spikein, -LibPrep, -ExtMethod) %>% 
  #filter(a == "12") %>% distinct(b) #10 is missing A in this table
  mutate(treatmt = case_when(str_detect(a, "10") ~"CU",
                             str_detect(a, "11") ~"CU+CsCl",
                             str_detect(a, "12") ~"Pellet+CsCl"),
         treatmt = factor(treatmt, levels = rev(c("CU", "CU+CsCl", "Pellet+CsCl")))) 

##Plot the average lenth of contig note one replicate
Contig_table %>%
  group_by(treatmt) %>% 
  summarise(MeanLength = mean(Len),
            serror = se(Len)) %>% 
  #pivot_longer(cols = -treatmt) %>% 
  ggplot(., aes(x = MeanLength, y = treatmt, fill = treatmt)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(width = .1, size = .3, aes(xmin=MeanLength-serror, xmax=MeanLength+serror)) +
  scale_fill_manual(name = "",
                    breaks = c("CU", "CU+CsCl", "Pellet+CsCl"),
                    values = c("grey", "steelblue", "grey50")) +
  #facet_wrap(~name, scales = "free") +
  labs(x = "Contig Length",
       y = "Treatment") +
  theme_minimal()

## Plot the number of contigs and the standard errors
Contig_table %>% 
  group_by(treatmt,b) %>% 
  summarise(NoContigs = n()) %>% 
  group_by(treatmt) %>% 
  summarise(MeanContigs = mean(NoContigs),
            serror = se(NoContigs)) %>% 
  ggplot(., aes(x = MeanContigs, y = treatmt, fill = treatmt)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(width = .1, size = .3, aes(xmin=MeanContigs-serror, xmax=MeanContigs+serror)) +
  scale_fill_manual(name = "",
                    breaks = c("CU", "CU+CsCl", "Pellet+CsCl"),
                    values = c("grey", "steelblue", "grey50")) +
  #facet_wrap(~name, scales = "free") +
  labs(x = "Number of Contigs",
       y = "Treatment") +
  theme_minimal()
ggsave(last_plot(), filename = "FinalReport/SpikednoContigs.pdf", width = 6, height = 4)  
  
### vOTU lengths Note descrepancy from Lab Journal 13Nov24
vOTUs_raw %>%
  filter(str_detect(Sample_ID, "(10[AB])|(11[AB])|(12[AB])"), LibPrep == "Ovation") %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, "RNA","RNA_"),
         Sample_ID = str_replace_all(Sample_ID, "(\\d{2})","\\1_")) %>% 
  separate(Sample_ID, into = c("c","a", "b"), sep = "_") %>% 
  #filter(a == "10") #%>% distinct(b) #10 is missing A in this table
  mutate(treatmt = case_when(str_detect(a, "10") ~"CU",
                             str_detect(a, "11") ~"CU+CsCl",
                             str_detect(a, "12") ~"Pellet+CsCl"),
         treatmt = factor(treatmt, levels = rev(c("CU", "CU+CsCl", "Pellet+CsCl")))) %>% 

  group_by(treatmt) %>% 
  summarise(MeanLength = mean(Len),
            serror = se(Len)) %>% 
  ggplot(., aes(x = MeanLength, y = treatmt, fill = treatmt)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(width = .1, size = .3, aes(xmin=MeanLength-serror, xmax=MeanLength+serror)) +
  scale_fill_manual(name = "",
                    breaks = c("CU", "CU+CsCl", "Pellet+CsCl"),
                    values = c("grey", "steelblue", "grey50")) +
  #facet_wrap(~name, scales = "free") +
  labs(x = "vOTU Seq Length",
       y = "Treatment") +
  theme_minimal()
ggsave(last_plot(), filename = "FinalReport/SpikedvOTUlength.pdf", width = 6, height = 4)  





######### Unspiked see lab journal 13Nov24 - One replicate, no pellet
Diversity_Table %>% 
  #filter(str_detect(Sample_ID, "7"))
  filter(str_detect(Sample_ID, "(7[AB])|(8[AB])|(9AB)"), LibPrep == "Ovation") %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, "RNA","RNA_"),
         Sample_ID = str_replace_all(Sample_ID, "(\\d{1})","\\1_")) %>% 
  separate(Sample_ID, into = c("c","a", "b"), sep = "_") %>% 
  select(-c, -Spikein, -LibPrep, -ExtMethod) %>% 
  pivot_longer(cols = c(-a,-b)) %>% 
  mutate(treatmt = case_when(str_detect(a, "7") ~"CU",
                             str_detect(a, "8") ~"CU+CsCl",
                             str_detect(a, "9") ~"Pellet+CsCl"),
         treatmt = factor(treatmt, levels = rev(c("CU", "CU+CsCl", "Pellet+CsCl")))) %>% 
  group_by(treatmt,name) %>% 
  summarise(avgvalue = mean(value),
            serror = se(value)) %>% 
  ggplot(., aes(x = avgvalue, y = treatmt, fill = treatmt)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(width = .1, size = .3, aes(xmin=avgvalue-serror, xmax=avgvalue+serror)) +
  scale_fill_manual(name = "",
                    breaks = c("CU", "CU+CsCl", "Pellet+CsCl"),
                    values = c("grey", "steelblue", "grey50")) +
  labs(x="Average Value",
       y="Treatment") +
  facet_wrap(~name, scales = "free") +
  theme_minimal()
ggsave(last_plot(), filename = "FinalReport/unSpikedDiversity.pdf", width = 7, height = 5)  

# Barplot illustrating community diversity metrics for spiked experiment. 
X %>% distinct(Sample_ID)
## assembly stuff note from Lab Journal 13Nov24, missing RNA10A
#

Contig_table = X %>% data.frame() %>% 
  filter(str_detect(Sample_ID, "(7[AB])|(8[AB])|(9[AB])"), LibPrep == "Ovation") %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, "RNA","RNA_"),
         Sample_ID = str_replace_all(Sample_ID, "(\\d{1})","\\1_")) %>% 
  separate(Sample_ID, into = c("c","a", "b"), sep = "_") %>% 
  select(-c, -Spikein, -LibPrep, -ExtMethod) %>% 
  #filter(a == "12") %>% distinct(b) #10 is missing A in this table
  mutate(treatmt = case_when(str_detect(a, "7") ~"CU",
                             str_detect(a, "8") ~"CU+CsCl",
                             str_detect(a, "9") ~"Pellet+CsCl"),
         treatmt = factor(treatmt, levels = rev(c("CU", "CU+CsCl", "Pellet+CsCl")))) 

##Plot the average lenth of contig
Contig_table %>%
  group_by(treatmt) %>% 
  summarise(MeanLength = mean(Len),
            serror = se(Len)) %>% 
  #pivot_longer(cols = -treatmt) %>% 
  ggplot(., aes(x = MeanLength, y = treatmt, fill = treatmt)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(width = .1, size = .3, aes(xmin=MeanLength-serror, xmax=MeanLength+serror)) +
  scale_fill_manual(name = "",
                    breaks = c("CU", "CU+CsCl", "Pellet+CsCl"),
                    values = c("grey", "steelblue", "grey50")) +
  #facet_wrap(~name, scales = "free") +
  labs(x = "Contig Length",
       y = "Treatment") +
  theme_minimal()
ggsave(last_plot(), filename = "FinalReport/unSpikedContigLength.pdf", width = 6, height = 4)  

## Plot the number of contigs and the standard errors
Contig_table %>% 
  group_by(treatmt,b) %>% 
  summarise(NoContigs = n()) %>% 
  group_by(treatmt) %>% 
  summarise(MeanContigs = mean(NoContigs),
            serror = se(NoContigs)) %>% 
  ggplot(., aes(x = MeanContigs, y = treatmt, fill = treatmt)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(width = .1, size = .3, aes(xmin=MeanContigs-serror, xmax=MeanContigs+serror)) +
  scale_fill_manual(name = "",
                    breaks = c("CU", "CU+CsCl", "Pellet+CsCl"),
                    values = c("grey", "steelblue", "grey50")) +
  #facet_wrap(~name, scales = "free") +
  labs(x = "Number of Contigs",
       y = "Treatment") +
  theme_minimal()
ggsave(last_plot(), filename = "FinalReport/unSpikednoContigs.pdf", width = 6, height = 4)  


### vOTU lengths unspiked
vOTUs_raw %>%
  filter(str_detect(Sample_ID, "(7[AB])|(8[AB])|(9[AB])"), LibPrep == "Ovation") %>% 
  mutate(Sample_ID = str_replace_all(Sample_ID, "RNA","RNA_"),
         Sample_ID = str_replace_all(Sample_ID, "(\\d{1})","\\1_")) %>% 
  separate(Sample_ID, into = c("c","a", "b"), sep = "_") %>% 
  #filter(a == "10") #%>% distinct(b) #10 is missing A in this table
  mutate(treatmt = case_when(str_detect(a, "7") ~"CU",
                             str_detect(a, "8") ~"CU+CsCl",
                             str_detect(a, "9") ~"Pellet+CsCl"),
         treatmt = factor(treatmt, levels = rev(c("CU", "CU+CsCl", "Pellet+CsCl")))) %>% 
  group_by(treatmt) %>% 
  summarise(MeanLength = mean(Len),
            serror = se(Len)) %>% 
  ggplot(., aes(x = MeanLength, y = treatmt, fill = treatmt)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(width = .1, size = .3, aes(xmin=MeanLength-serror, xmax=MeanLength+serror)) +
  scale_fill_manual(name = "",
                    breaks = c("CU", "CU+CsCl", "Pellet+CsCl"),
                    values = c("grey", "steelblue", "grey50")) +
  #facet_wrap(~name, scales = "free") +
  labs(x = "vOTU Seq Length",
       y = "Treatment") +
  theme_minimal()
ggsave(last_plot(), filename = "FinalReport/unSpikedvOTULength.pdf", width = 6, height = 3)  

 
#Make a function to plot whittaker plots for select samples
rankCurvefn = function(SAMPLE) {
  X= vOTU_table %>% 
    select(SAMPLE) %>% 
    arrange(desc(.)) %>% 
    filter(SAMPLE > 0) %>% 
    mutate(Rank = seq(1:nrow(.))) %>% 
    #rownames_to_column("Rank") %>% 
    #mutate(Rank = as.numeric(Rank))
    ggplot(., aes_string(x="Rank", y=SAMPLE)) +
    geom_line() +
    labs(y = "Read abundance",
         x = "vOTU rank",
         title = SAMPLE) +
    theme_minimal()
  return(X)
}

#Generate Rank abundance curves and concatenate in one plot
rankCurvefn("RNA7B") +
rankCurvefn("RNA8B") +
rankCurvefn("RNA10A") + rankCurvefn("RNA10B") +
rankCurvefn("RNA11A") + rankCurvefn("RNA11B") +
  patchwork::plot_layout(nrow = 2)
