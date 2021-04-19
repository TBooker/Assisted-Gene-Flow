#####################################################################
##Jan, 2021
## The Whitlock Lab.
##script to plot the main and supplementary figures for the AGF paper
## Plots by Andréa T. Thomaz (deatthomaz@gmail.com)
#####################################################################
#### loading packages ####
library(ggplot2)
library(stringr)
library(plyr)
library(dplyr)
library(cowplot)
library(reshape2)
#install.packages("devtools")
#library(devtools)
#devtools::install_github("teunbrand/ggh4x")
library(ggh4x) #used for facet_nested function

#### reading and making pretty labels for graphing ####
# setting working dir
setwd("~/Dropbox/Postdoc_UBC/LabProject_AGF/A.0.6_processed/")

#reading data.frame for 10K pop. size
df <- read.csv("10000.mean.A.0.6.output.csv", header = T, stringsAsFactors = F)
#transforming NA in 0
df[is.na(df)] <- 0

#digesting the migration for single vs. multiple pulses
df$mig <- str_count(df$migration_scenario, ',') #1 = single pulse, 9 = 5 pulses 
df$mig[df$mig == 9] <- 5
df$migGen <- ifelse(df$mig == 5, 
                    as.numeric(sapply(strsplit(df$migration_scenario, ","), tail, 2)[2,]) - as.numeric(sapply(strsplit(df$migration_scenario, ","), tail, 2)[1,]),0)
df$mig_split <- as.numeric(sapply(strsplit(df$migration_scenario, ","), "[[", 1))
df$migTOTAL <- df$mig_split*df$mig

#removing sM = 0.2 and 10, and 20 OD pairs from all analyses 
df <- df %>% filter(singleSToRescale_MAC != 0.2, singleSToRescale_MAC != 10, nbPairsLoci_OD != 20)

#renaming factors for showing what is it
##migration parameters
df$migGen <- factor(df$migGen, 
                    levels=levels(factor(df$migGen)),
                    labels=c("'Single'", "'Every 1'", "'Every 2'", "'Every 4'"))

df$mig <- factor(df$mig, levels=levels(factor(df$mig)),
                 labels=c("Migration Single", "Migration Pulse (5x)"))

df$migTOTAL <- factor(df$migTOTAL, 
                      levels=levels(factor(df$migTOTAL)),
                      labels=c(expression(italic(T[f])*" = 0.5%"),
                               expression(italic(T[f])*" = 5%"), 
                               expression(italic(T[f])*" = 50%")))

##MAC
df$singleSToRescale_MAC <- factor(df$singleSToRescale_MAC, 
                                  levels=levels(factor(df$singleSToRescale_MAC)),
                                  labels=c(expression(italic(Delta[MA])~"= 0"), 
                                           expression(italic(Delta[MA])%~~%~"9%"),
                                           expression(italic(Delta[MA])~"= 50%")))

df$nbLoci_MAC <- factor(df$nbLoci_MAC,
                        levels=levels(factor(df$nbLoci_MAC)),
                        labels=c("`0 MA`", "`1 MA`", "`5 MAs`", "`50 MAs`"))

##PAC
df$singleSToRescale_PAC <- factor(df$singleSToRescale_PAC, 
                                  levels=levels(factor(df$singleSToRescale_PAC)),
                                  labels=c(expression(italic(Delta[PA])~"= 0"), 
                                           expression(italic(Delta[PA])~"= 10%"),
                                           expression(italic(Delta[PA])~"= 50%")))
df$nbLoci_PAC <- factor(df$nbLoci_PAC,
                        levels=levels(factor(df$nbLoci_PAC)),
                        labels=c("`0 PA`", "`1 PA`", "`5 PAs`", "`50 PAs`"))

##OD
df$S_OD <- factor(df$S_OD,
                  levels=levels(factor(df$S_OD)),
                  labels=c("0%", "20%", "60%", "20%", "60%", "20%", "60%"))

df$nbPairsLoci_OD <- factor(df$nbPairsLoci_OD, 
                            levels=levels(factor(df$nbPairsLoci_OD)),
                            labels = c("0 OD pairs","2 OD pairs","10 OD pairs","100 OD pairs"))

#re-calculating the hybrid index to be introgression index (0 = no introgression, 1 = total introgression)
df$hybridIndexCOR <- 1 - df$hybridIndex


########## PLOTS MAIN PAPER ###############
#### Figure 1 ####
#Single event of migration and MAC dominance factor = 0.5
df_sub <- df %>% filter(migGen == "'Single'", 
                        dominanceCoef_MAC == 0.5, 
                        nbLoci_MAC == "`5 MAs`" | nbLoci_MAC == "`0 MA`", 
                        migTOTAL != 'italic(T[f]) * " = 0.5%"', 
                        nbLoci_PAC == "`5 PAs`" | nbLoci_PAC == "`0 PA`", 
                        nbPairsLoci_OD == "0 OD pairs" | nbPairsLoci_OD == "10 OD pairs", 
                        Generation <= 40)

#need to create a new dataframe with the max and min to use in geom_blank to set the range for the y-axis at each graph, got example from here: https://chrischizinski.github.io/rstats/using_geom_blank/
minmax_sub <- df_sub %>% group_by(migTOTAL,singleSToRescale_PAC, singleSToRescale_MAC) %>% 
  summarise(min = min(relFitness), max = max(relFitness))
minmax_sub$dev1 <- ifelse((1 - minmax_sub$min) >= (minmax_sub$max-1), 
                          1 - minmax_sub$min,
                          minmax_sub$max-1)
minmax_sub$min <- 1 - minmax_sub$dev1
minmax_sub$max <- minmax_sub$dev1 + 1
minmax_sub <- melt(minmax_sub[,-6])

### PLOT! ###
ts <- ggplot() +
  geom_hline(yintercept=1, linetype="longdash", color = "grey", size = 0.8) +
  geom_line(data = df_sub[df_sub$S_OD != "0%",], aes(Generation, relFitness, colour = S_OD), size = 0.7, alpha = 0.7) +
  geom_line(data = df_sub[df_sub$S_OD == "0%",], aes(Generation, relFitness, colour = S_OD), size = 1, alpha = 1) +
  geom_blank(data = minmax_sub, aes(y=value)) + #here we use the minmax to set the range
  labs(y = "Relative Fitness", colour = expression(italic(Delta[OD])), fill = expression(italic(O[d]))) +
  #geom_point(data = df_sub_fit1, 
  #           mapping = aes(x = GenMin, y = 1, fill = S_OD), shape = 21, color = "black") +
  scale_color_manual(values=c("black", 
                              "#2c7fb8", #"#41b6c4", "#2c7fb8","#253494", 
                              "#e6550d")) + #"#fd8d3c", "#e6550d", "#a63603")) +
  scale_x_continuous(breaks = c(1,20,40), limits = c(0,40), expand = c(0,0)) +
  scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  theme_bw() +
  #theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill="white"), axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_grid(migTOTAL + singleSToRescale_PAC ~ singleSToRescale_MAC, scale = "free");
  facet_nested(migTOTAL + singleSToRescale_PAC ~ singleSToRescale_MAC, scale = "free_y",
               labeller = label_parsed);

#manually setting the ticks on y-axis 
#got example from here: https://community.rstudio.com/t/separate-axis-breaks-for-facets/10352/2
my_breaks <- function(x) { if (max(x) > 0.6) seq(0, 0.8, 0.4) 
  else if (max(x) > 0.3) seq(0, 0.5, 0.25)
  else if (max(x) > 0.1) seq(0, 0.1, 0.05)
  else if (max(x) > 0.05) seq(0, 0.05, 0.025)}

th <- ggplot(df_sub, aes(Generation, hybridIndexCOR, colour = S_OD), 
             group = interaction(S_OD,migTOTAL,singleSToRescale_MAC,singleSToRescale_PAC)) +
  geom_line(data = df_sub[df_sub$S_OD != "0%",], aes(Generation, hybridIndexCOR, colour = S_OD), size = 0.7, alpha = 0.7) +
  geom_line(data = df_sub[df_sub$S_OD == "0%",], aes(Generation, hybridIndexCOR, colour = S_OD), size = 1, alpha = 1) +
  #geom_line(size = 0.8) +
  labs(y = "Local Genomic Replacement", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#2c7fb8", #"#41b6c4", "#2c7fb8","#253494", 
                              "#e6550d")) + #"#fd8d3c", "#e6550d", "#a63603")) +
  scale_x_continuous(breaks = c(1,20,40), limits = c(0,40), expand = c(0,0)) +
  scale_y_continuous(breaks = my_breaks, limits = c(0, NA)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), strip.background = element_rect(fill="white"), axis.text.x = element_text(angle = 45, hjust = 1)) +
  #facet_grid(migTOTAL + singleSToRescale_PAC ~ singleSToRescale_MAC, scales = "free");
  facet_nested(migTOTAL + singleSToRescale_PAC ~ singleSToRescale_MAC, scales = "free",
               labeller = label_parsed);

legend <- get_legend(
  th + theme(legend.position = "bottom")
)

g <- plot_grid(ts + theme(legend.position = "none"), 
               th + theme(legend.position = "none"), 
               labels = c("A", "B"),
               ncol = 2, nrow = 1,
               rel_widths = c(1,1))
all <- plot_grid(g, legend,
                 nrow = 2,
                 rel_heights = c(1, 0.05))

pdf(file='./SCRIPTandFIGSfinal_2021/Figure1.pdf', width=7,height=5.2)
all;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/Figure1.png', width=7,height=5.2, 
    units="in", res=600)
all;
dev.off()

#### Figure 2 ####

#Single event of migration for m = 0.05 and MAC dominance factor = 0.5
df_sub <- df %>% filter(migGen == "'Single'", 
                        migTOTAL == 'italic(T[f]) * " = 5%"', 
                        dominanceCoef_MAC == 0.5,
                        singleSToRescale_MAC == 'italic(Delta[MA]) ~ "= 50%"', 
                        nbLoci_MAC == "`1 MA`" | nbLoci_MAC == "`5 MAs`", 
                        singleSToRescale_PAC == 'italic(Delta[PA]) ~ "= 0"' | singleSToRescale_PAC == 'italic(Delta[PA]) ~ "= 50%"')

df_sub$S_OD <- factor(df_sub$S_OD,
                  levels=levels(factor(df_sub$S_OD)),
                  labels=c(expression(italic(Delta[OD])~"= 0%"),
                           expression(italic(Delta[OD])~"= 20%"),
                           expression(italic(Delta[OD])~"= 60%")))

df_sub$nbLoci_PAC <- factor(df_sub$nbLoci_PAC,
                            levels=levels(factor(df_sub$nbLoci_PAC)),
                            labels=c(expression("'0 PA'"),
                                     expression("1 PA;"~italic(s)~"= 0.5"),
                                     expression("5 PAs;"~italic(s)%~~%~"0.08"),
                                     expression("50 PAs;"~italic(s)%~~%~"0.008")))

#Single event of migration for m = 0.05 and MAC dominance factor = 0.5
ts <- ggplot(df_sub, aes(Generation, relFitness, 
                         colour = nbPairsLoci_OD), 
             group = interaction(nbPairsLoci_OD, factor(nbLoci_MAC), singleSToRescale_PAC, singleSToRescale_MAC, nbLoci_PAC, S_OD, factor(nbLoci_OD))) +
  geom_line(aes(linetype = nbLoci_MAC), size = 0.6)+
  scale_linetype_discrete(labels = c(expression("1 MA;"~italic(s)~"= -0.5"), 
                                     expression("5 MAs;"~italic(s)%~~%~"-0.13"))) +
  labs(y = "Relative Fitness", colour = "OD loci pairs", linetype = "MA loci") +
  scale_color_manual(values=c("gray20", 
                              "#fc9272", 
                              "#ef3b2c", 
                              "#a50f15"))+
  scale_x_continuous(breaks = c(1, 50, 100), limits = c(0,100), expand = c(0,0)) +
  scale_y_continuous(breaks = c(1, 1.25, 1.5)) +
  theme_bw() +
  geom_hline(yintercept=1, linetype="longdash", color = "grey", size = 0.6) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text.align = 0) +
  theme(strip.background =element_rect(fill="white")) +
  facet_nested(nbLoci_PAC ~ S_OD, 
               labeller = label_parsed); ts

pdf(file="./SCRIPTandFIGSfinal_2021/Figure2.pdf", width=7,height=5.3)
ts
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/Figure2.png', width=7,height=5.3, 
    units="in", res=600)
ts
dev.off()

#### Figure 3 ####
df_sub <- df %>% filter(dominanceCoef_MAC == 0.5, 
                        nbLoci_MAC == "`5 MAs`", 
                        migTOTAL == 'italic(T[f]) * " = 5%"', 
                        nbLoci_PAC == "`5 PAs`",
                        singleSToRescale_PAC == 'italic(Delta[PA]) ~ "= 10%"',
                        singleSToRescale_MAC == 'italic(Delta[MA]) %~~% ~"9%"', 
                        nbPairsLoci_OD == "0 OD pairs" | nbPairsLoci_OD == "10 OD pairs",
                        Generation <= 40)

df_min <- read.csv("N10000_IndividualReps.csv", header = T, stringsAsFactors = F)
#transforming NA in 0
df_min[is.na(df_min)] <- 0
df_min <- df_min %>% filter(singleSToRescale_MAC != 0.2, nbPairsLoci_OD != 20)

#renaming labels in df_min
df_min$S_OD <- factor(df_min$S_OD,
                      levels=levels(factor(df_min$S_OD)),
                      labels=c("0%", "20%", "60%", "20%", "60%", "20%", "60%"))

#digesting the migration for single vs. multiple pulses
df_min$mig <- str_count(df_min$migration_scenario, ',') #1 = single pulse, 9 = 5 pulses 
df_min$mig[df_min$mig == 9] <- 5
df_min$migGen <- ifelse(df_min$mig == 5, 
                        as.numeric(sapply(strsplit(df_min$migration_scenario, ","), tail, 2)[2,]) - as.numeric(sapply(strsplit(df_min$migration_scenario, ","), tail, 2)[1,]),0)
df_min$mig_split <- as.numeric(sapply(strsplit(df_min$migration_scenario, ","), "[[", 1))
df_min$migTOTAL <- df_min$mig_split*df_min$mig

df_min$migGen <- factor(df_min$migGen, 
                        levels=levels(factor(df_min$migGen)),
                        labels=c("Single", "Every 1", "Every 2", "Every 4"))
df_min$mig <- factor(df_min$mig, levels=levels(factor(df_min$mig)),
                     labels=c("Migration Single", "Migration Pulse (5x)"))

df_min_sub <- df_min %>% filter(dominanceCoef_MAC == 0.5, 
                                migTOTAL == 0.05, 
                                nbPairsLoci_OD == 0 | nbPairsLoci_OD == 10) %>% 
  group_by(rep, S_OD, migGen) %>% 
  summarise(min_Fit = min(relFitness))

#renaming labels necessary for correct format
df_sub$migGen <- factor(df_sub$migGen, 
                    levels=levels(factor(df_sub$migGen)),
                    labels=c("Single", "Every 1", "Every 2", "Every 4"))

ts <- ggplot(df_sub, aes(Generation, relFitness, 
                         colour = S_OD),
             groups = interaction(S_OD,migTOTAL,factor(nbPairsLoci_OD),migGen,singleSToRescale_MAC,singleSToRescale_PAC)) +
  geom_hline(yintercept=1, linetype="longdash", color = "grey", size = 0.8) +
  geom_line(data = df_sub[df_sub$S_OD != "0%",], aes(Generation, relFitness, colour = S_OD), size = 0.7, alpha = 0.7) +
  geom_line(data = df_sub[df_sub$S_OD == "0%",], aes(Generation, relFitness, colour = S_OD), size = 1, alpha = 1) +
  #geom_line(size = 0.8)+
  labs(y = "Relative Fitness", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#2c7fb8",
                              "#e6550d"))+ 
  scale_x_continuous(breaks = c(1, 20, 40), limits = c(0,40), expand = c(0,0)) +
  scale_y_continuous(limits = c(0.93,1.025)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))+ #center title
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        plot.margin = margin(t=5,r=5,b=0,l=20)) +
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(. ~ migGen);

th <- ggplot(df_sub, aes(Generation, hybridIndexCOR, 
                         colour = S_OD),
             groups = interaction(S_OD,migTOTAL,factor(nbPairsLoci_OD),migGen,singleSToRescale_MAC,singleSToRescale_PAC)) +
  geom_line(data = df_sub[df_sub$S_OD != "0%",], aes(Generation, hybridIndexCOR, colour = S_OD), size = 0.7, alpha = 0.7) +
  geom_line(data = df_sub[df_sub$S_OD == "0%",], aes(Generation, hybridIndexCOR, colour = S_OD), size = 1, alpha = 1) +
  #geom_line(size = 0.8)+
  labs(y = "Local Genomic    \nReplacement    ", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#2c7fb8", 
                              "#e6550d"))+ 
  scale_x_continuous(breaks = c(1, 20, 40), limits = c(0,40), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.025, 0.05)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5), plot.margin = margin(t=3,r=5,b=5,l=8))+ #center title
  theme(strip.background =element_rect(fill="white")) +
  facet_grid(. ~ migGen) +
  theme(strip.text.x = element_blank());

legend1 <- get_legend(
  th + theme(legend.position = "bottom")
)

AB <- plot_grid(ts + theme(legend.position = "none"), 
                th + theme(legend.position = "none"),
                #legend1, 
                labels = c("A","B"),
                ncol = 1, nrow = 2,
                rel_heights = c(0.9,1));

bp <- ggplot(df_min_sub, aes(x = S_OD, y = min_Fit)) +
  geom_hline(yintercept=1, linetype="longdash", color = "grey") +
  geom_boxplot(aes(fill = migGen),
               position = position_dodge(0.9)) +
  scale_fill_manual(values=c("#c994c7", 
                             "#ffffe5", "#d9f0a3", 
                             "#78c679")) +
  xlab(expression(italic(Delta[OD]))) +
  ylab("Minimum Fitness Experienced") +
  labs(fill = "Pulsed scenario") +
  theme_bw();

legend2 <- get_legend(
  bp + theme(legend.position = "bottom")
)

g <- plot_grid(AB, 
               bp + theme(legend.position = "none"), 
               labels = c("", "C"),
               ncol = 2, nrow = 1,
               rel_widths = c(1,0.7));

all <- plot_grid(g, 
                 legend1,
                 legend2,
                 ncol = 1, nrow = 3,
                 rel_heights = c(1, 0.09, 0.08)); all

pdf(file='./SCRIPTandFIGSfinal_2021/Figure3.pdf', width=7,height=4)
all;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/Figure3.png', width=7,height=4, 
    units="in", res=600)
all;
dev.off()

########## SUPPLEMENTARY PLOTS ############
#### Figure S1 ***TOM***####


#### Figure S2 #####
df_sub <- df %>% filter(migGen == "'Single'",
                        migTOTAL == 'italic(T[f]) * \" = 5%\"',
                        dominanceCoef_MAC == 0.5)

ts <- ggplot(df_sub, aes(Generation, relFitness, 
                         colour = factor(S_OD):factor(nbPairsLoci_OD)), 
             group = interaction(S_OD,factor(nbLoci_OD),factor(migTOTAL),factor(singleSToRescale_MAC))) +
  geom_line()+
  labs(y = "Relative Fitness", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#a1dab4","#41b6c4",
                              "#253494", 
                              "#fdbe85", "#fd8d3c",
                              "#a63603")) +
  scale_x_continuous(breaks = c(0,50,100), limits = c(0,100), expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") + #center title
  theme(strip.background =element_rect(fill="white")) +
  facet_nested(nbLoci_PAC + singleSToRescale_PAC ~ nbLoci_MAC + singleSToRescale_MAC, scale = 'free',
               labeller = label_parsed); ts

pdf(file='./SCRIPTandFIGSfinal_2021/SuppFigure2.pdf', width=9,height=7)
ts;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/SuppFigure2.png', width=9,height=7, 
    units="in", res=600)
ts;
dev.off()

#### Figure S4 ####
th <- ggplot(df_sub, aes(Generation, hybridIndexCOR, 
                         colour = factor(S_OD):factor(nbPairsLoci_OD)), 
             group = interaction(S_OD,factor(nbLoci_OD),factor(migTOTAL),factor(singleSToRescale_MAC))) +
  geom_line()+
  labs(y = "Local Genomic Replacement", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#a1dab4","#41b6c4",
                              "#253494", 
                              "#fdbe85", "#fd8d3c",
                              "#a63603")) +
  scale_x_continuous(breaks = c(0,50,100), limits = c(0,100), expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right") + #center title
  theme(strip.background =element_rect(fill="white")) +
  facet_nested(nbLoci_PAC + singleSToRescale_PAC ~ nbLoci_MAC + singleSToRescale_MAC, scale = 'free',
               labeller = label_parsed); th

pdf(file='./SCRIPTandFIGSfinal_2021/SuppFigure4.pdf', width=9,height=7)
th;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/SuppFigure4.png', width=9,height=7, 
    units="in", res=600)
th;
dev.off()

#### Figure S3 ####
#Single migration event with MAC dominance of 0.5 to compare different m rates
df_sub <- df %>% filter(migGen == "'Single'",
                        dominanceCoef_MAC == 0.5,
                        nbLoci_MAC == "`5 MAs`")

ts <- ggplot(df_sub, aes(Generation, relFitness, 
                         colour = factor(S_OD):factor(nbPairsLoci_OD),
                         shape = ), 
             group = interaction(S_OD,factor(nbLoci_OD),factor(migTOTAL),factor(singleSToRescale_MAC))) +
  geom_line()+
  labs(y = "Relative Fitness", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#a1dab4", "#41b6c4","#253494", 
                              "#fdbe85", "#fd8d3c","#a63603")) +
  scale_x_continuous(breaks = c(0,50,100), limits = c(0,100), expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background =element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_nested(factor(migTOTAL) + factor(singleSToRescale_MAC) ~ factor(nbLoci_PAC) + factor(singleSToRescale_PAC), 
               scale = 'free', labeller = label_parsed); ts

pdf(file='./SCRIPTandFIGSfinal_2021/SuppFigure3.pdf', width=9,height=7)
ts;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/SuppFigure3.png', width=9,height=7, 
    units="in", res=600)
ts;
dev.off()

#### Figure S5 ####
tr <- ggplot(df_sub, aes(Generation, hybridIndexCOR, 
                         colour = factor(S_OD):factor(nbPairsLoci_OD),
                         shape = ), 
             group = interaction(S_OD,factor(nbLoci_OD),factor(migTOTAL),factor(singleSToRescale_MAC))) +
  geom_line()+
  labs(y = "Local Genomic Replacement", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#a1dab4", "#41b6c4","#253494", 
                              "#fdbe85", "#fd8d3c","#a63603")) +
  scale_x_continuous(breaks = c(0,50,100), limits = c(0,100), expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background =element_rect(fill="white"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_nested(factor(migTOTAL) + factor(singleSToRescale_MAC) ~ factor(nbLoci_PAC) + factor(singleSToRescale_PAC), 
               scale = 'free', labeller = label_parsed); tr

pdf(file='./SCRIPTandFIGSfinal_2021/SuppFigure5.pdf', width=9,height=7)
tr;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/SuppFigure5.png', width=9,height=7, 
    units="in", res=600)
tr;
dev.off()

#### Figure S6 #### 
df_sub <- df %>% filter(migGen == "'Single'",
                        migTOTAL == 'italic(T[f]) * \" = 5%\"', 
                        nbLoci_MAC == "`5 MAs`")

df_sub$dominanceCoef_MAC <- factor(df_sub$dominanceCoef_MAC, 
                                   levels=levels(factor(df_sub$dominanceCoef_MAC)),
                                   labels = c("`MA dominance = 0`", "`MA dominance = 0.5`"))

dh <- ggplot(df_sub, aes(Generation, hybridIndexCOR, 
                         colour = factor(S_OD):factor(nbPairsLoci_OD)), 
             group = interaction(S_OD,factor(nbLoci_OD),factor(migTOTAL),factor(singleSToRescale_MAC))) +
  geom_line()+
  labs(y = "Local Genomic Replacement", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#a1dab4", "#41b6c4", "#253494", 
                              "#fdbe85", "#fd8d3c", "#a63603")) +
  scale_x_continuous(breaks = c(0,50,100), limits = c(0,100), expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.background =element_rect(fill="white")) +
  facet_nested(factor(nbLoci_PAC) + factor(singleSToRescale_PAC) ~ factor(dominanceCoef_MAC) + factor(singleSToRescale_MAC), 
               scale = 'free', labeller = label_parsed); dh

pdf(file='./SCRIPTandFIGSfinal_2021/SuppFigure6.pdf', width=9,height=7)
dh;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/SuppFigure6.png', width=9,height=7, 
    units="in", res=600)
dh;
dev.off()

#### Figure S7 ####
d <- ggplot(df_sub, aes(Generation, relFitness, 
                        colour = factor(S_OD):factor(nbPairsLoci_OD)), 
            group = interaction(S_OD,factor(nbLoci_OD),factor(migTOTAL),factor(singleSToRescale_MAC))) +
  geom_line()+
  labs(y = "Relative Fitness", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#a1dab4", "#41b6c4", "#253494", 
                              "#fdbe85", "#fd8d3c", "#a63603")) +
  scale_x_continuous(breaks = c(0,50,100), limits = c(0,100), expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.background =element_rect(fill="white")) +
  facet_nested(factor(nbLoci_PAC) + factor(singleSToRescale_PAC) ~ factor(dominanceCoef_MAC) + factor(singleSToRescale_MAC), 
               scale = 'free', labeller = label_parsed); d

pdf(file='./SCRIPTandFIGSfinal_2021/SuppFigure7.pdf', width=9,height=7)
d;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/SuppFigure7.png', width=9,height=7, 
    units="in", res=600)
d;
dev.off()

#### Figure S8 ####
#Single and pulses of migration event with 5 MAC loci with dominance of 0.05 
df_sub <-  df %>% filter(dominanceCoef_MAC == 0.5,
                         nbLoci_MAC == "`5 MAs`",
                         migTOTAL == 'italic(T[f]) * " = 5%"')

ts <- ggplot(df_sub, aes(Generation, relFitness, 
                         colour = factor(S_OD):factor(nbPairsLoci_OD)), 
             group = interaction(S_OD,factor(nbLoci_OD),factor(migTOTAL),factor(singleSToRescale_MAC))) +
  geom_line()+
  labs(y = "Relative Fitness", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#a1dab4", "#41b6c4", "#253494", 
                              "#fdbe85", "#fd8d3c", "#a63603")) +
  scale_x_continuous(breaks = c(0,50,100), limits = c(0,100), expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background =element_rect(fill="white")) +
  facet_nested(factor(nbLoci_PAC) + factor(singleSToRescale_PAC) ~ factor(migGen) + factor(singleSToRescale_MAC), 
               scale = 'free', labeller = label_parsed); ts

pdf(file='./SCRIPTandFIGSfinal_2021/SuppFigure8.pdf', width=10.5,height=7)
ts;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/SuppFigure8.png', width=9,height=7, 
    units="in", res=600)
ts;
dev.off()

#### Figure S9 ####
ts <- ggplot(df_sub, aes(Generation, hybridIndexCOR, 
                         colour = factor(S_OD):factor(nbPairsLoci_OD)), 
             group = interaction(S_OD,factor(nbLoci_OD),factor(migTOTAL),factor(singleSToRescale_MAC))) +
  geom_line()+
  labs(y = "Local Genomic Replacement", colour = expression(italic(Delta[OD]))) +
  scale_color_manual(values=c("black", 
                              "#a1dab4", "#41b6c4", "#253494", 
                              "#fdbe85", "#fd8d3c", "#a63603")) +
  scale_x_continuous(breaks = c(0,50,100), limits = c(0,100), expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background =element_rect(fill="white")) +
  facet_nested(factor(nbLoci_PAC) + factor(singleSToRescale_PAC) ~ factor(migGen) + factor(singleSToRescale_MAC), 
               scale = 'free', labeller = label_parsed); ts

pdf(file='./SCRIPTandFIGSfinal_2021/SuppFigure9.pdf', width=10.5,height=7)
ts;
dev.off()

png(file='./SCRIPTandFIGSfinal_2021/SuppFigure9.png', width=10.5,height=7, 
    units="in", res=600)
ts;
dev.off()
