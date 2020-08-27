rm(list = ls())
library(ggplot2)
library(ggthemes)


all <- read.csv('~/UBC/GeneticRescue/RemakeAitkenWhitlockPlot/mean.allReplicates.csv')

all$S_OD[is.na( all$S_OD ) ]<- 0

library(ggplot2)
library(ggthemes)

all$nbLoci_PAC <- factor(as.factor(all$nbLoci_PAC), levels = c(0, 5), labels = c("Outbreeding depression only (~20% in F1)", "Outbreeding depression + 5 preadapted climate alleles"))

max(all$relFitness)

ggplot( data = all, aes( x = Generation, y = relFitness, col = as.factor(S_OD)))+
  geom_hline(aes(yintercept = 1), lty = 2, alpha = 0.5)+
  geom_line(lwd = 1.3)+
  scale_color_colorblind("Strength of selection\non outbreeding depression")+
  facet_wrap( as.factor(nbLoci_PAC)~., scales = 'free_y')+
  scale_x_continuous(limits=c(1,100))+
  scale_y_continuous("Fitness\n(Relative to Ancestral Population)")+
    theme_bw()
?facet_wrap

all <- all[all$nbLoci_PAC == 0,]

mA1 <- all[all$nbPairsLoci_OD == 2,]
mA1$experiment <- 'Outbreeding Depression'
mA1$OD <- '2 pairs'
mA1$simulator = 'Remi'

mA2 <- all[all$nbPairsLoci_OD == 10,]
mA2$experiment <- 'Outbreeding Depression'
mA2$OD <- '10 pairs'
mA2$simulator = 'Remi'

mA3 <- all[all$nbPairsLoci_OD == 20,]
mA3$experiment <- 'Outbreeding Depression'
mA3$OD <- '20 pairs'
mA3$simulator = 'Remi'

mA4 <- all[all$nbPairsLoci_OD == 100,]
mA4$experiment <- 'Outbreeding Depression'
mA4$OD <- '100 pairs'
mA4$simulator = 'Remi'


mike_A1 <- c(1., 0.9891, 0.987644, 0.987167, 0.987333, 0.987444, 0.9886, 0.989778, 
             0.990767, 0.991256, 0.992222, 0.993256, 0.993978, 0.994533, 0.995378, 
             0.996167, 0.996311, 0.996822, 0.997222, 0.9976, 0.997767, 0.998067, 
             0.998333, 0.998456, 0.998644, 0.998922, 0.999, 0.999189, 0.999522, 
             0.999644, 0.999678, 0.999711, 0.999744, 0.9998, 0.999767, 0.999789, 
             0.999867, 0.999878, 0.9999, 0.9999, 0.9999, 0.999911, 0.999944, 
             0.9999, 0.999911, 0.999956, 0.999978, 0.999978, 0.999967, 0.999978, 
             0.999967, 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
             1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
             1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.)
A1 <- data.frame( relFitness = mike_A1, experiment <- 'Outbreeding Depression', OD = '2 loci - 1 pair', Generation = 0:100, simulator = 'Mike')

mike_A2 <- c(1., 0.989187, 0.985478, 0.984119, 0.983916, 0.9841, 0.984178, 
             0.984527, 0.985023, 0.985193, 0.985623, 0.98584, 0.986255, 0.986269, 
             0.986481, 0.986774, 0.986929, 0.987349, 0.987894, 0.987823, 0.988352, 
             0.988752, 0.988974, 0.98928, 0.989368, 0.98919, 0.989463, 0.989704, 
             0.990036, 0.990255, 0.990542, 0.990846, 0.991256, 0.991171, 0.991587, 
             0.992023, 0.992374, 0.992398, 0.992504, 0.99264, 0.992879, 0.993, 
             0.993284, 0.993362, 0.993633, 0.993818, 0.993817, 0.993866, 0.994179, 
             0.994369, 0.994346, 0.9946, 0.994672, 0.994849, 0.994838, 0.994855, 
             0.995201, 0.995356, 0.99537, 0.995578, 0.995731, 0.995886, 0.99605, 
             0.996236, 0.996337, 0.99627, 0.996385, 0.996621, 0.99664, 0.996746, 
             0.996644, 0.996741, 0.996997, 0.997013, 0.997058, 0.997152, 0.997245, 
             0.99729, 0.997397, 0.997542, 0.997612, 0.997622, 0.997703, 0.997786, 
             0.997781, 0.997787, 0.997812, 0.997841, 0.997875, 0.997814, 0.99799, 
             0.998039, 0.9981, 0.997965, 0.998111, 0.998282, 0.99832, 0.998325, 
             0.998354, 0.998413, 0.998501)
A2 <- data.frame( relFitness = mike_A2, experiment <- 'Outbreeding Depression', OD = '10 loci - 5 pairs', Generation = 0:100, simulator = 'Mike')
str(A2)






mike_A3 <- c(1., 0.989159, 0.985859, 0.984357, 0.984168, 0.983982, 0.984231, 
0.984525, 0.98467, 0.984985, 0.985105, 0.985242, 0.98543, 0.985728, 
0.986014, 0.985886, 0.986147, 0.985918, 0.985966, 0.985948, 0.986151, 
0.986313, 0.986296, 0.986568, 0.986666, 0.986792, 0.987004, 0.987141, 
0.987258, 0.987513, 0.987856, 0.987917, 0.988134, 0.988231, 0.988299, 
0.988268, 0.988499, 0.988617, 0.988771, 0.989044, 0.989112, 0.989173, 
0.989066, 0.989235, 0.98923, 0.989183, 0.989205, 0.989215, 0.989261, 
0.989539, 0.989685, 0.989668, 0.989729, 0.990027, 0.990287, 0.990398, 
0.990608, 0.990732, 0.990885, 0.990853, 0.99089, 0.99098, 0.991331, 
0.991603, 0.991716, 0.991822, 0.992086, 0.992095, 0.992269, 0.99243, 
0.9927, 0.992861, 0.992805, 0.992641, 0.992699, 0.992767, 0.992834, 
0.992915, 0.992885, 0.992885, 0.992994, 0.992942, 0.993074, 0.993112, 
0.993245, 0.993447, 0.993481, 0.993443, 0.993524, 0.993515, 0.993696, 
0.993553, 0.993455, 0.993653, 0.993816, 0.993895, 0.994027, 0.994045, 
0.994176, 0.99419, 0.994246)
A3 <- data.frame( relFitness= mike_A3, experiment <- 'Outbreeding Depression', OD = '20 loci - 10 pairs', Generation = 0:100, simulator = 'Mike')




A1_plot <- ggplot(data = A1, aes( x = Generation, y = relFitness, col = simulator))+
  geom_line( alpha = 0.75, lwd = 1.3)+
  geom_line( data = mA1, aes( x = Generation-1, y = relFitness, col = simulator), lwd = 2, alpha = 0.75)+
  scale_y_continuous('Fitness relative to ancestral population')+
  scale_color_colorblind('Simulator')+
  geom_hline(yintercept = 1, lty = 2, alpha = 0.5)+
  ggtitle(A1$OD)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    axis.text.y =  element_text(size = 10, color = 'black'),
    strip.text.y =  element_text(size = 12.5, color = 'black'),
    axis.text.x =  element_text(size = 10, color = 'black'),
    legend.text = element_text(size = 10, color = 'black'),
    axis.title.x = element_text(size = 12.5, color = 'black'),
    axis.title.y = element_text(size = 12.5,  vjust = 1.5, color = 'black'),
    panel.grid = element_line(),         # All grid lines
    panel.grid.major = element_line(),   # Major grid lines
    panel.grid.minor = element_line(),   # Minor grid lines
    panel.grid.major.y = element_blank(), # Horizontal major grid lines
    panel.grid.minor.y = element_blank(),  # Vertical major grid lines
    panel.grid.major.x = element_blank(), # Horizontal major grid lines
    panel.grid.minor.x = element_blank()  # Vertical major grid lines
  )


A2_plot <- ggplot(data = A2, aes( x = Generation, y = relFitness, col = simulator))+
  geom_line( alpha = 0.75, lwd = 2)+
  geom_line( data = mA2, aes( x = Generation-1, y = relFitness, col = simulator), lwd = 2, alpha = 0.75)+
  scale_y_continuous('Fitness relative to ancestral population')+
  scale_color_colorblind('Simulator')+
  geom_hline(yintercept = 1, lty = 2, alpha = 0.5)+
  ggtitle(A2$OD)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    axis.text.y =  element_text(size = 10, color = 'black'),
    strip.text.y =  element_text(size = 12.5, color = 'black'),
    axis.text.x =  element_text(size = 10, color = 'black'),
    legend.text = element_text(size = 10, color = 'black'),
    axis.title.x = element_text(size = 12.5, color = 'black'),
    axis.title.y = element_text(size = 12.5,  vjust = 1.5, color = 'black'),
    panel.grid = element_line(),         # All grid lines
    panel.grid.major = element_line(),   # Major grid lines
    panel.grid.minor = element_line(),   # Minor grid lines
    panel.grid.major.y = element_blank(), # Horizontal major grid lines
    panel.grid.minor.y = element_blank(),  # Vertical major grid lines
    panel.grid.major.x = element_blank(), # Horizontal major grid lines
    panel.grid.minor.x = element_blank()  # Vertical major grid lines
  )



A3_plot <- ggplot(data = A3, aes( x = Generation, y = relFitness, col = simulator))+
  geom_line( alpha = 0.75, lwd = 2)+
  geom_line( data = mA3, aes( x = Generation-1, y = relFitness, col = simulator), lwd = 2, alpha = 0.75)+
  scale_y_continuous('Fitness relative to ancestral population')+
  scale_color_colorblind('Simulator')+
  geom_hline(yintercept = 1, lty = 2, alpha = 0.5)+
  ggtitle(A3$OD)+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    axis.text.y =  element_text(size = 10, color = 'black'),
    strip.text.y =  element_text(size = 12.5, color = 'black'),
    axis.text.x =  element_text(size = 10, color = 'black'),
    legend.text = element_text(size = 10, color = 'black'),
    axis.title.x = element_text(size = 12.5, color = 'black'),
    axis.title.y = element_text(size = 12.5,  vjust = 1.5, color = 'black'),
    panel.grid = element_line(),         # All grid lines
    panel.grid.major = element_line(),   # Major grid lines
    panel.grid.minor = element_line(),   # Minor grid lines
    panel.grid.major.y = element_blank(), # Horizontal major grid lines
    panel.grid.minor.y = element_blank(),  # Vertical major grid lines
    panel.grid.major.x = element_blank(), # Horizontal major grid lines
    panel.grid.minor.x = element_blank()  # Vertical major grid lines
  )



library(ggpubr)
ggarrange(  A1_plot, A2_plot, A3_plot, nrow = 1, ncol = 3 , common.legend = T)





