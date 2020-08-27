rm(list = ls())

## Read in the data

mA1 <- read.csv('~/work/GeneticRescue/simulations/Simulations/mean.A.0.1.1.csv')
mA1$experiment <- 'Outbreeding Depression'
mA1$OD <- '2 loci'
aA1 <- read.csv('~/work/GeneticRescue/simulations/Simulations/all.A.0.1.1.csv')
aA1$experiment <- 'Outbreeding Depression'
aA1$OD <- '2 loci'

mA2 <- read.csv('~/work/GeneticRescue/simulations/Simulations/mean.A.0.1.2.csv')
mA2$experiment <- 'Outbreeding Depression'
mA2$OD <- '10 loci'
aA2 <- read.csv('~/work/GeneticRescue/simulations/Simulations/all.A.0.1.2.csv')
aA2$experiment <- 'Outbreeding Depression'
aA2$OD <- '10 loci'

mA3 <- read.csv('~/work/GeneticRescue/simulations/Simulations/mean.A.0.1.3.csv')
mA3$experiment <- 'Outbreeding Depression'
mA3$OD <- '20 loci'
aA3 <- read.csv('~/work/GeneticRescue/simulations/Simulations/all.A.0.1.3.csv')
aA3$experiment <- 'Outbreeding Depression'
aA3$OD <- '20 loci'

mA4 <- read.csv('~/work/GeneticRescue/simulations/Simulations/mean.A.0.1.4.csv')
mA4$experiment <- 'Outbreeding Depression'
mA4$OD <- '100 loci'
aA4 <- read.csv('~/work/GeneticRescue/simulations/Simulations/all.A.0.1.4.csv')
aA4$experiment <- 'Outbreeding Depression'
aA4$OD <- '100 loci'


### Outbreeding Depression  +  preadapted alleles

mA5 <- read.csv('~/work/GeneticRescue/simulations/Simulations/mean.A.0.1.5.csv')
mA5$experiment <- 'Outbreeding Depression + preadapted alleles'
mA5$OD <- '2 loci'
aA5 <- read.csv('~/work/GeneticRescue/simulations/Simulations/all.A.0.1.5.csv')
aA5$experiment <- 'Outbreeding Depression + preadapted alleles'
aA5$OD <- '2 loci'

mA6 <- read.csv('~/work/GeneticRescue/simulations/Simulations/mean.A.0.1.6.csv')
mA6$experiment <- 'Outbreeding Depression + preadapted alleles'
mA6$OD <- '10 loci'
aA6 <- read.csv('~/work/GeneticRescue/simulations/Simulations/all.A.0.1.6.csv')
aA6$experiment <- 'Outbreeding Depression + preadapted alleles'
aA6$OD <- '10 loci'

mA7 <- read.csv('~/work/GeneticRescue/simulations/Simulations/mean.A.0.1.7.csv')
mA7$experiment <- 'Outbreeding Depression + preadapted alleles'
mA7$OD <- '20 loci'
aA7 <- read.csv('~/work/GeneticRescue/simulations/Simulations/all.A.0.1.7.csv')
aA7$experiment <- 'Outbreeding Depression + preadapted alleles'
aA7$OD <- '20 loci'

mA8 <- read.csv('~/work/GeneticRescue/simulations/Simulations/mean.A.0.1.8.csv')
mA8$experiment <- 'Outbreeding Depression + preadapted alleles'
mA8$OD <- '100 loci'
aA8 <- read.csv('~/work/GeneticRescue/simulations/Simulations/all.A.0.1.8.csv')
aA8$experiment <- 'Outbreeding Depression + preadapted alleles'
aA8$OD <- '100 loci'


mA_OD <- rbind(mA1, mA2, mA3, mA4)
mA_ODa <- rbind(mA5, mA6, mA7, mA8)

aA <- rbind(aA1, aA2, aA3, aA4, aA5, aA6, aA7, aA8)

library(ggplot2)
library(ggthemes)

aA$OD <- factor(aA$OD, levels = c('2 loci', '10 loci', '20 loci', '100 loci'),  labels = c('2 Loci', '10 Loci', '20 Loci', '100 Loci'))
mA_OD$OD <- factor(mA_OD$OD, levels = c('2 loci', '10 loci', '20 loci', '100 loci'),  labels = c('2 Loci', '10 Loci', '20 Loci', '100 Loci'))
mA_ODa$OD <- factor(mA_ODa$OD, levels = c('2 loci', '10 loci', '20 loci', '100 loci'),  labels = c('2 Loci', '10 Loci', '20 Loci', '100 Loci'))


ggplot(data = aA, aes( x = Generation, y = relFitness, col = OD, group = rep))+
  geom_line( alpha = 0.5)+
  facet_grid( OD ~ experiment )+
  scale_color_colorblind('# Overdominant loci')+
  theme_bw()


ggplot(data = aA, aes( x = Generation, y = meanHI, col = OD, group = rep))+
  geom_line( alpha = 0.5)+
  facet_grid( OD ~ experiment )+
  scale_color_colorblind('# Overdominant loci')+
  theme_bw()

plot_A <- ggplot(data = mA_OD, aes( x = Generation, y = relFitness, col = OD))+
  geom_line( alpha = 0.75, lwd = 2)+
  scale_y_continuous('Fitness relative to ancestral population')+
  scale_color_colorblind('# Overdominant loci')+
  geom_hline(yintercept = 1, lty = 2, alpha = 0.5)+
  ggtitle('Outbreeding depression')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    axis.text.y =  element_text(size = 10, color = 'black'),
    strip.text.y =  element_text(size = 12.5, color = 'black'),
    axis.text.x =  element_text(size = 10, color = 'black'),
    legend.title = element_blank(),
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


plot_B <- ggplot(data = mA_ODa, aes( x = Generation, y = relFitness, col = OD))+
  geom_line( alpha = 0.75, lwd = 2)+
  scale_y_continuous('Fitness relative to ancestral population')+
  scale_color_colorblind('# Overdominant loci')+
  geom_hline(yintercept = 1, lty = 2, alpha = 0.5)+
  ggtitle('Outbreeding depression + preadapted climate alleles')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    axis.text.y =  element_text(size = 10, color = 'black'),
    strip.text.y =  element_text(size = 12.5, color = 'black'),
    axis.text.x =  element_text(size = 10, color = 'black'),
    legend.title = element_blank(),
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

output <- ggarrange(plot_A, plot_B, ncol = 2, nrow = 1, common.legend = T, legend = 'bottom', labels	= 'auto')

pdf('~/work/GeneticRescue/simulations/AitkenWhitlock_remake.pdf', height =  5, width = 10)
print(output)
dev.off()
