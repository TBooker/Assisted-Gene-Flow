##  COnverting N.mean.A.0.6.output.csv files into super set for the shiny app

library(tidyverse)

# Load each a.0.6 file 

#start100 = read.csv("/Users/whitlock/Desktop/AGF_Shiny/A.0.6_processed/100.mean.A.0.6.output.csv")
# NOte - I believe we decided to not include N=100 cases.  We don't have every migratin scenario for all N

start1000 = read.csv("/Users/whitlock/Desktop/AGF_Shiny/A.0.6_processed/1000.mean.A.0.6.output.csv")
start10000 = read.csv("/Users/whitlock/Desktop/AGF_Shiny/A.0.6_processed/10000.mean.A.0.6.output.csv")


#  Drop unnecessary columns
#start100= start100 %>% select(Generation, Fitness, relFitness, hybridIndex,N,S_OD, 
                              # dominanceCoef_MAC,dominanceCoef_PAC,migration_scenario,nbLoci,nbLoci_MAC,nbLoci_OD,nbLoci_PAC,
                              # nbPairsLoci_OD,singleSToRescale_MAC,singleSToRescale_PAC)
start1000= start1000 %>% select(Generation, Fitness, relFitness, hybridIndex,N,S_OD, 
                              dominanceCoef_MAC,dominanceCoef_PAC,migration_scenario,nbLoci,nbLoci_MAC,nbLoci_OD,nbLoci_PAC,
                              nbPairsLoci_OD,singleSToRescale_MAC,singleSToRescale_PAC)
start10000= start10000 %>% select(Generation, Fitness, relFitness, hybridIndex,N,S_OD, 
                              dominanceCoef_MAC,dominanceCoef_PAC,migration_scenario,nbLoci,nbLoci_MAC,nbLoci_OD,nbLoci_PAC,
                              nbPairsLoci_OD,singleSToRescale_MAC,singleSToRescale_PAC)

# Concat
start = rbind( start1000,start10000)
#rm(start100)
rm(start1000)
rm(start10000)

#Add new necessary columns (e.g. pulse type)

## Making a colum for Local Genomic Replacement (LGR) instead of hybridIndex

start$LGR= 1 - start$hybridIndex
start = start %>% select(-hybridIndex)

##  Make a column that contains the per pulse migration rate

start$migrationRate=rep(0,length(start$migration_scenario))

start$migrationRate[which(substr(start$migration_scenario,1,3)=="0.5")]=0.5
start$migrationRate[which(substr(start$migration_scenario,1,3)=="0.1")]=0.1
start$migrationRate[which(substr(start$migration_scenario,1,4)=="0.01")]=0.01
start$migrationRate[which(substr(start$migration_scenario,1,5)=="0.001")]=0.001
start$migrationRate[which(substr(start$migration_scenario,1,4)=="0.05")]=0.05
start$migrationRate[which(substr(start$migration_scenario,1,5)=="0.005")]=0.005

## Number of pulses
start$numPulses = rep(0,length(start$migration_scenario))
start$numPulses[which(start$migration_scenario=="0.001,0.001,0.001,0.001,0.001,1,2,3,4,5")]=5
start$numPulses[which(start$migration_scenario=="0.001,0.001,0.001,0.001,0.001,1,3,5,7,9")]=5
start$numPulses[which(start$migration_scenario=="0.001,0.001,0.001,0.001,0.001,1,5,9,13,17")]=5
start$numPulses[which(start$migration_scenario=="0.005,1")]=1

start$numPulses[which(start$migration_scenario=="0.01,0.01,0.01,0.01,0.01,1,2,3,4,5")]=5
start$numPulses[which(start$migration_scenario=="0.01,0.01,0.01,0.01,0.01,1,3,5,7,9")]=5
start$numPulses[which(start$migration_scenario=="0.01,0.01,0.01,0.01,0.01,1,5,9,13,17")]=5
start$numPulses[which(start$migration_scenario=="0.05,1")]=1

start$numPulses[which(start$migration_scenario=="0.1,0.1,0.1,0.1,0.1,1,2,3,4,5")]=5
start$numPulses[which(start$migration_scenario=="0.1,0.1,0.1,0.1,0.1,1,3,5,7,9")]=5
start$numPulses[which(start$migration_scenario=="0.1,0.1,0.1,0.1,0.1,1,5,9,13,17")]=5
start$numPulses[which(start$migration_scenario=="0.5,1")]=1

##  INterval between pulses
start$pulseInterval = rep(0,length(start$migration_scenario))
start$pulseInterval[which(start$migration_scenario=="0.001,0.001,0.001,0.001,0.001,1,2,3,4,5")]=1
start$pulseInterval[which(start$migration_scenario=="0.001,0.001,0.001,0.001,0.001,1,3,5,7,9")]=2
start$pulseInterval[which(start$migration_scenario=="0.001,0.001,0.001,0.001,0.001,1,5,9,13,17")]=4
start$pulseInterval[which(start$migration_scenario=="0.005,1")]=0

start$pulseInterval[which(start$migration_scenario=="0.01,0.01,0.01,0.01,0.01,1,2,3,4,5")]=1
start$pulseInterval[which(start$migration_scenario=="0.01,0.01,0.01,0.01,0.01,1,3,5,7,9")]=2
start$pulseInterval[which(start$migration_scenario=="0.01,0.01,0.01,0.01,0.01,1,5,9,13,17")]=4
start$pulseInterval[which(start$migration_scenario=="0.05,1")]=0

start$pulseInterval[which(start$migration_scenario=="0.1,0.1,0.1,0.1,0.1,1,2,3,4,5")]=1
start$pulseInterval[which(start$migration_scenario=="0.1,0.1,0.1,0.1,0.1,1,3,5,7,9")]=2
start$pulseInterval[which(start$migration_scenario=="0.1,0.1,0.1,0.1,0.1,1,5,9,13,17")]=4
start$pulseInterval[which(start$migration_scenario=="0.5,1")]=0

# DeltaOD
start$S_OD[which(is.na(start$S_OD))]=0
start$DeltaOD[start$S_OD == "0.009205015" | start$S_OD == "0.04688023" | start$S_OD == "0.09595823" | start$S_OD == "0.5811388"] <- "60%"
start$DeltaOD[start$S_OD == "0.002233927" | start$S_OD == "0.01121965" | start$S_OD == "0.02256518" | start$S_OD == "0.118034"] <- "20%"
start$DeltaOD[start$S_OD == "0"] <- "0%"

# Making column for total replacement fraction
start$totalReplacementFraction = start$migrationRate*start$numPulses

#setting singleSToRescale_MAC and singleSToRescale_PAC values to 0 if no loci of that type
start$singleSToRescale_MAC[which(start$nbLoci_MAC==0)]=0
start$singleSToRescale_PAC[which(start$nbLoci_PAC==0)]=0

#  reducing size of data file by trimming some choices

start = start %>% filter(nbPairsLoci_OD !=20)
start = start %>% filter(singleSToRescale_MAC != 0.2)
start = start %>% filter(Generation<21 | (Generation<31 & Generation>20 & Generation %% 2==0) | Generation %% 5==0)
start =start %>% filter(dominanceCoef_MAC==0.5)

# dropping columns not used in shiny app
start = start %>% select(-migration_scenario, -dominanceCoef_MAC,-dominanceCoef_PAC)

# Put into wide format  and save (probably easiest to have separate files for relative fitness and the LGR)
dataCharGen=start
dataCharGen$genChar= str_pad(as.character(dataCharGen$Generation),3,side="left", pad="0") 
dataCharGen=dataCharGen %>% select(-Generation) 
dataCharGen$genChar=paste("g",dataCharGen$genChar, sep="")
  
wideVersion_relFitness = dataCharGen %>% select(-LGR,-Fitness) %>% spread(key=genChar, value = relFitness)
wideVersion_LGR = dataCharGen %>% select(-relFitness,-Fitness) %>% spread(key=genChar, value = LGR)

write_csv(wideVersion_relFitness,"/Users/whitlock/Desktop/AGF_Shiny/supersetMeansRelFitnessA06.csv")
write_csv(wideVersion_LGR,"/Users/whitlock/Desktop/AGF_Shiny/supersetMeansReplacementA06.csv")


#  Let's try to reconstruct the long format:

test = gather(wideVersion_relFitness, key=genChar, value = relFitness, "g001":"g100", factor_key=TRUE)
test$genChar=substr(test$genChar,2,4)
test$Generation=as.numeric(test$genChar)
head(test)


