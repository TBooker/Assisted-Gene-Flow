#!/usr/bin/Rscript --vanilla

args <- commandArgs( trailingOnly = T)
## This little script converts Remi's simulation list file into something that is more handy to use in Python
# Read in Remi's List
remiList <- dget(args[1])

# Make a vector that is nrows long
newVec = array(dim = nrow(remiList))

# convert the migration scenarios into a comma-separated-format
for (i in seq(nrow(remiList))){
  newVec[i] = paste(unlist(remiList$migration_scenario[[i]] ), collapse = ',')
} 

remiList$migration_scenario <- newVec

# Save the dataframe as a CSV
write.csv(remiList, file  = args[2])

