# Whitlock Lab Assisted Gene Flow Project

![](newRescue/SuppFigure1.png)
This repository contains the code for running, analysing and plotting the SimBit simulations we performed lab for our Assisted Gene Flow project.

Remi Matthey-Doret ran a whole lot of simuations (hundreds of thousands) of short simulations modelling assisted gene flow as a means to acheive genetic rescue. Remi used SimBit to perform these simualtions. For this project, [SimBit](https://github.com/RemiMattheyDoret/SimBit) has several advantages over other simulation packages (e.g. SLiM or SFS_CODE). SimBit can model whole chromosomes relatively cheaply as it does not necessarily assume a discrete nucleotide model like SLiM does. Add to that, pairs of epistatic alleles are natively supported and that's a large part of what we were interested in. 

We did a semi-factorial simulation design, which resulted in lots and lots of parameter combinations (many 10s of thousands). This produced a whole heck of a lot of results files. We ran each parameter combination 50 times and each replicate produced 2 output files. This gave millions of individual files! This is not the ideal way to store the output as each time a computer has to look into output directory, unless you are careful, it enumerates all the files and this can waste a lot of time. This can be a huge slow down when you are trying to analyse the results. If we were to do this again, it would be far preferable to output all the results into a single (or small number) of files. If you want to repeat the simulations and want to avoid this issue, do get in touch and we'll figure out how to make this quicker for you. However, hindsight is 20:20, and below we outline how we analysed the results we had.

## Analysis pipeline

The simulation output files are all stored at WEBSITE, they are far too big to host on GitHub. Once you have them, it is assumed that they are stored in a directory called  ```A.0.6/outputs```. The name A.0.6 comes from the fact that this was the sixth version of the simulations that we ran. This file architechture is present in this directory, but we do not include the raw files here on GitHub.

```
cd A.0.6/

# The first thing that I did was convert Remi's list of parameter combinations into an easily navigable CSV file. 
Rscript --vanilla ../bin/convertRemiList.R paramGrid_A.0.6.txt paramGrid_A.0.6.csv

# Then, I moved all the output files into directories that corresponded to their parameter combination.
sh ../bin/mkdirAllFiles.sh paramGrid_A.0.6.csv

## this takes a while...
```

Once that's done, we run a short AWK script to summarise the results from each of the simulation sets. Note that the AWK script is relatively fragile. If you do more than 50 simulation reps and change the file architechture, you'll need to tweak things. 

**Disclaimer: I have not tested this on a naive computer, so I have not verified that relative paths are correct**
```
# Make a directory to store the mean fitness and hybrid index for all parameter sets
mkdir summary

## Store the number of parameter combinations as a variable
numCombos=$(< "paramGrid.A.0.6.csv" wc -l)

# Now get the means for each parameter set - THIS IS THE STEP THAT TAKES AGES!
parallel "sh ../bin/awk_mean.sh {}" ::: $(seq 1 $numCombos) 

## The above will dump ***A LOT*** of summary files into the summary/ dir - (numCombos x 2  files to be precise).

## I'm a bozo. 
## I neglected to check whether AWK is zero-indexed  when writing the script... it is.
## I needed to grab out the 0th generation from each parameter combination.
## This is because the fitnesses in the simulations are specified
## such that the maximum possible is 1.0, so a population can start with 
## a fitness <1 if positive selection is being modelled

## The following line of AWK-ish can get the 2nd line from each file, adding the file name at the beginning. I'll then read that into the parsing script and store it as a dictionary.
## Each parameter combination SHOULD have the exact same fitness in generation 0, as no selection has yet been imposed. I'll ponly store the first line  

awk 'FNR==2 {print FILENAME, $0}' outputs/A.0.6.*_1.fitStats > A.0.6_generation_0_fitnesses.txt
## This takes a few minutes


# Now take the means, and combine them into one big flipping file that can be used to explore parameter space... T
python ../bin/parseSimBitOutput_awkOutput.py --input summary/ --simbit paramGrid_A.0.6.csv --output A.0.6.output.csv --gen0 A.0.6_generation_0_fitnesses.txt

## The above will generate a file for each of the population sizes that were simulated.

## I sent those to Dea and she then worked her R magic to make beautiful plots.

```

The script that Dea wrote for making plots is at [here](bin/here) 
