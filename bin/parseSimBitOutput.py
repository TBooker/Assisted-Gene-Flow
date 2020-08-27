## A lil python script to combine the fitness data from across several of Remi's simulations

import pandas as pd, glob, argparse

def main():
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str,
			help = "The directory containing the input files")
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str,
			help = "The name of the output file. Provide a suffix.")
	parser.add_argument("--trim", 
			required = False,
			dest = "trim",
			type = int,
			help = "Specify a maximum generation. [100]",
			default = 100)
	args = parser.parse_args()

## Let's start with fitness
# Make an array to dump the processed dataframes into
	allFitness = []
	
# Loop over the Simbit output files for fitness
	for f in glob.glob(args.input + '/outputs/*.fitStats'):
		print(f)
		
# Identify the simulation replicate number
#		rep = f.split('.')[-2]
		rep = f.split('/')[-1].split('_')[-1].split('.')[0]
		

# Read in the file
		fitness = pd.read_csv(f, sep = '\t')
# Extract the fitness for the first generation. It'll be used to calc. rel. fitness
		initialFitness = fitness[fitness['Generation'] == 0].P0_meanFit

# To recreate Mike and Sally's plots, I trim off all generations beyond the first 100
		fitness = fitness[fitness['Generation'] <= args.trim]
# Calculate rel. Fitness and store in a new column
		fitness['relFitness'] = fitness['P0_meanFit']/float(initialFitness)
# Store the sim. replicate in a new column
		fitness['rep'] = int(rep)
# Add the processed file to the list
		allFitness.append(fitness)

	bigFitnessDF = pd.concat(allFitness)
	bigFitnessDF.to_csv('fitness.'+args.output, index = False)



## Now let's repeat the above process, but for hybrid index
## I don't comment the loop except where it is substantially different to the above
	allHybrid = []

	for h in glob.glob(args.input + '/outputs/*.T1HI'):
		rep = h.split('.')[-2]

		hybrid = pd.read_csv(h, sep = '\t')
		hybrid = hybrid[hybrid['Generation'] <= args.trim]
		hybrid['rep'] = int(rep)

# Get the mean hybrid index for all individuals in the population (for each row, take an average across columns)
		hybrid['meanHI'] = hybrid[list(hybrid)[1:]].mean(axis=1)
# Extract only the mean, generation and replicate data
		hybridSlice = hybrid[['Generation','meanHI','rep']]
		allHybrid.append(hybridSlice)

# Concatenate the lists of dataframes into full dataframes
	bigFitnessDF = pd.concat(allFitness)
	bigHybridDF = pd.concat(allHybrid)

# Merge the two big dataframes on Generation and replicate
	fullDF = bigFitnessDF.merge(bigHybridDF, left_on=['Generation','rep'], right_on = ['Generation','rep'], how='left')

# Write the big dataframe to a file
	fullDF.to_csv('all.'+args.output, index = False)

# Take the mean for each statistic across Generations
	meanDF = fullDF.groupby('Generation').mean().drop('rep', axis = 1)
# Write the means to a file
	meanDF.to_csv('mean.'+args.output, index = True)

main()



##
## I ran the script in parallel at the command line using:
##		 parallel "python ../bin/parseSimBitOutput.py -i A.0.0.{}/ -o A.0.0.{}.csv" ::: $(seq 1 8)
##
