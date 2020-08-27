## A lil python script to combine the fitness data from across several of Remi's simulations

import pandas as pd, glob, argparse


def parseSimBit( inputDir, trim ):
	print(inputDir)
## Let's start with fitness
# Make an array to dump the processed dataframes into
	allFitness = []
# Loop over the Simbit output files for fitness
	for f in glob.glob(inputDir + '_*.fitStats'):
#		print(f)
# Identify the simulation replicate number

		rep = f.split('.')[-2].split('_')[1]
# Read in the file
		fitness = pd.read_csv(f, sep = '\t')

# Extract the fitness for the first generation. It'll be used to calc. rel. fitness
		initialFitness = fitness[fitness['Generation'] == 0].P0_meanFit

# To recreate Mike and Sally's plots, I trim off all generations beyond the first 100
		fitness = fitness[fitness['Generation'] <= trim]
# Calculate rel. Fitness and store in a new column
		fitness['relFitness'] = fitness['P0_meanFit']/float(initialFitness)
# Store the sim. replicate in a new column
		fitness['rep'] = int(rep)
# Add the processed file to the list
		allFitness.append(fitness)

## Now let's repeat the above process, but for hybrid index
## I don't comment the loop except where it is substantially different to the above

	if len(allFitness) == 0:
		return None
	fitnessDF = pd.concat(allFitness)

	meanDF_fitness = fitnessDF.groupby('Generation').mean().drop('rep', axis = 1)

	allFitness = []

	allHybrid = []

	for h in glob.glob(inputDir + '_*ntrlLociOnly.T1AverageHI'):
# Identify the simulation replicate number

		rep = h.split('_')[-1].split('n')[0]
# Read in the file
#		print(h)
		hybridIndex = pd.read_csv(h, sep = '\t')
		
		hybridIndex['rep'] = int(rep)

# I trim off all generations beyond the first 100
		hybridIndex = hybridIndex[hybridIndex['Generation'] <= trim]

# Add the processed file to the list
		allHybrid.append(hybridIndex)

	if len(allHybrid) == 0:
		return None
		
	hybridDF = pd.concat(allHybrid)

	meanDF_hybrid = hybridDF.groupby('Generation').mean().drop('rep', axis = 1)

	meanDF_fitness['hybridIndex'] = meanDF_hybrid['P0']

	return meanDF_fitness


##
## I ran the script in parallel at the command line using:
##		 parallel "python ../bin/parseSimBitOutput.py -i A.0.0.{}/ -o A.0.0.{}.csv" ::: $(seq 1 8)
##
def main():

	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str,
			help = "The directory containing the SimBit output directories")
	parser.add_argument("--simbit","-s", 
			required = True,
			dest = "simbit",
			type =str,
			help = "The file describing the SimBit simulation parameters")
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


	simbit = pd.read_csv( args.simbit )
	if 'simulationBigPath' in list(simbit):
		simbit = simbit.drop(['simulationBigPath'], axis=1)

	simbit = simbit.set_index('simulationSetID')
	sim_dict = simbit.to_dict('index')
	
	population_sizes = set(list(simbit.N))

	count = 0

	all_sims = []

	for s in sim_dict.keys():
#	for d in glob.glob(args.input + '/A.*fitStats'):
		name = s.split('/')[-1]
		print(name)
		d_df = parseSimBit( args.input+s , args.trim)
	
		d_df['setID'] = s

		if d_df is not None:
			pass
		else:
			print("no files found for "+s)
			continue

		count += 1
		d_df = d_df.assign(**sim_dict[name]) 
		
		all_sims.append(d_df)
#		if count == 100: break

# concatenate the big list of simulations
	fullDF = pd.concat(all_sims)

# Write the big dataframe to one file per population size
	for N in  population_sizes:
		
		N_DF = fullDF[fullDF.N == N]
		N_DF.to_csv(str(N) + '.mean.'+args.output, index = True)






main()

