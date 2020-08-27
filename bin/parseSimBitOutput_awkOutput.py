## A lil python script to combine the fitness data from across several of Remi's simulations

import pandas as pd, glob, argparse

def get_gen0_fitness(gen0_file):
	fit_dict = {}
	with open(gen0_file, 'r') as ff:
		for l in ff:
			lineList = l.strip().split()
			sim_ID = lineList[0].split('/')[-1].split('_')[0]
			fitness = float( lineList[2] )
			
			fit_dict[ sim_ID ] = fitness
	return fit_dict

def main():

	parser = argparse.ArgumentParser(description="")
	parser.add_argument("--input","-i", 
			required = True,
			dest = "input",
			type =str,
			help = "The directory containing the Awk output files")
	parser.add_argument("--simbit","-s", 
			required = True,
			dest = "simbit",
			type =str,
			help = "The file describing the SimBit simulation parameters")
	parser.add_argument("--output","-o", 
			required = True,
			dest = "output",
			type =str,
			help = "The name of the output file. Provide a suffix - you'll get a file for each population size.")
	parser.add_argument("--gen0", "-g", 
			required = True,
			dest = "gen0",
			type =str,
			help = "The name of the file with the generation 0 fitnesses")
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

	gen0_fitnesses = get_gen0_fitness( args.gen0 )

	
	population_sizes = set(list(simbit.N))

	count = 0

	all_sims = []

	for s in sim_dict.keys():
## The way the files have been simulated, we have been working on batches of parameter combinations. The try/except conditions below make the script flexible to missing parameter combos.
		
		try:
			s_fitness = gen0_fitnesses[s] ## This is the fitness in generation 0
		except KeyError:
			continue
		mean_file = args.input + "/" + s + "_mean.fitStats"
		introgression_file = args.input + "/" + s + "_mean.AverageHI"

		try:
			means = pd.read_csv( mean_file , sep = ' ', header = None, names = ['Generation', 'Fitness'])
			introgression = pd.read_csv( introgression_file , sep = ' ', header = None, names = ['Generation', 'hybridIndex'])
			
		except FileNotFoundError:
			continue

		means['relFitness'] = means['Fitness']/ s_fitness
		means['hybridIndex'] = introgression['hybridIndex']

		means['setID'] = s

		if means is not None:
			pass
		else:
			print("no files found for "+s)
			continue

		count += 1
		means = means.assign(**sim_dict[s]) 
		
		all_sims.append(means)
		if count % 100 == 0:
			print( count )

# concatenate the big list of simulations
	fullDF = pd.concat(all_sims)

# Write the big dataframe to one file per population size
	for N in  population_sizes:
		N_DF = fullDF[fullDF.N == N]
		N_DF.to_csv(str(int(N)) + '.mean.'+args.output, index = True)







main()

