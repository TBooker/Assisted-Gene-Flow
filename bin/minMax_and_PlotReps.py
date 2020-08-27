
import pandas as pd
import sys, glob

def main():
	
	migrationScenarios = {"0.001,0.001,0.001,0.001,0.001,1,2,3,4,5":[0.005,"Every 1"],
	"0.001,0.001,0.001,0.001,0.001,1,3,5,7,9":[0.005,"Every 2"],
	"0.001,0.001,0.001,0.001,0.001,1,5,9,13,17":[0.005,"Every 4"],
	"0.005,1":[0.005,"Single"],
	"0.01,0.01,0.01,0.01,0.01,1,2,3,4,5":[0.05,"Every 1"],
	"0.01,0.01,0.01,0.01,0.01,1,3,5,7,9":[0.05,"Every 2"],
	"0.01,0.01,0.01,0.01,0.01,1,5,9,13,17":[0.05,"Every 4"],
	"0.05,1":[0.05,"Single"],
	"0.1,0.1,0.1,0.1,0.1,1,2,3,4,5":[0.5,"Every 1"],
	"0.1,0.1,0.1,0.1,0.1,1,3,5,7,9":[0.5,"Every 2"],
	"0.1,0.1,0.1,0.1,0.1,1,5,9,13,17":[0.5,"Every 4"],
	"0.5,1":[0.5,"Single"]}
	
	simSets = pd.read_csv("~/work/GeneticRescue/simulations/temp/A.0.6/paramGrid_A.0.6.csv")
	focalSet = simSets[simSets["dominanceCoef_MAC"] == 0.5]
	focalSet = focalSet[focalSet["dominanceCoef_PAC"] == 0.5]
	focalSet = focalSet[focalSet["nbLoci_MAC"] == 5]
	focalSet = focalSet[focalSet["nbLoci_PAC"] == 5]
	focalSet = focalSet[focalSet["N"] == 10000]
	focalSet = focalSet[focalSet["singleSToRescale_PAC"] == 0.1]
	focalSet = focalSet[focalSet["singleSToRescale_MAC"] == 0.1]

	repOutput = []
	for i in focalSet.simulationSetID:
		print(i)
		reps = []
		for f in glob.glob("outputs/"+i+"/*.fitStats"):
			temp =  pd.read_csv(f, sep = "\t") 
			temp["rep"] = f.split(".")[-2]
			temp["simulationSetID"] = i
			initialFitness = list(temp[temp["Generation"]==0].P0_meanFit)[0] 
			temp["simulationSetID"] = i
			temp["relFitness"] = temp["P0_meanFit"]/initialFitness

			reps.append(temp)
		repOut = pd.concat(reps)
		meta = focalSet[focalSet["simulationSetID"] == i].copy()
		migrationDetails = migrationScenarios[meta.migration_scenario.values[0]]
		meta['overallRate'] = migrationDetails[0]
		meta['pulse'] = migrationDetails[1]
		repOutput.append( pd.merge(repOut, meta, on = "simulationSetID") ) 
	output = pd.concat(repOutput).drop(["generalOutputPath","generalCommandsPath"], axis = 1)
	output.to_csv("N10000_IndividualReps.csv",index = False)


main()
