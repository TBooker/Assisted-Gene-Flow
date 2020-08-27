## Move all the information for each runs into individual directories

# Get the number of parameter sets - assuming that they are named sequentially
numFiles=$(cat $1 |wc -l )
for i in $(seq 1 $numFiles)
	do
	mkdir outputs/A.0.6.$i/
	mv outputs/A.0.6.${i}_*  outputs/A.0.6.$i/
	done
