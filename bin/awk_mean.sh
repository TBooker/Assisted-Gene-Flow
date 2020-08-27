

## Pseudo-code
# cat all files to pass them to AWK
# for each file
#	For i in (1 to 100)
#		make a vector, 'total' that is 100 elements long (all are set to 0)
#		Grab the entry for the ith row of the second column and make that total[i]
#	When that's done, divide it by 50, the total number of   
cat outputs/A.0.6.$1_*.fitStats | awk '
BEGIN {
	for (i =1; i <=100; i++)
		total[i] =0
	}
{ total[$1]+=$2 }

END {
	for (i = 1; i <=100; i++)
		print i, total[i]/50
		}
' > summary/A.0.6.$1_mean.fitStats 


cat outputs/A.0.6.$1_*ntrlLociOnly.T1AverageHI| awk '
BEGIN {
	for (i =1; i <=100; i++)
		total[i] =0
	}
{ total[$1]+=$2 }

END {
	for (i = 1; i <=100; i++)
		print i, total[i]/50
		}
' > summary/A.0.6.$1_mean.AverageHI 

