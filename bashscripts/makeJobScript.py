#import numpy as np
import os

start = 1547
end=3582

experiment='TRACK2'
Ascript = '/sps/nemo/scratch/Edinburgh_student/Katherine/analysisPipeline.cpp'

#rootArrays = np.arange(start,end+1)
for i in range(start,end):
	dataFile = f'/sps/nemo/scratch/Edinburgh_student/Katherine/data/reconstructed_runs_old/run-{i}.root'
	#print(dataFile)
	# if the data file doesn't exist, then move on
	if (not os.path.isfile(dataFile)):
		print(f"{i}.root non-existent")
		continue
	folder = f'/sps/nemo/scratch/Edinburgh_student/Katherine/{experiment}/run-{i}'
	os.mkdir(folder)
	script=f"""#!/bin/bash
#SBATCH --job-name={experiment}_{i}          # Job name
#SBATCH --output={folder}/output_%x_%j.log   # Standard output and error log

#SBATCH --cpus-per-task 1   # nb of cpus per tasks
#SBATCH --time 0-06:00      # time (D-HH:MM)
#SBATCH --partition htc
#SBATCH --ntasks 1
#SBATCH --mem 5500

set -euo pipefail
echo {dataFile}
#cd into folder
cd {folder}
echo $(pwd)

#need to load root module
module load root
root -l -b -q '{Ascript}("{dataFile}")'
echo "Finished Successfully"
	"""
	fd = open(f'{folder}/runScript.sh','w')
	fd.write(script)
	fd.close()
	#print(f"Wrote to {folder}")

