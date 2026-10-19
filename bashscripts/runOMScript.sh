#!/bin/bash

#SBATCH --job-name=htc_101          # Job name
#SBATCH --output=output_%x_%j.log   # Standard output and error log

#SBATCH --cpus-per-task 1   # nb of cpus per tasks
#SBATCH --time 0-06:00      # time (D-HH:MM)
#SBATCH --partition htc
#SBATCH --ntasks 1
#SBATCH --mem 5500

set -euo pipefail

JOBHOST=$(hostname)
WORKDIR=$(pwd)

# your command here
# root analysisScript.cpp datafile.root
module load root
root -l -b -q 'OMLookup.cpp'

echo ""
echo "Job submission ended OK \o/"

