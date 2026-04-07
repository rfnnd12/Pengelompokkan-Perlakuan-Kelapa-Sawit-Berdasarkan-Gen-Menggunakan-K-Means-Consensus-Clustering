#!/bin/bash
#running mapping
#SBATCH -p short      # Specifying queue
#SBATCH -n 16         # Specifying the number of CPUs (1-56)
#SBATCH --mem 27G     # Specifying memory usage (5-1500G)
#SBATCH -t 24:00:00   # Specifying elapsed time limits (1-24hours)

############################################################################
# SETTING VARS
############################################################################

# Now let's keep track of some information just in case anything goes wrong
start=$(date +%s%N)
echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "=========================================================="

module load bioinformatics/anaconda
conda init
conda activate skripsi

Rscript Consensus_clustering.R


echo "done running"
echo "=========================================================="


# print out some diagnostic stuff
end=$(date +%s%N)
duration=$(((end - start) / 1000000000 /60)) #in mins

echo "Stop time is $(date)"
echo "Duration: ${duration}"

