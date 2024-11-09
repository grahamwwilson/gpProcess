#!/bin/bash
#SBATCH --job-name=gpProcess          # Job name
#SBATCH --partition=bigjay            # Partition Name (Required)
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=gwwilson@ku.edu   # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --time=0-2:00:00              # Time limit days-hrs:min:sec
#SBATCH --output=gpProcess_%j.log     # Standard output and error log

# Note - originally I specified 4 GB, but recent runs have been close to 
# that limit and have taken more time, so have recently upped this to 8 GB.

date
echo 'Running lumiprep job as '
echo 'User '$USER

pwd
hostname
date

RUN=$1
echo '$RUN '${RUN}

./lumiprep.sh ${RUN}

date
exit
