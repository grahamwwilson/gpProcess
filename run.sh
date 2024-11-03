#!/bin/sh
#
# Submit a SLURM job to do post-processing of specified GP file
#
# Example ./run.sh Z-126
#
RUN=$1

echo '$RUN = '${RUN}

sbatch slurm_lumiprep.sh ${RUN}

exit
