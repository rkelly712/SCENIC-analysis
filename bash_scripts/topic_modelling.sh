
#!/bin/bash

#SBATCH --job-name Topic_modelling
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-type TIME_LIMIT_50,TIME_LIMIT_90,TIME_LIMIT
#SBATCH --cpus-per-task=15
#SBATCH --mem=250g
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:150

## Load the mamba/conda environment in the HPC space.
source /data/$USER/conda/etc/profile.d/conda.sh && source /data/$USER/conda/etc/profile.d/mamba.sh;

source myconda;

export MAMBA_NO_BANNER=1;

mamba activate scenicplus;

# Load python and run the script
python /data/kellyrc/SCENIC+/scripts/python_scripts/runModels_lda_cgs.py \
        -i /data/kellyrc/SCENIC+/TopicModelling/WT_cistopic_obj.pkl \
        -o /data/kellyrc/SCENIC+/TopicModelling/models/reik_WT_models_500iter_noDBL.pkl \
        -nt 2,5,10,15,20,25,30,35,40,45,50 \
        -c 11 \
        -it 500 \
        -a 50 \
        -abt True \
        -e 0.1 \
        -ebt False \
        -sp /lscratch/$SLURM_JOB_ID/ \
        -s 555 \
        -td /lscratch/$SLURM_JOB_ID/
        
