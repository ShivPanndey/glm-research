#!/bin/bash
#SBATCH --job-name=tnf_featurize
#SBATCH --partition=batch
#SBATCH --array=0-49              # 50 tasks; each handles every 50th genome
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/tnf_%A_%a.out
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sp74395@uga.edu

# Adjust to your environment (check `module spider Python` on Sapelo2, or use a conda env).
module load Python
source ~/envs/tnf/bin/activate      # pip install -r requirements.txt

DATA=/scratch/$USER/glm/family_1000      # extracted dataset
python scripts/run_baseline.py featurize \
    --genomes $DATA/genomes \
    --out /scratch/$USER/glm/tnf_cache_w5k \
    --window 5000 --stride 2500
