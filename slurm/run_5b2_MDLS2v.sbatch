#!/bin/bash
#SBATCH --account=def-ilafores
#SBATCH --job-name=SBC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=cpularge
#SBATCH --cpus-per-task=48
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=/net/nfs-bio/jbod2/def-ilafores/analysis/orchardPhyllosphere2/slurm/logs/%x-%j.out
#SBATCH --export=NONE

set -euo pipefail

export REPO=/net/nfs-bio/jbod2/def-ilafores/analysis/orchardPhyllosphere2
export SCRIPT=$REPO/src/hiermod/16S/5b_MDLS2v/5b.2_MDLS2v_16S_calibration.R
cd "$REPO"

module load StdEnv/2023 r/4.4.0 gcc/12.3

Rscript -e 'setwd(Sys.getenv("REPO")); source(Sys.getenv("SCRIPT"), print.eval = TRUE)'

# USAGE 
# sbatch $PWD/slurm/run_5b2_MDLS2v.sbatch