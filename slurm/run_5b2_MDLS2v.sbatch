#!/bin/bash
#SBATCH --account=def-ilafores
#SBATCH --job-name=SBC
#SBATCH --nodes=1
#SBATCH --partition=cpularge
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=240G
#SBATCH --time=24:00:00
#SBATCH --output=/net/nfs-bio/jbod2/def-ilafores/analysis/orchardPhyllosphere2/slurm/logs/%x-%j.out

set -euo pipefail

export REPO=/net/nfs-bio/jbod2/def-ilafores/analysis/orchardPhyllosphere2
export SCRIPT=$REPO/src/hiermod/16S/5b_MDLS2v/5b.2_MDLS2v_16S_calibration.R
cd "$REPO"
git pull

# >>> FILL IN: however you load R + cmdstanr's toolchain on this cluster <<<
# e.g. module load StdEnv/2023 gcc r/4.4   -- or `source activate <env>`
module load StdEnv/2023 r/4.4.0 gcc/12.3

# setwd() (not just the shell's cd) fixes the entire nested source() chain
# regardless of how this job's initial cwd behaves.
Rscript -e 'setwd(Sys.getenv("REPO")); source(Sys.getenv("SCRIPT"), print.eval = TRUE)'

# USAGE 
# sbatch $PWD/slurm/run_5b2_MDLS2v.sbatch