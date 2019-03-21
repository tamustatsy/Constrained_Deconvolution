#!/bin/bash
#SBATCH -N 1 
#SBATCH -n 10
#SBATCH --mem=30G 
#SBATCH --time=30:00:00 
#SBATCH --partition=stats_short

module load R

cd ~/Research/unimodal_density_measurement_error/Bayesian_codes/code_for_paper_update/Simulation_2/heteroscedastic/normal_error/
R CMD BATCH --no-save --no-restore run_Parallel_t_mix_HET.R run_Parallel_t_mix_HET.Rout
  
