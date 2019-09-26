#!/bin/bash
#SBATCH --partition=soap
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48 
#SBATCH --job-name=ud_grade_maps

time python ud_grade_maps.py ~/s/simonsobs/so_pysm_models_data/equatorial/websky/0.3/
