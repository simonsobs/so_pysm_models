#!/bin/bash
#SBATCH --partition=soap
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48 
#SBATCH --job-name=rotate_maps

time python rotate_maps.py ~/s/simonsobs/so_pysm_models_data/
