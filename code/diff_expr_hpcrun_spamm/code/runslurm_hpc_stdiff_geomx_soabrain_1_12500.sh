#!/bin/bash
#SBATCH --job-name=stidff_geomx_brain
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code/logs/output_geomx_brain_%A_%a.out
#SBATCH --error=/home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code/logs/error_geomx_brain_%A_%a.err
#SBATCH -a 1-12500%1000

cd /home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code

args=$(cat ../data/geomx_gene_meta_combo_soabrain.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Start analysis: "
echo $(date)

module load R
Rscript /home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code/spatial_model_tests_spamm.R $args

cd /home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code

echo "End analysis: "
echo $(date)
