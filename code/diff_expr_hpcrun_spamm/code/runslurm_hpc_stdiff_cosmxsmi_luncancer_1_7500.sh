#!/bin/bash
#SBATCH --job-name=stidff_cosmx_lung
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/4472525/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code/logs/output_cosmx_lung_%A_%a.out
#SBATCH --error=/home/4472525/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code/logs/error_cosmx_lung_%A_%a.err
#SBATCH -a 1-7500%1000

cd /home/4472525/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code

args=$(cat ../data/smi_gene_meta_combo_lungcancer.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Start analysis: "
echo $(date)

module load R
Rscript /home/4472525/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code/spatial_model_tests_spamm.R $args

cd /home/4472525/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm/code

echo "End analysis: "
echo $(date)
