#!/bin/bash
#SBATCH --job-name=stidff_cosmx_lung
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm_pairwise_tests/code/logs/output_cosmx_lung_%A_%a.out
#SBATCH --error=/home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm_pairwise_tests/code/logs/error_cosmx_lung_%A_%a.err
#SBATCH -a 1-1200%600

cd /home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm_pairwise_tests/code

args=$(cat ../data/smi_gene_meta_combo_lungcancer_pairwise_tests.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)

echo "Start analysis: "
echo $(date)

module load R
Rscript /home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm_pairwise_tests/code/spatial_model_tests_spamm_pairwise_tests.R $args

cd /home/hpc_user/spatialge_stdiff_manuscript/stdiff_hpcrun_spamm_pairwise_tests/code

echo "End analysis: "
echo $(date)
