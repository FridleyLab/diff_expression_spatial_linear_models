##
# Prepare list of genes in Hallmark pathway for STdiff tests
# NOTE: Requires first running "prepare_data_hpcrun_visium_raviglioblastoma.R"

library('tidyverse')
library('msigdbr')
library('spatialGE')

# Read STlist with data and annotations
fp = 'diff_expr_hpcrun_spamm/data/visium_stlist_w_clusters_raviglioblastoma.RDS'
visium = readRDS(fp)
rm(fp) # Clean env

# Get genes in Hallmark pathways
gene_sets = msigdbr(species="Homo sapiens", category='H')
gene_sets = filter(gene_sets, gs_subcat == "")
genes_hmk = unique(gene_sets[['human_gene_symbol']])

# Read randomly sampled genes to avoid running same gene twice
fp = 'diff_expr_hpcrun_spamm/data/visium_gene_meta_combo_raviglioblastoma.txt'
random_genes = read.table(fp, header=F, sep=' ')
random_genes = random_genes %>% select(c('V1', 'V4')) %>% distinct()

# Extract genes align a gradient of average expression
# Make combinations of gene x cluster
tmp_df = list()
for(i in c('UKF243', 'UKF275')){
  genes = genes_hmk[genes_hmk %in% rownames(visium@tr_counts[[i]])]
  genes = genes[!(genes %in% random_genes[['V1']][random_genes$V2 == i]) ]
  meta = unique(visium@spatial_meta[[i]][['annots']])
  tmp_df[[i]] = expand_grid(genes, meta, '../data/visium_stlist_w_clusters_raviglioblastoma.RDS', i)
  rm(meta, genes) # Clean env
}
tmp_df = do.call(bind_rows, tmp_df)
write.table(tmp_df, col.names=F, row.names=F, quote=F,
            file='diff_expr_hpcrun_spamm/data/visium_gene_meta_combo_raviglioblastoma_hallmark.txt')

rm(tmp_df) # Clean env

