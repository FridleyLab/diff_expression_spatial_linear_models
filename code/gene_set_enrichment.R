##
# Calculation of GSEA scores using STdiff p-values as ranking statistic

library('tidyverse')
library('fgsea')
library('msigdbr')

# Read list of genes selected for analysis
fp = 'diff_expr_hpcrun_spamm/data/visium_gene_meta_combo_raviglioblastoma_hallmark.txt'
hm_genes = read_delim(fp, show_col_types=F, col_names=F)

rm(fp) # Clean env

fp = 'diff_expr_hpcrun_spamm/code/spatial_test_all_compiled.csv'
df_res = read_delim(fp, delim=',', show_col_types=F) %>%
  filter(st_tech == 'visium') %>%
  filter(sample %in% c("UKF243", "UKF275")) %>%
  filter(gene %in% unique(hm_genes[['X1']]))

rm(fp) # Clean env

col_pal = list(
  visium=c(UKF243="#F4A736", UKF275="#D1BBD7")
)

# Get list of Hallmark pathways
gene_sets = msigdbr(species='Homo sapiens', category='H') %>%
  split(x=.$gene_symbol, f=.$gs_name)

# Get vector of genes in Hallmark pathways
genes_hmk = unique(unlist(gene_sets))

fgsea_res = list()
for(i in unique(df_res[['sample']])){
  df_tmp = df_res %>%
    filter(sample == i)

  fgsea_res[[i]] = list()
  for(cl in unique(df_tmp[['cluster']])){
    df_tmp2 = df_tmp %>%
      filter(cluster == cl) %>%
      add_column(adj_pval=p.adjust(.[['p_val']], method='BH')) %>%
      add_column(adj_pval_exp=p.adjust(.[['exp_p_val']], method='BH'))

    if(any(is.na(df_tmp2$p_val)) | any(is.na(df_tmp2$exp_p_val))){
      warning('There are NAs in p-values')
    }

    # Prepare named vector of p-values
    nonsp_pval = df_tmp2  %>% select(adj_pval) %>% unlist()
    names(nonsp_pval) = df_tmp2 %>% select(gene) %>% unlist()

    sp_pval = df_tmp2 %>% select(adj_pval_exp) %>% unlist()
    names(sp_pval) = df_tmp2 %>% select(gene) %>% unlist()

    cat(crayon::green(paste0('Non-spatial p-values, Sample: ', i, '; Cluster: ', cl, '\n')))
    nonsp_fgsea = fgseaMultilevel(pathways=gene_sets,
                                  stats=nonsp_pval,
                                  scoreType='pos',
                                  minSize=5,
                                  #maxSize=500, eps=0, sampleSize=50,
                                  nPermSimple=10000)

    cat(crayon::green(paste0('Spatial p-values, sample: ', i, '; Cluster: ', cl, '\n')))
    sp_fgsea = fgseaMultilevel(pathways=gene_sets,
                               stats=sp_pval,
                               scoreType='pos',
                               minSize=5,
                               #maxSize=500, eps=0, sampleSize=50,
                               nPermSimple=10000)

    # Rename columns to differentiate them during merging
    colnames(nonsp_fgsea) = paste0('nonsp_', colnames(nonsp_fgsea))
    colnames(sp_fgsea) = paste0('sp_', colnames(sp_fgsea))

    fgsea_res[[i]][[cl]] = nonsp_fgsea %>% full_join(., sp_fgsea, by=c('nonsp_pathway'='sp_pathway'))

    rm(df_tmp2, nonsp_fgsea, sp_fgsea, nonsp_pval, sp_pval) # Clean env
  }
  rm(df_tmp) # Clean env
}

# Save GSEA results
saveRDS(fgsea_res, '../data/fgsea_results.RDS')

