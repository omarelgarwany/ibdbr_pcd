library(tidyverse)
library(patchwork)
library(ggplotify)
source('/nfs/users/nfs_o/oe2/ibdbr/scripts/gwas/extract_pheno.R')

phenos <- c('pcd')
output_dir <- '/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/output/pheno/'
about_f <- '/lustre/scratch123/hgi/projects/ibdgwas_bioresource/phenotype_data/release_20220404/raw/IBD_BioRes_phenotypes_20220404.txt'

for (pheno in phenos) {
  output_f <- paste0(output_dir,pheno,'.pheno')
  dat <- read.csv(about_f,header=T,sep='~')
  
  pheno_dat <- extract_pheno_func_list[[pheno]](dat)
  
  #IID,phenotype
  ctrl_samples <- pheno_dat[['0']] %>% distinct(SPid_1) %>% dplyr::rename(c('V2'='SPid_1')) %>% mutate(pheno=0)
  case_samples <- pheno_dat[['1']] %>% distinct(SPid_1)%>% dplyr::rename(c('V2'='SPid_1')) %>% mutate(pheno=1)
  
  #Loading fam data
  chroms <- 1:22
  fam_f_list <- setNames(paste0('/lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/post_imputation_b04b06b15b19b20/imp_QC/HRC/nodup/hg38/chr',chroms,'_hg38.fam'),
                         chroms)
  fam_dat <- map_dfr(fam_f_list,read.csv,sep=' ',.id='chrom',header=F) %>% distinct(V1,V2) 
  
  #Joining ctrl/case data with fam data
  ctrl_dat <- ctrl_samples %>% inner_join(fam_dat) %>% select(V1,V2,pheno)
  case_dat <- case_samples %>% inner_join(fam_dat) %>% select(V1,V2,pheno)
  
  case_ctrl_dat <- rbind(case_dat,ctrl_dat)
  
  write.table(case_ctrl_dat,output_f,sep='\t',col.names=F,row.names=F,quote=F)
  
}
