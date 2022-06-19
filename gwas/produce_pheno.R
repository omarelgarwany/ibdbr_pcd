library(tidyverse)
library(patchwork)
library(ggplotify)

about_f <- '/lustre/scratch123/hgi/projects/ibdgwas_bioresource/phenotype_data/release_20220404/raw/IBD_BioRes_phenotypes_20220404.txt'
output_f <- '/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/output/pheno/pcd.pheno'
dat <- read.csv(about_f,header=T,sep='~')

pheno_dat <- list()
####################################################################
###########This part is customizable for different phenotypes#######
####################################################################
#Setting up data of patients who responded to phenotype question
pheno_dat <- dat %>% drop_na(perianal) 

#Identifying 1 and 0
pheno_yes_dat <- pheno_dat %>% filter(perianal==1) %>% filter(peri_type___2==1 | peri_type___3==1 | peri_type___4==1 | peri_type___5==1)
pheno_no_dat <- pheno_dat %>% filter(perianal==2) %>% filter(peri_type___2==0 & peri_type___3==0 & peri_type___4==0 & peri_type___5==0)

####################################################################
####################################################################
####################################################################

#IID,phenotype
ctrl_samples <- pheno_no_dat %>% distinct(SPid_1) %>% dplyr::rename(c('V2'='SPid_1')) %>% mutate(pheno=0)
case_samples <- pheno_yes_dat %>% distinct(SPid_1)%>% dplyr::rename(c('V2'='SPid_1')) %>% mutate(pheno=1)

#Loading fam data
chroms <- 1:22
fam_f_list <- setNames(paste0('/lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/post_imputation_b04b06b15b19b20/imp_QC/HRC/nodup/hg38/chr',chroms,'_hg38.fam'),
                       chroms)
fam_dat <- map_dfr(fam_f_list,read.csv,sep=' ',.id='chrom',header=F) %>% distinct(V1,V2) 

#Joining ctrl/case data with fam data
ctrl_dat <- ctrl_samples %>% inner_join(fam_dat) %>% select(V1,V2,pheno)
case_dat <- ctrl_samples %>% inner_join(fam_dat) %>% select(V1,V2,pheno)

case_ctrl_dat <- rbind(case_dat,ctrl_dat)

write.table(case_ctrl_dat,output_f,sep='\t',col.names=F,row.names=F,quote=F)
