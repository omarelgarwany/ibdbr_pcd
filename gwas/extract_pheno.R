extract_pcd_pheno_dat <- function(dat) {
  #Setting up data of patients who responded to phenotype question
  pheno_dat_list <- list()
  pheno_dat <- dat %>% drop_na(perianal) 
  
  #Identifying 1 and 0
  pheno_yes_dat <- pheno_dat %>% filter(perianal==1) %>% filter(peri_type___2==1 | peri_type___3==1 | peri_type___4==1 | peri_type___5==1)
  pheno_no_dat <- pheno_dat %>% filter(perianal==2) %>% filter(peri_type___2==0 & peri_type___3==0 & peri_type___4==0 & peri_type___5==0)
  
  pheno_dat_list <- list(
    '1'=pheno_yes_dat,
    '0'=pheno_no_dat
  )
  return(pheno_dat_list)
}

extract_surgery_pheno_dat <- function(dat) {
  #Setting up data of patients who responded to phenotype question
  pheno_dat_list <- list()
  pheno_dat <- dat %>% drop_na(perianal) 
  
  #Identifying 1 and 0
  pheno_yes_dat <- pheno_dat %>% filter(surgery==1) 
  pheno_no_dat <- pheno_dat %>% filter(surgery==2)
  
  pheno_dat_list <- list(
    '1'=pheno_yes_dat,
    '0'=pheno_no_dat
  )
  return(pheno_dat_list)
}


extract_pheno_func_list <- list(
  'pcd'=extract_pcd_pheno_dat,
  'surgery'= extract_surgery_pheno_dat
)