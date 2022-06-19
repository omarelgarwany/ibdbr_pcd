from pathlib import Path
import glob
import subprocess as sp
import pandas as pd
import numpy as np




PATHS_FILE=config['PATHS_FILE']
CONFIGS_FILE=config['CONFIGS_FILE']
BDIR=config['BDIR']
IDPSNPS_DIR=config['IDPSNPS_DIR']

GRM_GWAS_DIR=config['GRM_GWAS_DIR']

phenos=['pcd','surgery']

chr_nums=22
chrs=[str(i) for i in range(1,chr_nums+1)]
# localrules: create_bed_paths_file,collect_idp_snps,make_sparse_grm, run_gwas

#RULES
wildcard_constraints:
        chr=r"|".join(set(chrs))

rule all:
        input:
                expand("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/output/gwas/{pheno}.fastGWA",pheno=phenos)

#/STARTSTEP/
#1 use independent SNPs to create GRM...
rule get_idp_snps:
        input:
                BDIR+"chr{chr}_hg38.bed"
        output:
                IDPSNPS_DIR+"chr{chr}.prune.in",IDPSNPS_DIR+"chr{chr}.prune.out"
        params:
                PATHS_FILE=PATHS_FILE,
                CONFIGS_FILE=CONFIGS_FILE
        resources:
                mem_mb=4000
        shell:
                'source {params.PATHS_FILE}; source {params.CONFIGS_FILE}; "${{PLINK_BIN}}" --bfile "${{BDIR}}"chr{wildcards.chr}_hg38 --exclude "${{exclude_regions_f}}" --maf ${{IDPSNPS_MAF}} --indep-pairwise "${{IDPSNPS_WINDOW}}" "${{IDPSNPS_ITERATION}}" "${{IDPSNPS_R2THRESH}}" --out "${{IDPSNPS_DIR}}"chr{wildcards.chr}'
#... and collecting all independent SNPs
rule collect_idp_snps:
    input:
        expand(IDPSNPS_DIR+"chr{chr}.prune.in",chr=chrs)
    output:
        IDPSNPS_DIR+"gwas.idpSNP"
    shell:
        "cat {input} > {output}"
#/ENDSTEP/

#2 Compute GRMs
#/STARTSTEP/
#First we collect all bed paths into a list...
rule create_bed_paths_file:
    output:
        GRM_GWAS_DIR+"geno_chrs.txt"
    params:
        PATHS_FILE=PATHS_FILE,
        CONFIGS_FILE=CONFIGS_FILE,
        chr_nums=chr_nums

    shell:
        'source {params.PATHS_FILE}; source {params.CONFIGS_FILE}; for i in $(seq 1 {params.chr_nums}); do echo "${{BDIR}}"chr${{i}}_hg38; done > {output}'

#...Second
rule make_dense_grm:
    input:
        mbfile=GRM_GWAS_DIR+"geno_chrs.txt",idpsnp_list=IDPSNPS_DIR+"gwas.idpSNP"
    output:
        GRM_GWAS_DIR+'geno_grm.grm.bin'
    params:
        PATHS_FILE=PATHS_FILE,
        CONFIGS_FILE=CONFIGS_FILE,
        chr_nums=chr_nums,
        grm_output_prefix=GRM_GWAS_DIR+'geno_grm'
    resources:
        threads=10,
        mem_mb=100000
    shell:
        'source {params.PATHS_FILE}; source {params.CONFIGS_FILE}; "${{GCTA_BIN}}" --mbfile {input.mbfile} --extract {input.idpsnp_list} --make-grm --threads {resources.threads} --out {grm_output_prefix}'

#...Third...
rule make_sparse_grm:
    input:
        GRM_GWAS_DIR+'geno_grm.grm.bin'
    output:
        GRM_GWAS_DIR+'sparse_grm.grm.sp'
    params:
        PATHS_FILE=PATHS_FILE,
        CONFIGS_FILE=CONFIGS_FILE,
        dense_grm_prefix=GRM_GWAS_DIR+'geno_grm',
        sparse_grm_prefix=GRM_GWAS_DIR+'sparse_grm'
    resources:
        mem_mb=4000
    shell:
        'source {params.PATHS_FILE}; source {params.CONFIGS_FILE}; "${{GCTA_BIN}}" --grm {params.dense_grm_prefix} --make-bK-sparse 0.05 --out {params.sparse_grm_prefix}'
#3 calculate PCs

#4 Running GWAS...
rule run_gwas_pheno:
    input:
        mbfile=GRM_GWAS_DIR+"geno_chrs.txt",
        pheno_file='/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/output/pheno/{pheno}.pheno'
    output:
        "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/output/gwas/{pheno}.fastGWA"
    params:
        PATHS_FILE=PATHS_FILE,
        CONFIGS_FILE=CONFIGS_FILE,
        sparse_grm_prefix=GRM_GWAS_DIR+'sparse_grm',
        gwas_output_prefix="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/output/gwas/{pheno}",
        wgs_samples="/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/analysis_22_03_2022/pheno/WGS.sample",
        cd_samples="/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/analysis_22_03_2022/pheno/CD.sample",
        pca_f="/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/analysis_22_03_2022/GRM_sparse/GWAS/geno_pc.eigenvec"
    resources:
        mem_mb=10000,
        threads=10
    shell:
        'source {params.PATHS_FILE}; source {params.CONFIGS_FILE}; /software/team152/oe2/bin/gcta64 --mbfile {input.mbfile} --out {params.gwas_output_prefix} --fastGWA-mlm-binary --grm-sparse {params.sparse_grm_prefix} --joint-covar --pheno {input.pheno_file} --remove {params.wgs_samples} --threads {resources.threads} --keep {params.cd_samples} --qcovar {params.pca_f}'

rule run_gwas_surgery:
    input:
        mbfile=GRM_GWAS_DIR+"geno_chrs.txt",
        pheno_file='/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/scripts/pheno/{pheno}_pheno.txt'
    output:
        "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/output/gwas/{pheno}.fastGWA"
    params:
        PATHS_FILE=PATHS_FILE,
        CONFIGS_FILE=CONFIGS_FILE,
        sparse_grm_prefix=GRM_GWAS_DIR+'sparse_grm',
        gwas_output_prefix="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/output/gwas/",
        wgs_samples="/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/analysis_22_03_2022/pheno/WGS.sample",
        cd_samples="/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/analysis_22_03_2022/pheno/CD.sample",
        pca_f="/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/analysis_22_03_2022/GRM_sparse/GWAS/geno_pc.eigenvec"
    resources:
        mem_mb=10000,
        threads=10
    shell:
        'source {params.PATHS_FILE}; source {params.CONFIGS_FILE}; /software/team152/oe2/bin/gcta64 --mbfile {input.mbfile} --out {params.gwas_output_prefix}{wildcards.pheno} --fastGWA-mlm-binary --grm-sparse {params.sparse_grm_prefix} --joint-covar --pheno {input.pheno_file} --remove {params.wgs_samples} --threads {resources.threads} --keep {params.cd_samples} --qcovar {params.pca_f}'
