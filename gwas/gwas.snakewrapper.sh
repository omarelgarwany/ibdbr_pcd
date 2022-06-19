#!/bin/bash

source /lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/scripts/paths.sh


PATHS_FILE=/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/scripts/paths.sh
CONFIGS_FILE=/lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/scripts/configs.sh


group="team152"
workdir="$IBDBR_SNAKEMAKE_WORKDIR"gwas/
worfklow_prefix="ibdbr"

mkdir -p "$GRM_GWAS_DIR"
mkdir -p "$IDPSNPS_DIR"

snakemake -j 100 --until all --latency-wait 90 -p --keep-going --default-resources threads=1 mem_mb=2000 --config PATHS_FILE="${PATHS_FILE}" CONFIGS_FILE="${CONFIGS_FILE}" BDIR="${BDIR}" IDPSNPS_DIR="${IDPSNPS_DIR}" GRM_GWAS_DIR="${GRM_GWAS_DIR}" --directory "${workdir}" --use-conda --conda-prefix "${IBDBR_CONDA_ENV}" --cluster "rm logs/cluster/${worfklow_prefix}_{rule}/{rule}.{wildcards}.out logs/cluster/${worfklow_prefix}_{rule}/{rule}.{wildcards}.err ; mkdir -p logs/cluster/${worfklow_prefix}_{rule}; bsub -R \"rusage[mem={resources.mem_mb}] select[model==Intel_Platinum && mem>{resources.mem_mb}] span[hosts=1]\" -M {resources.mem_mb} -m \"modern_hardware\" -n {resources.threads} -J \"${worfklow_prefix}_{rule}.{wildcards}\" -G ${group} -o logs/cluster/${worfklow_prefix}_{rule}/{rule}.{wildcards}.out -e logs/cluster/${worfklow_prefix}_{rule}/{rule}.{wildcards}.err" -s gwas.Snakefile
