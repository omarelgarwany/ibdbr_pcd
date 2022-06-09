#!/bin/bash

source /lustre/scratch123/hgi/projects/ibdgwas_bioresource/oe2/scripts/paths.sh

for((i=1;i<23;i++))
do
       plink --bfile "$BDIR"chr${i}_hg38 --exclude "$exclude_regions_f" --maf "${IDPSNPS_MAF}" --indep-pairwise "$IDPSNPS_WINDOW" "$IDPSNPS_ITERATION" "$IDPSNPS_R2THRESH" --out "$IDPSNPS_DIR"chr${i}
done
