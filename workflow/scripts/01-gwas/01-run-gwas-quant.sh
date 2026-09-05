#!/usr/bin/env bash

# desc: Run GWAS for quantitative trait (assumes covars have been regressed out)


PGEN=$1
PVAR=$2
PSAM=$3
SAMPLES=$4
FREQS=$5
VARIANTS=$6
TRAIT=$7
OUT=$8

PLINK2="/path/to/plink2"


start=`date +%s`

$PLINK2 --pgen $PGEN \
        --pvar $PVAR \
        --psam $PSAM \
        --keep $SAMPLES \
        --read-freq $FREQS \
        --extract $VARIANTS \
        --pheno $TRAIT \
        --threads 2 \
        --memory 32000 \
        --glm log10 allow-no-covars \
        --vif 999 \
        --out $OUT

end=`date +%s`
runtime=$((end-start))

echo $runtime >> $OUT.runtime
