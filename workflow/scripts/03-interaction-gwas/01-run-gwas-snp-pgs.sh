#!/usr/bin/env bash

# desc: Run GWAS for model: phen ~ PGS + PGS*SNP + SNP
#       For LOCO PGS we run: phen_res_full ~ PGS_LOCO + PGS_LOCO*SNP + SNP (to avoid spurious assocations)


PHEN=$1
CHR=$2
PGEN=$3
PVAR=$4
PSAM=$5
SAMPLES=$6
FREQS=$7
VARIANTS=$8

PLINK2="/path/to/plink2"



# full PGS
$PLINK2 --pgen $PGEN \
        --pvar $PVAR \
        --psam $PSAM \
        --keep ../data/sample-ids/filtered/${SAMPLES}-ids.tab \
        --read-freq $FREQS \
        --extract $VARIANTS \
        --pheno ../results/02-pgs/ancestry-interactions/${PHEN}/r2_0.9-kb_500/${PHEN}-${SAMPLES}-resid-covars-pgs-full_mc-ac-centre.tab \
        --pheno-name ${PHEN}.res.cov.pgs_full_mc.ac.ctr \
        --covar ../results/02-pgs/mean-corrected/${PHEN}/r2_0.9-kb_500/${PHEN}-${SAMPLES}-pgs-full_mc.tab \
        --covar-name pgs_full_mc \
        --threads 2 \
        --memory 32000 \
        --glm log10 interaction \
        --vif 999 \
        --out ../results/03-interaction-gwas/plink-output/${PHEN}/snp-pgs/${PHEN}.${SAMPLES}.full.chr${CHR}


# LOCO PGS
$PLINK2 --pgen $PGEN \
        --pvar $PVAR \
        --psam $PSAM \
        --keep ../data/sample-ids/filtered/${SAMPLES}-ids.tab \
        --read-freq $FREQS \
        --extract $VARIANTS \
        --pheno ../results/02-pgs/ancestry-interactions/${PHEN}/r2_0.9-kb_500/${PHEN}-${SAMPLES}-resid-covars-pgs-full_mc-ac-centre.tab \
        --pheno-name ${PHEN}.res.cov.pgs_full_mc.ac.ctr \
        --covar ../results/02-pgs/mean-corrected/${PHEN}/r2_0.9-kb_500/${PHEN}-${SAMPLES}-pgs-loco${CHR}_mc.tab \
        --covar-name pgs_loco${CHR}_mc \
        --threads 2 \
        --memory 32000 \
        --glm log10 interaction \
        --vif 999 \
        --out ../results/03-interaction-gwas/plink-output/${PHEN}/snp-pgs/${PHEN}.${SAMPLES}.loco.chr${CHR}
