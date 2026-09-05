#!/usr/bin/env bash

# Run interaction GWAS for LOCO independent hits


PHEN=$1
CHR=$2
PGEN=$3
PVAR=$4
PSAM=$5
SAMPLES=$6
FREQS=$7
VARIANTS=$8
SNP=$9

PLINK2="/path/to/plink2"



# interaction GWAS
## make temporary file with A1 allele for this SNP
line=$(grep "${SNP//_/:}" ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/hits/${PHEN}-${SAMPLES}-loco-indep-hits.tab)
echo $line | awk '{print $3, $7}' > ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}-A1-${CHR}.txt

## extract genotypes
SNP_CHR="${SNP%%_*}"
$PLINK2 --pgen ../data/imputed-genotypes/ukb_imp_chr${SNP_CHR}_v3.pgen \
        --pvar ../data/imputed-genotypes/alternative_pvar_files/ukb_imp_POSID_INFO_chr${SNP_CHR}_v3.pvar \
        --psam $PSAM \
        --keep ../data/sample-ids/filtered/${SAMPLES}-ids.tab \
        --snp "${SNP//_/:}" \
        --threads 2 \
        --memory 32000 \
        --export A \
        --export-allele ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}-A1-${CHR}.txt \
        --out ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}

## transform them into form required for input into plink
awk '{print $1, $2, $7}' ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw > ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.tmp
mv -f ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.tmp ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw
echo "#$(cat ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw)" > ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw
## add column with square of genotype
awk 'NR==1 {print $0, $3"_SQ"} NR>1 {print $0, ($3)^2}' ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw > ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.tmp
mv -f ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.tmp ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw

## check that correlation between SNP and SNP^2 is not >=0.999
cat > ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.R <<- EOF
library(data.table)

geno <- fread("../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw",
              data.table = FALSE)

cat(cor(geno[, 3], geno[, 4]))
EOF

dom_cor=$(Rscript ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.R)

## log correlation
if (( $CHR == 1 ))
then
    echo "Correlation between SNP genotype and its square is ($dom_cor)." > \
         ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${PHEN}.${SAMPLES}.${SNP}.int.loco.snp_snp2_cor.log
fi

if (( $(echo "$dom_cor < 0.999" | bc -l) ))
then

    ## run GWAS with SNP^2 variable
    $PLINK2 --pgen $PGEN \
            --pvar $PVAR \
            --psam $PSAM \
            --keep ../data/sample-ids/filtered/${SAMPLES}-ids.tab \
            --read-freq $FREQS \
            --extract $VARIANTS \
            --pheno ../results/01-gwas/residuals-covars/${PHEN}/${PHEN}-${SAMPLES}-resid-covars.tab \
            --pheno-name ${PHEN}.res.cov \
            --covar ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw \
            --covar-col-nums 3, 4 \
            --threads 2 \
            --memory 32000 \
            --glm log10 interaction \
            --parameters 1-4 \
            --vif 999 \
            --out ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${PHEN}.${SAMPLES}.${SNP}.int.loco.chr${CHR}

else

    ## run GWAS without SNP^2 variable
    $PLINK2 --pgen $PGEN \
            --pvar $PVAR \
            --psam $PSAM \
            --keep ../data/sample-ids/filtered/${SAMPLES}-ids.tab \
            --read-freq $FREQS \
            --extract $VARIANTS \
            --pheno ../results/01-gwas/residuals-covars/${PHEN}/${PHEN}-${SAMPLES}-resid-covars.tab \
            --pheno-name ${PHEN}.res.cov \
            --covar ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw \
            --covar-col-nums 3 \
            --threads 2 \
            --memory 32000 \
            --glm log10 interaction \
            --vif 999 \
            --out ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${PHEN}.${SAMPLES}.${SNP}.int.loco.chr${CHR}

fi

gzip -f ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${PHEN}.${SAMPLES}.${SNP}.int.loco.chr${CHR}.${PHEN}.res.cov.glm.linear

## remove temporary files
rm ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}-A1-${CHR}.txt
rm ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.raw
rm ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.log
rm ../results/03-interaction-gwas/indep-hits/${PHEN}/snp-pgs/gwas/${SNP}.geno.${CHR}.R
