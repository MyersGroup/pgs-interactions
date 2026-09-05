# desc: Compute allele counts/freqs, genotype missing rates and HWE p-values


PGEN=$1
PVAR=$2
PSAM=$3
SAMPLES=$4  # list of sample IDs to consider
OUT_PATH=$5
OUT_PREF=$6

PLINK2="/path/to/plink2"



# minor allele frequencies
$PLINK2 --pgen $PGEN \
        --pvar $PVAR \
        --psam $PSAM \
        --keep $SAMPLES \
        --freq \
        --threads 4 \
        --out "${OUT_PATH}/allele-freq/${OUT_PREF}"


# minor allele counts
$PLINK2 --pgen $PGEN \
        --pvar $PVAR \
        --psam $PSAM \
        --keep $SAMPLES \
        --freq counts \
        --threads 4 \
        --out "${OUT_PATH}/allele-counts/${OUT_PREF}"


# genotype missing rates
$PLINK2 --pgen $PGEN \
        --pvar $PVAR \
        --psam $PSAM \
        --keep $SAMPLES \
        --missing variant-only \
        --threads 4 \
        --out "${OUT_PATH}/missing-rates/${OUT_PREF}"


# Hardy-Weinberg eq p-values
$PLINK2 --pgen $PGEN \
        --pvar $PVAR \
        --psam $PSAM \
        --keep $SAMPLES \
        --hardy midp \
        --threads 4 \
        --out "${OUT_PATH}/hardy-eq/${OUT_PREF}"
