#!/usr/bin/env bash

# desc: Compute PGS with PLINK2 --score


PGEN=$1
PVAR=$2
PSAM=$3
SAMPLES=$4
FREQS=$5
COEFF=$6
OUT=$7

PLINK2="/path/to/plink2"


start=`date +%s`

# get number of columns in coefficients file
COLS=$(awk '{print NF}' $COEFF | sort -nu | tail -n 1)

if [[ $COLS -gt 5 ]]
then
    $PLINK2 --pgen $PGEN \
            --pvar $PVAR \
            --psam $PSAM \
            --keep $SAMPLES \
            --read-freq $FREQS \
            --score $COEFF header-read \
            --score-col-nums 5-$COLS \
            --threads 2 \
            --memory 32000 \
            --out $OUT
elif [[ $COLS -eq 5 ]]
then
    $PLINK2 --pgen $PGEN \
            --pvar $PVAR \
            --psam $PSAM \
            --keep $SAMPLES \
            --read-freq $FREQS \
            --score $COEFF header-read \
            --score-col-nums 5 \
            --threads 1 \
            --memory 16000 \
            --out $OUT
else
    echo "Error: coefficients file has fewer than 5 columns"
    exit 1
fi

end=`date +%s`
runtime=$((end-start))

echo $runtime >> $OUT-runtime.tab
