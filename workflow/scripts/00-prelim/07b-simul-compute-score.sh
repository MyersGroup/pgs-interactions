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

$PLINK2 --pgen $PGEN \
        --pvar $PVAR \
        --psam $PSAM \
        --keep $SAMPLES \
        --read-freq $FREQS \
        --score $COEFF header \
        --score-col-nums 3-4 \
        --threads 1 \
        --memory 16000 \
        --out $OUT

end=`date +%s`
runtime=$((end-start))

echo $runtime >> $OUT-runtime.tab
