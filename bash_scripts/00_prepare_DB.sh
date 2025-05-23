#!/bin/bash

# This script prepares a database for taxonomic assignation by trimming primers from the input sequences using Cutadapt.
# It requires the following software: Cutadapt 5.0.
# It also requires the following parameters: raw database input fasta file (can be compressed), forward primer, reverse primer.

module purge
module load bioinfo/Cutadapt/5.0

DIR="$(pwd)"
cd $DIR

while getopts "hi:o:f:r:" option
do
        case $option in
                h)
                    echo "Usage: $0 -i <input_file> -f <forward_primer> -r <reverse_primer> -o <output_directory>"
                    echo "  -i  Input raw database file (FASTA format, can be compressed)"
                    echo "  -f  Forward primer sequence"
                    echo "  -r  Reverse primer sequence"
                    echo "  -o  Output directory (default: current directory)"
                    exit 0
                    ;;
                i)
                    INPUT="$OPTARG"
                    ;;
                f)
                    PRIMER_F="$OPTARG"
                    ;;
                r)
                    PRIMER_R="$OPTARG"
                    ;;
                o)
                    DIR="$OPTARG"
                    ;;
        esac
done

if [ -z "$INPUT" ]; then
    echo "Error: An input file must be specified."
    exit 1
fi

if (file $INPUT | grep -q "compressed data");
then
    command="zcat"
else
    command="cat"
fi

PRIMER_F_SPECIFIED=false
PRIMER_R_SPECIFIED=false

if [ -n "$PRIMER_F" ]; then
    PRIMER_F_SPECIFIED=true
fi

if [ -n "$PRIMER_R" ]; then
    PRIMER_R_SPECIFIED=true
fi

if [ "$PRIMER_F_SPECIFIED" = false ] && [ "$PRIMER_R_SPECIFIED" = false ]; then
    echo "Error: At least one primer must be specified."
    exit 1
fi

# Define variables and output files primers 
BASENAME=$(basename "$INPUT")
OUTPUT=$DIR/${BASENAME%.*}_database.fasta
LOG=$DIR/${BASENAME%.*}_database.log
OUTPUT_R=$DIR/${BASENAME%.*}_database_R.fasta
LOG_R=$DIR/${BASENAME%.*}_database_R.log

# Reverse complement of the reverse primer
if [ "$PRIMER_R_SPECIFIED" = false ]; then
    ANTI_PRIMER_R=$( echo "${PRIMER_R}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
fi

# Paramerers of cutadapt
MIN_LENGTH=32 # minimum sequence length
MIN_F=$(( ${#PRIMER_F} / 2 )) # should match at least 50% of the forward primer
MIN_R=$(( ${#PRIMER_R} / 2 )) # should match at least 50% of the reverse primer
ERROR_RATE=0.2
CUTADAPT="cutadapt --cores 0 --discard-untrimmed --minimum-length ${MIN_LENGTH} --no-indels -e ${ERROR_RATE}"


# Trim forward
if [ "$PRIMER_F_SPECIFIED" = true ]; then
    ${command} "${INPUT}" | sed '/^>/ ! s/U/T/g' | \
    ${CUTADAPT} -g "${PRIMER_F}" -O "${MIN_F}" - 2> "${LOG}" | \
    sed '/^>/ s/;/|/g ; /^>/ s/ /_/g' > "${OUTPUT}"
else 
    OUTPUT=$INPUT
fi

# Trim reverse
if [ "$PRIMER_R_SPECIFIED" = true ]; then
    ${command} "${OUTPUT}" | sed '/^>/ ! s/U/T/g' | \
    ${CUTADAPT} -a "${PRIMER_R}" -O "${MIN_R}" - 2> "${LOG_R}" | \
    sed '/^>/ s/;/|/g ; /^>/ s/ /_/g' > "${OUTPUT_R}"
fi