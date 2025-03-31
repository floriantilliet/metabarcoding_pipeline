#!/bin/bash

module purge
module load bioinfo/Cutadapt/5.0
module load bioinfo/VSEARCH/2.29.3
module load bioinfo/FastQC/0.12.1

OUTPUT_DIR=OTU_CLUSTERING
nb_cores=2

while getopts "o:i:d:m:f:r:" option
do
        case $option in
                h)
			        echo "help message"
                    exit 0
                    ;;
                o)
                    OUTPUT_DIR="$OPTARG"
                    ;;
                i)
			        INPUT_FILES="$OPTARG"
                    ;;
                f)
                    PRIMER_F="$OPTARG"
                    ;;
                r)
                    PRIMER_R="$OPTARG"
                    ;;
                m)
                    MAPPING="$OPTARG"
        esac
done

DIR="$(pwd)"
cd $DIR
mkdir $OUTPUT_DIR

# Reverse complement the reverse primer
ANTI_PRIMER_R=$( echo "${PRIMER_R}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )

# Step 1-A: Prepare the data

# List all samples, remove the extension, and store the names
ls $INPUT_FILES/*.f* | sed 's#.*/##' | sed 's/.f.*//' | sed 's/_R[12].*//' | sort -u > $OUTPUT_DIR/list_sample.txt

# Check the quality encoding (33 or 64?)
OUTPUT="$OUTPUT_DIR/Quality.encoding.log"
FIRST_SAMPLE=$(head -n 1 $OUTPUT_DIR/list_sample.txt)
INPUT_SAMPLE=$(ls $INPUT_FILES | grep -E "^${FIRST_SAMPLE}_R1\.(fastq\.gz|fastq)$" | sed "s|^|$INPUT_FILES/|")

if file "$INPUT_SAMPLE" | grep -q "compressed"; then
    zcat "$INPUT_SAMPLE" > $OUTPUT_DIR/temp_decompressed.fastq
fi
vsearch \
    --threads 0 \
    --fastq_chars $OUTPUT_DIR/temp_decompressed.fastq 2> ${OUTPUT}
# rm -f temp_decompressed.fastq
QUALITY_ENCODING=$(grep "Guess: Original" $OUTPUT | awk -F'[+]' '{print $2}' | awk -F')' '{print $1}')

# Step 1-B: Merge the paired-end reads

mkdir $OUTPUT_DIR/merged_reads/

for sample in $(cat $OUTPUT_DIR/list_sample.txt); do
    
    INPUT_R1=$(ls $INPUT_FILES | grep -E "^${sample}_R1\.(fastq\.gz|fastq)$" | sed "s|^|$INPUT_FILES/|")
    INPUT_R2=$(ls $INPUT_FILES | grep -E "^${sample}_R2\.(fastq\.gz|fastq)$" | sed "s|^|$INPUT_FILES/|")
    
    OUTPUT="$OUTPUT_DIR/merged_reads/"$sample".fastq"

    vsearch \
        --threads 0 \
        --fastq_mergepairs ${INPUT_R1} \
        --reverse ${INPUT_R2} \
        --fastq_ascii $QUALITY_ENCODING \
        --fastqout $OUTPUT \
        --quiet 2>> ${OUTPUT/.fastq/.log}
    
done

# Step 1-C: Checking the quality with FASTQC

# FastQC is available on https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip)

mkdir $OUTPUT_DIR/FastQC/

for sample in $(cat $OUTPUT_DIR/list_sample.txt); do
    
    INPUT="$OUTPUT_DIR/merged_reads/"$sample".fastq"
    
    fastqc -t $nb_cores ${INPUT} -o $OUTPUT_DIR/FastQC/
    
done

#quality filtering
for sample in $(cat $OUTPUT_DIR/list_sample.txt); do

    echo $sample

    INPUT="$OUTPUT_DIR/merged_reads/"$sample".fastq"   # fastq
    OUTPUT="$OUTPUT_DIR/merged_reads/"$sample".fasta"  # fasta

    vsearch \
    --threads 0 \
    --fastq_filter $INPUT \
    --fastq_maxns 0 \
    --fastq_maxee 2 \
    --fastaout "${OUTPUT}"

done

# Check if the mapping file is specified
if [ -z "$MAPPING" ]; then
    mkdir $OUTPUT_DIR/merged_derep_reads/
    for sample in $(cat $OUTPUT_DIR/list_sample.txt); do

    echo $sample

    INPUT="$OUTPUT_DIR/merged_reads/"$sample".fasta"  # fasta
    OUTPUT="$OUTPUT_DIR/merged_derep_reads/"$sample"_dereplicated.fasta"  # fasta

     # Dereplicate at the sample level
    vsearch --quiet \
        --threads 0 \
        --derep_fulllength "temp_"$sample".fasta" \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel_sha1 \
        --output "${OUTPUT}"

    done
    exit 0
fi

# Step Demultiplex the sequences

barcodeFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'barcodeFw' | cut -d: -f1) - 1 ))" # getting the index of column with forward barcode etc. "-1" is applied as array below counts from 0
primerFwColumnIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'primerFw' | cut -d: -f1) -1 ))"
barcodeRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'barcodeRev' | cut -d: -f1) -1 ))"
primerRevColumnsIdx="$(( $(sed -n $'1s/\t/\\\n/gp' $MAPPING | grep -nx 'primerRev' | cut -d: -f1) -1 ))"
 
# Define binaries, temporary files and output files
mkdir "$OUTPUT_DIR/Demultiplexed_data/"
MIN_LENGTH=200

for sample in $(cat $OUTPUT_DIR/list_sample.txt); do
    INPUT="$OUTPUT_DIR/merged_reads/${sample}.fasta"
    INPUT_REVCOMP="${INPUT/.fasta/_RC.fasta}"

    # Reverse complement fastq file
    vsearch --quiet \
        --threads 0 \
        --fastx_revcomp "${INPUT}" \
        --fastaout "${INPUT_REVCOMP}"

    while read -r line; do
        if [[ ! "$line" =~ ^#.* && ! "$line" =~ ^Sample.* ]]; then
            IFS=$'\t' read -r -a array <<< "$line"

            # Get sequences
            FwBarcode="${array[${barcodeFwColumnIdx}]}"
            FwPrimer="${array[${primerFwColumnIdx}]}"
            RevBarcode="${array[${barcodeRevColumnsIdx}]}"
            RevBarcodeRC=$( echo "${RevBarcode}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
            RevPrimer="${array[${primerRevColumnsIdx}]}"
            RevPrimerRC=$( echo "${RevPrimer}" | tr ACGTacgtYyMmRrKkBbVvDdHh TGCAtgcaRrKkYyMmVvBbHhDd | rev )
            SAMPLE_NAME="${array[0]}"
            
            # Output file names
            LOG="$OUTPUT_DIR/Demultiplexed_data/${SAMPLE_NAME}_${sample}.log"
            FINAL_FASTA="$OUTPUT_DIR/Demultiplexed_data/${SAMPLE_NAME}_${sample}.fas"
            
            # Some information
            echo "${SAMPLE_NAME} (${sample}) is being processed.."
            echo "Barcode Fw: ${FwBarcode}"
            echo "Primer Fw: ${FwPrimer}"
            echo "Primer Rev (RC): ${RevPrimerRC}"
            echo "Barcode Rev (RC): ${RevBarcodeRC}"

            if [ -f "${FINAL_FASTA}" ]; then
                echo "${SAMPLE_NAME} (${sample}) has already been processed. Skipping..."
                continue
            fi
            
            function trim_without_ambiguity {

                SEQTOT="${FwBarcode}${FwPrimer}${RevPrimerRC}${RevBarcodeRC}"
                MIN_MATCHED=${#SEQTOT}
                ERROR_RATE=0

                cat "${INPUT}" "${INPUT_REVCOMP}" | cutadapt -j 0 -g "${FwBarcode}${FwPrimer}...${RevPrimerRC}${RevBarcodeRC}" --discard-untrimmed --minimum-length "${MIN_LENGTH}" -O ${MIN_MATCHED} -e "${ERROR_RATE}" - 2> "${LOG}" > "$OUTPUT_DIR/temp_${sample}.fasta"
            }
            

            trim_without_ambiguity

            # Dereplicate at the study level
            vsearch --quiet \
                --threads 0 \
                --derep_fulllength "$OUTPUT_DIR/temp_${sample}.fasta" \
                --sizein \
                --sizeout \
                --fasta_width 0 \
                --relabel_sha1 \
                --output "${FINAL_FASTA}" 2>> "${LOG}"

        fi
    done < "${MAPPING}"
done

# Note: the option sha1 (encoding system) is giving the same names to the identical amplicons across samples

# Remove the files that are no longer useful
# rm -f "${INPUT}" "${INPUT_REVCOMP}" "temp.fastq" "temp.fasta"

