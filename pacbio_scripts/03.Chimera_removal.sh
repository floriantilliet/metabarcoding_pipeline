#!/bin/bash

module load bioinfo/VSEARCH/2.29.3
module load tools/ripgrep/12.1.1

## Script for sample-wise chimera removal
##   1. reference-based chimera removal
##   2. de novo chimera search. Sequences identified as putative chimeras will have lower priority in clustering

## Input = FASTQ, Output = FASTQ

# $INPUT_FILE = input file name
# $OUTPUT_PREFIX = output prefix
# $LOG_FILE = log file
# $REFERENCE_DB = reference database path

INPUT_FILE="../02extractITS/Extracted_ITS_full.fq.gz"
OUTPUT_PREFIX="chimera_removal"
LOG_FILE="chimera_removal.log"
REFERENCE_DB="UNITE_9.1_beta_reference.fasta"

echo -e "Start " "$INPUT_FILE" > "$LOG_FILE"
TMPDIR="chimera_tmp"

mkdir -p "$TMPDIR"

echo -e "\n### Dereplication\n" >> "$LOG_FILE"

## Dereplicate
vsearch \
  --fastx_uniques "$INPUT_FILE" \
  --fasta_width 0 \
  --threads 1 \
  --sizein --sizeout \
  --fastaout "tmp.fa" \
  2>> "$LOG_FILE" \

gzip -c tmp.fa > "$TMPDIR"/derep.fa.gz

echo -e "\n### Reference-based chimera removal\n" >> "$LOG_FILE"

## Reference-based
vsearch \
  --uchime_ref "$TMPDIR"/derep.fa.gz \
  --db "$REFERENCE_DB" \
  --fasta_width 0 \
  --sizein --sizeout \
  --threads 1 \
  --chimeras "$TMPDIR"/Ref_Chimeras.fasta \
  --nonchimeras "$TMPDIR"/NoChimera_Ref.fasta \
  2>> "$LOG_FILE" \

gzip -c "$TMPDIR"/Ref_Chimeras.fasta > "$TMPDIR"/Ref_Chimeras.fasta.gz

echo -e "\n### De novo chimera removal\n" >> "$LOG_FILE"

## De-novo
numseqs=$(grep -c "^>" "$TMPDIR"/derep.fa.gz)
if [[ $numseqs -gt 3 ]] ; then

  vsearch \
    --uchime_denovo "$TMPDIR"/derep.fa.gz \
    --chimeras "$TMPDIR"/Chimeras_Denovo.fasta \
    --nonchimeras "$TMPDIR"/NoChimera_DeNovo.fasta \
    --fasta_width 0 \
    --sizein --xsize \
    2>> "$LOG_FILE" \

gzip -c "$TMPDIR"/Chimeras_Denovo.fasta > "$TMPDIR"/Chimeras_Denovo.fasta.gz

else
  echo "The number of sequences is too small. Do not perform de novo search." >> "$LOG_FILE"
fi


echo -e "\n### Sequence extraction\n" >> "$LOG_FILE"

## Get sequences
if [ -f "$TMPDIR"/NoChimera_DeNovo.fasta.gz ]; then
  zcat "$TMPDIR"/NoChimera_DeNovo.fasta.gz | awk '$0 !~ /^>/ {print toupper($0)}' |  sort > "$TMPDIR"/tmp_ok.txt
elif [ -f "$TMPDIR"/NoChimera_Ref.fasta ]; then
  awk '$0 !~ /^>/ {print toupper($0)}' "$TMPDIR"/NoChimera_Ref.fasta | sort > "$TMPDIR"/tmp_ok.txt
else
  touch "$TMPDIR"/tmp_ok.txt
fi

if [ -f "$TMPDIR"/Chimeras_Denovo.fasta ]; then
  awk '$0 !~ /^>/ {print toupper($0)}' "$TMPDIR"/Chimeras_Denovo.fasta | sort > "$TMPDIR"/tmp_chim_denovo.txt
else
  touch "$TMPDIR"/tmp_chim_denovo.txt
fi

if [ -f "$TMPDIR"/Ref_Chimeras.fasta.gz ]; then
  zcat "$TMPDIR"/Ref_Chimeras.fasta.gz | awk '$0 !~ /^>/ {print toupper($0)}' |  sort > "$TMPDIR"/chim_ref.txt
else
  touch "$TMPDIR"/chim_ref.txt
fi


## Remove de-novo chimeras identified with reference database
comm -13 "$TMPDIR"/chim_ref.txt "$TMPDIR"/tmp_ok.txt > "$TMPDIR"/ok.txt

## Remove seqs identified as reference chimeras
comm -13 "$TMPDIR"/chim_ref.txt "$TMPDIR"/tmp_chim_denovo.txt > "$TMPDIR"/chim_denovo.txt

echo -e "Number of unique non-chimeric sequences: " $(wc -l < "$TMPDIR"/ok.txt) >> "$LOG_FILE"
echo -e "Number of unique reference-based chimeric sequences: " $(wc -l < "$TMPDIR"/chim_ref.txt) >> "$LOG_FILE"
echo -e "Number of unique de novo chimeric sequences: " $(wc -l < "$TMPDIR"/chim_denovo.txt) >> "$LOG_FILE"

#### To improve speed, split patterns into chunks (with 100 lines each), then ripgrep

## Good sequences
if [[ -s "$TMPDIR"/ok.txt ]]; then

    split -l 100 "$TMPDIR"/ok.txt "$TMPDIR"/patt.split.
    for CHUNK in "$TMPDIR"/patt.split.* ; do
      rg -z -B 1 -A 2 -x -f "$CHUNK" "$INPUT_FILE" >> "$TMPDIR"/tmp_NoChimera.fq
    done
    rm "$TMPDIR"/patt.split.*

    sed '/^--$/d' "$TMPDIR"/tmp_NoChimera.fq | gzip -6 > "$OUTPUT_PREFIX"_NoChimera.fq.gz
fi

## De novo chimeras
if [[ -s "$TMPDIR"/chim_denovo.txt ]]; then

  split -l 100 "$TMPDIR"/chim_denovo.txt "$TMPDIR"/patt.split.
  for CHUNK in "$TMPDIR"/patt.split.* ; do
    rg -z -B 1 -A 2 -x -f "$CHUNK" "$INPUT_FILE" >> "$TMPDIR"/tmp_ChimDeNov.fq
  done
  rm "$TMPDIR"/patt.split.*

  sed '/^--$/d' "$TMPDIR"/tmp_ChimDeNov.fq | gzip -6 > "$OUTPUT_PREFIX"_ChimDeNov.fq.gz
fi

## Reference-based chimeras
if [[ -s "$TMPDIR"/chim_ref.txt ]]; then
  split -l 100 "$TMPDIR"/chim_ref.txt "$TMPDIR"/patt.split.
  for CHUNK in "$TMPDIR"/patt.split.* ; do
    rg -z -B 1 -A 2 -x -f "$CHUNK" "$INPUT_FILE" >> "$TMPDIR"/tmp_ChimRef.fq
  done
  rm "$TMPDIR"/patt.split.*

  sed '/^--$/d' "$TMPDIR"/tmp_ChimRef.fq  | gzip -6 > "$OUTPUT_PREFIX"_ChimRef.fq.gz
fi

## Clean up
#rm -r "$TMPDIR"

echo -e "\nDone " "$INPUT_FILE" >> "$LOG_FILE"

## Example:
# ./chimera_rm.sh Sample.fq.gz Samp Samp.log Reference_chimera_DB.fasta.gz
