#!/bin/bash

# This script performs OTU clustering using VSEARCH and generates a contingency table.
# It requires the following software: VSEARCH 2.29.3, Python 3.11.1.
# It also requires the following parameters: path to the directory containing input fasta files, path to the database, path to the scripts, output directory, abundance threshold, and identity percentage.

module purge
module load bioinfo/VSEARCH/2.29.3
module load devel/python/Python-3.11.1

nb_cores=2
OUTPUT_DIR=OTU_CLUSTERING
IDENTITY_PERCENTAGE=0.97
path_scripts=$0/../python_scripts

while getopts "i:s:d:o:a:c:" option
do
        case $option in
                h)
                    echo "Usage: $0 -i <input_files> -s <path_scripts> -d <path_db> -o <output_directory> -c <identity_percentage>"
                    echo "  -i  Input directory containing fasta files (can be compressed)"
                    echo "  -d  Path to the database"
                    echo "  -o  Output directory (default: OTU_CLUSTERING)"
                    echo "  -c  Identity percentage for clustering (default: 0.97)"
                    echo "  -h  Show this help message"
                    exit 0
                    ;;
                i)
			        INPUT_FILES="$OPTARG"
                    ;;
                d)
                    path_to_db="$OPTARG"
                    ;;
                o)
                    OUTPUT_DIR="$OPTARG"
                    ;;    
                c)
                    IDENTITY_PERCENTAGE="$OPTARG"
                    ;;

        esac
done

DIR="$(pwd)"
cd $DIR
mkdir $OUTPUT_DIR

#STEP 2
#OTU Clustering

# indicate the name of the database for taxonomic assignation (generated in Step 0):

cat $INPUT_FILES/*.f* > $OUTPUT_DIR/reads_amplicon.fasta

vsearch \
    --derep_fulllength $OUTPUT_DIR/reads_amplicon.fasta \
    --sizein \
    --sizeout \
    --relabel_sha1 \
    --fasta_width 0 \
    --output $OUTPUT_DIR/reads_amplicon_derep.fasta

rm $OUTPUT_DIR/reads_amplicon.fasta

vsearch -sortbysize $OUTPUT_DIR/reads_amplicon_derep.fasta -output $OUTPUT_DIR/reads_amplicon_sorted.fasta -minsize 1

rm $OUTPUT_DIR/reads_amplicon_derep.fasta

vsearch -cluster_size  $OUTPUT_DIR/reads_amplicon_sorted.fasta \
    --threads $nb_cores \
    --id $IDENTITY_PERCENTAGE --centroids $OUTPUT_DIR/reads_OTU97.fasta \
    --uc $OUTPUT_DIR/clusters_OTU97.uc \
    --sizein --sizeout

vsearch --fasta_width 0 --sortbysize $OUTPUT_DIR/reads_OTU97.fasta --output $OUTPUT_DIR/reads_OTU97_final.fasta

python3 $path_scripts/map2qiime.py $OUTPUT_DIR/clusters_OTU97.uc > $OUTPUT_DIR/reads_mapped_OTU97.txt

python3 $path_scripts/make_stats.py $OUTPUT_DIR/reads_OTU97_final.fasta > $OUTPUT_DIR/stats_file_OTU97.txt

vsearch --uchime_denovo $OUTPUT_DIR/reads_OTU97_final.fasta \
    --uchimeout $OUTPUT_DIR/reads_OTU97.uchime \
    --nonchimeras $OUTPUT_DIR/reads_OTU97_nonchimeras.fasta

vsearch --usearch_global $OUTPUT_DIR/reads_OTU97_nonchimeras.fasta \
    --threads $nb_cores \
    --dbmask none \
    --qmask none \
    --rowlen 0 \
    --notrunclabels \
    --userfields query+id1+target \
    --maxaccepts 0 \
    --maxrejects 32 \
    --top_hits_only \
    --output_no_hits \
    --db $path_to_db \
    --id 0.5 \
    --iddef 4 \
    --userout $OUTPUT_DIR/taxonomy_OTU97.txt

STATS="$OUTPUT_DIR/stats_file_OTU97.txt"
OTUS="$OUTPUT_DIR/reads_mapped_OTU97.txt"
REPRESENTATIVES="$OUTPUT_DIR/reads_OTU97_final.fasta"
UCHIME="$OUTPUT_DIR/reads_OTU97.uchime"
ASSIGNMENTS="$OUTPUT_DIR/taxonomy_OTU97.txt"
OTU_TABLE="$OUTPUT_DIR/OTU_table_OTU97.txt"

SCRIPT=$path_scripts"/OTU_contingency_table.py" 

python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${OTUS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    $INPUT_FILES/*.f* > "${OTU_TABLE}"
   
FILTERED="${OTU_TABLE/.txt/_filtered.txt}"
head -n 1 "${OTU_TABLE}" > "${FILTERED}"
cat "${OTU_TABLE}" | awk '$5 == "N" && $4 >= 200 && $2 >= 10' >> "${FILTERED}"

cd $OUTPUT_DIR

python3 $DIR/$path_scripts/identity_distribution.py ../$OTU_TABLE

cd $DIR
