#!/bin/bash

module purge

module load devel/java/17.0.6
module load containers/singularity/3.9.9
module load bioinfo/Nextflow/24.10.0
module load bioinfo/VSEARCH/2.29.3
module load devel/python/Python-3.11.1

nextflow pull vmikk/NextITS

#download the chimera db from the example: https://owncloud.ut.ee/owncloud/s/iaQ3i862pjwYgdy
#mkdir chimeraDB
#cd chimeraDB
#wget https://owncloud.ut.ee/owncloud/s/iaQ3i862pjwYgdy
#cd ..

nextflow run vmikk/NextITS -r main \
  -profile singularity \
  -resume \
  --step "Step1" \
  --input          "data/" \
  --demultiplexed true \
  --chimera_db "$PWD/chimeraDB/UN95_chimera.udb" \ 
  --primer_forward "TACACACCGCCCGTCG" \
  --primer_reverse "CCTSCSCTTANTDATATGC" \
  --outdir         "Step1_Results/Run01"

  nextflow run vmikk/NextITS -r main \
  -resume \
  -profile     singularity \
  --step       "Step2" \
  --data_path  "$PWD/Step1_Results/" \
  --outdir     "Step2_Results" \
  --clustering_method "vsearch" \
  --otu_id 0.98


mkdir Step3_Results
path_scripts="../scripts/python_scripts"
path_to_db=../data/databases/UNITE_DBs\(ITS\)/ITS1/sh_general_release_dynamic_s_all_04.04.2024.fasta
vsearch --fasta_width 0 --sortbysize Step2_Results/03.Clustered_VSEARCH/Clustered.fa.gz --output Step3_Results/reads_OTU97_final.fasta

python3 $path_scripts/map2qiime.py <(zcat Step2_Results/03.Clustered_VSEARCH/Clustered.uc.gz) > Step3_Results/reads_mapped_OTU97.txt

python3 $path_scripts/make_stats.py Step3_Results/reads_OTU97_final.fasta > Step3_Results/stats_file_OTU97.txt

vsearch --uchime_denovo Step3_Results/reads_OTU97_final.fasta \
  --uchimeout Step3_Results/reads_OTU97.uchime \
  --nonchimeras Step3_Results/reads_OTU97_nonchimeras.fasta

vsearch --usearch_global Step3_Results/reads_OTU97_final.fasta \
    --threads 0 \
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
    --userout Step3_Results/taxonomy_OTU97.txt

STATS="Step3_Results/stats_file_OTU97.txt"
OTUS="Step3_Results/reads_mapped_OTU97.txt"
REPRESENTATIVES="Step3_Results/reads_OTU97_final.fasta"
UCHIME="Step3_Results/reads_OTU97.uchime"
ASSIGNMENTS="Step3_Results/taxonomy_OTU97.txt"
OTU_TABLE="Step3_Results/OTU_table_OTU97.txt"

SCRIPT=$path_scripts"/OTU_contingency_table.py" 

python3 \
    "${SCRIPT}" \
    "${REPRESENTATIVES}" \
    "${STATS}" \
    "${OTUS}" \
    "${UCHIME}" \
    "${ASSIGNMENTS}" \
    <(zcat Step2_Results/01.Dereplicated/Dereplicated.fa.gz) > "${OTU_TABLE}"

FILTERED="${OTU_TABLE/.txt/_filtered.txt}"
head -n 1 "${OTU_TABLE}" > "${FILTERED}"
cat "${OTU_TABLE}" | awk '$5 == "N" && $4 >= 200 && $2 >= $ABUNDANCE' >> "${FILTERED}"

python3 $path_scripts/identity_distribution.py $OTU_TABLE
