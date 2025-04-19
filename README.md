# Fungi metabarcoding analysis pipeline

This repository contains the code and resources used during my internship at CRBE. The pipeline is designed to process and analyze fungi metabarcoding data by using the following process:
- Step 0: Create a database for taxonomic assignation.
- Step 1: Merge and demultiplex the raw reads from metabarcoding.
- Step 2: OTU Clustering a 97% and performs a taxonomic assignation of the OTUs.

## Prerequisites
These scripts require a shell/Linux computing environment and the following software:
- Cutadapt 5.0
- VSEARCH 2.29.3
- FastQC 0.12.1
- Python 3.11.1
- R 4.4.3

- Nextflow 24.10.0

## Usage
You can either run the whole script at once, or each of the steps individually.

Running the whole script:
```bash
./Illumina_Pipeline.sh 
```
- `-h` Display the help message.

Running each step individually:

Step 0: Creating the database for taxonomic assignation.
``` bash
./00_prepare_DB.sh -i path_to_raw_database [-o output_directory] [-f forward_primer_sequence] [-r reverse_primer_sequence] [-h]
```
- `-h` Display the help message.
- `-i` Provide the path to input FASTA file of the raw database. Compressed formats are accepted.
- `-o` Specify an output directory (default: current directory).
- `-f` `-r` Specify the forward/reverse primer sequence. Ex: "ACAGTCGTCGAT"  

Step 1: Merging/demultiplexing the data.
```bash
./01_merging_demultiplexing.sh -i path_to_directory_containing_raw_reads [-o output_directory] [-m path_to_mapping_file] [-h]
```
- `-h` Display the help message.
- `-i` Provide the path to the directory containing the FASTQ files of the raw reads. Compressed formats are accepted.
- `-o` Specify an output directory (default: current directory).
- `-m` Specify a mapping file with the following columns: SampleID, barcodeFw, primerFw, primerRev, barcodeRev. Demultiplexing will only be done if a mapping file is specified. Otherwise only the merging will be done.

Step 2: Clustering at 97% and taxonomic assignation.
```bash
./02_OTU_Clustering.sh -i path_to_directory_containing_reads -d path_to_database [-o output_directory] [-c identity_percentage]
```
- `-h` Display the help message.
- `-i` Provide the path to the directory containing the FASTA files of the reads. Compressed formats are accepted.
- `-o` Specify an output directory (default: current directory).
- `-d` Provide the path to the database FASTA file for taxonomic assignation.
- `-c` Specify an identity percentage for the clustering (default: 0.97).

## Repository Structure
```
/home/florian/work/CRBE_STAGE2A/gitrep/
├── data/           # Input and output data
├── scripts/        # Core pipeline scripts
├── tests/          # Unit tests
├── README.md       # Project documentation
└── requirements.txt # Python dependencies
```

## Contact
For questions or support, please contact me at .