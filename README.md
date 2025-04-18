# Pipeline Repository

This repository contains the code and resources used during my internship at CRBE. The pipeline is designed to process and analyze fungi metabarcoding data.

## Process
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

## Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/yourusername/your-repo.git
    ```
2. Navigate to the directory:
    ```bash
    cd CRBE_STAGE2A
    ```
3. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```

## Usage
You can either run the whole script at once, or each of the steps individually.

- Running the whole script:

- Running each step individually:
-

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