#!/bin/bash

#

while getopts "i:f:r:" option
do
        case $option in
                i)
                    INPUT="$OPTARG"
                    ;;
                f)
                    PRIMER_F="$OPTARG"
                    ;;
                r)
                    PRIMER_R="$OPTARG"
                    ;;
        esac
done

./00_prepare_DB.sh

./01_run_cutadapt.sh

./02_run_dada2.sh

Rscript -e "rmarkdown::render('$DIR/../report.rmd', params = list(dir='$DIR', output='$OUTPUT', depth='$DEPTH'), output_dir = '$DIR/', output_file = 'report_$OUTPUT.html')"
