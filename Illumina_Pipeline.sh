#!/bin/bash

#
DIR=$(pwd)
OUTPUT=OUTPUT

while getopts "hi:o:f:r:m:d:c:" option
do
        case $option in
                h)
                    echo "help"
                    exit 0
                    ;;
                i)
                    INPUT="$OPTARG"
                    ;;
                o)
                    OUTPUT="$OPTARG"
                    ;;
                f)
                    PRIMER_F="$OPTARG"
                    ;;
                r)
                    PRIMER_R="$OPTARG"
                    ;;
                m)
                    MAPPING="$OPTARG"
                    ;;
                d)
                    DATABASE="$OPTARG"
                    ;;
                c)
                    IDENTITY_PERCENTAGE="$OPTARG"
                    ;;

        esac
done

mkdir $OUTPUT

./$0/bash_scripts/00_prepare_DB.sh -i $DATABASE -f $PRIMER_F -r $PRIMER_R

./$0/bash_scripts/01_merging_demultiplexing.sh -i $INPUT -o $OUTPUT -m $MAPPING

./$0/bash_scripts/02_OTU_clustering.sh -c $IDENTITY_PERCENTAGE -o $OUTPUT

Rscript -e "rmarkdown::render('$0/R_scripts/report.rmd', params = list(dir='$DIR', output='$OUTPUT', depth='$DEPTH'), output_dir = '$DIR/', output_file = 'report_$OUTPUT.html')"
