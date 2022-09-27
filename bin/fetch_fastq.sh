#!/usr/bin/sh

OUTDIR=$1
INFO=$2

module load sratoolkit

if [ ! -d "$OUTDIR" ]
then 
mkdir $OUTDIR
fi

while read SRR NAME
do 

    if [ ! -f "$OUTDIR/$NAME.fastq.gz" ]; then
        echo -e "\nDownloading $SRR.sra\n"
        prefetch $SRR

        fastq-dump --outdir $OUTDIR --gzip --split-files ~/ncbi/public/sra/$SRR.sra

        mv $OUTDIR/$SRR\_1.fastq.gz $OUTDIR/$NAME.fastq.gz
    else
        echo -e "$NAME\.fastq.gz already exists"
    fi

done < $INFO
