#! /bin/bash
#$ -S /bin/bash
#$ -N samtools
#$ -q transabyss.q
#$ -l mem_token=2G,mem_free=2G,h_vmem=2G
#$ -V

echo "Please make sure the correct genome reference is used."
echo "The default reference is hg19!"
/home/rcorbett/aligners/samtools/samtools-0.1.17/samtools mpileup \
    -f /home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa \
    -l {{mpileup_positions}} \
    {{bam_file}} \
    > {{mpileup_output}}

if [ $? -eq 0 ]
    then
        touch {{mpileup_stamp}}
        echo "mpileup job finished successfully at:" `date`
else
    echo "ERROR: samtools mpileup did not finish correctly!" \
         >> {{log_file}}
    echo "Failed patient is: {{patient_status}} " \
         >> {{log_file}}
    exit
fi
/gsc/software/linux-x86_64/python-2.7.2/bin/python \
    /home/szong/python/PileupParser_CountSnpBases.py -f \
    -i {{mpileup_output}}    \
    -o {{mpileup_AFcounts}}
if [ $? -eq 0 ]
    then
    echo "mpileup and post-processing job finished successfully at:" `date`
fi


