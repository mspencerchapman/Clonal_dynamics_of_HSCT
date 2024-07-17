PROJECT=$1
SAMPLES=$(ls /nfs/cancer_ref01/nst_links/live/${PROJECT})

for SAMPLE in $SAMPLES; do
    echo $SAMPLE
    alleleCounter \
    -b /nfs/cancer_ref01/nst_links/live/${PROJECT}/${SAMPLE}/${SAMPLE}.sample.dupmarked.bam \
    -l /lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/targeted_sequencing/baitset1/baitset1_SNVs.bed \
    -o allele_counter/${SAMPLE}.allelecounts.txt
done

