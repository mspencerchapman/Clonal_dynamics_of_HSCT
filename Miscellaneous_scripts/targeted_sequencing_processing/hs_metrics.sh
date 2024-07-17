SAMPLE=$1
PROJECT=$2
DATA=/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/targeted_sequencing/baitset1/bams
OUT=/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/targeted_sequencing/baitset1/metrics_files
REF=/lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/genome.fa
java -jar $PICARD CollectHsMetrics \
	I=$DATA/$PROJECT/$SAMPLE/mapped_sample/$SAMPLE.sample.dupmarked.bam \
	O=$OUT/$SAMPLE.hs_metrics.txt \
	R=$REF \
	BAIT_INTERVALS=/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/targeted_sequencing/baitset1/Covered_bed_files/Covered_bed_files_combined.interval_list \
	TARGET_INTERVALS=/lustre/scratch119/casm/team154pc/ms56/Zur_HSCT/targeted_sequencing/baitset1/Target_files/Targets_combined.interval_list
