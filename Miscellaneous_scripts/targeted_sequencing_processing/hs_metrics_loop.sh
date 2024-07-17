SAMPLES=$(cat fail_samples.txt)
PROJECT=2737
# Found IDs

for SAMPLE in $SAMPLES; do
                bsub -o $PWD/log.%J -e $PWD/err.%J -q normal -R'select[mem>1000] rusage[mem=1000]' -M1000 -J $SAMPLE ./hs_metrics.sh $SAMPLE $PROJECT
done
