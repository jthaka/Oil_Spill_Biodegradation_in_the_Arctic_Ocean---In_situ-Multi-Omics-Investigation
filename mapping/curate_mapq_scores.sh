#This script will pull mapQ scores from your bam file into a text file.
#The .txt output can then be used in the plot_mapq_scores.R script in R.
#Flags should be edited to define the type of alignments you want to analyze.

module load samtools   # or: conda activate <env_name>
indir=path/to/bams
outdir=path/to/outdir
for bam in ${indir}/*.bam; do
  echo "== $bam =="
  samtools view -@ 3 -F 4 "$bam" | cut -f5 | sort -n | uniq -c
done > ${outdir}/mapq_raw_single.txt
