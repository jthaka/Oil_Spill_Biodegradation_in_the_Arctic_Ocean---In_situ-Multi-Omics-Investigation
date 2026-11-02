#This script will extract and count the PRIMARY and PRIMARY PAIRED alignments from a bam file.
#Alignments were only kept if >10% of the read was aligned OR/AND if the read covered >50% of the gene.
#This script was used for estimating the abundance and expression levels of genes in metagenome
#and metatranscriptome samples, respectively.
#Metagenome and metatranscriptome sample reads were mapped against the metagenome coassembly.
#This script is written as an array job for a clusters from the Digital Research Alliance of Canada (DRAC).

#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=2
#SBATCH --job-name=bedtools_intersect
#SBATCH --array=1-10

set -euo pipefail

module load bedtools
module load samtools

pool='EΧP'
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /path/to/sample/names.txt)
bam=/path/to/bam/dir/${sample}.sorted.bam
gff=/path/to/reference/assembly/gene/calls/.gff3
fileprefix=${sample}_mtcomg   #Will be the prefix for output files
outdir=path/to/output/directory
mkdir -p "${outdir}"

#Double check that files exist
if [ ! -f "${bam}" ]
then
echo "BAM file not found"
exit
fi

if [ ! -f "${gff}" ]
then
echo "gff file not found"
exit
fi

#Count primary alignments and primary paired alignments.
#Also remove alignments that constitute <10% of the read length OR/AND <50% of the coding sequence in the gff.

#Count primary alignments
samtools view --threads $SLURM_CPUS_PER_TASK -b -F 0x100,0x4,0x800 ${bam}\
| bedtools intersect -a ${gff} -b stdin -c -F 0.10 -f 0.50 -e > ${outdir}/${fileprefix}_prim.tsv

#Count primary paired alignments (a subset of primary alignments where both mates map)
samtools view --threads $SLURM_CPUS_PER_TASK -b -f 0x2 -F 0x100,0x4 ${bam}\
| bedtools intersect -a ${gff} -b stdin -c -F 0.10 -f 0.50 -e > ${outdir}/${fileprefix}_pair.tsv

#Optional: Fix formatting of output to have a column with just the protein calling softwares unique ID
#This is a step specific to the output from the Bakta gene calling software. 
gawk -i inplace -F'\t' -v OFS='\t' \
'NR==1{printf "contig\tprogram\tregion\tsequence_type\tstart\tstop\textra\tstring\tprotein_id\tread_count\tbakta_id\n";next} {$(NF+1)= substr($9, 4, 15) ; print}' \
${outdir}/${fileprefix}_prim.tsv

gawk -i inplace -F'\t' -v OFS='\t' \
'NR==1{printf "contig\tprogram\tregion\tsequence_type\tstart\tstop\textra\tstring\tprotein_id\tread_count\tbakta_id\n";next} {$(NF+1)= substr($9, 4, 15) ; print}' \
${outdir}/${fileprefix}_pair.tsv

#Combine the two read count files with 2 read count fields based on unique gene ID from your gene calling software.
#You may have to join field number depending on the setup of your output file.
#If your files don't have headers, take out "tail -n +2".
join -t $'\t' -1 11 -2 11 -a 1 -o auto \
 <(tail -n +2 "${outdir}/${fileprefix}_prim.tsv" | sort -t $'\t' -k11,11) \
 <(tail -n +2 "${outdir}/${fileprefix}_pair.tsv" | sort -t $'\t' -k11,11) \
  > ${outdir}/${fileprefix}_readcounts.tsv

#Optional: Reapply the header and remove unecessary columns.
awk 'BEGIN{FS=OFS="\t"} NR==1{print "bakta_id","contig","program","sequence_type","start","stop","strand","protein_id","read_count_all_primary","read_count_paired_only"; next} {print $1,$2,$3,$4,$5,$6,$8,$10,$11,$21}' \
${outdir}/${fileprefix}_readcounts.tsv > ${outdir}/${fileprefix}_readcounts_trim.tsv

#Clean up intermediate files
mv ${outdir}/${fileprefix}_readcounts_trim.tsv ${outdir}/${fileprefix}_readcounts.tsv
rm ${outdir}/${fileprefix}_prim.tsv
rm ${outdir}/${fileprefix}_pair.tsv
