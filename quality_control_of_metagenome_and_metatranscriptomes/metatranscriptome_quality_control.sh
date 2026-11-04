#___________________________

# trimmomatic

#___________________________

module load trimmomatic
for R1 in *_1.fq.gz ; do
  R2="${R1%%_1.*}_2.fq.gz"
  java -jar ../../20250826-CMO-tx/clean_reads/databases/trimmomatic.jar PE -threads 8 \
  "$R1" "$R2" \
  "${R1%%.*}_out_paired.fq.gz" "${R1%%.*}_out_unpaired.fq.gz" \
  "${R2%%.*}_out_paired.fq.gz" "${R2%%.*}_out_unpaired.fq.gz" \
  ILLUMINACLIP:../../20250826-CMO-tx/clean_reads/databases/TruSeq3-PE.fa:2:30:10:8:True SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:36
done

#___________________________

# rcorrector

#___________________________

module load apptainer
module load jellyfish
module load perl
apptainer exec -B /cvmfs/soft.computecanada.ca -B /home/jwinter/scratch/ -B /cvmfs --env PATH="${JELLY_PATH}:/usr/bin:/bin:/usr/local/bin" /home/jwinter/scratch/20250826-CMO-tx/clean_reads/rcorrector_all/rcorrector_latest.sif perl /rcorrector/run_rcorrector.pl -1 ../trim_out/paired/${sample}_1_out_paired.fq.gz -2 ../trim_out/paired/${sample}_2_out_paired.fq.gz -t 10 -od rcorrector_out

#___________________________

# deduplication pipeline

#___________________________

module load StdEnv/2020
module load bbmap

# Making all the directories we will use
mkdir -p ./clean_fastq/
mkdir -p ./clean_fastq/normalized/

# Initialize read count log
readcount_log=./clean_fastq/readcounts_${sample}.txt
echo "Step,Forward_reads,Reverse_reads" > $readcount_log

#Set the “current” clean files to the initial phix output
current_1=./clean_fastq/${sample}_clean_1.fastq.gz
current_2=./clean_fastq/${sample}_clean_2.fastq.gz
norm_1=./clean_fastq/normalized/${sample}_norm_1.fastq.gz
norm_2=./clean_fastq/normalized/${sample}_norm_2.fastq.gz

# Remove phix contamination with bbduk
tmp_1=./clean_fastq/${sample}_tmp_1.fastq.gz
tmp_2=./clean_fastq/${sample}_tmp_2.fastq.gz
bbduk.sh \
in1=./rcorrector/rcorrector_out/${sample}_1_out_paired.cor.fq.gz \
in2=./rcorrector/rcorrector_out/${sample}_2_out_paired.cor.fq.gz \
out=$tmp_1 \
out2=$tmp_2 \
ref=$(dirname $(which bbduk.sh))/resources/phix174_ill.ref.fa.gz \
k=31 \
hdist=1 \
stats=./CMO2024_fastq_cleaning_temp/${sample}_phix_stats.txt \
ow=t \
threads=$SLURM_CPUS_PER_TASK \
-Xmx40g
mv $tmp_1 $current_1
mv $tmp_2 $current_2

# Log read counts after phix removal
echo "phix_removal,$(zcat $current_1 | wc -l | awk '{print $1/4}'),$(zcat $current_2 | wc -l | awk '{print $1/4}')" >> $readcount_log

# Remove human contamination with bbduk
tmp_1=./clean_fastq/${sample}_tmp_1.fastq.gz
tmp_2=./clean_fastq/${sample}_tmp_2.fastq.gz
removehuman.sh \
path=../../20250826-CMO-tx/clean_reads/databases/hg19_index/ \
in1=$current_1 \
in2=$current_2 \
outu1=$tmp_1 \
outu2=$tmp_2 \
outm1=./human_contam/${sample}_humancon_1.fastq.gz \
outm2=./human_contam/${sample}_humancon_2.fastq.gz \
ow=t \
threads=$SLURM_CPUS_PER_TASK \
-Xmx40g
mv $tmp_1 $current_1
mv $tmp_2 $current_2

# Log read counts after human removal
echo "human_removal,$(zcat $current_1 | wc -l | awk '{print $1/4}'),$(zcat $current_2 | wc -l | awk '{print $1/4}')" >> $readcount_log

#Deduplication with Clumpify.
#This will look at the Illumina flow cell information in the heading of each read and identify
#poor-quality areas of the flowcell and remove replicate reads that arise from poor signals. ecc=t (Error-correct reads. Requires multiple passes for complete correction).
#passes=[int] (Use this many error-correction passes. 6 passes are suggested, though more will be more through).
#optical=t (If true, mark or remove optical duplicates only. This means they are Illumina reads within a certain distance on the flowcell. Normal Illumina names needed.
#dedupe=t (Remove duplicate reads. For pairs, both must match. By default, deduplication does not occur.)
#dupedist=[int] (Max distance to consider for optical duplicates. Higher removes more duplicates but is more likely to remove PCR rather than optical duplicates.
#This is platform-specific; recommendations: NextSeq=40 (and spany=t); HiSeq-1T=40; HiSeq-2500=40; HiSeq-3k/4k=2500; Novaseq6000=12000; NovaseqX+=50)
#Transcripts are NextSeq so using 40
#Also for tile-edge and well duplicates).

clumpify.sh \
in=$current_1 \
in2=$current_2 \
out=$norm_1 \
out2=$norm_2 \
ecc=t passes=4 dedupe=t dupedist=40 optical=f \
t=$SLURM_CPUS_PER_TASK

# Log read counts after clumpify
echo "clumpify,$(zcat $current_1 | wc -l | awk '{print $1/4}'),$(zcat $current_2 | wc -l | awk '{print $1/4}')" >> $readcount_log

#Normalization with BBNorm.
#This tool is useful for samples with too much data or uneven coverage (including metagenomes).
#This is helpful for things like high copy number plasmids or just really high abundant organisms. The result is minimizing the redundancy in decision making for assemblers. NOTE: we will be mapping NON-normalized reads to our assemblies to generate differential coverage information. prefilter=t (True is slower, but generally more accurate; filters out low-depth kmers from the main hashtable. The prefilter is more memory-efficient because it uses 2-bit cells). mindepth=[int] (Kmers with depth below this number will not be included when calculating the depth of a read). target=[int] (Target normalization depth. NOTE: All depth parameters control kmer depth, not read depth. For kmer depth Dk, read depth Dr, read length R, and kmer size K: Dr=Dk*(R/(R-K+1)). WE ARE NOT OVERWRITING THE CLEAN FILES.

bbnorm.sh \
in=$norm_1 \
in2=$norm_2 \
out=$tmp_1 \
out2=$tmp_2 \
target=100 mindepth=2 \
prefilter=t \
hist=./clean_fastq/normalized/${sample}_bbnorm_hitin.txt \
histout=./clean_fastq/normalized/${sample}_bbnorm_hitout.txt \
peaks=./clean_fastq/normalized/${sample}_bbnorm_peaks.txt \
ow=t \
t=$SLURM_CPUS_PER_TASK \
-Xmx40g

mv $tmp_1 $norm_1
mv $tmp_2 $norm_2

# Log read counts after normalization
echo "bbnorm,$(zcat $norm_1 | wc -l | awk '{print $1/4}'),$(zcat $norm_2 | wc -l | awk '{print $1/4}')" >> $readcount_log
