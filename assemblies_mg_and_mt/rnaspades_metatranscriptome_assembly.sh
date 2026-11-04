#rnaSPAdes was used to assembly metatranscriptome samples from seawater and sea-ice. Both single assemblies and coassemblies
#were performed.

#!/bin/bash
## Job Name
#SBATCH --time=05:00:00
#SBATCH --array=1-66%30
#SBATCH --mem=250G
#SBATCH --cpus-per-task=32

#single assembly

module load spades
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../samplenames.txt)
spades.py --rna -1 ../clean_fastq/normalized/${sample}_norm_1.fastq.gz -2 ../clean_fastq/normalized/${sample}_norm_2.fastq.gz -o spades_single_mt/${sample} -m 250 -t 32

#coassembly

spades.py --rna -1 ../clean_fastq/exp_pooled_samples_R1.fq.gz -2 ../clean_fastq/exp_pooled_samples_R2.fq.gz -o spades_co_mt/exp -m 250 -t 32
