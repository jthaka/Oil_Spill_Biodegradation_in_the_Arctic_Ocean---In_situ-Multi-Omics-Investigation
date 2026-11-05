#This script will assign taxonomy at the gene-level
#Genes on contigs were classified as of eukaryotic versus 
#prokaryotic origin using Tiara v1.0.3 (Karlicki et al., 2022). 
#Eukaryotic genes were taxonomically annotated with the MarFERReT database 
#with DIAMOND (Buchfink et al., 2015; Groussman et al., 2023), 
#while prokaryotic genes were annotated with MarMicroDB with Kaiju (Hogle, 2019; Menzel et al., 2016). 
#Both MarFERReT and MarMicroDB were chosen as these databases are specific to marine metagenomes and 
#include some representation from Arctic environments. 

#SBATCH --time=02:00:00
#SBATCH --array=1-66%30
#SBATCH --cpus-per-task=4
#SBATCH --mem=62G

sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samplenames.txt)

module load apptainer
module load diamond
module load kaiju

OUTPUT_DIR="taxonomy_single_mt/${sample}"
mkdir -p "$OUTPUT_DIR"

# run tiara

apptainer run -C -B "$PWD:$PWD" ../20250826-CMO-tx/db_temp/tiara/tiara_1.0.2.sif tiara \
    -t "$SLURM_CPUS_PER_TASK" \
    -i "$PWD/rnaspades/spades_single_mt/${sample}/transcripts.fasta" \
    -o "$PWD/${OUTPUT_DIR}/tiara_result.txt" \
    -m 1000 \
    -p 0.65 0.65

# run diamond

diamond blastp \
    --threads "$SLURM_CPUS_PER_TASK" \
    --db ../20250826-CMO-tx/db_temp/marferret_db/MarFERReT.v1.1.1.dmnd \
    --query ./danaseq_single_mt/${sample}_single_mt/annotation/bakta/extra/annotation.faa \
    --no-self-hits \
    --outfmt '102' \
    --out "${OUTPUT_DIR}/marferret_all.tsv" \
    --include-lineage \
    --compress '0' \
    --matrix 'BLOSUM62' \
    --comp-based-stats '1' \
    --masking 'tantan' \
    --max-target-seqs '25' \
    --evalue '0.001' \
    --id '0.0' \
    --approx-id '0.0' \
    --query-cover '0.0' \
    --subject-cover '0.0' \
    --block-size '2.0' \
    --motif-masking 0 \
    --soft-masking 0 \
    --algo 0 \
    --file-buffer-size "${DIAMOND_FILE_BUFFER_SIZE:-67108864}"

# run kaiju

kaiju \
    -t ../20250826-CMO-tx/db_temp/marmicrodb_db/nodes.dmp \
    -z "$SLURM_CPUS_PER_TASK" \
    -f ../20250826-CMO-tx/db_temp/marmicrodb_db/MARMICRODB.fmi \
    -i ./danaseq_single_mt/${sample}_single_mt/annotation/bakta/extra/annotation.faa \
    -o "${OUTPUT_DIR}/${sample}.kaiju.out" \
    -p \
    -v

kaiju-addTaxonNames \
    -t ../20250826-CMO-tx/db_temp/marmicrodb_db/nodes.dmp \
    -n ../20250826-CMO-tx/db_temp/marmicrodb_db/names.dmp \
    -i "${OUTPUT_DIR}/${sample}.kaiju.out" \
    -r superkingdom,phylum,class,order,family,genus,species | \
    awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2, $3, $NF}' > "${OUTPUT_DIR}/${sample}_all_levels_names.tsv"

OUTPUT_DIR="taxonomy_single_mt/${sample}"

FINAL_TAXONOMY="${OUTPUT_DIR}/${sample}_final_taxonomy.tsv"
BAKTA_TSV="./danaseq_single_mt/${sample}_single_mt/annotation/bakta/extra/annotation.tsv"

echo -e "Query_ID\tTaxID\tLineage\tSource_DB" > "$FINAL_TAXONOMY"

awk -F'\t' '
    NR==FNR { if(NR>1) tiara[$1]=$2; next }
    /^#Sequence Id/ {
        sub(/^#/, "", $1)
        for(i=1; i<=NF; i++) {
            if($i=="Sequence Id") seq_col=i
            if($i=="Locus Tag") locus_col=i
        }
        next
    }
    /^#/ { next }
    $seq_col in tiara { print $locus_col, tiara[$seq_col] }
' "${OUTPUT_DIR}/tiara_result.txt" "$BAKTA_TSV" > "${OUTPUT_DIR}/locus_map.txt"

awk '$2=="eukarya" || $2=="organelle" {print $1}' "${OUTPUT_DIR}/locus_map.txt" > "${OUTPUT_DIR}/euk_ids.txt"
awk '$2=="bacteria" || $2=="archaea" || $2=="prokarya" {print $1}' "${OUTPUT_DIR}/locus_map.txt" > "${OUTPUT_DIR}/prok_ids.txt"

# filter marferret
RAW_MARFERRET="${OUTPUT_DIR}/marferret_all.tsv" 

if [ -s "${OUTPUT_DIR}/euk_ids.txt" ] && [ -s "$RAW_MARFERRET" ]; then
    awk -F'\t' 'BEGIN{OFS="\t"}
        NR==FNR { ids[$1]; next }                     
        $1 in ids && $2!="0" && $4!="N/A" {         
            print $1, $2, $4, "MarFERReT" 
        }
    ' "${OUTPUT_DIR}/euk_ids.txt" "$RAW_MARFERRET" >> "$FINAL_TAXONOMY"
fi

# subset kaiju
RAW_KAIJU="${OUTPUT_DIR}/${sample}_all_levels_names.tsv" 

if [ -s "${OUTPUT_DIR}/prok_ids.txt" ] && [ -s "$RAW_KAIJU" ]; then
    awk -F'\t' 'BEGIN{OFS="\t"}
        NR==FNR { ids[$1]; next }                     
        $2 in ids && $3!="0" {                        
            lineage = $4
            sub(/;[ \t]*$/, "", lineage)
            print $2, $3, lineage, "MarMicroDB"
        }
    ' "${OUTPUT_DIR}/prok_ids.txt" "$RAW_KAIJU" >> "$FINAL_TAXONOMY"
fi

rm "${OUTPUT_DIR}/locus_map.txt" "${OUTPUT_DIR}/euk_ids.txt" "${OUTPUT_DIR}/prok_ids.txt"
