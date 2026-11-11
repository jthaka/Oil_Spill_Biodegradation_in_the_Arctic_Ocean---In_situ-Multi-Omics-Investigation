#This script will classify taxonomy of gene calls based on the most popular classification
#for all genes on one contig.
#CAT/BAT was used as inspiration for this task: https://github.com/lbl-cbg/CAT

taxa <- read.csv("path/to/gene/level/tax.tsv", sep = "\t")

meta <- read_tsv("path/to/gene/level/annotations.tsv", comment = "#", col_names = FALSE, show_col_types = FALSE)
colnames(meta) <- c("Sequence_Id", "Type", "Start", "Stop", "Strand", "Locus", "Tag", "Gene", "Product")
meta$node <- str_split_fixed(meta$Sequence_Id, "_", 3)[,2]
meta <- meta[,c("node", "Locus")]
colnames(meta) <- c("node", "Query_ID")

taxa <- merge(taxa, meta)

taxa$TaxID <- as.character(taxa$TaxID)

catbat <- taxa %>% group_by(node) %>% count(TaxID) %>% dplyr::slice(which.max(n))

meta_names <- unique(taxa[,c("TaxID", "Lineage")])

catbat <- merge(catbat, meta_names)
colnames(catbat) <- c("Overall_TaxID", "node", "n", "Overall_Lineage")
taxa <- merge(taxa, catbat)

write.csv(taxa, "~/Downloads/CMO/mg_taxonomy_24/mg_catbat_co_exp.csv")
