library(tidyverse)
library(edgeR)
library(dplyr)
library(stringr)
library(KEGGREST)
library(tidytext)
library(pheatmap)

counts_folder <- "~/Downloads/CMO/diffexp_JT/read_counts_mt_tomg/EXP/"

# merge read counts
#file_list <- list.files(counts_folder, pattern = "_readcounts\\.txt$")
file_list <- list.files(counts_folder, pattern = "readcounts")

master_counts <- NULL

# for (f in file_list) {
#   print(f)
#   sample_name <- str_split_fixed(f, "_readcounts", 2)[, 1]
#   tmp <- read.delim(paste0(counts_folder, f), comment.char = "#", stringsAsFactors = FALSE)
#   tmp_counts <- tmp[, c(1, 7)]
#   colnames(tmp_counts) <- c("Geneid", sample_name)
#   if (is.null(master_counts)) {
#     master_counts <- tmp_counts
#   } else {
#     master_counts <- full_join(master_counts, tmp_counts, by = "Geneid")
#   }
# }

for (f in file_list) {
  print(f)
  sample_name <- str_split_fixed(f, "_readcounts", 2)[, 1]
  sample_name <- str_split_fixed(sample_name, "short_version_", 2)[, 2]
  tmp <- read.delim(paste0(counts_folder, f), comment.char = "#", stringsAsFactors = FALSE)
  colnames(tmp) <- c("Geneid", sample_name)
  if (is.null(master_counts)) {
    master_counts <- tmp
  } else {
    master_counts <- full_join(master_counts, tmp, by = "Geneid")
  }
}

master_counts[is.na(master_counts)] <- 0
rownames(master_counts) <- master_counts$Geneid
master_counts$Geneid <- NULL
master_matrix <- as.matrix(master_counts)

# add metadata
# sample_meta <- read.csv("/Users/jordini/Downloads/CMO/sample_meta.csv") %>%
#   filter(pool == "EXP" & type == "si")
# sample_meta <- sample_meta %>%
#   mutate(
#     day_num = as.numeric(str_extract(time_since_oil, "\\d+")),
#     week_since_oil = case_when(
#       time_since_oil == "pre" ~ "pre",
#       !is.na(day_num)         ~ paste0("week", ceiling(day_num / 7)),
#       TRUE                    ~ NA_character_
#     )
#   ) %>%
#   select(-day_num)

# select meta data

sample_meta <- read.csv("/Users/jordini/Downloads/CMO/sample_meta.csv") %>%
  filter(pool == "EXP" & type == "si")
sample_meta <- sample_meta %>%
  mutate(day_num = as.numeric(str_extract(time_since_oil, "\\d+")),
    week_since_oil = case_when(
      time_since_oil == "pre" ~ "pre",
      !is.na(day_num) ~ paste0("week", ceiling(day_num / 7)),
      TRUE ~ NA_character_)) %>%
  select(-day_num)

# use this for filtering out singletons
# sample_meta <- sample_meta %>% add_count(week_since_oil, name = "week_count") %>%
#   filter(week_count > 1) %>%
#   select(-week_count)

# or this
#sample_meta <- sample_meta[,c("code","sample_location_no_replicate","time_since_oil","days_relative_to_oil_spill","pool","type","medium","week_since_oil")]
#sample_meta <- unique(sample_meta)

# combine with meta data

colnames(sample_meta) <- c("sample_mt_long","genohub_ID","sample_name","pool","type","time","ice_location","time_since_oil","week_since_oil")
sample_meta$mt_small_code_pt1 <- str_split_fixed(sample_meta$sample_mt_long, "_", 3)[,2]
sample_meta$mt_small_code_pt1 <- as.numeric(str_extract(sample_meta$mt_small_code_pt1, "\\d+"))
sample_meta$mt_small_code <- str_split_fixed(sample_meta$sample_mt_long, "_", 3)[,3]
sample_meta$sample <- paste0(sample_meta$mt_small_code_pt1, "_", sample_meta$mt_small_code, "_mtcomg_", sample_meta$pool)

sample_meta <- sample_meta %>% add_count(time_since_oil, name = "time_count") %>%
  filter(time_count > 1) %>%
  select(-time_count)

# sample_meta <- subset(sample_meta, sample_name != "CMO_EXP_DAY63")
# sample_meta <- sample_meta %>% add_count(week_since_oil, name = "week_count") %>%
#   filter(week_count > 1) %>%
#   select(-week_count)

pool_meta <- data.frame(sample = colnames(master_matrix)) %>%
  inner_join(sample_meta, by = "sample")
time_groups <- factor(pool_meta$time_since_oil)

master_matrix_exp <- master_matrix[, pool_meta$sample, drop = FALSE]

# start edgeR

d <- DGEList(counts = master_matrix_exp, group = time_groups)
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]
d <- normLibSizes(d)

# choose levels

d$samples$group <- factor(d$samples$group, levels = c("day02", "day03", "day05", "day07", "day09", "day11", "day50"))
#d$samples$group <- factor(d$samples$group, levels = c("week1","week2","week8"))
#d$samples$group <- factor(d$samples$group, levels = c("day01","day03","day06","day11","day15","day49"))
#d$samples$group <- factor(d$samples$group, levels = c("pre","week1","week2"))
num_groups <- length(levels(d$samples$group))

# check data using MDS

rainbow_palette <- hcl(h = seq(15, 375, length.out = num_groups + 1)[1:num_groups], c = 100, l = 50)

#png("mt_ctl_exp_mds_plot.png", width=1050, height=800, res=120)
set.seed(666)
par(mar=c(5.1, 4.1, 4.1, 10.1))
plotMDS(d, method="logFC", 
        col=rainbow_palette[as.numeric(d$samples$group)],
        xlim=c(-3.5, 10))
legend("topleft", inset=c(1.02, 0), 
       legend=levels(d$samples$group), 
       col=rainbow_palette, 
       pch=20, xpd=TRUE, bty="n")

dev.off()

# contine edgeR

design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d <- estimateDisp(d, design.mat, tagwise=T)
fit <- glmQLFit(d, design.mat)

# choose baseline to compare other days/weeks to

group_names  <- colnames(design.mat)
baseline_day <- "day02"
#baseline_day <- "week1"
#baseline_day <- "pre"
other_days   <- group_names[group_names != baseline_day]

master_edgeR_table <- tibble()

for (target_day in other_days) {
  comparison_name <- paste0(target_day, "_vs_", baseline_day)
  
  my_contrast <- makeContrasts(contrasts = paste0(target_day, "-", baseline_day), levels = design.mat)
  qlf <- glmQLFTest(fit, contrast = my_contrast[, 1])
  
  df_data <- topTags(qlf, n = Inf)$table
  df_data$Geneid <- rownames(df_data)
  
  df_data$Comparison <- comparison_name
  df_data$Baseline <- baseline_day
  df_data$TargetDay <- target_day
  
  master_edgeR_table <- bind_rows(master_edgeR_table, df_data)
}

master_edgeR_table <- master_edgeR_table %>%
  mutate(Global_FDR = p.adjust(PValue, method = "BH")) %>%
  mutate(Significance = case_when(
    Global_FDR < 0.05 & logFC > 2 ~ "Up-regulated",
    Global_FDR < 0.05 & logFC < -2 ~ "Down-regulated",
    TRUE ~ "Not Significant"))

all_canthyd_genes <- tibble()
lines <- read_lines("~/Downloads/CMO/diffexp_JT/exp_danaseq/canthyd_results.tblout")
data_lines <- lines[!grepl("^#", lines)]

if (length(data_lines) > 0) {
  all_canthyd_genes <- tibble(text = data_lines) %>%
    separate_wider_regex(text, 
      patterns = c(target_name = "^\\S+", "\\s+", target_accession = "\\S+", "\\s+", query_name = "\\S+"),
      too_few = "align_start") %>%
    filter(!is.na(target_name) & !is.na(query_name)) %>%
    select(Geneid = target_name, CANTHYD_Query_Name = query_name)
}

all_canthyd_genes_clean <- all_canthyd_genes %>%
  distinct(Geneid, .keep_all = TRUE)

annotated_canthyd <- master_edgeR_table %>%
  inner_join(all_canthyd_genes_clean, by = "Geneid")

#write.csv(annotated_canthyd, file = "~/Downloads/CMO/edgeR/exp/sw/canthyd-lenient-logfc-all.csv", row.names = FALSE)
