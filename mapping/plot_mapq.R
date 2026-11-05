#This script will plot the mapQ scores as fractions of total aligned reads.

library(tidyverse)
library(ggplot2)

### read in raw mapq file output firect from the bam files
lines <- readLines("path/to/mapq_scores.txt")
lines <- lines[nzchar(trimws(lines))]

is_hdr <- grepl("^== ", lines)
samples <- sub("^== (.*) ==$", "\\1", lines[is_hdr])
samples <- sub("\\.bam$", "", samples)
body <- trimws(lines[!is_hdr])


### create mapq long table from the line object
mapq_long <- tibble(
  sample = samples[cumsum(is_hdr)[!is_hdr]],
  count  = as.numeric(sub(" .*", "", body)),
  mapq   = as.integer(sub("^[0-9]+ +", "", body))
)


### pivot to one row per sample and one column per mapq value
mapq <- mapq_long |>
  arrange(mapq) |>                       # keeps columns in numeric order
  pivot_wider(names_from = mapq, values_from = count,
              names_prefix = "MAPQ_", values_fill = 0) |>
  arrange(sample)


### make line chart for all samples
mapq_long |>
  group_by(sample) |>
  mutate(prop = count / sum(count)) |>
  ggplot(aes(x = mapq, y = prop, group = sample, color = sample)) +
    geom_line(alpha = 0.4) +
  labs(title = "Line chart for ALL samples", 
       y = "Proportion of aligned reads with that score", 
       x = "MapQ Score (integer values)") +
    theme_bw()

### make histogram for one sample
mapq_filtered <- mapq_long |> 
  filter(sample == "insert sample name") |>
  mutate(prop = count / sum(count))
mapq_col <- ggplot(mapq_filtered, aes(x = mapq, y = prop)) +
  geom_col() +
  labs(title = "Bar chart for ONE sample", 
       y = "Proportion of aligned reads with that score", 
       x = "MapQ Score (integer values)")
print(mapq_col)
