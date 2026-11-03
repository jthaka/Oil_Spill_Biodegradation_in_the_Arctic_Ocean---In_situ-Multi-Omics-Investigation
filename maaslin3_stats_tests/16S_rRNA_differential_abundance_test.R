#This R script was used to perform differential abundance testing on 16S rRNA amplicon data
#From the CMO experiments. Time is accounted for in this stats test and the goal is to
#see if there are any significant differences in the abundance and/or prevelence of ASVs
#with the diesel treatment, while accounting for other variables (i.e., seawater vs. sea-ice, along with time).

#Load packages
source("C:/Users/thaka/OneDrive/Desktop/required_micro_packages.R")
if (!require('openxlsx')) install.packages('openxlsx'); library('openxlsx')
if (!require('magrittr')) install.packages('magrittr'); library('magrittr')
if (!require('microeco')) install.packages('microeco'); library('microeco')
if (!require('picante')) install.packages('picante'); library('picante')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('tibble')) install.packages('tibble'); library('tibble')
if (!require('tidyr')) install.packages('tidyr'); library('tidyr')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('fuzzyjoin')) install.packages('fuzzyjoin'); library('fuzzyjoin')
if (!require('Rsubread')) BiocManager::install("Rsubread"); library('Rsubread')
if (!require('reshape')) install.packages("reshape"); library('reshape')
if (!require('scales')) install.packages("scales"); library('scales')
if (!require('stringr')) install.packages("stringr"); library('stringr')
if (!require('svglite')) install.packages("svglite"); library('svglite')
if (!require('readr')) install.packages("readr"); library('readr')
if (!require('maaslin3')) install.packages("maaslin3"); library('maaslin3')

#Load phyloseq object
ps=load("path/to/phyloseq")

#Filter out rare ASVs
gGlommed= phyloseq::prune_taxa(phyloseq::taxa_sums(ps) > 0, ps)
gGlommed= phyloseq::prune_samples(phyloseq::sample_sums(gGlommed) > 0, gGlommed)
gGlommed=gGlommed %>% tax_filter(min_prevalence = 0.1) #Cut taxa present in less than 10% of samples
gGlommed=tax_transform(gGlommed,rank="Genus",trans="identity")
gotu=as(otu_table(gGlommed), "matrix")
gotu=data.frame(gotu)

#Load metadata
meta=read.csv("path/to/metadata.csv")
meta_df <- data.frame(meta)

#Factor groups for maaslin3 test (for this study, the variables are: pool and medium)
meta_df$pool=factor(meta_df$pool, levels = c('CTL', 'EXP'))
meta_df$medium_simple=factor(meta_df$medium_simple, levels = c('seawater', 'ice'))
meta_df=column_to_rownames(meta_df, var="genohub") #Ensure rownames match the OTU table from your phyloseq

# Fit models
fit_out_w_medium <- maaslin3(input_data = gotu,
                    input_metadata = meta_df,
                    output = 'maaso_w_medium_pool_as_grouping',
                    formula = 'group(pool) + ordered(time)',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = FALSE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    save_models = TRUE)
