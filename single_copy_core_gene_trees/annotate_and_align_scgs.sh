#Goal: annotate single copy core genes (scgs) and construct a maximum-likelihood tree with IQ-tree/model-finder
#Note: genomes/MAGs that contain >1 copy of an scg are removed from the alignment of that particular scg and
#are representated as gaps in the concatenated alignment for that particular region of the scg.
#Inputs: genome fasta files & gene calls
#Outputs: partitioned maximum-likelihood tree
