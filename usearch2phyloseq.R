#!/usr/bin/R

library("tidyverse")
library("ape")
library("phyloseq")

args <- commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Two argument must be supplied!!!", call.=FALSE)
} else if (length(args)==1) {
  grpbase <- args[2]
  sample_file <- args[1]
}

### read uc to otu_table
otutab_file <- paste0(grpbase, "_otutab.txt")
otutab <- read_tsv(otutab_file) %>%
  column_to_rownames("#OTU ID") %>%
  as.matrix()
print("finish step 1")

### samtab
samdf <- read.csv(sample_file, sep="\t", header=TRUE)
rownames(samdf) <- samdf$SampleID
#SampleID   FecusID MouseID PrimerF PrimerR BarcodeF    BarcodeR    Feed    TreatTime   SamplingDate    BatchPlate
#S0001  F0001   M001    ACTCCTACGGGNGGCWGCAG    GGACTACHVGGGTWTCTAAT    CGATGT  ACTGAT  common_control  0M  2016.12.29  P1
#S0002  F0002   M002    ACTCCTACGGGNGGCWGCAG    GGACTACHVGGGTWTCTAAT    CGATGT  ATGAGC  common_control  0M  2016.12.29  P1
meta.cols <- c("FecusID", "MouseID", "TreatTime", "SamplingDate", "BatchPlate")
samtab <- samdf[colnames(otutab), meta.cols]
print("finish step 2")

### taxonomy
sintax_file <- paste0(grpbase, "_otus.sintax")
taxonomy.cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
taxtab <- read_tsv(sintax_file, col_names=c("OTUID", "TAXPLIST", "SIMBOL", "TAXLIST")) %>%
  separate(col = TAXLIST, taxonomy.cols, ",") %>%
  filter(OTUID %in% rownames(otutab)) %>%
  column_to_rownames("OTUID") %>%
  select(one_of(taxonomy.cols)) %>%
  as.matrix()

print("finish step 3")

### phylogeny
tree_file <- paste0(grpbase, "_otus.nwk")
tree <- read.tree(tree_file)
print("finish step 4")

#### The full suite of data for this study – the sample-by-sequence feature table, the sample metadata, sequence taxonomies, and the phylogenetic tree – can now be combined into a single object.
ps <- phyloseq(tax_table(taxtab),
               sample_data(samtab),
               otu_table(otutab, taxa_are_rows = TRUE),
               phy_tree(tree))
saveRDS(ps, file=paste0(grpbase, "_ps.rds"))
print("finish all...")
