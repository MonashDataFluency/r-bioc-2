
library(Biostrings)      # Provides DNAString, DNAStringSet, etc
library(BSgenome)        # Provides getSeq()
library(GenomicRanges)   # Provides GRanges containing genomic ranges, etc
library(GenomicFeatures) # Provides TxDb objects containing genes/transcripts/exons
library(rtracklayer)     # Provides import() and export()
library(Gviz)            # Provides plotting of genomic features

readDNAStringSet("Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz")

sessionInfo()

# [Home](index.html)
