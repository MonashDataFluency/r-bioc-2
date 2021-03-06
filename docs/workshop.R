


#////////////////////////////////////////
# 1 Installing Bioconductor packages ----
#
# Bioconductor packages are installed using the CRAN package
# "BiocManager".

## 1. Install BiocManager to your hard disk:
#
# install.packages("BiocManager")
#
## 2. Use BiocManager to install a BioConductor package to your hard disk:
#
# BiocManager::install("Biostrings")
#
## 3. Load the package into memory in your current R session:
#
# library(Biostrings)
#

# To update any out-of-date packages use install() with no arguments.
# Restart your R session after doing this!

# BiocManager::install()

# Bioconductor packages usually have useful documentation in the form of
# "vignettes". These are readable on the Bioconductor website
# https://bioconductor.org/packages/release/BiocViews.html#___Software ,
# from the RStudio "help" pane, or within R:

help(package="Biostrings")
# -> select "User guides, package vignettes, and other documentation"

vignette("BiostringsQuickOverview", package="Biostrings")



#///////////////////////////////////////////
# 2 Main packages used in this workshop ----

library(Biostrings)      # Provides DNAString, DNAStringSet, etc
library(BSgenome)        # Provides getSeq()
library(GenomicRanges)   # Provides GRanges containing genomic ranges
library(rtracklayer)     # Provides import() and export()
library(Gviz)            # Provides plotting of genomic features



#/////////////////
# 3 Sequences ----

# 3.1 DNAString ----
#
# Package Biostrings offers classes for storing DNA strings, DNAString,
# amino acid sequences, AAString, or anything else in a BString. These
# are like character strings, but a variety of biologically meaningful
# functions can be applied to them.

myseq <- DNAString("CCGCGCACCAAC")
myseq

class(myseq)

reverseComplement(myseq)
translate(myseq)

subseq(myseq, 3,5)
myseq[3:5]

as.character(myseq)

# You can see a complete set of functions that work with DNAString with:

methods(class="DNAString")

# You can see internal details of a class with:

getClass("DNAString")

# The "slots", used to store data internally in these objects, are
# accessible with, eg, myseq@length. However, wherever possible, we
# should leave the slots alone and use "methods" to access information
# in the object, eg length(myseq)! The classes that this class "extends"
# tell us some further ways we can interact with it.
#
# You can get help on the DNAString class as below. Try this with
# classes it extends too!

?"DNAString-class"

# 3.2 DNAStringSet ----
#
# Often we want to work with a list of sequences, such as chromosomes.

myset <- DNAStringSet( list(chrI=myseq, chrII=DNAString("ACGTACGT")) )
myset

# A DNAStringSet is list-like
myset$chrII
# or myset[["chrII"]]
# or myset[[2]]

# 3.3 Challenge: sequences ----
#
# 1. Reverse complement the following DNA sequence. Translate the result
# to an amino acid sequence.

"GTAGAGTAATATGGA"

# 2. Also translate bases 3 to 11 of the reverse complement.
#


#//////////////////////
# 4 Genomic ranges ----
#
# We may then wish to refer to regions of these sequences, often with an
# associated strand. This is done with the GRanges type. GRanges builds
# on IRanges, "integer ranges". An IRanges has starts and ends. A
# GRanges additionally has sequence names and strand information.

range1 <- GRanges("chrI", IRanges(start=3,end=5), "+")
range1
getSeq(myset, range1)

range2 <- GRanges("chrI", IRanges(start=3,end=5), "-")
getSeq(myset, range2)

# Accessing GRanges data:

seqnames(range1)
start(range1)
end(range1)
width(range1)
strand(range1)
as.data.frame(range1)

# Internals of an "S4" object such as a GRanges can be accessed using @,
# but this is discouraged. It is better to use the accessor functions
# above. Observe the completions RStudio offers when you type range1@.

# Look at completions for
# range1@

# Further manipulations:

# GRanges are like vectors:
c(range1, range2)

# GRanges can have metadata columns, and are often used like data frames:
mcols(range1)$wobble <- 10
range1
mcols(range1)
range1$wobble

# A handy way to create a GRanges
as("chrI:3-5:+", "GRanges")

# 4.1 Question ----
#
# Based on what we saw for DNAString, how can we learn more about using
# GRanges and IRanges objects?
#
# 4.2 Visualization ----
#
# A couple of options for showing ranges are the GViz and ggbio
# packages. ggbio is a ggplot2 extension. ggbio is a great idea, but (as
# at late 2021) the package we could use some love to polish and update
# it. GViz takes a grid graphics approach.
#
# Here we'll use the GViz package, which is somewhat confusing to use
# but produces nice graphics out of the box.

options(ucscChromosomeNames=FALSE)
look <- function(features, labels=seq_along(features)) {
    axis_track <- GenomeAxisTrack()
    feature_track <- AnnotationTrack(
        features, id=labels, name="",
        showFeatureId=TRUE, fontcolor.feature="black")
    plotTracks(
        list(axis_track, feature_track),
        extend.left=1, extend.right=2)
}

look( c(range1, range2) )



#/////////////////////
# 5 Loading files ----

# 5.1 Loading sequences ----
#
# DNA sequences are generally stored in FASTA format, a simple text
# format. These can be loaded with readDNAStringSet from Biostrings.
# Let's load the genome of the *Caeonrhabditis elegans* worm, obtained
# from the Ensembl FTP site. This is the "soft masked" (dna_sm) version,
# with repetitive regions given in lower case, but we won't be making
# use of this today.

### After "gunzip"ing, the start of the .fa file looks like this:
# >I dna_sm:chromosome chromosome:WBcel235:I:1:15072434:1 REF
# gcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaa
# gcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaa
# gcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaagcctaa
# ...
# ttttttagaaaaattatttttaagaatttttcattttAGGAATATTGTTATTTCAGAAAA
# TAGCTAAATGTGATTTCTGTAATTTTGCCTGCCAAATTCGTGAAATGCAATAAAAATCTA
# ATATCCCTCATCAGTGCGATTTCCGAATCAGTATATTTTTACGTAATAGCTTCTTTGACA
# TCAATAAGTATTTGCCTATATGACTTTAGACTTGAAATTGGCTATTAATGCCAATTTCAT
# ...

seqs <- readDNAStringSet("Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz")
seqs

# Our chromosome name is too verbose.
# Remove everything from the name after the first space.
names(seqs)
names(seqs) <- sub(" .*","",names(seqs))
names(seqs)

# Conversely, a DNAStringSet can be written to a file with
# writeXStringSet.

# 5.2 Loading features ----
#
# Genome annotations are available in a variety of text formats such as
# GFF3 and GTF. They can be loaded with the import function from
# rtracklayer. This GTF file is also from Ensembl, and gives the
# locations of the genes in the genome, and features within them.

### TSV with columns: seqname, source, start, end, score, strand, frame, attributes
### The start of the .gtf file looks like this:
# #!genome-build WBcel235
# #!genome-version WBcel235
# #!genome-date 2012-12
# #!genome-build-accession GCA_000002985.3
# #!genebuild-last-updated 2019-01
# V       WormBase        gene    9244402 9246360 .       -       .       gene_id "WBGene00000003"; gene_name "aat-2"; gene_source "WormBase"; gene_biotype "protein_coding";
# V       WormBase        transcript      9244402 9246360 .       -       .       gene_id "WBGene00000003"; transcript_id "F07C3.7.1"; gene_name "aat-2"; gene_source "WormBase"; gene_biotype "protein_coding"; transcript_source "WormBase"; transcript_biotype "protein_coding";
# V       WormBase        exon    9246080 9246360 .       -       .       gene_id "WBGene00000003"; transcript_id "F07C3.7.1"; exon_number "1"; gene_name "aat-2"; gene_source "WormBase"; gene_biotype "protein_coding"; transcript_source "WormBase"; transcript_biotype "protein_coding"; exon_id "F07C3.7.1.e1";
# V       WormBase        CDS     9246080 9246352 .       -       0       gene_id "WBGene00000003"; transcript_id "F07C3.7.1"; exon_number "1"; gene_name "aat-2"; gene_source "WormBase"; gene_biotype "protein_coding"; transcript_source "WormBase"; transcript_biotype "protein_coding"; protein_id "F07C3.7.1";
# V       WormBase        start_codon     9246350 9246352 .       -       0       gene_id "WBGene00000003"; transcript_id "F07C3.7.1"; exon_number "1"; gene_name "aat-2"; gene_source "WormBase"; gene_biotype "protein_coding"; transcript_source "WormBase"; transcript_biotype "protein_coding";
# ...

features <- import("Caenorhabditis_elegans.WBcel235.104.gtf.gz")

# Optional: just retain the columns of metadata we need
mcols(features) <- mcols(features)[,c("type","gene_name","gene_id","transcript_biotype","transcript_id")]

features



#///////////////////////////////////
# 6 Seqinfo and genome versions ----
#
# Various objects have associated Seqinfo information, listing the
# chromosomes they refer to and their lengths. This allows some extra
# sanity checking, and is also necessary for some tasks.

seqinfo(features)
seqinfo(seqs)

myseqinfo <- seqinfo(seqs)
isCircular(myseqinfo) <- c(rep(FALSE,6),TRUE)

seqinfo(features) <- myseqinfo

# One thing you need to care about is the genome version. We often refer
# to genome versions by their name in the [UCSC genome
# browser](https://genome.ucsc.edu/). For example the two most recent
# versions of the worm genome are "ce10" from 2010 and "ce11" from 2013.
# Similarly with human data you may find data for either the "hg19" or
# "hg38" versions of the human genome.

Seqinfo(genome="ce10")
Seqinfo(genome="ce11")

# Notice the slightly different chromosome lengths! Features from one
# genome version will not be correctly located on a different genome
# version.
#
# Our ENSEMBL sequences match "ce11", also called "WBcel235". However
# the chromosome names are slightly different, a common source of pain.
# See ?seqlevelsStyle.



#////////////////////////////////
# 7 Querying GRanges objects ----
#
# The metadata columns let us query the GRanges, for example for a
# feature type:

subset(features, type == "CDS")

# Equivalently:
#   features[features$type == "CDS"]

# Note: subset is a generic R function. It is similar to dplyr's filter.
# The second argument is special, in it you can refer to columns of the
# GRanges directly.



#/////////////////////////////////////////////
# 8 Representation of genes in a GTF file ----
#
# Let's have a look at the different feature types in a particular gene.

trx1_features <- subset(features, gene_name == "trx-1")
# Equivalently:
#   features[features$gene_name == "trx-1" & !is.na(features$gene_name)]

trx1_features
look(trx1_features, trx1_features$type)

# Look at the different types in the "type" column. Each "gene" may have
# multiple "transcript" features (isoforms). Each transcript in turn has
# a set of "exon" features, and if it is a protein coding gene, a set of
# "CDS" (coding sequence) features. The CDS features cover a subset of
# the bases covered by the exon features. The "gene_id" and
# "transcript_id" columns let us know which transcript and gene each
# exon belongs to.

# --------------------------------------------------> gene
#
# -------------------------------------------->       transcript
# ---------->         --->    ---------------->       exon
#       ---->         --->    ---------->             CDS
#
#
#                -----------------------------------> transcript
#                -------->       ---->    ----------> exon
#                     --->       ---->    -->         CDS

# Let's look at this in the Integrative Genome Browser (IGV):
#
# * Open IGV.
# * Select the genome "C. elegans (ce11)" (top left drop-down box).
# * Load the GTF file ("File" menu, "Load from file...").
# * Search for "trx-1" in the location box (next to the "Go" button).
#
# IGV turns all these features into a nice diagram of the gene's
# transcripts. The ENSEMBL gene annotation contains 5'UTR and 3'UTR
# regions, which are missing from the default gene annotations shown by
# IGV (probably from NCBI's RefSeq). Gviz or ggbio can be used to
# produce diagrams like this from within R (with a little work).



#/////////////////////////////////////
# 9 Further operations on GRanges ----
#
# 💡 Using accessor methods like "start()" and "strand()" frees the
# authors of GRanges to store data internally using whatever efficiency
# tricks they want, and even to change this in a new version of
# Bioconductor. This is good, but it can be better:
#
# 💡 Double-stranded DNA has a rotational symmetry. None of the physics
# of DNA is changed if we look at it rotated end-to-end 180 degrees.
# Labelling of the strands as "+" and "-" is arbitrary. Also, it is
# usually not the exact position of features in a genome that is
# important so much as their relative position. If we restrict outselves
# to using functions that work just as well with either labelling of
# strands, and that are not affected by absolute position, our code will
# more directly express our intention and be less prone to bugs.

# 9.1 Intra-range ----
#
# Various useful manipulations of individual ranges are defined.

?"intra-range-methods"

# Note: How these make use of the strand is a little haphazard. For
# example flank(), resize(), and promoters() respect strand but shift()
# does not.
#
# For example, suppose we are interested in the first 100 bases of a
# transcript:

feat <- features[2]
feat_initial <- resize(feat, 100, fix="start")

look( c(feat, feat_initial) )
getSeq(seqs, feat_initial)

# resize can fix either the fix="start" or fix="end" of the sequence.
#
# flank can be either flank the start (start=TRUE) or end (start=FALSE).

#                                  5'     3'
# input                            ------->
#                                  .      .
# resize(input, n, fix="start")    -->    .
# resize(input, n, fix="end")      .    -->
#                                  .      .
# flank(input, n, start=TRUE)   -->.      .
# flank(input, n, start=FALSE)     .      .-->
#                                  .      .
# promoters(input, n1,n2)     -------->   .
#

# 9.2 Inter-range ----

?"inter-range-methods"
?"nearest-methods"
?"setops-methods"

# One compelling feature of GenomicRanges is that it is able to find
# overlapping ranges very quickly.

query <- as("II:7752600-7754000:-", "GRanges")
hits <- findOverlaps(query, features, ignore.strand=TRUE)
hits
subjectHits(hits)
features[subjectHits(hits)]

findOverlaps(query, features, ignore.strand=FALSE)

# With findOverlaps, we can use genomic location as the key when
# *joining* disparate types of data, so this is an important tool for
# *integrative* analysis. See also the related functions nearest,
# precede, follow, and distance.
#
# GenomicRanges also provides:
#
# * range - get a feature that spans from the start to the end of all
# features in a GRanges.
# * reduce - merge overlapping features, so that the same bases are
# covered by a reduced collection of features.
# * disjoin - as with reduce, but broken at each start and end of the
# input features.
# * setdiff - subtracts one set of features from another, could be used
# with range on a set of exons to get introns. Might need to use
# GenomicRanges::setdiff if also using dplyr.

# input
#         --------->
#              ----------------->
#                      ------>
#                                     ---->
#
# range   -------------------------------->
#
# reduce  ----------------------->    ---->
#
# disjoin ---->---->-->------>--->    ---->
#
# setdiff(range(input),input)     --->

trx1_exons <- subset(trx1_features, type == "exon")
look( trx1_exons )
look( range(trx1_exons) )
look( reduce(trx1_exons) )
look( disjoin(trx1_exons) )

# 9.3 Challenge: transcript initiation ----

transcripts <- subset(features, type == "transcript")

# 1. Make a GRanges of all transcripts with a biotype of
# "protein_coding".
#
# 2. Make a GRanges of the two bases immediately preceding each of these
# transcripts.
#
# 3. Get the DNA sequences of the ranges from part 2. Use table() to see
# if any sequences are particularly common.
#
# Advanced: Get the four bases spanning the start of the transcript (two
# upstrand of the start and two downstrand of the start).
#
# Note: In C. elegans transcript initiation is complicated by trans-
# splicing. I wonder what exactly the 5' ends of these transcript
# annotations represent?
#


#////////////////////////////////////////////////////////////
# 10 Example: Examining the PolyAdenylation Signal (PAS) ----
#
# We'll now look at an example where the location of the thing we're
# looking for is more uncertain.
#
# We suspect there is some motif that causes cleavage and
# polyadenylation of RNA transcripts, close to the 3' end of each
# transcript. Bioconductor will let us explore sequence around the ends
# of transcripts. It can also help prepare data for command-line
# software such as the MEME Suite.
#
# We'll limit our attention to protein coding transcripts.

transcripts <- subset(features,
    type == "transcript" & transcript_biotype == "protein_coding")

ends <- resize(transcripts, width=30, fix="end")
end_seqs <- getSeq(seqs, ends)
names(end_seqs) <- ends$transcript_id

# We can check for general biasses in base composition:

library(seqLogo)

letter_counts <- consensusMatrix(end_seqs)

props <- prop.table(letter_counts[1:4,], 2)
seqLogo(props, ic.scale=FALSE)
seqLogo(props)

# Saving the sequences to a file lets us apply command-line tools:

writeXStringSet(end_seqs, "end_seqs.fasta")

# This is suitable for use with, for example, the streme tool from MEME.

# streme --rna --minw 4 --maxw 15 --p end_seqs.fasta

# Ok ok, so actually there is known PolyAdenylation Signal (PAS) motif,
# canonically "AAUAAA". Let's now look for this known motif.

counts  <- vcountPattern("AATAAA", end_seqs)
table(counts)

matches <- vmatchPattern("AATAAA", end_seqs)
plot( tabulate(start(unlist(matches)),30), xlab="start position",ylab="matches")

# We could also look for this pattern in the genome in general. Here is
# some code to find matches on either strand, using vmatchPattern().

query <- DNAString("AATAAA")
max.mismatch <- 0

fwd <- vmatchPattern(query, seqs, max.mismatch=max.mismatch)
fwd <- as(fwd, "GRanges")
strand(fwd) <- "+"

rev <- vmatchPattern(reverseComplement(query), seqs, max.mismatch=max.mismatch)
rev <- as(rev, "GRanges")
strand(rev) <- "-"

matches <- c(fwd, rev)
matches

# Once we have a set of ranges, they can be related to other sets of
# ranges. For example, to reproduce the earlier result:

overlaps <- findOverlaps(matches, ends, type="within")
table(table( subjectHits(overlaps) ))

# The matches can be examined in IGV if saved in a GFF file. With a
# "tabix" index, IGV avoids the need to load the whole file.

export(matches, "motif-matches.gff", index=TRUE)

# The pairwiseAlignment() function provides more flexibility than
# vmatchPattern(), if needed. Command-line tools such as BLAST will be
# faster for longer query sequences.
#
# Similar exploration can be performed around other regions of interest,
# such as peaks identified from ChIP. R may either provide a complete
# solution or serve as the glue plugging command-line tools together.



#/////////////////////////
# 11 BAMs and bigWigs ----
#
# import() can read various other file types. This section will
# demonstrate reading read alignments from a BAM file and producing a
# bigWig depth of coverage file. The reads in the BAM file used here is
# a small sample from GSE57993 sample "N2 Rep1".
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57993
#
# Both BAM and bigWig files are designed so that specific regions can be
# accessed efficiently. Genome browsers can use this to only load data
# as needed. BigWig files also contain versions of their data at
# different resolutions, so a zoomed-out view can be quickly loaded.
# Bioconductor packages will also usually offer some way to load
# specific regions. This is useful because these files can be very
# large.

alignments <- import("example.bam")
alignments

depth <- coverage(alignments)
depth

# Depth of coverage for gene rps-12
plot(depth$III[5679200:5680200], type="l")

export(depth, "depth.bw")

# Explore the documentation for "GAlignments", "Rle", and "RleList".
# "Rle" and "RleList" objects are memory-efficient vectors and lists of
# vectors. They support arithmetic operation such as scaling with * and
# adding with +. Rather neat! Sometimes a function won't support "Rle"s.
# When this happens, an Rle can be converted to a conventional vector
# with as.numeric().
#
# Filtered import of a BAM file can be performed with a suitable use of
# ScanBamParam from the Rsamtools package, import("example.bam",
# param=ScanBamParam( ... )).
#
# * Specific strand.
# * Specific region of reference genome.
# * Specific set of single cells, based on CB attribute.
#
# Try loading the "motif-matches.gff", "example.bam", and "depth.bw"
# files into IGV. Look for example at gene "rps-12".
#
# (return to slideshow)



#////////////////////////////////////////
# 12 Genome and annotation resources ----
#
# Besides software, Bioconductor includes various types of reference
# data.

# 12.1 Example packages ----
#
# * *BSgenome.Hsapiens.UCSC.hg38* Biostrings genome, Homo sapiens, from
# the UCSC browser, version hg38. DNA for chromosomes, usable in the
# same way as the DNAStringSet used above.
#
# * *TxDb.Hsapiens.UCSC.hg38.knownGene* Transcript database, Homo
# sapiens, from UCSC browser, genome verison hg38, "knownGene" gene
# annotations. GRanges information for genes and transcripts, much as we
# loaded from a GTF file above.
#
# * *org.Hs.eg.db* Organism Homo sapiens, primary key is Entrez Gene,
# database. Translation of gene ids from various databases, assignment
# to GO terms, KEGG pathways, etc. Entrez Gene ids are used as the
# primary key.
#
# * *GO.db* GO Gene Ontology term descriptions and relationships.
# However to find which terms apply to specific genes, use the GOALL
# column in the relevant species' organism database.
#
# Loading any of these packages gives you an AnnotationDb object
# (defined in the package AnnotationDbi). The main way to query these
# objects is with the select() function. mapIds() can be used if you
# need to enforce a one-to-one mapping.

library(org.Ce.eg.db)

keytypes(org.Ce.eg.db)
columns(org.Ce.eg.db)

AnnotationDbi::select(org.Ce.eg.db,
    keys="trx-1", keytype="SYMBOL", columns=c("ENSEMBL","ENTREZID"))

head( AnnotationDbi::select(org.Ce.eg.db,
    keys="WBGene00015062", keytype="ENSEMBL", columns="GOALL") )


library(GO.db)

AnnotationDbi::select(GO.db,
    keys="GO:0005622", keytype="GOID", columns="TERM")

# 12.2 Transcript database objects ----
#
# We've been using our genomic features as one big unstructured GRanges.
# This is messy. TxDb objects offer a more structured representation of
# genes, transcripts, exons, and coding sequences.

library(GenomicFeatures)

txdb <- makeTxDbFromGRanges(features)
txdb

# txdb is a "TxDb" object. Some ways to extract data from a TxDb are:

genes(txdb)
transcriptsBy(txdb, by="gene")
exonsBy(txdb, use.names=TRUE)
cdsBy(txdb, use.names=TRUE)

cds_ranges <- cdsBy(txdb, use.names=TRUE)
cds_ranges$B0228.5a.1
cds_ranges[["B0228.5a.1"]]
unlist(cds_ranges)

# cds_ranges here is a "GRangesList". That is, a list containing GRanges
# objects. To get the transcript sequence or the coding sequence, these
# each need to be retrieved and then concatenated together.
# extractTranscriptSeqs() can do this for us.

extractTranscriptSeqs(seqs, cds_ranges)

# A TxDb can also be accessed like an AnnotationDb.

keytypes(txdb)
columns(txdb)

# Under the hood these are all SQLite databases. select() joins tables
# together as needed to answer a query.

DBI::dbListTables( dbconn(txdb) )
DBI::dbListFields( dbconn(txdb), "transcript" )

# There's much more to explore here. Have a look at the documentation
# for the GenomicFeatures and ensembldb packages. To annotate a set of
# genomic features such as peaks from a ChIP-seq experiment, see for
# example the ChIPseeker package.

# 12.3 biomaRt ----
#
# [BioMart](http://www.biomart.org/) servers accessed using the biomaRt
# package are another way to get information such as translation of gene
# ids, gene sets, and gene information. For reproducibility, if using
# the ENSEMBL servers, make sure to specify a specific version of
# Ensembl to use. Your scripts may fail in future if the server you are
# using goes down or changes how it is accessed.

library(biomaRt)

mart <- useEnsembl(biomart = "genes", dataset = "celegans_gene_ensembl", version=104)

getBM(
    attributes=c("external_gene_name", "ensembl_gene_id", "entrezgene_id"),
    filters="external_gene_name", values="trx-1", mart=mart)

# 12.4 AnnotationHub ----
#
# AnnotationHub is a way to retrieve data from a more comprehensive set
# of organisms and data providers than are provided as individual
# Bioconductor packages. The retrieved data is returned in an
# appropriate Bioconductor type. If data is being updated over time (eg
# improved annotation of a genome), each version receives a unique ID in
# AnnotationHub, making it possible to write reproducible analyses.
#
# Files are cached, so they will only be downloaded once.
#
# In the example below, the yeast genome and Ensembl annotations are
# retrieved:

library(AnnotationHub)
ah <- AnnotationHub()

# ah contains a large collection of records that can be retrieved
ah
length(ah)
colnames( mcols(ah) )
table( ah$rdataclass )

# There is an interactive Shiny search interface
display(ah)

# query() searches for terms in an unstructured way
query(ah, c("Ensembl", "96", "Saccharomyces cerevisiae"))

# Having located records of interest,
# your R script can refer to the specific AH... record,
# so it always uses the same version of the data.
sc_genome  <- ah[["AH70449"]]
sc_granges <- ah[["AH69700"]]
sc_txdb    <- ah[["AH69265"]]

# sc_genome is a TwoBitFile
# Can use getSeq on it without loading everything into memory
seqinfo(sc_genome)
getSeq(sc_genome, as("IV:1-10","GRanges"))
import(sc_genome)



#//////////////////////////////
# 13 Package versions used ----
#
# Use sessionInfo() to check the versions of packages used.

sessionInfo()

# [Home](index.html)
