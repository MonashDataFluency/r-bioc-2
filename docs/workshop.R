


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

vignette()
vignette(package="Biostrings")
vignette("BiostringsQuickOverview", package="Biostrings")



#//////////////////////////////////////
# 2 Packages used in this workshop ----

library(Biostrings)      # Provides DNAString, DNAStringSet, etc
library(BSgenome)        # Provides getSeq()
library(GenomicRanges)   # Provides GRanges containing genomic ranges, etc
library(GenomicFeatures) # Provides TxDb objects containing genes/transcripts/exons
library(rtracklayer)     # Provides import() and export()
library(Gviz)            # Provides plotting of genomic features
library(seqLogo)         # Provides sequence logo plots for motifs, etc



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

# 3.3 Challenge ----
#
# Reverse complement the following DNA sequence and then translate to an
# amino acid sequence:

GCTTTCGTTTTCGCC
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
        from=min(start(features))-2,
        to=max(end(features))+2)
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

# Conversely, a DNAStringSet can also be written to a file with
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
# See ?seqlevelsStyle().



#////////////////////////////////
# 7 Querying GRanges objects ----
#
# The metadata columns let us query the GRanges, for example for a
# feature type:

cds <- subset(features, type == "CDS")
cds
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

# 9.1 Intra-range ----
#
# Various useful manipulations of individual ranges are defined.

?"intra-range-methods"

# Note: How these make use of the strand is a little haphazard. For
# example flank() and resize() respect strand but shift() does not.
#
# Earlier we translated a coding sequence. Coding sequences are
# terminated by a stop codon. Let's extend the CDS feature to include
# this.

feat <- features[4]
feat_stop <- resize(feat, width(feat)+3)
seq_stop <- getSeq(seqs, feat_stop)
translate(seq_stop)

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



#////////////////////////////////////////////////////////////
# 10 Example: Examining the PolyAdenylation Signal (PAS) ----
#
# Say we suspect there is some motif that causes cleavage and
# polyadenylation of RNA transcripts, close to the 3' end of each
# transcript. Bioconductor will let us explore sequence around the ends
# of transcripts. It can also serve as to prepare data for command-line
# software such as the MEME Suite.
#
# We'll limit our attention to protein coding transcripts.

ce_transcripts <- subset(features, type == "transcript" & transcript_biotype == "protein_coding")

ce_ends <- resize(ce_transcripts, width=30, fix="end")
ce_end_seqs <- getSeq(seqs, ce_ends)
names(ce_end_seqs) <- ce_ends$transcript_id

# We can check for general biasses in base composition:

letter_counts <- consensusMatrix(ce_end_seqs)
probs <- prop.table(letter_counts[1:4,], 2)
seqLogo(probs, ic.scale=FALSE)
seqLogo(probs)

# Saving the sequences to a file lets us apply command-line tools:

writeXStringSet(ce_end_seqs, "ce_end_seqs.fasta")

# This is suitable for use with, for example, the streme tool from MEME.

# streme --rna --minw 4 --maxw 15 --p ce_end_seqs.fasta

# Ok, so actually it is well known that there is a PolyAdenylation
# Signal (PAS) motif, canonically "AAUAAA". Let's look for this known
# motif.

# Which contain the canonical PAS?
counts <- vcountPattern("AATAAA", ce_end_seqs)
table(counts)

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

overlaps <- findOverlaps(matches, ce_ends, type="within")
table(table( subjectHits(overlaps) ))

# The matches can be examined in IGV if save in a GFF file.

export(matches, "motif-matches.gff")

# The pairwiseAlignment() function provides more flexibility if needed.
# Command-line tools such as BLAST may be faster.
#
# Similar exploration can be performed around other regions of interest,
# such as peaks identified from ChIP. R will either provide a complete
# solution, or serve as the glue plugging command-line tools together.



#/////////////////////////
# 11 BAMs and Bigwigs ----
#
# import() can read various other file types. This section will
# demonstrate reading read alignments from a BAM file and producing a
# Bigwig depth of coverage file. The reads in the BAM file used here are
# a small sample from GSE57993 sample "N2 Rep1".
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57993
#
# Both BAM and Bigwig files are designed so that specific regions can be
# accessed efficiently. Genome browsers can use this to only load data
# as needed. Bioconductor packages will also usually offer some way to
# load specific regions. This is useful because these can be very large
# files.

alignments <- import("example.bam")
alignments

depth <- coverage(alignments)
depth

# Depth of coverage for gene rps-12
plot(depth$III[5679200:5680200], type="l")

export(depth, "depth.bw")

# Explore the documentation for "GenomicAlignments", "Rle", and
# "RleList". "Rle" and "RleList" objects are memory-efficient vectors
# and lists of vectors. They support arithmetic operation such as
# scaling with * and adding with +. Rather neat!
#
# Sometimes a function won't support "Rle"s. When this happens, an Rle
# can be converted to a conventional vector with as.numeric().
#
# Filtered import can be performed with a suitable use of ScanBamParam,
# import(bam, param=ScanBamParam( ... )):
#
# * specific strand
# * specific region of reference genome
# * specific set of cells, based on BAM attributes
#
# Try loading the "motif-matches.gff", "example.bam", and "depth.bw"
# files into IGV. Look for example at gene "rps-12".
#
# ...
#
# ...
#
# ...
#
# At this point, your brain should be full. We now want to finish be
# touching lightly on some further facilities of Bioconductor without
# going into too much detail.
#
# ...
#
# ...
#
# ...



#////////////////////////////////////
# 12 Transcript database objects ----
#
# We've been using our genomic features as one big unstructured GRanges.
# This is messy. Furthermore, eukaryote genes contain exons and introns.

txdb <- makeTxDbFromGRanges(features)
txdb

# txdb is a TxDb object.

genes(txdb)
transcriptsBy(txdb, by="gene")
exonsBy(txdb, use.names=TRUE)
cdsBy(txdb, use.names=TRUE)

cds_ranges <- cdsBy(txdb, use.names=TRUE)
cds_ranges$B0228.5a.1
cds_ranges[["B0228.5a.1"]]
unlist(cds_ranges)

# cds_ranges is a GRangesList. That is, a list containing GRanges
# objects. To get the transcript sequence or the coding sequence, these
# each need to be retrieved and then concatenated together.
# extractTranscriptSeqs can do this for us.

cds_seqs <- extractTranscriptSeqs(seqs, cds_ranges)
cds_seqs

# There's much more to explore here. Have a look at the documentation
# for the GenomicFeatures and ensembldb packages. To annotate a set of
# genomic features such as peaks from a ChIP-seq experiment, see for
# example the  ChIPseeker package.



#////////////////////////////////////////
# 13 Genome and annotation resources ----

# (return to slideshow)

# Besides software, Bioconductor includes various types of data. An
# important type of data is data describing model organisms. This is
# either supplied through data packages or through the newer
# AnnotationHub system. It is generally derived from these central
# repositories:
#
# * The NCBI's [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/)
# * The EBI's [Ensembl genome browser](https://ensembl.org/)
# * The [UCSC genome browser](https://genome.ucsc.edu/cgi-bin/hgGateway)
#
# UCSC was the original web-based genome browser. UCSC's "KnownGene"
# gene annotations used to be the cutting edge gene annotation source,
# but UCSC now relies on other sources for gene annotations. Many [file
# types](https://genome.ucsc.edu/FAQ/FAQformat.html) that remain
# important today were developed for the UCSC genome browser, such as
# "bed", "bigWig", and "2bit".
#
# These organizations will generally obtain genome assemblies from the
# same ultimate sources. For example, all of the above use the Genome
# Reference Consortium's GRCh38 DNA sequence for homo sapiens. UCSC
# calls this "hg38" but it is the same DNA sequence. These DNA sequences
# serve as a common frame of reference. However the three organizations
# above will differ on their exact set of gene and transcript
# annotations, and use different gene and transcript ID systems.
#
# Genome assemblies are released infrequently. GRCh38 (hg38) was
# released in 2013. The previous assembly, GRCh37 (hg19) was released in
# 2009. Some people haven't updated yet, you will find plenty of data
# using "hg19" positions! Gene and transcript annotations are updated
# far more frequently.
#
# As well as the chromosomes in the "primary assembly" a genome assembly
# may have further sequences, which may have been added after the
# initial release:
#
# * patch sequences: fixes that would change the sizes of chromosomes
# * alt loci: a way to represent alleles, genetic diversity in the
# species

# 13.1 Example packages ----
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

# 13.2 biomaRt ----
#
# [BioMart](http://www.biomart.org/) servers accessed using the biomaRt
# package are another way to get information such as translation of gene
# ids, gene sets, and gene information.
#
# Beware replicability. BioMart servers are not guaranteed to be
# available forever, and BioMart doesn't have a clear story around
# versioning.

# 13.3 AnnotationHub ----
#
# AnnotationHub is a way to retrieve data from a more comprehensive set
# of organisms and data providers than are provided as individual
# Bioconductor packages. The retrieved data is returned in an
# appropriate Bioconductor type. If data is being updated over time (eg
# improved annotation of a genome), each version receives a unique ID in
# AnnotationHub, making it possible to write reproducable analyses.
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
records <- query(ah, c("Ensembl", "96", "Saccharomyces cerevisiae"))
records

mcols(records)
mcols(records)[,c("title","rdataclass")]

# Having located records of interest,
# your R script can refer to the specific AH... record,
# so it always uses the same version of the data.
sc_genome <- ah[["AH70449"]]
sc_granges <- ah[["AH69700"]]
sc_txdb <- ah[["AH69265"]]

# sc_genome is a TwoBitFile
# Can use getSeq on it without loading everything into memory
seqinfo(sc_genome)
getSeq(sc_genome, as("IV:1-10","GRanges"))
import(sc_genome)

# An OrgDb contains information about genes in an organism, and lets you map between different identifiers
query(ah, c("OrgDb", "Saccharomyces cerevisiae"))
sc_orgdb <- ah[["AH84129"]]
columns(sc_orgdb)
head( keys(sc_orgdb, "ENSEMBL") )
select(sc_orgdb, "YFL039C", keytype="ENSEMBL", columns=c("GENENAME", "DESCRIPTION"))

# As well as IDs, genes have short, easy to remember "symbols" (also often called "names")
# We can use the OrgDb to look up gene IDs from symbols
# Notice a problem here!
select(sc_orgdb, c("ACT1", "COS2"), keytype="GENENAME",  columns=c("ENSEMBL"))



#////////////////////
# 14xxxxxxxxxxxx ----

# 14.1 Under the hood ----
#
# A TxDb is a subclass of AnnotationDb. Under the hood it is an SQLite
# database.

# Poking around inside
DBI::dbListTables( dbconn(txdb) )
DBI::dbListFields( dbconn(txdb), "gene" )
DBI::dbListFields( dbconn(txdb), "transcript" )
DBI::dbListFields( dbconn(txdb), "exon" )
DBI::dbListFields( dbconn(txdb), "cds" )
DBI::dbListFields( dbconn(txdb), "splicing" )

# AnnotationDb objects contain a set of tables. The select method will
# join tables together as necessary to output the requested set of
# columns. Depending on the columns, you may get less or more rows.

keytypes(txdb)
columns(txdb)

AnnotationDbi::select(txdb, key="WBGene00015062", keytype="GENEID", columns="TXID")
AnnotationDbi::select(txdb, key="WBGene00015062", keytype="GENEID",
                      columns=c("TXID","EXONID","EXONSTART","EXONEND","EXONSTRAND"))

# Various other reference information is given by Bioconductor in the
# form of AnnotationDb objects. We will see some of these later.

# 14.2 Challenge ----
#
# Get the 30 bases upstrand of the ends of transcripts which have the
# "protein_coding" biotype.
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
#
# Use sessionInfo() to check the versions of packages you have used.

sessionInfo()

# [Home](index.html)
