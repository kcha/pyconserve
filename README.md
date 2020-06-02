Python package that helps intersect BED files with sequence conservation tracks.

This package was inspired by the need to extract conservation scores for a given
set of BED intervals. It is based on a solution described by Dave Tang in this
[blog post](https://davetang.org/muse/2012/08/07/sequence-conservation-in-vertebrates/).

`pyconserve` takes as input the following:

  1. A BED file containing the intervals of interest
  2. Phastcons or PhyloP bedGraph files (one per chromosome)

It then uses [pybedtools] to perform genomic intersections between the input BED
file and bedGraph files. As each chromosome is queried independently, this can
be run using parallelization to speed up runtime.

**Note: this software is in beta and may contain bugs.**

# Getting started

## Installation

  1. Install bedtools
  2. Clone the pyconserve GitHub repo

          git clone https://www.github.com/kcha/pyconserve.git
          cd pyconserve

  3. Install the package

          python setup.py install
 

## Preparing conservation files

The following steps describe how to download conservation tracks and convert
them to bedGraph format using the UCSC tools `wigToBigWig` and
`bigWigToBedGraph`. These steps are adapted from this [blog post](https://davetang.org/muse/2012/08/07/sequence-conservation-in-vertebrates/).

  1. Download Phastcons or PhyloP wiggle files (one chromosome per file) from UCSC 
  1. Convert wiggle to bigWig using [`wigToBigWig`](https://anaconda.org/bioconda/ucsc-wigtobigwig). 

          wigToBigWig chr1.phastCons100way.wigFix.gz hg19.chrom_sizes.txt \
              chr1.phastCons100way.bigWig
            
     This step requires a file containing the chromosome sizes of your species.

  1. Convert bigWig to bedGraph using [`bigWigToBedGraph`](https://anaconda.org/bioconda/ucsc-bigwigtobedgraph)

          bigWigToBedGraph chr1.phastCons100way.bigWig chr1.phastCons100way.bedGraph

  1. (optional) To save space, compress the bedGraph file and remove the wiggle
     and bigWig files

          rm -v chr1.phastCons100way.wigFix.gz chr1.phastCons100way.bigWig
          gzip -vf chr1.phastCons100way.bedGraph

## Run `pyconserve`

Given a BED file and a set of bedGraph conservation files, `pyconserve` can be
run as follows:

    pyconserve a.bed chr*.phastCons100way.bedGraph.gz > a_conservation.bed

This will perform `bedtools intersect` and save the output to
`a_conservation.bed`.

## Aggregate conservation scores

To aggregate these results by computing the mean conservation score for each
interval in the BED file, use the command `summarize_conserve`:

    summarize_conserve a_conservation.bed > summarized.bed

Alternatively, the above two commands can be chained together as follows:

    pyconserve a.bed chr*.phastCons100way.bedGraph.gz | \
        summarize_conserve - > summarized.bed

## Acknowledgements

This software is inspired by previous work from others:

  - Extracting sequence conservation: https://davetang.org/muse/2012/08/07/sequence-conservation-in-vertebrates/
  - Multiprocessing using pybedtools: https://daler.github.io/pybedtools/3-brief-examples.html#example-3-count-reads-in-introns-and-exons-in-parallel

