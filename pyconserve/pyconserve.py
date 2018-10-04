import pdb
import sys
import os
import os.path
import argparse
import multiprocessing
from functools import partial
import pybedtools
import re
import pandas as pd
import numpy as np

from . import reduce_df
from . import summarize as sc
from .version import __version__


def getoptions():
    desc = "Intersect BED file with conservation bedGraph files " + \
           "and return conservation scores"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('bedfile', metavar='BED', nargs=1,
                        help="Input BED file")
    parser.add_argument('consfiles', metavar='BEDGRAPH', nargs='+',
                        help="Conservation bedGraph files to intersect. "
                             "e.g. could be all chr*.bedGraph.gz files")
    parser.add_argument('-s', '--summarize', action="store_true",
                        help="Summarize conservation scores by taking the "
                             "average per BED interval [%(default)s]")
    parser.add_argument('-c', '--cores', type=int, default=4,
                        help="Number of processing cores [%(default)s]")
    parser.add_argument('-d', '--cores2', type=int, default=2,
                        help="Number of processing cores for summary step. "
                        "[%(default)s]")
    parser.add_argument('-S', '--splitdir', type=str, default=None,
                        help="Directory to keep intersections separate for "
                             "each conservation file. e.g. chromosome-specific. "
                             "Output won't be written to stdout. [%(default)s]")
    parser.add_argument('-t', '--temp', type=str,
                        help="set temp directory [{}]".
                        format(pybedtools.get_tempdir()))
    parser.add_argument('-v', '--version', action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()

    if args.splitdir and not os.path.exists(args.splitdir):
        os.makedirs(args.splitdir)

    return args


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_chrom_from_file(filename):
    """
    Extract chromosome from filename
    """
    m = re.match(r"^(chr[^.]+)\.(.+).bedGraph.gz", os.path.basename(filename))
    if m:
        return m.groups()
    else:
        raise Exception("Unable to determine chromosome prefix")


def subset_chrom(bd, chrom):
    """
    Filter BED file for specified chromosome
    """
    result = bd.filter(lambda x: x.chrom == chrom).saveas()
    #return pybedtools.BedTool(result)
    return result


def remove_unused_bedgraphs(files, bed_chroms):
    """
    Go through the list of inputted bedGraphs and keep only the ones that
    have matching chromosomes. This is to ensure that we don't perform
    unncessary intersects
    """
    new_consfiles = []
    removed = []
    c = 0
    for bg in files:
        bg_chrom = get_chrom_from_file(bg)
        if bg_chrom in bed_chroms:
            new_consfiles.append(bg)
        else:
            removed.append(bg_chrom)
            c += 1
    eprint("Removed %d unused chromosomes: %s" % (c, ','.join(removed)))
    return new_consfiles


def subset_conservation(bg, bd):
    """
    Perform intersectBed between BED file and chromosome-specific
    bedGraph file
    """
    # Load bedGraph
    cons = pybedtools.BedTool(bg)
    chrom, track = get_chrom_from_file(bg)

    # Filter bed file
    chrom_bd = subset_chrom(bd, chrom)
    if chrom_bd.file_type == 'empty':
        eprint("Skipping %s" % chrom)
        fn = None
    else:
        assert len(chrom_bd) > 0
        outfile = os.path.join(pybedtools.get_tempdir(),
                               'pybedtools.%s.%s.tmp' % (chrom, track))
        inter = chrom_bd.intersect(cons, wo=True, sorted=True, output=outfile)
        fn = inter.fn

    os.remove(chrom_bd.fn)
    return fn


def summarize(filename, delete_input=False):
    if filename is None:
        return None
    df = sc.load_intersect(filename)
    summed = sc.groupby_cons(df)
    if delete_input:
        os.remove(filename)
    return summed


def process(bg, bd, summarize_scores=False):
    """
    This function will be called by parallel processing.
    """
    fn = subset_conservation(bg, bd)

    if summarize_scores:
        results = summarize(fn, delete_input=True)
        return results

    return fn
    

def main():
    """
    See https://daler.github.io/pybedtools/3-brief-examples.html#example-3-count-reads-in-introns-and-exons-in-parallel
    as example of multithreading using pybedtools
    """
    args = getoptions()

    # Load BED file
    eprint("Loading BED file")
    bed = pybedtools.BedTool(args.bedfile[0])

    # Multipool
    eprint("Intersecting...")
    pool = multiprocessing.Pool(args.cores)
    func = partial(subset_conservation, bd=bed)
    results = pool.map(func, args.consfiles)
    pool.close()

    if args.summarize:
        pool = multiprocessing.Pool(args.cores2)
        func = partial(summarize, delete_input=True)
        results = pool.map(summarize, results)
        pool.close()

    for item in results:
        if item is not None:
            if args.summarize:
                item.to_csv(sys.stdout, sep="\t", index=False, header=False)
            else:
                if args.splitdir:
                    os.rename(item, 
                              os.path.join(args.splitdir, 
                                           os.path.basename(item)
                                          )
                             )
                else:
                    with open(item, 'r') as fin:
                        for line in fin:
                            print(line.strip())
                    os.remove(item)


if __name__ == '__main__':
    main()
