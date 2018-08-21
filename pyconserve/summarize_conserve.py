import sys
import argparse
import pandas as pd

from . import reduce_df
from .version import __version__

def getoptions():
    desc = "Summarize conservation scores by BED interval"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('bedfile', metavar='BED', nargs=1,
                        help="Input BED file")
    parser.add_argument('-v', '--version', action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    return args


def load_intersect(filename):
    """
    Load intersected BED file as DataFrame
    """
    colnames = ['chr', 'start', 'end', 'name', 'score', 'strand',
                'chr2', 'start2', 'end2', 'score2', 'length']
    usecolnames = ['chr', 'start', 'end', 'name', 'score', 'strand',
                   'score2', 'length']
    df = pd.read_table(filename, names=colnames)
    return reduce_df.reduce_mem_usage(df)[0]
   

def groupby_cons(df):
    """
    Group by BED interval and report the average score
    """
    df['sum'] = df['score2'] * df['length']
    result = df.groupby(df.columns[0:6].tolist())\
               .apply(compute_average_cons)\
               .reset_index(name="avg_cons")
    return result


def compute_average_cons(data):
    """
    Compute averge conservation 
    """
    return sum(data.score2 * data.length) / sum(data.length)


def main():
    """
    Summarize the collected conservation scores by grouping the BED intervals
    """
    # Convert BedTool to DataFrame
    args = getoptions()

    if args.bedfile[0] == '-':
        df = load_intersect(sys.stdin)
    else:
        df = load_intersect(args.bedfile[0])

    df = groupby_cons(df)

    df.to_csv(sys.stdout, sep="\t", index=False, header=False)


if __name__ == '__main__':
    main()
