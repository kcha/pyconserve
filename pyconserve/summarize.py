import sys
import argparse
import pandas as pd
import multiprocessing

from . import reduce_df
from .version import __version__

def getoptions():
    desc = "Summarize conservation scores by BED interval"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('bedfile', metavar='BED', nargs=1,
                        help="Input BED file")
    parser.add_argument('-c', '--cores', type=int, default=4,
                        help="Number of processing cores [%(default)s]")
    parser.add_argument('-b', '--batchsize', type=int,
                        default=2000000,
                        help="Chunk size [%(default)s]")
    parser.add_argument('-v', '--version', action='version',
        version='%(prog)s ' + __version__)
    args = parser.parse_args()
    return args


def load_intersect(filename, chunksize=2000000, **kwargs):
    """
    Load intersected BED file as DataFrame
    """
    colnames = ['chr', 'start', 'end', 'name', 'score', 'strand',
                'chr2', 'start2', 'end2', 'score2', 'length']
    usecolnames = ['name', 
                   'score2', 'start2', 'end2']
    column_types = {   
        'chr': 'category',
        'end': 'uint32',
        'length': 'uint16',
        'name': 'object',
        'score2': 'float32',
        'start': 'uint32',
        'start2': 'uint32',
        'end2': 'uint32'
    }

    reader = pd.read_table(filename, names=colnames, 
                               usecols=usecolnames,
                               dtype=column_types,
                               chunksize=chunksize)
    return reader
   

def sum_scores(data):
    le = data.end2 - data.start2
    N = sum(le)
    total_score = sum(data.score2 * le)
    return pd.DataFrame({'N': [N], 'total_score': [total_score]})


def count_scores(df):
    result = df.groupby('name')\
               .apply(sum_scores)
    return result


# def groupby_cons(df):
#     """
#     Group by BED interval and report the average score
#     """
#     result = df.groupby('name')\
#                .apply(compute_average_cons)\
#                .reset_index(name="avg_cons")
#     return result


def compute_average_cons(data):
    """
    Compute average conservation 
    """
    #le = data.end2 - data.start2
    #N = sum(le)
    #return sum(data.score2 * le) / N
    data = data.groupby('name').agg('sum')
    data['avg_cons'] = data.total_score / data.N
    return data


def main():
    """
    Summarize the collected conservation scores by grouping the BED intervals
    """
    # Convert BedTool to DataFrame
    args = getoptions()

    if args.bedfile[0] == '-':
        reader = load_intersect(sys.stdin, args.batchsize)
    else:
        reader = load_intersect(args.bedfile[0], args.batchsize)

    pool = multiprocessing.Pool(args.cores)
    df = pd.concat(pool.map(count_scores, reader))
    pool.close()

    df = compute_average_cons(df)

    df.loc[:,'avg_cons'].to_csv(sys.stdout, sep="\t", index=True, header=False)


if __name__ == '__main__':
    main()
