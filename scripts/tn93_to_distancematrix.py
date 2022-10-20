#! /usr/bin/env python3
'''
Convert the CSV output from tn93 into the skbio.stats.distance.DistanceMatrix format Qiita supports for distance matrices
'''
from os.path import isfile
from skbio.stats.distance import DistanceMatrix
from sys import argv
import argparse

# main functionality
if __name__ == "__main__":
    # parse user args
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', required=False, type=str, default='stdin', help="Input TN93 Distances (CSV format)")
    parser.add_argument('-o', '--output', required=True, type=str, help="Output DistanceMatrix File")
    args = parser.parse_args()
    if isfile(args.output):
        raise ValueError("Output file exists: %s" % args.output)
    if args.input.lower() == 'stdin':
        from sys import stdin as infile
    elif not isfile(args.input):
        raise ValueError("Input file not found: %s" % args.input)
    elif args.input.lower().endswith('.gz'):
        raise TypeError("Gzipped TN93 CSVs are not supported. Use stdin:\nzcat <tn93.csv.gz> | %s -o <output_file>" % argv[0])
    elif not args.input.lower().endswith('.csv'):
        raise TypeError("Input file must be CSV: %s" % args.input)
    else:
        infile = open(args.input)

    # load distance matrix
    dist = dict() # dist[u][v] = distance between u and v
    for line in infile:
        if line.strip() == 'ID1,ID2,Distance':
            continue
        u, v, d = [x.strip() for x in line.replace(', ',' ').split(',')]; d = float(d)
        if u not in dist:
            dist[u] = dict()
        if v not in dist:
            dist[v] = dict()
        dist[u][v] = d
    infile.close()
    IDs = sorted(dist.keys())

    # create DistanceMatrix object and save to file
    def get_dist(u, v):
        if u == v:
            return 0.
        elif v in dist[u]:
            return dist[u][v]
        elif u in dist[v]:
            return dist[v][u]
        else:
            return 1.
    dm = DistanceMatrix.from_iterable(IDs, get_dist, keys=IDs)
    dm.write(args.output)
