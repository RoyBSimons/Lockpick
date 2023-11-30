#!/usr/bin/python3
# Change targeted coverage file into genomic coverage file
import os
import csv
import gzip
import matplotlib.pyplot as plt
from argparse import ArgumentParser

def main():
    parser = get_arg_parser()  # Parse input to variables
    args = vars(parser.parse_args())
    targeted_cov_path = args["targeted_cov_path"]
    genome_cov_path = args["genome_cov_path"]
    output_dir = args["output_dir"]

    outrows = []
    target_dict = {}
    with gzip.open(targeted_cov_path,'rt') as infile:
        reader = csv.reader(infile,delimiter= '\t')
        for row in reader:
            target = row[0].split('_')[0]
            if target[6:] in target_dict:
                target_dict[target[6:]] += 1
            else:
                target_dict[target[6:]] = 1
            chrom = row[0].split('_')[1].split(':')[0]
            start = str(int(row[0].split(':')[1].split('-')[0]) + int(row[1]) - 5) #4N padding
            end = str(int(start) + 2)
            methylation_percentage = row[3]
            methylated_cytosines = row[4]
            unmethylated_cytosines = row[5]
            outrows.append([chrom,start,end,methylation_percentage, methylated_cytosines, unmethylated_cytosines])

    probe_nr_list = [int(probenr) for probenr in target_dict.keys()]
    coverage_list = [int(coverage) for coverage in target_dict.values()]
    plt.bar(probe_nr_list,coverage_list)
    plt.savefig(output_dir + 'Probe_distribution.png')

    with open(genome_cov_path, 'w') as outfile:
        writer = csv.writer(outfile, delimiter = '\t')
        writer.writerows(outrows)

def get_arg_parser():
    # parse all files
    parser = ArgumentParser()
    parser.add_argument("-t", "--targeted_cov_path", dest="targeted_cov_path",
                        help="Input cov file", metavar="targeted_cov_path")
    parser.add_argument("-g", "--genome_cov_path", dest="genome_cov_path",
                        help="Output cov file", metavar="genome_cov_path"),
    parser.add_argument("-o", "--output_dir", dest="output_dir",
                        help="Output directory", metavar="output_dir")
    return parser

if __name__ == '__main__':
    main()
