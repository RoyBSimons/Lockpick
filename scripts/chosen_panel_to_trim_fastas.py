#!/usr/bin/python3
# Create primer files for adapter and probe arm removal with cutadapt
import csv
from argparse import ArgumentParser
from Bio.Seq import Seq

def main():
    parser = get_arg_parser()  # Parse input to variables
    args = vars(parser.parse_args())
    chosen_panel_path = args["chosen_panel_path"]
    primer_a_fasta_path = args["primer_a_fasta_path"]
    primer_A_fasta_path = args["primer_A_fasta_path"]

    with open(chosen_panel_path) as input:
        reader = csv.reader(input)
        upstream_arm_list = [row[3] for row in reader]
    with open(chosen_panel_path) as input:
        reader = csv.reader(input)
        downstream_arm_list = [row[5] for row in reader]

    with open(primer_A_fasta_path, 'w') as output:
        for i,row in enumerate(downstream_arm_list[1:]):
            output.write('>primers_'+str(i)+'\n')
            output.write(row+'...'+upstream_arm_list[i+1]+'\n')

    with open(primer_a_fasta_path, 'w') as output:
        for i,row in enumerate(downstream_arm_list[1:]):
            output.write('>primers_'+str(i)+'\n')
            output.write(str(Seq(upstream_arm_list[i+1]).reverse_complement())+'...'+str(Seq(row).reverse_complement())+'\n')

def get_arg_parser():
    # parse all files
    parser = ArgumentParser()
    parser.add_argument("-c", "--chosen_panel", dest="chosen_panel_path",
                        help="The chosen panel csv file created by Locksmith", metavar="chosen_panel_path")
    parser.add_argument("-a", "--primer_a", dest="primer_a_fasta_path",
                        help="Output filename for primer_a as input for cutadapt", metavar="primer_a_fasta_path")
    parser.add_argument("-A", "--primer_A", dest="primer_A_fasta_path",
                        help="Output filename for primer_A as input for cutadapt", metavar="primer_A_fasta_path")
    return parser
if __name__ == '__main__':
    main()
