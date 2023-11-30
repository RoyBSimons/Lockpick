#!/usr/bin/python3
# Create reference genome (amplicons) from target sequences
import os
import csv
from argparse import ArgumentParser
from Bio.Seq import Seq
from Bio import SeqIO

def main():
    parser = get_arg_parser()  # Parse input to variables
    args = vars(parser.parse_args())
    chosen_panel_path = args["chosen_panel_path"]
    targets_fasta_path = args["targets_fasta_path"]
    targets_bed_path = args["targets_bed_path"]
    reference_genome_path = args["reference_genome_path"]
    reference_amplicons_path = args["reference_amplicons_path"]

    with open(chosen_panel_path) as input:
        reader = csv.reader(input)
        probe_target_regions = [(row[20],row[9]) for row in reader]

    with open(targets_bed_path, 'w', newline='') as output:
        writer = csv.writer(output,delimiter= '\t')
        for row in probe_target_regions[1:]:
            chrom = row[0].split(':')[0]
            start = row[0].split(':')[1].split('-')[0]
            end = str(int(row[0].split(':')[1].split('-')[1])+1)
            strand = row[1]
            writer.writerow([chrom,start,end,'.','.', strand])
    
    os.system('bedtools getfasta -fi ' + reference_genome_path  +  ' -bed ' + targets_bed_path + ' -fo ' + targets_fasta_path)

    records = list(SeqIO.parse(targets_fasta_path, "fasta"))
    with open(reference_amplicons_path, 'w', newline='') as output:
        for i,record in enumerate(records):
            amplicon = 'NNNN' + str(record.seq) + 'NNNN' # Add N-padding for bismark
            output.write('>Target' + str(i) + '_' +record.id + '\n')
            output.write(amplicon + '\n')

def G_to_A_convert(sequence):
    converted_sequence = sequence.replace('CG','__').replace('G','A').replace('__','CG')  # As the gDNA converts the C to T, the gap which is filled by extending the probe (which is the complement) must change from G to A, while keeping CGs intact.
    return converted_sequence

def get_arg_parser():
    # parse all files
    parser = ArgumentParser()
    parser.add_argument("-c", "--chosen_panel", dest="chosen_panel_path",
                        help="The chosen panel csv file created by Locksmith", metavar="chosen_panel_path")
    parser.add_argument("-t", "--targets_fasta", dest="targets_fasta_path",
                        help="Output filename for targets fasta file", metavar="targets_fasta_path")
    parser.add_argument("-b", "--targets_bed", dest="targets_bed_path",
                        help="Output filename for probe targets bed file", metavar="targets_bed_path")
    parser.add_argument("-r", "--reference_genome", dest="reference_genome_path",
                        help="Input filename for reference genome", metavar="reference_genome_path")
    parser.add_argument("-a", "--reference_amplicons", dest="reference_amplicons_path",
                        help="Output filename for reference amplicons", metavar="reference_amplicons_path")
    return parser

if __name__ == '__main__':
    main()
