#!/usr/bin/python3
from argparse import ArgumentParser
import csv
from itertools import compress

def main():
    parser = get_arg_parser()  # Parse input to variables
    args = vars(parser.parse_args())
    read_id_mapped_path = args["read_id_mapped_path"]
    origin_path = args["origin_path"]
    outputfile_path = args["outputfile_path"]
    read_ids_mapped = set()
    def off_target(origin): #returns True when the read is off_target
        fil = ['' != value for value in origin]
        origin_fil = list(compress(origin,fil))
        if len(origin_fil) > 0:
            boolean = [origin_fil[0]]*len(origin_fil) == origin_fil #Check whether the filtered list only contains the same primer origin (This means the read contains only arms from the same probe)
        else:
            boolean = True #Returns False when no arms can be found in read
        return not boolean



    with open(read_id_mapped_path, 'r') as infile_mapped:
        reader = csv.reader(infile_mapped,delimiter='\t')
        for row in reader:
            read_id = row[0].split(' ')[0]
            read_ids_mapped.add(read_id)

    mapped_read_count = 0
    total_read_count = 0
    chimeric_read_count = 0 
    single_read_count = 0
    no_probe_read_count = 0

    chimeric_read_insert_count = 0
    chimeric_read_no_insert_count = 0
    single_read_insert_count = 0
    single_read_no_insert_count = 0
    origin_dict_chim = {}
    origin_dict_chim_no_insert = {}
    origin_dict_single = {}
    origin_dict_single_no_insert = {}
    origin_id_dict_single = {}
    with open(origin_path, 'r') as infile:
        reader = csv.reader(infile,delimiter = '\t')
        for row in reader:
            total_read_count += 1
            read_id = row[0].split(' ')[0]
            a1_origin = row[1]
            g1_origin = row[2]
            A2_origin = row[3]
            G2_origin = row[4]
            R1_insert = row[5]
            R2_insert = row[6]
            if R1_insert == 'no_insert':
                R1_insert = None
            else:
                R1_insert = int(R1_insert)
            if R2_insert == 'no_insert':
                R2_insert = None
            else:
                R2_insert = int(R2_insert)
            #Count reads that are mapped
            if read_id in read_ids_mapped:
                mapped_read_count += 1
            else:
                origins = [a1_origin,g1_origin,A2_origin,G2_origin]
                key = ':'.join(row[1:7])
                if off_target(origins): # When true, the read originates from multiple probes.: Chimeric read.
                    chimeric_read_count += 1 #Count reads that originate from multiple probes
                    if R1_insert !=None and R1_insert > 5 and R2_insert!=None and R2_insert > 5:
                        chimeric_read_insert_count += 1
                        if key in origin_dict_chim.keys():
                            origin_dict_chim[key] = origin_dict_chim[key] + 1
                        else:
                            origin_dict_chim[key] = 1
                    else:
                        chimeric_read_no_insert_count += 1
                        if key in origin_dict_chim_no_insert.keys():
                            origin_dict_chim_no_insert[key] = origin_dict_chim_no_insert[key] + 1
                        else:
                            origin_dict_chim_no_insert[key] = 1

                elif origins == ['','','','']:
                    no_probe_read_count += 1 # count reads that do not seem to originate from a probe.
                else:
                    single_read_count += 1 #count reads that originate from a single probe
                    if R1_insert!=None and R1_insert > 5 and R2_insert!=None and R2_insert > 5:
                        single_read_insert_count += 1
                        if key in origin_dict_single.keys():
                            origin_dict_single[key] = origin_dict_single[key] + 1
                            origin_id_dict_single[key].append(read_id)
                        else:
                            origin_dict_single[key] = 1
                            origin_id_dict_single[key] =  [read_id]
                    else:
                        single_read_no_insert_count += 1
                        if key in origin_dict_single_no_insert.keys():
                            origin_dict_single_no_insert[key] = origin_dict_single_no_insert[key] + 1
                        else:
                            origin_dict_single_no_insert[key] = 1

    # Per read group output the read_ids 
    # Per read group list the most common reads and their parameters:
            # which probe arms
            # length of insert
            # sequence of insert

    print('total_read_count: ' + str(total_read_count)+'\n')
    print('mapped_read count: ' + str(mapped_read_count)+ '\n')
    print('Chimeric read count: ' + str(chimeric_read_count))
    print('\tChimeric read insert count: ' + str(chimeric_read_insert_count))
    print(list(sorted(origin_dict_chim.items(), key = lambda item: item[1]))[-10:-1])
    print('\n')
    print('\tChimeric read no insert count: ' + str(chimeric_read_no_insert_count))
    print(list(sorted(origin_dict_chim_no_insert.items(), key = lambda item: item[1]))[-10:-1])
    print('\n')
    print('no probe read count: ' + str(no_probe_read_count))
    print('\n')
    print('single_read_count: ' + str(single_read_count))
    print('\tsingle_read_insert_count: ' + str(single_read_insert_count))
    print(list(sorted(origin_dict_single.items(), key = lambda item: item[1]))[-10:-1])
    print(sum([item for item in list(origin_dict_single.values()) if item > 10]))
    print(sum([item for item in list(origin_dict_single.values()) if item <= 10]))

    print('\n')
    print('\tsingle_read_no_insert_count: ' + str(single_read_no_insert_count))
    print(list(sorted(origin_dict_single_no_insert.items(), key = lambda item: item[1]))[-10:-1])
    print('\n')
    with open(outputfile_path,'w') as outputfile:
        outputfile.write('Total,Mapped,Chimeric,Chimeric_insert,Chimeric_no_insert,Single,Single_insert,Single_no_insert,No_probe' + '\n')
        outputfile.write(str(total_read_count) + ',' + str(mapped_read_count) + ',' + str(chimeric_read_count) + ',' + str(chimeric_read_insert_count) + ',' + str(chimeric_read_no_insert_count) + ',' + str(single_read_count) + ',' + str(single_read_insert_count) + ',' + str(single_read_no_insert_count) + ',' + str(no_probe_read_count))

def get_arg_parser():
    # parse all files
    parser = ArgumentParser()
    parser.add_argument("-m", "--read_id_mapped_path", dest="read_id_mapped_path",
                        help="tsv file containing the list of mapped read ids", metavar="read_id_mapped_path")
    parser.add_argument("-r", "--origin_path", dest="origin_path",
                        help="Input tab-delimited file containing origin info for each read", metavar="origin_path")
    parser.add_argument("-o", "--outputfile_path", dest="outputfile_path",
                        help= "Output csv file to summarize read groups", metavar="outputfile_path")
                
    return parser
if __name__ == '__main__':
    main()
