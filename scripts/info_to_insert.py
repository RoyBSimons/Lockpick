import csv
from argparse import ArgumentParser
def main():
    parser = get_arg_parser()  # Parse input to variables
    args = vars(parser.parse_args())
    outfile_path = args["outfile_path"]
    a1_info_path = args["a1_info_path"]
    A2_info_path = args["A2_info_path"]
    g1_info_path = args["g1_info_path"]
    G2_info_path = args["G2_info_path"]

    insert_list = []

    a1_loc_list = []
    A2_loc_list = []
    g1_loc_list = []
    G2_loc_list = []
    with open(a1_info_path,'r') as infile:
        reader = csv.reader(infile, delimiter = '\t')
        for row in reader:
            if row[1] == '-1': #no arm found
                a1_loc_list.append('no arm')
            else:
                a1_loc_list.append(row[2])
    with open(A2_info_path,'r') as infile:
        reader = csv.reader(infile, delimiter = '\t')
        for row in reader:
            if row[1] == '-1': #no arm found
                A2_loc_list.append('no arm')
            else:
                A2_loc_list.append(row[2])
    with open(g1_info_path,'r') as infile:
        reader = csv.reader(infile, delimiter = '\t')
        for row in reader:
            if row[1] == '-1': #no arm found
                g1_loc_list.append('no arm')
            else:
                g1_loc_list.append(row[3])
    with open(G2_info_path,'r') as infile:
        reader = csv.reader(infile, delimiter = '\t')
        for row in reader:
            if row[1] == '-1': #no arm found
                G2_loc_list.append('no arm')
            else:
                G2_loc_list.append(row[3])

    for i in range(len(a1_loc_list)):
        a1_loc = a1_loc_list[i]
        A2_loc = A2_loc_list[i]
        g1_loc = g1_loc_list[i]
        G2_loc = G2_loc_list[i]
        if a1_loc == 'no arm' or g1_loc == 'no arm':
            R1_len = 'no_insert'
        else:
            R1_len = int(a1_loc) - int(g1_loc)
        if A2_loc == 'no arm' or G2_loc == 'no arm':
            R2_len = 'no_insert'
        else:
            R2_len = int(A2_loc) - int(G2_loc)
        insert_list.append([R1_len,R2_len])

    with open(outfile_path,'w') as outfile:
        writer = csv.writer(outfile,delimiter = '\t')
        for row in insert_list:
            writer.writerow(row)
def get_arg_parser():
    # parse all files
    parser = ArgumentParser()
    parser.add_argument("-a", "--a1_info_path", dest="a1_info_path",
                        help="Input info file for cutadapt removal of 3 prime adapter of read 1", metavar="a1_info_path")
    parser.add_argument("-g", "--g1_info_path", dest="g1_info_path",
                        help="Input info file for cutadapt removal of 5 prime adapter of read 1", metavar="g1_info_path")
    parser.add_argument("-A", "--A2_info_path", dest="A2_info_path",
                        help="Input info file for cutadapt removal of 3 prime adapter of read 1", metavar="A2_info_path")
    parser.add_argument("-G", "--G2_info_path", dest="G2_info_path",
                        help="Input info file for cutadapt removal of 5 prime adapter of read 1", metavar="G2_info_path")
    parser.add_argument("-o", "--outfile_path", dest="outfile_path",
                        help="Output file for tab-delimited insert lengths", metavar="outfile_path")
    return parser
if __name__ == '__main__':
    main()
