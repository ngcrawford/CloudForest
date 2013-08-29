#!/usr/bin/env python
# encoding: utf-8

import re
import os
import sys
import glob
import argparse
from Bio import AlignIO

filename = '/Users/ngcrawford/Desktop/UCEs/birds/after-trimming-changes/nexus-rename/uce-914.nex'


def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()

    parser.add_argument('input', 
                        type=str,
                        default=sys.stdin,
                        help='Path to input. (default is STDIN)')

    parser.add_argument('output',
                        type=str,
                        default=sys.stdout,
                        help='Path to output. (default is STDOUT)')

    args = parser.parse_args()
    return args


def truncate_sample_ID(ID):

    genus, species = ID.split("_") 
    trunc_name = "{0}{1}{2}"
    
    if len(genus) > 3:
        genus = genus[:3]

    if len(species) > 3:

        num_list = re.findall('\d', species)

        if len(num_list) is not 0:

            num = str(num_list[0])
            idx = 7 - len(num)

            species = species[:idx].title() + num
        
        else:
            species = species[:7].title()


    trunc_name = trunc_name.format(genus, species, "")

    return trunc_name

def process_file(input, output):
    sample_ids = {}

    with open(input,'rU') as fin:

        alignment = AlignIO.read(fin, "nexus")
   
        for record in alignment:
            t_id = truncate_sample_ID(record.id)
            
            sample_ids[record.id] = t_id
            record.id = t_id

        if len(sample_ids.keys()) != len(set(sample_ids.values())):
            print 'Sample truncation produces duplicate IDs.'
            sys.exit()

        else:
            fout = filename.strip(".nex") + ".phylip"
            AlignIO.write(alignment, open(output, 'w'), "phylip")

    return sample_ids


def make_sample_ids_2_trunc_output(d, args):
    ids = d.keys()
    ids.sort()

    if os.path.isdir(args.output) is True:
        fout = os.path.join(args.output,'ids_2_trunc_ids.txt')
    else:
        fout = os.path.split(args.output)[0]
        fout = os.path.join(fout,'ids_2_trunc_ids.txt')
    print fout

    fout = open(fout,'w')
    fout.write('Name\tID\n')
    for i in ids:
        line = "{}\t{}\n".format(i,d[i])
        fout.write(line)


def main(args):

    if os.path.isdir(args.input) is True:
        
        for i in  glob.glob(args.input + "/*.nex"):
            
            o = os.path.split(i)[1].replace('nex', 'phylip')
            o = os.path.join(args.output, o)
            process_file(i, o)            

    else:
        p = process_file(args.input, args.output)
        make_sample_ids_2_trunc_output(p, args)




if __name__ == '__main__':
    args = get_args()
    main(args)







