#!/usr/bin/env python
# encoding: utf-8

import os
import sys 
sys.path.append('.')
sys.path.append('bin')
import glob
import shlex
import tempfile
import random
import getpass
from tree import Tree
from subprocess import *

def oneliner2phylip(line):
    seqs = line.split(',')
    label_seqs = zip(seqs[:-1:2],seqs[1::2])
    taxa_count = len(label_seqs)
    seq_length = len(label_seqs[0][1])
    alignment = "%s %s\n" % (taxa_count, seq_length) # add header
    for taxa_name, seq in label_seqs:
        taxa_name = taxa_name.strip()
        alignment += '%-*s%s\n' % (10, taxa_name, seq)
    return alignment


def phyml(args):
    """sends the individual
    cat practice_alignments/3.align.oneliners.txt | ./seqcap.py fasttree
    """
    
    for line in sys.stdin:
        user = getpass.getuser()
        phylip = oneliner2phylip(line) # convert line to phylip
        
        # SETUP TEMP FILE
        temp_in = tempfile.NamedTemporaryFile(suffix='.out', dir="/mnt/var/lib/hadoop/tmp")
        for line in phylip:
            temp_in.write(line)
        temp_in.seek(0)     # move pointer to beginning of file
        
        # RUN PHYML
        # /user/nick/
        #cli = 'bin/./PhyML3linux64 --input %s' % (temp_in.name)
        cli = 'bin/./PhyML3 --input %s' % (temp_in.name)
        cli_parts = shlex.split(cli)
        ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE)
        ft.communicate()[0]
        
        # EXTRACT RESULTS
        temp_string = os.path.split(temp_in.name)[1].split('.')[0]
        treefile =  '/mnt/var/lib/hadoop/tmp/%s.out_phyml_tree.txt' % (temp_string)
        newick = open(treefile,'r').readlines()[0].strip()
        
        # ROOT TREE
        tree = Tree(newick)
        tree.set_outgroup('anoCar2')
        newick = tree.write(format=5)
        print newick.strip()
        
        # CLEAN UP TEMPFILES
        temp_in.close()
        for filename in glob.glob('/mnt/var/lib/hadoop/tmp/%s.*' % (temp_string)) :
            os.remove( filename )

def main():
    phyml(None)

if __name__ == '__main__':
    main()
        