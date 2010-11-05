#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Nick Crawford on 2010-10-28.
Copyright (c) 2010

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses

The author may be contacted at ngcrawford@gmail.com
"""

import os
import sys 
import glob
import shlex
import tempfile
import random
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
        phylip = oneliner2phylip(line) # convert line to phylip
        
        # SETUP TEMP FILE
        temp_in = tempfile.NamedTemporaryFile(suffix='.out',dir='tmp/') #c c
        for line in phylip:
            temp_in.write(line)
        temp_in.seek(0)     # move pointer to beginning of file
        
        # RUN PHYML
        cli = './PhyML_3.0 --input {0}'.format(temp_in.name)
        cli_parts = shlex.split(cli)
        ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE)
        ft.communicate()[0]
        
        # EXTRACT RESULTS
        temp_string = os.path.split(temp_in.name)[1].split('.')[0]
        treefile =  'tmp/%s.out_phyml_tree.txt' % temp_string
        newick = open(treefile,'r').readlines()[0].strip()
        
        # ROOT TREE
        tree = Tree(newick)
        tree.set_outgroup('anoCar2')
        newick = tree.write(format=5)
        print newick.strip()
        
        # CLEAN UP TEMPFILES
        temp_in.close()
        for filename in glob.glob('tmp/%s.*' % (temp_string)) :
            os.remove( filename )
        

def fasttree(args):
    """sends the individual
    cat practice_alignments/3.align.oneliners.txt | ./seqcap.py fasttree
    """
    for line in sys.stdin:
        ft = Popen(['./FastTree','-nt','-quiet'],stdin=PIPE,stderr=PIPE,stdout=PIPE)
        ft.stdin.write(oneliners2phylip(line))
        tree_str = ft.communicate()[0]
        tree_str = tree_str.strip()
        print tree_str
        tree = Tree(tree_str)
        tree.set_outgroup('anoCar2')
        newick = tree.write(format=5)
        print newick.strip()

def muscle(args):
    
    # [('seq1','aaatttc'),('seq2','ggaaaatc')] to be aligned, and
    # aln_from_stdout will be a string containing the fasta alignment, ripe
    # for the parsing.
    # 
    # Note that as of muscle 3.8 -stable is no longer available (stupid Bob
    # Edgar independent scientist :) so you'll need to drop the
    # ever-so-handy -stable flag if you're using a muscle after 3.7
    # 
    # mh = Popen(['muscle','-stable'],stdin=PIPE,stderr=PIPE,stdout=PIPE)
    # mh.stdin.write('\n'.join(['>%s\n%s' % (l,s) for l,s in label_seqs]))
    # aln_from_stdout,output_from_stderr = mh.communicate()
    
    
    mh = Popen(['muscle','-stable'],stdin=PIPE,stderr=PIPE,stdout=PIPE)
    mh.stdin.write('\n'.join(['>%s\n%s' % (l,s) for l,s in label_seqs]))
    alnstr = mh.communicate()[0]
    print ','.join(alnstr.strip().replace('>','').split('\n'))

def aligns2oneliners(args):
    """This code converts phylip alignments to a single tab-delimited line
        of taxa ids and sequence pairs"""
        
    # Storage Dictionaries
    taxa_id_dict = {}
    taxa_seq_dict = {}
    
    filenames = glob.glob(os.path.join('practice_alignments/',"*.phylip"))
    for fin in filenames:
        for count, line in enumerate(open(fin,'r')):    
        
            # Process taxa count and alignment length line
            if count == 0:
                taxa_count, align_len = line.strip().split(' ')
                taxa_count = int(taxa_count.strip()) 
                align_len = int(align_len.strip())
        
            count = count -1    # update count to skip header
        
            if len(line.strip()) == 0: 
                continue # skip blank lines
        
            # Do initial dictionary data population
            if 0 <= count <= taxa_count:
                taxa_name = line[:11].strip()
                sequence = line[10:].strip()
                taxa_id_dict[count] = taxa_name
                taxa_seq_dict[taxa_name] = sequence.replace(' ','')
        
            # Add extra seqeunces to taxa/seq dictionary
            if count > taxa_count:
                dict_id = count % tcount_plus_extra
                taxa_name = taxa_id_dict[dict_id]
                taxa_seq_dict[taxa_name] = taxa_seq_dict[taxa_name] + line.replace(' ','').strip()
        
            tcount_plus_extra = taxa_count + 1
    
        final_line = ""
        for key, value in taxa_seq_dict.iteritems():
            final_line += '%s,%s, ' % (key, value)
        
        final_line = final_line[:-2]    # trim space and comma from end of sequence
        final_line + "\n"
        print final_line




def mpest(args):
    
    def create_control_files(path):
        """explaination"""
        base_path, tree_file = os.path.split(path)

        control_file = '.'.join([tree_file.split('.')[0],'control'])
        fout = open(os.path.join(base_path, control_file),'w')
        # open up the input tree file and get the first tree
        tree_file_contents = open(path,'r').readlines()
        tree = tree_file_contents[0].strip()
        tree = Tree(tree)
        taxa_names =  tree.get_leaf_names()
        #pdb.set_trace()
        taxa_string  = ['{0}\t1\t{0}'.format(taxa) for taxa in taxa_names]
        template_info = {'filename':path,
                        'tree_count':len(tree_file_contents),
                        'numb_taxa':len(taxa_names),
                        'taxa_details':'\n'.join(taxa_string)}
        # needs to go to format()
        template = "%(filename)s\n0\n%(tree_count)s %(numb_taxa)s\n%(taxa_details)s\n0\n" % template_info
        fout.write(template)
        fout.close()
    
    def execute_mpest():
        seed = random.randint(0,1000000)
        cli = "./mpest-bf-v1 rooted.control %s mpest.species.tree" % (seed)
        cli_parts = shlex.split(cli)
        mp = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE)
        mp.communicate()[0]
        
        for line in open('mpest.species.tree'):
            print line
    
    fout = open('rooted.trees','w')
    for line in sys.stdin:
        fout.write(line)
    fout.close()
    
    create_control_files('rooted.trees')
    execute_mpest()
    
    # CLEAN UP TREE/STAT FILES
    for filename in glob.glob('rooted.*') :
        os.remove( filename )
    
    
def main():
    
    # if sys.argv[1] == "fasttree":
    #     fasttree(sys.argv[2:])
    
    if sys.argv[1] == "mpest":
        mpest(sys.argv[2:])

    if sys.argv[1] == "phyml":
        phyml(sys.argv[2:])
            
if __name__ == '__main__':
    main()

