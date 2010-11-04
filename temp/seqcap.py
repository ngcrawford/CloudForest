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
from subprocess import *

def oneliners2phylip(line):
    seqs = line.split(',')
    label_seqs = zip(seqs[:-1:2],seqs[1::2])
    taxa_count = len(label_seqs)
    seq_length = len(label_seqs[0][1])
    alignment = "%s %s\n" % (taxa_count, seq_length) # add header
    for taxa_name, seq in label_seqs:
        taxa_name = taxa_name.strip()
        alignment += '%-*s%s\n' % (10, taxa_name, seq)
    return alignment

def fasttree(args):
    """sends the individual 
    
    cat practice_alignments/1.phylip | ./FastTree -nt -quiet
    
    """
    for line in sys.stdin:
        
        
        ft = Popen(['./FastTree','-nt','-quiet'],stdin=PIPE,stderr=PIPE,stdout=PIPE)
        ft.stdin.write(oneliners2phylip(line))
        alnstr = ft.communicate()[0]
        print ','.join(alnstr.strip().replace('>','').split('\n'))
        
    
    # for line in sys.stdin:
    #     path = line.strip()
    #     # pathname = os.path.dirname(sys.argv[0])
    #     # pathname = os.path.abspath(pathname)
    #     # print pathname
    #     # this will probably need to be changed for aws
    #     command = './PhyML -nt -quiet %s' % (path)
    #     command = shlex.split(command)
    #     print subprocess.Popen(command, stderr=subprocess.STDOUT, \
    #                     stdout=subprocess.PIPE).communicate()# [0]
        

        
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
    
    for count, line in enumerate(sys.stdin):    
        
        # Process taxa count and alignment length line
        if count == 0:
            taxa_count, align_len = line.strip().split(' ')
            taxa_count = int(taxa_count.strip()) 
            align_len = int(align_len.strip())
        
        count = count -1    # update count to skip header
        
        if len(line.strip()) == 0: continue # skip blank lines
        
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

def fastas2oneliners(args):
    pass
    
def main():
    """
    Can be run as follows:
    
    
    cat practice_alignments/1.phylip | ./seqcap.py convert-file | ./seqcap.py fasttree
    
    
    
    """
    
    if sys.argv[1] == "fasttree":
        fasttree(sys.argv[2:])
    if sys.argv[1] == 'convert-file':
        aligns2oneliners(sys.argv[2:])

if __name__ == '__main__':
    main()

