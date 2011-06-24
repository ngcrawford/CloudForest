#!/usr/bin/env python
# encoding: utf-8

"""
This script reads in 

Command line Usage: cat practice_alignments/*.phylip | ./prep_aligns.py > oneliners.txt
 
 Note: you must provide multiple alignments with wildcard. 

Things to check for...

1.) missing taxa, read in first 5 alignments
2.) """

import sys

def aligns2oneliners():

    def printAlignment(count, taxa_seq_dict):
        final_line = ""
        for key, value in taxa_seq_dict.iteritems():
            final_line += '%s,%s,' % (key, value)
        final_line = final_line[:-1]
        final_line = str(count) + "\t" + final_line + "\n"
        sys.stdout.write(final_line)
    
    # STORAGE DICTIONARIES AND VARIABLES 
    taxa_id_dict = {}
    taxa_seq_dict = {}
    count = 0
    alignment_count = 0
    for line in sys.stdin:    
        line = line.strip()         # remove extra whitespace 
        
        # IDENTIFY START OF ALIGNMENT
        if len(line.split()) == 2:
            taxa_count, align_len = line.split(' ')
            
            # DOUBLE CHECK THAT LINE IS START OF ALIGNMENT
            if taxa_count.isdigit() == False: continue
            if align_len.isdigit() == False: continue
            
            # PRINT COMPLETE ALIGNMENT AS LINE  
            if len(taxa_seq_dict) != 0:
                printAlignment(alignment_count, taxa_seq_dict)
            
            # UPDATE VALUES
            taxa_count = int(taxa_count.strip())
            align_len = int(align_len.strip())
            taxa_id_dict = {}
            taxa_seq_dict = {}
            count = 0
            alignment_count += 1
        
        # SKIP BLANK LINES
        if len(line.strip()) == 0: 
            continue
            
        # INITIALIZE DICTS WITH TAXA ID'S AND INITIAL SEQS
        if 0 < count <= taxa_count:
            taxa_name = line[:11].strip()
            sequence = line[10:].strip()
            taxa_id_dict[count] = taxa_name
            taxa_seq_dict[taxa_name] = sequence.replace(' ','')
        
        # ADD ADDITIONAL LINES TO ALIGNMENT 
        if count > taxa_count:
            dict_id = (count % taxa_count)
            if (count % taxa_count) == 0:
                dict_id = max(taxa_id_dict.keys())
            taxa_name = taxa_id_dict[dict_id]
            taxa_seq_dict[taxa_name] = taxa_seq_dict[taxa_name] + line.replace(' ','').strip()
             
        count += 1
        
    # PRINT COMPLETE ALIGNMENT AS LINE  
    if len(taxa_seq_dict) != 0:
        printAlignment(alignment_count, taxa_seq_dict)

aligns2oneliners()