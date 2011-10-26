#!/usr/bin/env python
# encoding: utf-8

import os
import sys
import glob
import argparse

def get_args():
    """Parse sys.argv"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input-dir', required=True, help='The input directory containing the nexus files.')
    args = parser.parse_args()
    return args

def isSeqOK(seq):
    numb_of_qs = 0
    for char in seq:
        if char == "?" or char == "-":
            numb_of_qs +=1
    
    if numb_of_qs == len(seq):
        return False
    if numb_of_qs < len(seq):
        return True

def removeAmbiguousBases(seq):
    """ Converts Ambiguous bases to ?"""
    new_seq = ""
    for char in seq:
        if char not in  ["A", "T", "C", "G"]:
            char = "-"
        new_seq += char
    return new_seq

def nexus2oneliner(fin, name=None):
    
    as_one_liner = False
    in_data = False
    final_line = ''
    align_count = 0
    for count, line in enumerate(fin): 
        
        if "#NEXUS" in line:
            in_data = False
            if as_one_liner == True and len(final_line) != 0:
                sys.stdout.write(final_line + ";")
            if as_one_liner == False and len(final_line) != 0:
                sys.stdout.write(final_line + ";\n")
            final_line = ''
            align_count += 1
        
        if in_data == True:
            
            if ';' not in line:
   
                if len(line.strip()) == 0: continue
                taxon, seq = line.strip().split()
                #seq = removeAmbiguousBases(seq)     # remove any ambiguous bases (e.g,. not A,T,C, or G)
                if isSeqOK(seq) == False: continue  # skip taxa with not data (= all ?????)
    
                if '_' in taxon:
                    # taxon = taxon.split('_')[-1]
                    genus, species = taxon.split('_')[-2:]
                    genus_trunc = '%s%s' % (genus[0].upper(), genus[1:3])
                    species_cap =  '%s%s' % (species[0].upper(), species[1:])
                    taxon = genus_trunc + species_cap
                    taxon =  taxon[:8]
                else:
                    taxon = taxon[0].upper() + taxon[1:8]
                
                if final_line == '':
                    final_line = taxon +','+ seq
                else:
                    final_line += ',' + taxon +','+ seq
                        
        if 'matrix' in line:
            in_data = True
    
    if as_one_liner == True:
        sys.stdout.write(final_line + ";")
    else:
        sys.stdout.write(name + "|" + final_line + ";\n")

def processNexusFiles():
    args = get_args()
    in_dir = os.path.join(args.input_dir, "*.nex*")
    for count, nexus_file in enumerate(glob.glob(in_dir)):
        filename = os.path.split(nexus_file)[-1]
        fileID = os.path.splitext(filename)[0]
        fin = open(nexus_file,'rU')

        nexus2oneliner(fin, name=fileID)
        fin.close()

processNexusFiles()