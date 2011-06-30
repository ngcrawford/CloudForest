#!/usr/bin/env python
# encoding: utf-8
import sys

def nexus2oneliner():
    
    as_one_liner = True
    
    in_data = False
    final_line = ''
    align_count = 0
    for count, line in enumerate(sys.stdin): 
        
        if "#NEXUS" in line:
            in_data = False
            if as_one_liner == True and len(final_line) != 0:
                sys.stdout.write(final_line + ";")
            if as_one_liner == False and len(final_line) != 0:
                sys.stdout.write(final_line + "\n")
            final_line = ''
            align_count += 1
        
        if in_data == True:
            if ';' not in line:
                
                taxon, seq = line.strip().split()
                if '_' in taxon:
                    genus, species = taxon.split('_')[-2:]
                    genus_trunc = '%s%s' % (genus[0].upper(), genus[1:3])
                    species_cap =  '%s%s' % (species[0].upper(), species[1:])
                    taxon = genus_trunc + species_cap
                    taxon =  taxon[:8]
                else:
                    taxon = taxon[0].upper() + taxon[1:8]
                if final_line == '':
                    if as_one_liner == True:
                        final_line = taxon +','+ seq
                    else:
                        final_line = str(align_count) + "\t"+ taxon +','+ seq
                else:
                    final_line += ',' + taxon +','+ seq
                        
        if 'matrix' in line:
            in_data = True
    
    if as_one_liner == True:
        sys.stdout.write(final_line + ";")
    else:
        sys.stdout.write(final_line + "\n")
            
nexus2oneliner()