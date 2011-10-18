#!/usr/bin/env python
# encoding: utf-8
"""
Phybase.py

Created by Nicholas Crawford and Brant C. Faircloth Copyright (c) 2010 Nicholas Crawford and
Brant C. Faircloth. All rights reserved.

Parses gene trees from directories of tree files. Submits the trees to R and runs phybase to
generate species trees. Both Steac and Star trees are generated. Seperate functions are included
to parse NJ trees from Paup and ML trees from PhyMl as the tree formats are slightly different.

Dependencies:

   dendropy - http://packages.python.org/DendroPy/tutorial/index.html
   rpy2 - http://rpy.sourceforge.net/rpy2.html 
   phybase - http://cars.desu.edu/faculty/lliu/research/phybase.html

Future directions:

   - add commandline (optparse)
   - add NJtrees fuction
   - Parsing of trees from Garli 
   - Bootstrapping with Phybase
   - Nexus output of Steac and Star trees 
      - also png/svg output? 
"""

import os
import sys
import gzip
import glob
from pylab import *
import argparse
import dendropy
import rpy2.robjects as robjects

def interface():
    """ create commandline interface for script"""
    
    description="""
    
Phybase.py calculates species trees from gene trees.
Phybase.py is actually a wrapper script that runs an R package
of the same name (Liu & Yu 2010).

Dependancies:
-------------

Phybase R package: basic functions for phylogenetic analysis
    Installation: 
        At the R prompt type 'install.packages("Phybase")'
    Website: http://cran.r-project.org/web/packages/phybase/index.html

DendroPy: phylogenetic computing library:
    Installation: sudo easy_install -U dendropy
    Website: http://packages.python.org/DendroPy/

Rpy2: simple and efficient access to R from Python
    Installation: sudo easy_install -U rpy2
    Website: http://rpy.sourceforge.net/rpy2.html
    
Argparse: present in python 2.7 and later (I think)

References:
-----------
Liu, L., & Yu, L. (2010). Phybase: an R package for species tree 
analysis. Bioinformatics (Oxford, England). 
doi:10.1093/bioinformatics/btq062 
"""
     
    p = argparse.ArgumentParser(description,)
    
    p.add_argument('--input-file','-i',
        help='Path to input file.')
    p.add_argument('--outgroup','-o',
        help='Name of outgroup.')
    p.add_argument('--genetrees','-g', action='store_true',
        help='Set this flag if the input is genetrees and you want the species tree.')
    p.add_argument('--bootstraps','-b', action='store_true',
        help='Set this flag if the input is bootstraps and you want multiple species trees.')
    p.add_argument('--print-taxa','-p', action='store_true',
        help='Print all the taxon names in the first tree.')
    p.add_argument('--sorted','-s', action='store_true',
        help='Set this flag if your bootstraps are sorted by key.')

    args = p.parse_args()
    
    # check options for errors, etc.
    if args.input_file == None:
        print "Input directory required."
        print "Type 'python phybase.py -h' for details" 
        sys.exit()
    
    if args.genetrees == True and args.bootstraps == True:
        print "You must pick either genetrees or bootstraps,"
        print "but not both."
        print "Type 'python phybase.py -h' for details" 
        sys.exit()

    if args.genetrees == False and args.bootstraps == False and args.sorted == False and args.print_taxa == False:
        print "You must select either --genetrees, --bootstraps, --sorted, or --print-taxa"
        print "Type 'python phybase.py -h' for details" 
        sys.exit()
    
    if args.print_taxa != True and args.outgroup == None:
        print 'You much define an outgroup'
        print "Type 'python phybase.py -h' for details" 
        sys.exit()

    return args

def cleanPhybaseTree(tree):
    tree.strip("\"")
    tree = tree.split("\"")
    return tree[1]

def remove_branch_lengths(str_newick_tree):
    tree = dendropy.Tree()
    tree.read_from_string(str_newick_tree,'newick')
    for nd in tree:
        nd.label = None
    tree = tree.as_newick_string()
    return tree
    
    # para_pos = None
    # in_para = False
    # in_colon = False
    # final_tree = ""
    # slices = []
    # for count, char in enumerate(str_newick_tree):
    #     if char == ")":
    #         para_pos = count + 1
    #         in_para = True
    #         in_colon = False
    #     
    #     if char == "(":
    #         in_para = False
    #             
    #     if char == ":" and in_para == True and in_colon == False:
    #         slices.append([para_pos, count])
    #     
    #     if char == ":":
    #         in_colon = True  
    # 
    # final_tree = ''
    # slices = array(slices)
    # slices = slices.flatten()
    # for count, item in enumerate(slices):
    #     next_count = count + 1
    #     if count == 0:
    #         final_tree += (str_newick_tree[:item])
    #     
    #     if next_count + 1 > slices.shape[0]:
    #         start = slices[count]
    #         final_tree += str_newick_tree[item:]
    #         break
    #         
    #     if count % 2 != 0:
    #         start = slices[count]
    #         stop = slices[next_count]          
    #         final_tree += str_newick_tree[start:stop]
    # 
    # return final_tree
                

def branch_lengths_2_decimals(str_newick_tree):
    """replaces branch lengths in scientific notation with decimals"""
    colon_s = 0
    comma_back_paren_s = 0
    num = ''
    new_tree = ''
    for count, char in enumerate(str_newick_tree):
        if char == ':': 
            colon_s = count
            continue

        if char in (')',','): 
            comma_back_paren_s = 1
            num = '%f' % float(num)
            new_tree += ":" + num 
            colon_s = 0
            num = ''

        if colon_s != 0:
            num = num + char

        if colon_s == 0:
            new_tree += char
    new_tree = new_tree.strip('\'').strip('\"').strip('\'') + ";"
           
    return new_tree

def cleanPhyMLTree(tree):
    tree = tree.strip()
    tree = remove_branch_lengths(tree)
    tree = dendropy.Tree.get_from_string(tree, 'newick') # removes support values (=total hack)
    tree = tree.as_newick_string()
    tree = branch_lengths_2_decimals(tree)    # converts numbers in sci. notation
                                              # to decimals (e.g., 1e-22 = 0.000000)
    return tree
    
def phybase(trees, outgroup, all_taxa):
    """ generate Steac and Star trees from a list of trees. Requires Phybase and rpy2."""
    robjects.r['library']('phybase')
    trees = robjects.StrVector(trees)
    species_taxaname = robjects.StrVector(all_taxa)
    species_spname = species_taxaname                               # list of species in current tree
    matrix_size = len(species_taxaname)
    species_structure = robjects.r['diag'](1,matrix_size,matrix_size)
    star_sptree = robjects.r['star.sptree'](trees, species_spname, species_taxaname,\
                                            species_structure,outgroup,'nj')
    steac_sptree = robjects.r['steac.sptree'](trees, species_spname, species_taxaname,\
                                            species_structure,outgroup,'nj')
    
    star_sptree = cleanPhybaseTree(str(star_sptree))
    steac_sptree = cleanPhybaseTree(str(steac_sptree))
    
    return (star_sptree, steac_sptree)
    
def phyMLTrees(directory):
    """ get and format for phybase() PhyML 3.0 trees in a directory."""

    # this 'taxa_labels portion' corrects a bug where some trees are 
    # missing expected taxa. This should be fixed such that it simply
    # checks for an expected number of taxa and no duplicates

    tree_list = []
    for tree_file in glob.glob(os.path.join(directory,'*tree.txt')):
        for tree in open(tree_file,'r'):
            tree = cleanTree(tree)
            tree_list.append(tree)
    return tree_list
    
def consensus(tree_list, min_freq=0.5):
    trees = dendropy.TreeList()
    for tree in tree_list:
        t = dendropy.Tree()
        t.read_from_string(tree,'newick')
        trees.append(t)    
    con_tree = trees.consensus(min_freq)
    return con_tree.as_string('newick')

def getTaxa(tree):
    t = dendropy.Tree()
    t.read_from_string(tree, 'newick', preserve_underscores=True)
    taxa = t.taxon_set.labels()
    return taxa
    
def parseBootreps(args):  
    bootreps = {}
    fin = open(args.input_file,'rU')
    taxa = []
    for count, line in enumerate(fin):
        key, tree = line.split("\t")
        tree = cleanPhyMLTree(tree)
        if count == 0: 
            taxa = getTaxa(tree)
        if bootreps.has_key(key) != True: bootreps[key] = [tree]
        else: bootreps[key].append(tree)  
    
    star_fout = os.path.splitext(args.input_file)[0]
    star_fout += '.star.trees'
    
    steac_fout = os.path.splitext(args.input_file)[0]
    steac_fout += '.steac.trees'
    
    steac_trees = []
    star_trees = []
    print 'starting bootreps'
    for count, key in enumerate(bootreps.keys()):
        trees = bootreps[key]
        star_tree, steac_tree = phybase(trees, args.outgroup, taxa)
        steac_trees.append(steac_tree)
        star_trees.append(star_tree)
        star_fout.write(star_tree)
        steac_fout.write(steac_tree)
        print 'processed', count

    steac_consensus = consensus(steac_trees)
    star_consensus = consensus(star_trees)

    template =  """#NEXUS\n
begin trees;\n
tree 'STARConsensus' = %s\n
tree 'STEACConsensus' = %s\n
end;\n""" % (star_consensus, steac_consensus)

    steac_star_cons_out = os.path.splitext(args.input_file)[0]
    steac_star_cons_out = os.path.join(steac_star_cons_out,'steac_star.consensus.trees')
    steac_star_cons_out.write(template)
    

def parseSortedBootreps(args):
    """docstring for parseBootreps"""
    
    taxa = None
    line_id = None
    trees = []
    steac_trees = []
    star_trees = []
    
    # SETUP OUTPUT FILES
    star_file = os.path.splitext(args.input_file)[0]
    star_file += '.star.trees'
    star_fout = open(star_file,'w')
    
    steac_file = os.path.splitext(args.input_file)[0]
    steac_file += '.steac.trees'
    steac_fout = open(steac_file,'w')
    
    # LOOP THROUGH SORTED FILE
    if os.path.splitext(args.input_file)[-1] == '.gz':
        fin = gzip.open(args.input_file, 'r')
    else:
        fin = open(args.input_file,'rU')

    print 'processing bootreps'
    for count, line in enumerate(fin):
        if count == 0 and os.path.splitext(args.input_file)[-1] == '.gz': continue # skip first line in gzip file
        bootrep, tree = line.split("\t")
        bootrep = int(bootrep)
        tree = tree.strip(";")
        tree = cleanPhyMLTree(tree)

        if count == 1: 
            taxa = getTaxa(tree.strip())

        if line_id == None:
            line_id = bootrep
            trees.append(tree)
            continue
        
        if bootrep != line_id:
            star_tree, steac_tree = phybase(trees, args.outgroup, taxa)
            
            print 'processed', len(trees), \
            'trees of bootstrap replicate', line_id
            
            star_fout.write(star_tree.strip()+"\n")
            steac_fout.write(steac_tree.strip()+"\n")
            
            steac_trees.append(steac_tree)
            star_trees.append(star_tree)
            
            line_id = bootrep
            trees = []
            
        trees.append(tree)
    
    # Process final trees:
    star_tree, steac_tree = phybase(trees, args.outgroup, taxa) 
    print 'processed', len(trees), \
    'trees of bootstrap replicate', line_id
    
    steac_trees.append(steac_tree)
    star_trees.append(star_tree)
    
    # Clean up files 
    star_fout.close()
    steac_fout.close()
    
    # WRITE CONSENSUS TREES
    steac_consensus = consensus(steac_trees)
    star_consensus = consensus(star_trees)

    template =  """#NEXUS
begin trees;
tree 'STARConsensus' = %s
tree 'STEACConsensus' = %s
end;""" % (star_consensus, steac_consensus)

    steac_star_cons_out = os.path.splitext(args.input_file)[0]
    steac_star_cons_out += '.steac_star.consensus.trees'
    steac_star_cons_out = open(steac_star_cons_out,'w')
    steac_star_cons_out.write(template)
    steac_star_cons_out.close()
        
def parseGenetrees(args):
    """docstring for parse"""
    is_nexus = False
    if args.input_file.endswith('.nex') == True or args.input_file.endswith('.nexus') == True:
        is_nexus = True
    
    fin = open(args.input_file,'rU')
    taxa = []
    trees = []

    for count, line in enumerate(fin):
        tree = 'tree'
        if is_nexus == True:  # this doesn't work properly. 
            if len(line.strip().split("=")) == 2:
                tree = line.strip().split("=")[-1]
                tree = tree.strip()
            else: continue
        else: tree = line.strip()
        tree = tree.strip(";")
        tree = cleanPhyMLTree(tree)
        trees.append(tree)
        if count == 0: 
            taxa = getTaxa(tree)
            if args.outgroup not in taxa:
                print args.outgroup, 'not in taxa:', taxa
                sys.exit() 

    star_tree, steac_tree = phybase(trees, args.outgroup, taxa)

    template =  """#NEXUS
begin trees;
tree 'STAR' = %s
tree 'STEAC' = %s
end;""" % (star_tree, steac_tree)
    print template

def print_taxa(args):
    if os.path.splitext(args.input_file)[-1] == '.gz':
        fin = gzip.open(args.input_file, 'r')
    else:
        fin = open(args.input_file,'rU')
    line = fin.readline()
    
    for count, line in enumerate(fin): 
        if count == 1:
            tree = line.split('\t')[-1]
            taxa = getTaxa(tree)
            taxa.sort()
            taxa.pop(0) # remove
            for taxon_count, taxon in enumerate(taxa, 1):
                print taxon
            print "---------\n", taxon_count, 'total taxa.'
            break
    
def main():
    args = interface()
    
    if args.print_taxa == True:
        print_taxa(args)
        sys.exit()
    
    if args.bootstraps == True:
        parseBootreps(args)
    
    if args.genetrees == True:
        parseGenetrees(args)
    
    if args.sorted == True:
        parseSortedBootreps(args)
    
if __name__ == '__main__':
    main()
    pass