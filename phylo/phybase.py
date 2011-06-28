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
import glob
import argparse
import dendropy
import rpy2.robjects as robjects

def interface():
    """ create commandline interface for script"""
    p = argparse.ArgumentParser()
    
    p.add_argument('--input-file','-i',
        help='Path to input file.')
    p.add_argument('--outgroup','-o',
        help='Name of outgroup.')
    p.add_argument('--genetrees','-g', action='store_true',
        help='Name of outgroup.')
    p.add_argument('--bootstraps','-b', action='store_true',
        help='Name of outgroup.')
    p.add_argument('--print-taxa','-p', action='store_true',
        help='Print all the taxon names in the first tree')

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

    if args.genetrees == False and args.bootstraps == False:
        print "You must select either --genetrees or --bootstraps"
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
    new_tree += ";"
    return new_tree

def cleanPhyMLTree(tree):
    tree = tree.strip()
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
    star_sptree = robjects.r['star.sptree'](trees, species_spname, species_taxaname, \
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
    
def consensus(tree_list):
    trees = dendropy.TreeList()
    for tree in tree_list:
        t = dendropy.Tree()
        t.read_from_string(tree,'newick')
        tree.append(tree)

    con_tree = trees.consensus(min_freq=0.5)
    print(con_tree.as_string('newick'))

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
        if count == 0: taxa = getTaxa(tree.strip())
        if bootreps.has_key(key) != True: bootreps[key] = [tree]
        else: bootreps[key].append(tree)  
    
    steac_trees = []
    star_trees = []
    for count, key in enumerate(bootreps.keys()):
        trees = bootreps[key]
        star_tree, steac_tree = phybase(trees, 'HomSapie', taxa)
        steac_trees.append(steac_tree)
        star_trees.append(star_tree)
        print 'processed', count
    
    print 'steac_trees'
    for tree in steac_trees:
        print tree
    
    print 'star_trees'
    for tree in star_trees:
       print tree

def parseGenetrees(args):
    """docstring for parse"""
    fin = open(args.input_file,'rU')
    taxa = []
    trees = []
    for count, line in enumerate(fin):
        tree = line.strip()
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
    fin = open(args.input_file,'rU')
    line = fin.readline() 
    if len(line.split('\t')) == 2:
        line = line.split('\t')[-1]
    taxa = getTaxa(line)
    taxa.sort()
    for count, taxon in enumerate(taxa):
        print taxon
    print count, 'total taxa.'
    
def main():
    args = interface()
    
    if args.print_taxa == True:
        print_taxa(args)
        sys.exit()
    
    if args.bootstraps == True:
        parseBootreps(args)
    
    if args.genetrees == True:
        parseGenetrees(args)
    
if __name__ == '__main__':
    main()
    pass