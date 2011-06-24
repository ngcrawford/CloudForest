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
import optparse
import dendropy
import rpy2.robjects as robjects

def interface():
    """ create commandline interface for script"""
    p = optparse.OptionParser()
    
    p.add_option('--input-file','-i', type='string',
        help='Path to input file.')
    p.add_option('--output-dir','-o', type='string',
        help='Path to output directory.')
    p.add_option('--outgroup','-g', type='string',
        help='Name of outgroup.')
    p.add_option('--taxa','-t',type='string',
        help='List of all possible taxa in phylogy.\n \
              Allows for trees with missing taxa to be included.')
    
    options, arg = p.parse_args()
    
    # check options for errors, etc.
    if options.input_file == None:
        print "Input directory required."
        print "Type 'python phybase.py -h' for details" 
        sys.exit()
    
    if options.output_dir == None:
        print "No output directory supplied."
        print "Type 'python phybase.py -h' for details" 
        sys.exit()
    
    if options.taxa == None:
        print "No list of taxa supplied."
        print "Type 'python phybase.py -h' for details"
        
    return options

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
    
def consensus():

    trees = dendropy.TreeList()
    for tree_file in ['pythonidae.mb.run1.t',
            'pythonidae.mb.run2.t',
            'pythonidae.mb.run3.t',
            'pythonidae.mb.run4.t']:
        trees.read_from_path(
                tree_file,
                'nexus',
                tree_offset=20)
    con_tree = trees.consensus(min_freq=0.95)
    print(con_tree.as_string('newick'))


def parseEMR(file):
    bootreps = {}
    fin = open(file,'rU')
    for count, line in enumerate(fin):
        key, tree = line.split("\t")
        tree = cleanPhyMLTree(tree)
        if bootreps.has_key(key) != True: bootreps[key] = [tree]
        else: bootreps[key].append(tree)  
    
    taxa = ['mus_mus', 'cal_jac', 'oto_gar', 'pon_abe', 'pap_ham', \
            'rhe_mac', 'hom_sap', 'gor_gor', 'nom_leu', 'pan_tro']
    steac_trees = []
    star_trees = []
    for count, key in enumerate(bootreps.keys()):
        trees = bootreps[key]
        star_tree, steac_tree = phybase(trees, 'mus_mus', taxa)
        steac_trees.append(steac_tree)
        star_trees.append(star_tree)
        print 'processed', count
    
    print 'steac_trees'
    for tree in steac_trees:
        print tree
    
    print 'star_trees'
    for tree in star_trees:
        print tree
    
def main():
    
    options = interface()
    taxa = options.taxa.strip().split()
    
    trees = []
    for count, tree in enumerate(open(options.input_file, 'r')):
        tree = cleanPhyMLTree(tree)
        trees.append(tree)
    
    star_tree, steac_tree = phybase(trees, options.outgroup, taxa)
    
    basename = os.path.basename(options.input_file)
    star_out =  os.path.join(options.output_dir, basename + '.star.tre')
    steac_out = os.path.join(options.output_dir, basename + '.steac.tre')
    
    star_out = open(star_out,'w')
    steac_out = open(steac_out,'w')
    
    star_out.write(star_tree)
    steac_out.write(steac_tree)
        
    star_out.close()
    steac_out.close()
    
    
if __name__ == '__main__':
    #parseEMR("real_data/primate_bootreps_50")
    pass