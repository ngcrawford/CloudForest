#!/usr/bin/env python
# encoding: utf-8

import os
import sys 
import glob
import shlex
import tempfile
import random
from tree import Tree
from subprocess import *

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
    mpest(None)
    
if __name__ == '__main__':
    main()