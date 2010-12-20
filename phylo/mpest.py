#!/usr/bin/env python
# encoding: utf-8


import os
import sys 
sys.path.append('.')
sys.path.append('bin')
sys.path.append('osx.phylo.tar')
import glob
import shlex
import random
import tempfile
import platform
import ConfigParser
from tree import Tree
from subprocess import *

def setup_configuration(setupfile):
    # globals settings (e.g, config paths..)
    config = ConfigParser.RawConfigParser()
    config.read(setupfile)
    return config

def mpest(config, system):
    
    def create_control_files(rooted_trees_path):
        """path = directory to rooted trees"""
        control_file = os.path.join(config.get(system,'tempdir'),'rooted.control')
        fout = open(os.path.join(control_file),'w')
        
        # open up the input tree file and get the first tree
        tree_file_contents = open(rooted_trees_path,'r').readlines()
        tree = tree_file_contents[0].strip()
        tree = Tree(tree)
        taxa_names =  tree.get_leaf_names()
        taxa_string  = ['%s\t1\t%s' % (taxa, taxa) for taxa in taxa_names]
        template_info = {'filename':rooted_trees_path,
                        'tree_count':len(tree_file_contents),
                        'numb_taxa':len(taxa_names),
                        'taxa_details':'\n'.join(taxa_string)}
        template = "%(filename)s\n0\n%(tree_count)s %(numb_taxa)s\n%(taxa_details)s\n0\n" % template_info
        fout.write(template)
        fout.close()
    
    def execute_mpest(rooted_trees):
        
        # SETUP FILE PATHS AND SEED VALUE
        control_file = os.path.join(config.get(system, 'tempdir'), 'rooted.control')
        seed = random.randint(0,1000000)
        mpest_trees = os.path.join(config.get(system, 'tempdir'), 'mpest.species.tree')
       
        # RUN MPEST
        # ./mpest control.file int_seed output_file
        cli = "%s/./%s %s %s %s" \
            % (config.get(system, 'binarydir'), config.get(system, 'mpest_exe'), control_file, seed, mpest_trees)
        cli_parts = shlex.split(cli)
        mp = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE)
        mp.communicate()[0]
        
        if platform.system() == 'Darwin':
            for line in open('tmp/rooted.trees.mpest','r'):
                print line
        else:
            for line in open('tmp/rooted.trees.mpestout','r'):
                print line
    
    # MAIN CODE 
    rooted_trees_path = os.path.join(config.get(system,'tempdir'),'rooted.trees')
    rooted_trees = open(rooted_trees_path,'w')
    
    for line in sys.stdin:
        rooted_trees.write(line)
    rooted_trees.close()
    
    create_control_files(rooted_trees_path)
    execute_mpest(rooted_trees)
    
    # CLEAN UP TREE/STAT FILES
    tempfiles = os.path.join(config.get(system,'tempdir'),'rooted.*')
    for filename in glob.glob(tempfiles):
        os.remove(filename)

def main():
    
    if platform.system() == 'Darwin':
        config = setup_configuration('bin/setup.cfg')
        mpest(config, system = 'OSX_setup')
    else:
        config = setup_configuration('bin/setup.cfg')
        mpest(config, system = 'AWS_setup')

    
if __name__ == '__main__':
    main()