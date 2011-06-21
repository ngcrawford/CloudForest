#!/usr/bin/env python
# encoding: utf-8
"""
test.py

Created by Nick Crawford on 2010-12-23.
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

def mapper(key, line):
    import os
    import sys 
    sys.path.append('.')
    sys.path.append('bin')
    import glob
    import shlex
    import random
    import getpass
    import tempfile
    import platform
    import ConfigParser
    from tree import Tree
    from subprocess import Popen, PIPE
    
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
    
    if platform.system() == 'Darwin':
        system = 'OSX_setup'
        phyml_exe = 'PhyML3OSX'
        mpest_exe = 'mpestOSX'
    else:
        system = 'AWS_setup'
        phyml_exe = 'PhyML3linux32'
        mpest_exe = 'mpestEC2'
        
        
    user = getpass.getuser()
    phylip = oneliner2phylip(line) # convert line to phylip
    
    # SETUP TEMP FILE
    temp_in = tempfile.NamedTemporaryFile(suffix='.out', dir='tmp/')
    for line in phylip:
        temp_in.write(line)
    temp_in.seek(0)     # move pointer to beginning of file
    
    # RUN PHYML
    cli = '%s/%s/./%s --input %s' % (os.getcwd(), \
        'bin/', phyml_exe, temp_in.name)
        
    cli_parts = shlex.split(cli)
    ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE)
    ft.communicate()[0]
    
    # EXTRACT RESULTS
    temp_string = os.path.split(temp_in.name)[1].split('.')[0]
    treefile =  os.path.join('tmp','%s.out_phyml_tree.txt' % (temp_string))
    newick = open(treefile,'r').readlines()[0].strip()
    
    # ROOT TREE
    tree = Tree(newick)
    tree.set_outgroup('anoCar2')
    newick = tree.write(format=5)
    
    #CLEAN UP TEMPFILES
    temp_in.close()
    for filename in glob.glob('%s/%s.*' % ('tmp', temp_string)) :
        os.remove(filename)
    
    # PRINT TREE TO STOUT
    tree = newick.strip()
    yield 1, tree

def reducer(key, trees):
    import os
    import sys 
    sys.path.append('.')
    sys.path.append('bin')
    import glob
    import shlex
    import random
    import getpass
    import tempfile
    import platform
    import ConfigParser
    from tree import Tree
    from subprocess import Popen, PIPE
    
    if platform.system() == 'Darwin':
        system = 'OSX_setup'
        phyml_exe = 'PhyML3OSX'
        mpest_exe = 'mpestOSX'
    else:
        system = 'AWS_setup'
        phyml_exe = 'PhyML3linux32'
        mpest_exe = 'mpestEC2'
    
    
    def create_control_files(rooted_trees_path):
        """path = directory to rooted trees"""
        control_file = os.path.join('tmp','rooted.control')
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
        
        def get_mpest_tree(trees):
            mpest_tree = ''
            for line in reversed(trees):
                if 'mpestTree' in line:
                    pre, post = line.strip('\n').split('=')
                    post = post.strip(' ')
                    mpest_tree += post
            return mpest_tree
    
        def get_species_dict(control):
            """docstring for get_species_dict"""
            control_file = open(control, 'rU')
            [control_file.readline() for i in xrange(3)]
            species = []
            for line in control_file.readlines():
                ls = line.strip().split('\t')
                if ls[0] == '0':
                    break
                else:
                    species.append(ls[0])
            return dict(zip(xrange(1,len(species)+1),species))

        def convert_taxa_labels(species_tree, control_file):
            species_dict = get_species_dict(control_file)
            species_tree = Tree(species_tree)
            for lf in species_tree.iter_leaves():
                current = lf.name
                new = species_dict[int(current)]
                lf.name = new
            
            species_tree = species_tree.write(format=3)
            return species_tree
        
        # SETUP FILE PATHS AND SEED VALUE
        control_file = os.path.join('tmp', 'rooted.control')
        seed = random.randint(0,1000000)
        mpest_trees = os.path.join('tmp', 'mpest.species.trees')
       
        # RUN MPEST
        # ./mpest control.file int_seed output_file
        cli = "%s/./%s %s %s %s" % ('bin', mpest_exe, control_file, seed, mpest_trees)
        cli_parts = shlex.split(cli)
        mp = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE)
        mp.communicate()[0]
        
        trees = []
        if platform.system() == 'Darwin':
            for line in open('tmp/mpest.species.trees','r'):
                 trees.append(line)
        else:
            for line in open('tmp/mpest.species.trees','r'):
                trees.append(line)
                
        species_tree = get_mpest_tree(trees)
        species_tree = convert_taxa_labels(species_tree, control_file)
        return species_tree
    
    # MAIN CODE 
    rooted_trees_path = ('tmp/rooted.trees.mpest')
    rooted_trees = open(rooted_trees_path,'w')
    
    for line in trees:
        rooted_trees.write(line+'\n')
    rooted_trees.close()
    
    create_control_files(rooted_trees_path)
    species_tree = execute_mpest(rooted_trees)
    
    # CLEAN UP TREE/STAT FILES
    tempfiles = os.path.join('tmp','rooted.*')
    for filename in glob.glob(tempfiles):
         os.remove(filename)
    
    yield key, species_tree

if __name__ == "__main__":
    import dumbo
    dumbo.run(mapper, reducer)


