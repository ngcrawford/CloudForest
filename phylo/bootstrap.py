#!/usr/bin/env python
# encoding: utf-8
"""
bootstrap.py

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





class Bootstrap:
    def __init__(self):
        import os
        import sys 
        sys.path.append('.')
        sys.path.append('bin')
        import glob
        import shlex
        import random
        import getpass
        import itertools
        import tempfile
        import platform
        import ConfigParser
        import numpy as np
        # from tree import Tree
        from subprocess import Popen, PIPE
    
    def bootstrap(self, sample, replicates, start):
        import numpy as np
        from copy import deepcopy
        """bootstrap(initial_sample, replicates=10)

        Create boostrapped replicates of an a numpy array.  Results
        are returned in a python list.

        Parameters
        ----------
        sample : array_like
            An array, any object exposing the array interface, an
            object whose __array__ method returns an array, or any
            (nested) sequence. 

        replicates : int
            The number of bootstrap replicates to produce.

        Examples
        --------
        >>> bootstrap(array([1,2,3,4,5]),5)
        [array([1, 5, 1, 2, 2]),
         array([1, 5, 1, 4, 2]),
         array([2, 2, 2, 4, 4]),
         array([2, 2, 4, 4, 2]),
         array([2, 4, 3, 3, 3])]
        """
        replicates = int(replicates)
        sample_size = len(sample)

        bootstrap_replicates = []
        for boot_rep_num in range(start,start+replicates):
            choices = np.random.random_integers(0, sample_size-1, sample_size)  # generate index array of random choices

            if type(sample) == list:
                boot_rep = []
                for choice in choices:
                    element = deepcopy(sample[choice])
                    boot_rep.append(sample[choice])
            else:
                boot_rep = sample[choices]

            bootstrap_replicates.append(boot_rep)
        return bootstrap_replicates

    def onelinerAlignment2Array(self, line):
        import numpy as np
        seqs = line.split(",")
        label_seqs = zip(seqs[:-1:2],seqs[1::2])
        taxa = []
        bases = []
        for taxon, seq in label_seqs:
            bases.append([seq[count] for count, item in enumerate(seq)])
            taxon = taxon.strip()
            taxa.append(taxon)
        bases = np.array(bases)
        return (taxa, bases)
    
    def array2OnelinerAlignment(self, taxa, bases):
        import itertools
        oneliner = ''
        for count, seq in enumerate(bases):
            oneliner += taxa[count]+","+''.join(itertools.chain(bases[0])) + ","
        return oneliner

    
    def __call__(self, key, line):      
        loci = line.strip().split(';')
        loci = loci[:-1]
        bootreps = int(self.params["bootreps"])
        bootstapped_loci = self.bootstrap(loci, bootreps, 0)                                  # first bootstrap the loci
        for bcount, bootrep in enumerate(bootstapped_loci):                 
            for lcount, locus in enumerate(bootrep):
                taxa, numpy_alignment = self.onelinerAlignment2Array(locus)      # convert loci to 2d arrays
                bases_by_col = numpy_alignment.transpose()                  # transpose so columns are bootstapped   
                bootrep = self.bootstrap(bases_by_col, 1, 0)             # generate one replicate of bootstrapped columns
                bootrep = bootrep[0].transpose()                            # tranpose back to rows of sequences
                oneliner = self.array2OnelinerAlignment(taxa, bootrep)                      # back to oneliner
                yield bcount, oneliner
        
    def phyml(bootrep, line):
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

        # PRINT TREE TO STOUT
        tree = newick.strip()
        yield bootrep, tree


    def reducer(bootstrap_replicate, alignment):
        for item in alignment:
            yield bootstrap_replicate, item  

    def reducer2(bootstrap_replicate, tree):
        yield bootstrap_replicate, tree
    
if __name__ == '__main__':
    # fin = open('/Users/nick/Desktop/Code/BioAWS/phylo/practice_alignments/3.align.oneliners.txt','r')
    # line = fin.readlines()[0]
    # key = 1
    # Boot = Bootstrap()
    # Boot(key,line)
    import dumbo
    job = dumbo.Job()
    job.additer(Bootstrap, reducer, combiner=reducer)
    job.additer(phyml, reducer, combiner=reducer)
    job.run()
    
    