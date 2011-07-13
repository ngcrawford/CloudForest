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


import os
import sys 
sys.path.append('.')
sys.path.append('bin')
import glob
import tempfile
import platform
import itertools
import numpy as np
from copy import copy, deepcopy
from mrjob.job import MRJob
from subprocess import Popen, PIPE

class BootstrapAWS(MRJob):
    
    def configure_options(self):
        super(BootstrapAWS, self).configure_options()

        self.add_passthrough_option(
            '--bootreps', dest='bootreps2run', default=None, type='int',
            help='number of bootstrap replicates to generate')    
        
        self.add_passthrough_option(
            '--gene-trees', action="store_true", dest='gene_trees',
            help='generate gene trees from alignments')
        
        self.add_passthrough_option(
            '--full-analysis', action="store_true",
            help='Run full analysis')
    
    def oneliner2phylip(self, line):
        """Convert one-liner to phylip format."""
        seqs = line.strip(";").split(',')
        label_seqs = zip(seqs[:-1:2],seqs[1::2])
        taxa_count = len(label_seqs)
        seq_length = len(label_seqs[0][1])
        alignment = "%s %s\n" % (taxa_count, seq_length) # add header
        for taxa_name, seq in label_seqs:
            taxa_name = taxa_name.strip()
            alignment += '%-10s%s\n' % (taxa_name, seq)
        return alignment

    def bootstrap(self, sample, replicates, start):
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
        """Convert oneliner to 2d numpy array."""
        seqs = line.split(",")
        label_seqs = zip(seqs[:-1:2],seqs[1::2])
        taxa = []
        bases = []
        for taxon, seq in label_seqs:
            seq_list = copy([seq[count] for count, item in enumerate(seq)])
            seq_array = np.array(seq_list)
            bases.append(seq_array)
            taxon = taxon.strip()
            taxa.append(taxon)
        bases = np.array(bases)
        return (taxa, bases)

    def array2OnelinerAlignment(self, taxa, bases):
        """Convert array of array of taxa and an array of bases to one-liner."""
        oneliner = ''
        for count, seq in enumerate(bases):
            oneliner += taxa[count]+","+''.join(itertools.chain(bases[count])) + ","
        oneliner = oneliner[:-1]
        return oneliner

    def makeReps(self, key, line):
    
        loci = line.strip().split(';')            
        loci = loci[:-1]                                                    # remove empty cell due to trailing ";" 
        bootstapped_loci = bootstrap(loci, self.options.bootreps2run, 1)    # first bootstrap the loci
        for bcount, bootrep in enumerate(bootstapped_loci):                 
            for lcount, locus in enumerate(bootrep):
                taxa, numpy_alignment = onelinerAlignment2Array(locus)      # convert loci to 2d arrays
                bases_by_col = np.column_stack(numpy_alignment)             # transpose so columns are bootstapped   
                shuffled = bootstrap(bases_by_col, 1, 0)                    # generate one replicate of bootstrapped columns
                shuffled = np.column_stack(shuffled[0])                     # tranpose back to rows of sequences
                oneliner = array2OnelinerAlignment(taxa, shuffled)          # back to oneliner
                yield bcount, oneliner
    
    def bootstrapReplicates(self, key, line):
        """ This fuction is slightly different than the standard bootstrapping fuction:
        
        Rather that producing all the bootstrap replicates in one process, it takes a file
        of duplicated datasets equal to the number of replicates. One boostrap replicate is
        produced from the loci within each dataset, and one bootstap replicate is made on
        each locus of the bootstrapped loci. This speeds up the MapReduce algorithm because
        by parallelizing the bootstrapping opperation.
        
        This fuction also keeps the appropriate model of molecular evolution with each
        locus.
        """
        
        loci = line.strip().split(';')
        loci = loci[:-1]                                 # remove empty cell due to trailing ";" 
        bootstapped_loci = self.bootstrap(loci, 1, 0)    # first bootstrap the loci
        for bcount, bootrep in enumerate(bootstapped_loci):                 
            for lcount, locus in enumerate(bootrep):
                model, locus = locus.split("|")
                taxa, numpy_alignment = self.onelinerAlignment2Array(locus)     # convert loci to 2d arrays
                bases_by_col = np.column_stack(numpy_alignment)                 # transpose so columns are bootstapped   
                shuffled = self.bootstrap(bases_by_col, 1, 0)                   # generate one replicate of bootstrapped columns
                shuffled = np.column_stack(shuffled[0])                         # tranpose back to rows of sequences
                oneliner = self.array2OnelinerAlignment(taxa, shuffled)         # back to oneliner
                oneliner = "%s|%s" % (model, locus)
                yield key, oneliner
    
    
    def phyml(self, key, line):
        """Run PhyML on each line in the input file. 
        Parses out evolutionary model if it provided as first word in file
        and """
        
        # USE CORRECT BINARYS
        if platform.system() == 'Darwin':
            system = 'OSX_setup'
            phyml_exe = 'PhyML3OSX'
        else:
            system = 'AWS_setup'
            phyml_exe = 'PhyML3linux32'
        
        if len(line.split("\t")) == 2:
            key, line = line.split("\t")
        
        model = 'HKY85'                     # set default model
        if "|" in line:
            model, line = line.split("|")   # use provided model otherwise default is HKY85  
        phylip = self.oneliner2phylip(line) # convert line to phylip
        
        # SETUP TEMP FILE
        if system == "OSX_setup" and os.path.exists('tmp/') != True:
            os.mkdir('tmp/')
        temp_in = tempfile.NamedTemporaryFile(suffix='.out', dir='tmp/')
        for line in phylip:
            temp_in.write(line)
        temp_in.seek(0)     # move pointer to beginning of file
             
        # RUN PHYML
        cli = '%s/%s/./%s --input=%s --model=%s >/dev/null 2>&1' % (os.getcwd(),\
                'bin', phyml_exe, temp_in.name, model) 
        cli_parts = cli.split()
        ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()
        
        # EXTRACT RESULTS
        temp_string = os.path.split(temp_in.name)[1].split('.')[0]
        treefile =  os.path.join('tmp','%s.out_phyml_tree.txt' % (temp_string))
        newick = open(treefile,'r').readlines()[0].strip()
        
        yield key, newick
    
    def mrAIC(self, key, alignment):
        """ Run mr-aic.pl on the each one-liner in the input file."""
        
        # USE CORRECT BINARYS
        if platform.system() == 'Darwin':
            system = 'OSX_setup'
            phyml_exe = 'PhyML3OSX'
            mpest_exe = 'mpestOSX'
        else:
            system = 'AWS_setup'
            phyml_exe = 'PhyML3linux32'
            mpest_exe = 'mpestEC2'
        
        alignment = alignment.split("\t")[-1]     # weirdness parsing key, value
        phylip = self.oneliner2phylip(alignment)  # convert line to phylip
        
        # SETUP TEMPFILE
        if system == "OSX_setup" and os.path.exists('tmp/') != True:
            os.mkdir('tmp/')
            
        temp_in = tempfile.NamedTemporaryFile(suffix='.out', dir='tmp/')
        for line in phylip:
            temp_in.write(line)
        temp_in.seek(0) 
        
        temp_dir = os.path.dirname(temp_in.name)
        
        # EXECUTE MR-AIC (AIC output only)
        cli = "%s/%s/./mraic_mod.pl --infile=%s --output_dir=%s >/dev/null 2>&1" % \
            (os.getcwd(), 'bin', temp_in.name, temp_dir)
        cli_parts = cli.split()
        ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()
        
        # PARSE FILE NAMES IN TMP/ TO GET MODEL
        aic_dir_string = "tmp/%s.AICc-*tre*" % os.path.basename(temp_in.name)
        aic_file = glob.glob("tmp/*AICc-*tre*")[0]
        aic_model = aic_file.split('.')[-2].split("-")[-1]  
        oneliner = "%s|%s" % (aic_model, alignment)
        
        yield 1, oneliner  # give everything the same key so it can be reduced to a 
                           # 'oneliner' suitable for bootstrapping
      
    def duplicateOneliners(self, key, line):
        """Take lines and duplicate them the number of times
        specified by the --bootreps flag."""
        # line = line.split("\t")[1].strip("\"") 
        reps = self.options.bootreps2run
        while reps != 0: 
            yield reps, line
            reps -= 1
    
    # REDUCER FUNCTIONS     
    def basicReducer(self, key, line):
        """Yield lines from a generator."""
        for item in line:
            yield key, item   
        
    def lines2Oneliner(self, key, line):
        """Convert multiple alignments with the same key
        to a concatenated oneliner with the same key"""
        concatenated_line = "".join(line)
        yield 1, concatenated_line
        
    def lines2Genetrees(self, key, line):
        """Convert multiple alignments with the same key
        to a concatenated oneliner with the same key"""
        for item in line:
            yield key, line
        
    def steps(self):
        
        if self.options.full_analysis == True:
            return [self.mr(self.mrAIC, self.lines2Oneliner),
                    self.mr(self.duplicateOneliners, self.basicReducer),
                    self.mr(self.bootstrapReplicates, self.basicReducer),
                    self.mr(self.phyml, self.basicReducer)] 
                            
        if self.options.gene_trees == True:
            return [self.mr(self.mrAIC, self.lines2Genetrees)]
            
        # else:
        #     return [self.mr(self.makeReps, self.boot_reducer), 
        #             self.mr(self.phyml, self.boot_reducer),]
    

if __name__ == '__main__':
    BootstrapAWS.run()