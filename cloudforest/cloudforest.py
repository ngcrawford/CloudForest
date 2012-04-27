import os
import re
import glob
import platform
import tempfile
import numpy as np
from subprocess import Popen, PIPE

import pdb


class Process():
    """ """
    def __init__(self):
        pass

    def oneliner_to_phylip(self, line):
        """Convert one-liner to phylip format."""
        seqs = line.strip(";").split(',')
        label_seqs = zip(seqs[:-1:2], seqs[1::2])
        taxa_count = len(label_seqs)
        seq_length = len(label_seqs[0][1])
        # add header
        alignment = "%s %s\n" % (taxa_count, seq_length)
        for taxa_name, seq in label_seqs:
            taxa_name = taxa_name.strip()
            alignment += '%-10s%s\n' % (taxa_name, seq)
        return alignment

    def get_bootstraps(self, sample, replicates=1, return_choices=False):
        """Create boostrapped replicates of an a numpy array.

        Returns list

        """
        if not isinstance(sample, np.ndarray):
            try:
                sample = np.array(sample)
            except:
                raise TypeError("bootstrap() input must be list or numpy.ndarray")
        #replicates = int(replicates)
        size = len(sample)
        if replicates == 1:
            choices = np.random.random_integers(0, size - 1, size)
            if not return_choices:
                return sample[choices].tolist()
            else:
                return sample[choices].tolist(), choices

        else:
            return [sample[np.random.random_integers(0, size - 1, size)].tolist()
                            for i in xrange(replicates)]

    def oneliner_to_array(self, line):
        """Convert oneliner to 2d numpy array.

        Returns tuple(list, array)

        """
        seqs = line.split(",")
        label_seqs = zip(seqs[:-1:2], seqs[1::2])
        taxa, bases = [], []
        for taxon, seq in label_seqs:
            bases.append(list(seq))
            taxa.append(taxon.strip())
        return taxa, np.array(bases)

    def array_to_oneliner(self, taxa, bases):
        """Convert array of array of taxa and an array of bases to one-liner.

        Returns string

        """
        # join is faster and cleaner than concatenate
        # list comp is faster than for
        # enumerate gets index of taxon from taxa
        oneliner = ','.join([','.join([taxa[count], seq.tostring()]) for count, seq in enumerate(bases)])
        return oneliner

    def make_bootstrap_replicates(self, key, line):
        loci = line.strip().split(';')
        loci = loci[:-1]
        # first, bootstrap across loci
        multilocus_bstrap = self.get_bootstraps(loci)
        # second, bootstrap bases within loci
        for locus in multilocus_bstrap:
            # split keys from alignments
            args_dict, locus = self.split_one_liner(locus, return_locus=True)
            # convert alignments to arrays
            taxa, numpy_alignment = self.oneliner_to_array(locus)
            # transpose so we bootstrap by columns
            bases_by_col = np.column_stack(numpy_alignment)
            # bootstrap by columns
            shuffled = self.get_bootstraps(bases_by_col)
            # transpose back to rows of sequences
            shuffled = np.column_stack(shuffled)
            # convert back to oneliner
            oneliner = self.array_to_oneliner(taxa, shuffled)
            oneliner = "%s:%s" % (self.make_tree_name(args_dict), oneliner)
            yield key, oneliner

    def make_tree_name(self, args_dict):
        """Converts dictionary of arguement into a sorted string with the
        format:  """
        # join is faster and cleaner than concatenate
        # list comp is faster than for
        name = ','.join(['='.join([pair[0], pair[1]]) for pair in sorted(args_dict.items())])
        return name

    def process_stats_file(self, fin):
        """"given an input phyml stats file, return the log-likelihood of the tree"""
        result = None
        regex = re.compile("Log-likelihood:\s+(.+)")
        for line in fin:
            result = regex.search(line)
            if result:
                break
        if result is None:
            raise ValueError("No Log-likelihood found")
        return float(result.groups()[0])

    def split_oneliner(self, line, args_dict={}, return_locus=False):
        """From a oneline, split locus/taxon info into dict and locus into string

        Returns dict or tuple(dict, locus)

        """
        args, locus = line.split(":")
        args = args.split(',')
        for item in args:
            dict_key, value = item.split("=")
            args_dict[dict_key] = value
        if not return_locus:
            return args_dict
        else:
            return args_dict, locus

        
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
        
        # GET MODEL INFO FROM FILE IF AVAILABLE
        args_dict = {}
        args_dict['model'] = 'HKY85'                     # set default model
        if ":" in line:
            args_dict = self.processOnelinerData(line, args_dict)
            phylip = self.oneliner2phylip(line.split(":")[-1]) # convert line to phylip
        
        else:
            phylip = self.oneliner2phylip(line) # convert line to phylip
        
        # SETUP TEMP FILE
        temp_in = tempfile.NamedTemporaryFile(suffix='.out', dir='tmp/')
        for line in phylip:
            temp_in.write(line)
        temp_in.seek(0)     # move pointer to beginning of file
             
        # RUN PHYML
        cli = '%s --input=%s --model=%s >/dev/null 2>&1' % \
            (os.path.join("bin", phyml_exe), temp_in.name, args_dict['model']) 
            
        cli_parts = cli.split()
        ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()
        
        # EXTRACT RESULTS AND FORMAT AS NEXUS TREES
        temp_string = os.path.split(temp_in.name)[1].split('.')[0]
        
        treefile =  os.path.join('tmp','%s.out_phyml_tree.txt' % (temp_string))
        tree = open(treefile,'r').readlines()[0].strip().strip("\"")
        
        statsfile = os.path.join('tmp','%s.out_phyml_stats.txt' % (temp_string))        
        lnL = self.processStatsFile(open(statsfile,'r'))
        args_dict['lnL'] = lnL
        
        if self.options.gene_trees == True and self.options.mraic_opt == None:
            tree = "tree '" + self.makeTreeName(args_dict) + "' = [&U] " + tree
            yield key, tree
        
        else:
            yield key, tree
    
    def mrAIC(self, key, line, bin='bin'):
        """ Run mr-aic.pl on the each one-liner in the input file."""     
        oneliner = line.split("\t")[-1]     # weirdness parsing key, value

        # PARSE EXTRA ALIGNMENT INFO
        args_dict = {}
        if ":" in line:
            args_dict = self.processOnelinerData(oneliner, args_dict)
            phylip = self.oneliner2phylip(oneliner.split(":")[-1]) # convert line to phylip
                
        else:        
            phylip = self.oneliner2phylip(oneliner)  # convert line to phylip
                        
        temp_in = tempfile.NamedTemporaryFile(suffix='.out', dir='tmp/')
        for line in phylip:
            temp_in.write(line)
        temp_in.seek(0)          
        temp_dir = os.path.dirname(temp_in.name)
                
        # EXECUTE MR-AIC (AIC output only)            
        cli = "%s --infile=%s --output_dir=%s >/dev/null 2>&1" % \
            (os.path.join(bin, 'mraic_mod.pl'), temp_in.name, temp_dir)
        #pdb.set_trace()        
        cli_parts = cli.split()
        ft = Popen(cli_parts, stdin=PIPE, stderr=PIPE, stdout=PIPE).communicate()
                
        # PARSE FILE NAMES IN TMP/ TO GET MODEL          
        aic_file = glob.glob("tmp/%s.AICc-*tre*" % os.path.basename(temp_in.name))[0]
        aic_model = aic_file.split('.')[-2].split("-")[-1]
        args_dict['model'] = aic_model
        oneliner = "%s:%s" % (self.makeTreeName(args_dict), oneliner.split(":")[-1])       
        '''
        if self.options.gene_trees == True:
            aic_fin = open(aic_file,'rU')
            for line in aic_fin:
                yield key, "tree '%s' = [&U] %s" % (self.makeTreeName(args_dict), line.strip())
        else:
            yield 1, oneliner  # give everything the same key so it can be reduced to a 
                               # 'oneliner' suitable for bootstrapping
        '''
      
    def duplicate_oneliner(self, key, line, bootreps=500):
        """Take lines and duplicate them the number of times
        specified by the --bootreps flag."""
        # line = line.split("\t")[1].strip("\"") 
        try:
            reps = self.options.bootreps2run
        except AttributeError:
            reps = bootreps
        for i in reversed(xrange(1, reps + 1)):
            yield i, line

          
    def lines2Oneliner(self, key, line):
        """Convert multiple alignments with the same key
        to a concatenated oneliner with the same key"""
        concatenated_line = "".join(line)
        yield 1, concatenated_line
