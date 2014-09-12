#!/usr/bin/env python
# encoding: utf-8

"""
File: cloudforest/core.py

Created by Nicholas G. Crawford and Brant C. Faircloth
Copyright (c) 2011-2012 Nicholas G. Crawford and Brant C. Faircloth. All rights reserved.

Description: Core functions for cloudforest

"""

import os
import re
import sys
import shutil
import argparse
import platform
import tempfile
import subprocess
import numpy as np

#import pdb


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def format_oneliner_from_dict(taxa_seq_dict, name, model=None):
    assert len(taxa_seq_dict) > 0, ValueError("Phylip file has no content")
    seq = ['%s,%s' % (k, v) for k, v in taxa_seq_dict.iteritems()]
    if not model:
        return "chrm=%s:%s;" % (name.strip(), ','.join(seq))
    else:
        return "chrm=%s,model=%s:%s;" % (name.strip(), model, ','.join(seq))


def phylip_to_oneliner(phylip, locus, model=None):
    taxa_id_dict, taxa_seq_dict = {}, {}
    count = 0
    for lineno, line in enumerate(phylip.split('\n')):
        # skip blank lines
        if len(line.strip()) == 0:
            continue
        if lineno == 0:
            taxa_count, align_len = [int(i) for i in line.strip().split()]
        # INITIALIZE DICTS WITH TAXA ID'S AND INITIAL SEQS
        if 0 < count <= taxa_count:
            # clean up errant spaces
            sp = [val for val in line.split(' ') if val != '']
            name = sp[0]
            sequence = ''.join(sp[1:])
            taxa_id_dict[count] = name
            taxa_seq_dict[name] = sequence.replace(' ', '')
        # ADD ADDITIONAL LINES TO ALIGNMENT
        if count > taxa_count:
            num = count % taxa_count
            if num == 0:
                num = max(taxa_id_dict.keys())
            name = taxa_id_dict[num]
            taxa_seq_dict[name] += line.strip().replace(' ', '')
        count += 1
    oneliner = format_oneliner_from_dict(taxa_seq_dict, locus, model)
    return oneliner


def oneliner_to_phylip(line):
    """Convert one-liner to phylip format."""

    # Problem: Phylip has a 10 char limit that is not accounted for. (NGC)
    # "Strict Phylip expects the first character state to appear on Column 11
    #  for each and every sequence, no ifs, and, or buts."
    # http://www.phylo.org/tools/phylip.html

    # Remove additional info (e.g., chrm,model,etc.) from line (fixed above problem)
    if ':' in line:
        line = line.split(":")[-1]
    seqs = line.strip(";\n").split(',')
    label_seqs = zip(seqs[:-1:2], seqs[1::2])
    #print label_seqs
    taxa_count = len(label_seqs)
    seq_length = len(label_seqs[0][1])
    # pad all names to length of longest name + 1 space
    max_name_length = max([len(val) for val in seqs[:-1:2]]) + 1
    # add header
    header = "%s %s\n" % (taxa_count, seq_length)
    alignment = '\n'.join(['%s%s' % (i[0].ljust(max_name_length), i[1]) for i in label_seqs])
    return header + alignment


def get_bootstraps(sample, replicates=1, return_choices=False):
    """Create boostrapped replicates of an a numpy array.

    Returns list

    """
    # make sure we instantiate a Random number generator that doesn't
    # fall victim to the forking issues w/ multiprocessing
    #
    # Create generator w/ seed.  None as input draws from /dev/urandom
    rng = np.random.RandomState()
    if not isinstance(sample, np.ndarray):
        try:
            sample = np.array(sample)
        except:
            raise TypeError("bootstrap() input must be list or numpy.ndarray")
    size = len(sample)
    if replicates == 1:
        choices = rng.randint(0, size, size)
        # for debugging
        #sys.stdout.write(str(choices))
        if not return_choices:
            return sample[choices].tolist()
        else:
            return sample[choices].tolist(), choices

    else:
        return [sample[rng.randint(0, size, size)].tolist()
                        for i in xrange(replicates)]


def oneliner_to_array(line):
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


def array_to_oneliner(taxa, bases):
    """Convert array of array of taxa and an array of bases to one-liner.

    Returns string

    """
    # join is faster and cleaner than concatenate
    # list comp is faster than for
    # enumerate gets index of taxon from taxa
    oneliner = ','.join([','.join([taxa[count], seq.tostring()]) for count, seq in enumerate(bases)])
    return oneliner


def make_tree_name(args_dict):
    """Converts dictionary of arguement into a sorted string with the
    format:  """
    # join is faster and cleaner than concatenate
    # list comp is faster than for
    name = ','.join(['='.join([pair[0], pair[1]]) for pair in sorted(args_dict.items())])
    return name


def split_oneliner(line, args_dict=None, default_model=False):
    """From a oneliner, split locus/taxon info into dict and locus into string

    Returns dict or tuple(dict, locus)

    """
    if not args_dict:
        args_dict = {}
    if ':' in line:
        args, locus = line.split(":")
        args = args.split(',')
        for item in args:
            dict_key, value = item.split("=")
            args_dict[dict_key] = value
    else:
        locus = line
    if default_model:
        if 'model' not in args_dict.keys():
            args_dict['model'] = 'GTR'
    return args_dict, locus


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


class Process():
    """ """
    def __init__(self):
        pass

    def get_bootstrap_replicates(self, key, line):
        loci = line.strip().split(';')
        loci = loci[:-1]
        # first, bootstrap across loci
        multilocus_bstrap = get_bootstraps(loci)
        # second, bootstrap bases within loci
        for locus in multilocus_bstrap:
            # split keys from alignments
            args_dict, locus = split_oneliner(locus)
            # convert alignments to arrays
            taxa, numpy_alignment = oneliner_to_array(locus)
            # transpose so we bootstrap by columns
            bases_by_col = np.column_stack(numpy_alignment)
            # bootstrap by columns
            shuffled = get_bootstraps(bases_by_col)
            # transpose back to rows of sequences
            shuffled = np.column_stack(shuffled)
            # convert back to oneliner
            oneliner = array_to_oneliner(taxa, shuffled)
            oneliner = "%s:%s" % (make_tree_name(args_dict), oneliner)
            yield key, oneliner

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

    def get_genetrees(self, key, line, pth='bin', genetrees=True, no_model=False):
        """Compute genetrees using model in oneliner. Parses out evolutionary
         model if provided as first word in file.  Otherwise, runs GTR."""
        # TODO:  Why is this here?
        #    EMR often adds extra tab chars in front of lines.
        #    this fixes that particular edge case. (NGC)

        if len(line.split("\t")) == 2:
            key, line = line.split("\t")
        args_dict, locus = split_oneliner(line, default_model=True)

        phylip = oneliner_to_phylip(locus)
        phyml = Phyml(phylip, pth)
        # run phyml.  if no model, defaults to GTR
        # TOOD: Why do we need LnL?
        #   For comparing the quality of topologies
        args_dict['lnL'], tree = phyml.run(args_dict['model'])
        try:
            gtrees = self.options.gene_trees
            no_model = self.options.mraic_opt
        except:
            gtrees = genetrees
        if gtrees == True and no_model == False:
            yield key, "tree '%s' = [&U] %s" % (make_tree_name(args_dict), tree)
        else:
            yield key, tree

    def get_genetrees_and_models(self, key, line, pth='bin', genetrees=True):
        """Compute genetrees from the best fitting substitution model and return
        generator of genetrees and/or oneliner with models integrated"""
        # TODO: move into separate line_cleaner function?
        oneliner = line.split("\t")[-1].strip('\n')
        args_dict, locus = split_oneliner(line)
        phylip = oneliner_to_phylip(locus)
        phyml = Phyml(phylip, pth)
        model, tree = phyml.best_aicc_model_and_tree()
        args_dict['model'] = model
        oneliner = "%s:%s" % (make_tree_name(args_dict), oneliner.split(":")[-1])
        try:
            gtrees = self.options.gene_trees
        except AttributeError:
            gtrees = genetrees
        if gtrees == True:
            yield key, "tree '%s' = [&U] %s" % (make_tree_name(args_dict), tree)
        else:
            yield 1, oneliner

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

    def concatenate_oneliners(self, key, line):
        """Convert multiple alignments with the same key
        to a concatenated oneliner with the same key"""
        concatenated_line = "".join(line)
        yield 1, concatenated_line


class Phyml:
    """Use phyml to generate trees or help select models"""
    def __init__(self, phylip, pth='bin', temp_dir=None, exe=None,
                       starting_tree=None, constraint_tree=None):
        self.cwd = os.getcwd()
        if not temp_dir:
            working = tempfile.mkdtemp()
        else:
            # generate a tempdir in which we'll work
            working = tempfile.mkdtemp(dir=temp_dir)
        self.working = os.path.abspath(working)
        # if we get a file for phylip var, put in tempdir
        if os.path.exists(phylip):
            phylip = os.path.abspath(os.path.expanduser(phylip))
            self.phylip = os.path.join(self.working, os.path.basename(phylip))
            shutil.copyfile(
                    phylip,
                    self.phylip
                    )
        # if we get a string for a file, write to tempdir/tempfile
        elif type(phylip) == str:
            self.phylip = self._string_2_tempfile(string=phylip, suffix='phylip')
        else:
            raise TypeError("Input must be a phylip file or a phylip-formatted string")

        if starting_tree != None:
           self.starting_tree = self._string_2_tempfile(string=starting_tree, suffix='txt')

        if constraint_tree != None and starting_tree == None:
            raise Exception("A constraint tree also requires a starting tree.")

        if constraint_tree != None:
           self.constraint_tree = self._string_2_tempfile(string=constraint_tree, suffix='txt')


        self.phyml3 = self._get_phyml_pth(pth, exe)
        # container for model selection results
        self.aicc_results = None
        # model parameters taken from John Nylander's excellent mr_aic.pl
        # http://www.abc.se/~nylander/
        self.models = {
                'JC69': "+\nM\nM\nM\nM\nM\nR\nY\n",
                'JC69I': "+\nM\nM\nM\nM\nM\nV\nY\nR\nY\n",
                'JC69G': "+\nM\nM\nM\nM\nM\nY\n",
                'JC69IG': "+\nM\nM\nM\nM\nM\nV\nY\nY\n",
                'F81': "+\nM\nM\nM\nM\nM\nM\nM\nF\nR\nY\n",
                'F81I': "+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nR\nY\n",
                'F81G': "+\nM\nM\nM\nM\nM\nM\nM\nF\nY\n",
                'F81IG': "+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nY\n",
                'K2P': "+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nY\n",
                'K2PI': "+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nV\nY\nY\n",
                'K2PG': "+\nM\nM\nM\nM\nM\nM\nT\nY\nY\n",
                'K2PIG': "+\nM\nM\nM\nM\nM\nM\nT\nY\nV\nY\nY\n",
                'HKY': "+\nF\nT\nY\nR\nY\n",
                'HKYI': "+\nF\nT\nY\nR\nV\nY\nY\n",
                'HKYG': "+\nF\nT\nY\nY\n",
                'HKYIG': "+\nF\nT\nY\nV\nY\nY\n",
                'SYM': "+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nY\n",
                'SYMI': "+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
                'SYMG': "+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nY\n",
                'SYMIG': "+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n",
                'GTR': "+\nM\nM\nM\nF\nR\nY\n",
                'GTRI': "+\nM\nM\nM\nF\nR\nV\nY\nY\n",
                'GTRG': "+\nM\nM\nM\nF\nY\n",
                'GTRIG': "+\nM\nM\nM\nF\nV\nY\nY\n"
            }
        self.numparams = {
                'JC69': 0,
                'JC69I': 1,
                'JC69G': 1,
                'JC69IG': 2,
                'F81': 3,
                'F81I': 4,
                'F81G': 4,
                'F81IG': 5,
                'K2P': 1,
                'K2PI': 2,
                'K2PG': 2,
                'K2PIG': 3,
                'HKY': 4,
                'HKYI': 5,
                'HKYG': 5,
                'HKYIG': 6,
                'SYM': 5,
                'SYMI': 6,
                'SYMG': 6,
                'SYMIG': 7,
                'GTR': 8,
                'GTRI': 9,
                'GTRG': 9,
                'GTRIG': 10
            }
        # compile regex for LnL once
        self.ll = re.compile("Log-likelihood:\s+(.+)")
        self.dim = re.compile("\s*(\d+)\s+(\d+)")

    def __del__(self):
        """Cleanup empty dir - not guaranteed to be called during execution,
        only at close"""
        shutil.rmtree(self.working)

    def __str__(self):
        return "Phyml object of %s" % os.path.basename(self.phylip)

    def __repr__(self):
        return "<Phyml object of %s at %s>" % (os.path.basename(self.phylip), hex(id(self)))

    def _string_2_tempfile(self,string,suffix):
        """[Private] Create tempfile with data from string."""
        fd, file_path = tempfile.mkstemp(dir=self.working, suffix=suffix, text=True)
        os.write(fd, string)
        os.close(fd)
        return file_path

    def _get_phyml_pth(self, pth, exe):
        """[Private] Get path to phyml"""
        if not exe:
            # USE CORRECT BINARYS
            if platform.system() == 'Darwin':
                #system = 'OSX_setup'
                phyml3 = os.path.join(pth, 'PhyML3OSX')
            else:
                #system = 'AWS_setup'
                phyml3 = os.path.join(pth, 'PhyML3linux32')
        else:
            phyml3 = os.path.join(pth, exe)
        return os.path.abspath(os.path.expanduser(phyml3))

    def _get_taxon_and_char_data(self):
        """[Private] Parse the first line of a phylip file and return nchar and ntax"""
        # get taxon and character data for file
        first_line = open(self.phylip, 'rU').readline()
        self.taxa, self.nchar = [int(val) for val in self.dim.search(first_line).groups()]
        # calculate to keep results comparable to mr_aic.pl
        self.nbranch = (2 * self.taxa) - 3

    def _get_log_like(self, statfile, regex, phylip):
        """"[Private] Given an input phyml stats file, return the log-likelihood of the tree"""
        result = None
        for line in open(statfile, 'rU'):
            result = regex.search(line)
            if result:
                break
        if result is None:
            raise ValueError("No Log-likelihood found")
        return float(result.groups()[0])

    def _get_tree(self, treefile):
        """[Private] Return the tree produced for a given subs. model"""
        tree = None
        tree = open(treefile, 'rU').read().strip()
        if tree is None or tree == '':
            raise ValueError("No tree found")
        return tree

    def _compute_aicc(self, model, loglik, count_branches=True):
        """[Private] Compute aicc given loglik

        AICc: -2lnL + 2K + 2K(K+1)/n-K-1.

        we're not worried about AIC, since AICc > AIC with larger samples, and
        BIC doesn't seem as sensible beacause it mixes ML and Bayesian paradigms
        """
        if count_branches:
            params = self.numparams[model] + self.nbranch
        else:
            params = self.numparams[model]
        return -2. * loglik + 2. * params + ((2. * params * (params + 1.)) / (self.nchar - params - 1.))

    def _runner(self, phylip, model):
        """[Private] Given alignment and model, run phyml"""
        statfile, treefile = [''.join([phylip, ext]) for ext in ['_phyml_stats.txt', '_phyml_tree.txt']]
        # Delete existing files from previous runs
        try:
            [os.remove(f) for f in [statfile, treefile]]
        except OSError, e:
            # if no files, skip
            if e.errno == 2:
                pass
        template = "%s\n%s" % (phylip, model)
        cli = [self.phyml3]
        subprocess.Popen(
                cli,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE
            ).communicate(input=template)
        return statfile, treefile

    def _best_model_runner(self, return_aicc=False):
        """Compute the best model for an alignment using AICc"""
        self.lnl_results = {}
        self.aicc_results = {}
        # because of phyml, move to working dir
        os.chdir(self.working)
        self._get_taxon_and_char_data()
        for model_name, model in self.models.iteritems():
            phylip = os.path.basename(self.phylip)
            # self._runner generates out phyml locus and model template
            statfile, treefile = self._runner(phylip, model)
            lnl = self._get_log_like(statfile, self.ll, phylip)
            self.lnl_results[model_name] = lnl
            try:
                aicc = self._compute_aicc(model_name, lnl)
            except ZeroDivisionError:
                txt = "An alignment is shorter than necessary to compute AICc. " + \
                    "Ensure alignments are longer than (# taxa + 10)."
                raise IOError(txt)
            tree = self._get_tree(treefile)
            self.aicc_results[aicc] = [model_name, tree]
        # move back to cwd
        os.chdir(self.cwd)

    def best_aicc_model(self):
        """Return best model; Do not recompute if models have been run"""
        if not self.aicc_results:
            self._best_model_runner()
        best = min(self.aicc_results.keys())
        return self.aicc_results[best][0]

    def best_aicc_tree(self):
        """Return best lnl and tree; Do not recompute if models have been run"""
        if not self.aicc_results:
            self._best_model_runner()
        best = min(self.aicc_results.keys())
        name = self.aicc_results[best][0]
        return str(self.lnl_results[name]), self.aicc_results[best][1]

    def best_aicc_model_and_tree(self):
        """Return best model and tree; Do not recompute if models have been run"""
        if not self.aicc_results:
            self._best_model_runner()
        best = min(self.aicc_results.keys())
        return self.aicc_results[best][0], self.aicc_results[best][1]

    def aicc_model_results(self):
        """Return all model results; Do not recompute if models have been run"""
        if not self.aicc_results:
            self._best_model_runner()
        return self.aicc_results

    def run(self, model='GTR'):
        """
        Compute a phyml tree for the alignment given and subs model.  Tree output
        here is identical to that from best_aicc_tree, except that we jump right to
        optimum subs model here, without iterating across models.  Input to phyml
        is model structure from self.models, rather than model test, since phyml
        requires custom input for anything other than several standard models.
        """
        # because of phyml, move to working dir
        os.chdir(self.working)
        phylip = os.path.basename(self.phylip)
        # run phyml
        try:
            statfile, treefile = self._runner(phylip, self.models[model])
        except KeyError:
            raise KeyError("You must use a valid model: %s" % (','.join(sorted(self.models.keys()))))
        # get tree
        tree = self._get_tree(treefile)
        # get LnL
        lnl = self._get_log_like(statfile, self.ll, phylip)
        # move back to cwd
        os.chdir(self.cwd)
        return str(lnl), tree

