#!/usr/bin/env python
# encoding: utf-8
"""
File: model_selection.py
Author: Brant Faircloth

Created by Brant Faircloth on 26 April 2012 20:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import shutil
import platform
import tempfile
import subprocess

#import pdb


class Phyml:
    """Use phyml to generate trees or help select models"""
    def __init__(self, phylip, pth='bin', exe=None):
        self.phylip = os.path.abspath(os.path.expanduser(phylip))
        self.phyml3 = self._get_phyml_pth(pth, exe)
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

    def _get_phyml_pth(self, pth, exe):
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
        return phyml3

    def _get_taxon_and_char_data(self, regex):
        """parse the first line of a phylip file and return nchar and ntax"""
        # get taxon and character data for file
        first_line = open(self.phylip, 'rU').readline()
        self.taxa, self.nchar = [int(val) for val in regex.search(first_line).groups()]
        # calculate to keep results comparable to mr_aic.pl
        self.nbranch = (2 * self.taxa) - 3

    def _get_log_like(self, statfile, regex, phylip):
        """"given an input phyml stats file, return the log-likelihood of the tree"""
        result = None
        #stats = ''.join([phylip, '_phyml_stats.txt'])
        for line in open(statfile, 'rU'):
            result = regex.search(line)
            if result:
                break
        if result is None:
            raise ValueError("No Log-likelihood found")
        return float(result.groups()[0])

    def _get_tree(self, treefile):
        """return the tree produced for a given subs. model"""
        tree = None
        #treefile = ''.join([phylip, '_phyml_tree.txt'])
        tree = open(treefile, 'rU').read().strip()
        if tree is None or tree == '':
            raise ValueError("No tree found")
        return tree

    def _compute_aicc(self, model, loglik, count_branches=True):
        """
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
        """given alignment and model, run phyml"""
        template = "%s\n%s" % (phylip, model)
        cli = [self.phyml3]
        subprocess.Popen(
                cli,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE
            ).communicate(input=template)
        statfile, treefile = [''.join([phylip, ext]) for ext in ['_phyml_stats.txt', '_phyml_tree.txt']]
        return statfile, treefile

    def best_model(self, aicc=False):
        """compute the best model for an alignment using AICc"""
        # compile this once
        ll_regex = re.compile("Log-likelihood:\s+(.+)")
        dim_regex = re.compile("\s*(\d+)\s+(\d+)")
        results = {}
        # get current dir
        cwd = os.getcwd()
        # create tempdir to hold phyml output on
        # per locus level
        working = tempfile.mkdtemp(dir='tmp')
        os.chdir(working)
        # copy phylip to working - phyml does like long paths
        # copy file to tempdir once - work in tempdir - delete temp
        # out files after computing aicc.
        shutil.copyfile(self.phylip, os.path.basename(self.phylip))
        self._get_taxon_and_char_data(dim_regex)
        for model_name, model in self.models.iteritems():
            phylip = os.path.basename(self.phylip)
            # self._runner generates out phyml locus and model template
            statfile, treefile = self._runner(phylip, model)
            lnl = self._get_log_like(statfile, ll_regex, phylip)
            aicc = self._compute_aicc(model_name, lnl)
            tree = self._get_tree(treefile)
            results[aicc] = [model_name, tree]
            # need to delete files individually, because names are same on next iter
            [os.remove(f) for f in [statfile, treefile]]
        os.chdir(cwd)
        shutil.rmtree(working)
        # best is min(AICc)
        best = min(results.keys())
        if not aicc:
            return results[best][0], results[best][1]
        else:
            return results

    def run(self, model='GTR'):
        """
        compute the tree for an alignment given a subs model.  tree output here is identical
        to that from best_model(return_tree=True), except that we jump right to optimum here
        """
        cwd = os.getcwd()
        # create tempdir to hold phyml output and move there
        working = tempfile.mkdtemp(dir='tmp')
        os.chdir(working)
        # copy over phylip file for locus
        shutil.copyfile(self.phylip, os.path.basename(self.phylip))
        phylip = os.path.basename(self.phylip)
        # run phyml
        statfile, treefile = self._runner(phylip, self.models[model])
        # get tree
        tree = self._get_tree(treefile)
        # move back to cwd and delete temp dir
        os.chdir(cwd)
        shutil.rmtree(working)
        return tree


if __name__ == '__main__':
    import pdb
    phyml = Phyml('tests/alignments/phylip_primates/chr1_1036.phylip', '../../binaries')
    best = phyml.best_model()
    print "best model is ", best
    tree = phyml.run(best)
    print "best tree is ", tree
    pdb.set_trace()
