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
import sys
import shutil
import platform
import tempfile
import subprocess

import pdb


class Phyml:
    """ """
    def __init__(self, phylip, pth='bin', exe=None):
        self.phylip = os.path.abspath(os.path.expanduser(phylip))
        self.phyml3 = self._get_phyml_pth(pth, exe)
        self.models = {
                'JC69':"+\nM\nM\nM\nM\nM\nR\nY\n",
                'JC69I':"+\nM\nM\nM\nM\nM\nV\nY\nR\nY\n",
                'JC69G':"+\nM\nM\nM\nM\nM\nY\n",
                'JC69IG':"+\nM\nM\nM\nM\nM\nV\nY\nY\n",
                'F81':"+\nM\nM\nM\nM\nM\nM\nM\nF\nR\nY\n",
                'F81I':"+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nR\nY\n",
                'F81G':"+\nM\nM\nM\nM\nM\nM\nM\nF\nY\n",
                'F81IG':"+\nM\nM\nM\nM\nM\nM\nM\nF\nV\nY\nY\n",
                'K2P':"+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nY\n",
                'K2PI':"+\nM\nM\nM\nM\nM\nM\nT\nY\nR\nV\nY\nY\n",
                'K2PG':"+\nM\nM\nM\nM\nM\nM\nT\nY\nY\n",
                'K2PIG':"+\nM\nM\nM\nM\nM\nM\nT\nY\nV\nY\nY\n",
                'HKY':"+\nF\nT\nY\nR\nY\n",
                'HKYI':"+\nF\nT\nY\nR\nV\nY\nY\n",
                'HKYG':"+\nF\nT\nY\nY\n",
                'HKYIG':"+\nF\nT\nY\nV\nY\nY\n",
                'SYM':"+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nY\n",
                'SYMI':"+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nR\nV\nY\nY\n",
                'SYMG':"+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nY\n",
                'SYMIG':"+\nM\nM\nM\nM\nE\n0.25\n0.25\n0.25\n0.25\nK\n012345\n1.00\n1.00\n1.00\n1.00\n1.00\n1.00\nV\nY\nY\n",
                'GTR':"+\nM\nM\nM\nF\nR\nY\n",
                'GTRI':"+\nM\nM\nM\nF\nR\nV\nY\nY\n",
                'GTRG':"+\nM\nM\nM\nF\nY\n",
                'GTRIG':"+\nM\nM\nM\nF\nV\nY\nY\n"
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

    def run(self):
        pass

    def _get_log_like(self, regex, phylip):
        result = None
        stats = ''.join([phylip, '_phyml_stats.txt'])
        for line in open(stats, 'rU'):
            result = regex.search(line)
            if result:
                break
        if result is None:
            raise ValueError("No Log-likelihood found")
        return float(result.groups()[0])

    def _get_aic_tree(self, phylip):
        tree = None
        treefile = ''.join([phylip, '_phyml_tree.txt'])
        tree = open(treefile, 'rU').read().strip()
        if tree is None or tree == '':
            raise ValueError("No tree found")
        return tree

    def select_model(self):
        # compile this once
        regex = re.compile("Log-likelihood:\s+(.+)")
        results = {}
        # get current dir
        cwd = os.getcwd()
        for model, template in self.models.iteritems():
            sys.stdout.write('.')
            sys.stdout.flush()
            # create tempfile to hold phyml commands
            working = tempfile.mkdtemp(dir='tmp')
            phylip = os.path.basename(self.phylip)
            template = "%s\n%s" % (phylip, template)
            #pdb.set_trace()
            # copy phylip to working
            shutil.copyfile(
                self.phylip,
                os.path.join(working, os.path.basename(self.phylip))
                )
            os.chdir(working)
            cli = [self.phyml3]
            subprocess.Popen(
                    cli,
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE
                ).communicate(input=template)
            loglik = self._get_log_like(regex, phylip)
            tree = self._get_aic_tree(phylip)
            results[model] = [loglik, tree]
            os.chdir(cwd)
            shutil.rmtree(working)
        pdb.set_trace()


def main():
    phyml = Phyml('tests/alignments/phylip_primates/chr1_1036.phylip', '../../binaries')
    phyml.select_model()


if __name__ == '__main__':
    main()
