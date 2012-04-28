#!/usr/bin/env python
# encoding: utf-8
"""
File: test_cloudforest.py
Author: Brant Faircloth

Created by Brant Faircloth on 26 April 2012 16:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

"""

import os
import cPickle
import unittest
import numpy as np
from context import cloudforest

import pdb


class TestProcess(unittest.TestCase):

    def setUp(self):
        self.p = cloudforest.Process()
        self.one = open('alignments/3.oneliners', 'rU').readline()

    def test_oneliner_to_phylip(self):
        expected = cPickle.load(open('pickles/expected_phylip.pickle'))
        observed = self.p.oneliner_to_phylip(self.one)
        assert observed == expected

    def test_bootstrap(self):
        # we are assuming sampling with replacement works as advertised
        # and just checking to make sure we re-order according to sample
        # here.
        expected = self.prep_oneliner_array()
        taxa, align = self.p.oneliner_to_array(expected)
        bs, choices = self.p.get_bootstraps(align, return_choices=True)
        for k, aln in enumerate(bs):
            assert (bs[k] == align[choices[k]]).all()

    def prep_oneliner_array(self):
        locus = self.one.strip().split(';')
        locus = locus[0].split(':')[1]
        return locus

    def test_oneliner_to_array(self):
        exp_taxa = ['MusMuscu', 'GorGoril', 'PanTrogl']
        exp_align = cPickle.load(open('pickles/expected_align_to_array.pickle'))
        locus = self.prep_oneliner_array()
        obs_taxa, obs_align = self.p.oneliner_to_array(locus)
        assert obs_taxa == exp_taxa
        # comparing arrays
        assert (obs_align == exp_align).all()

    def test_array_to_oneliner(self):
        expected = self.prep_oneliner_array()
        taxa, align = self.p.oneliner_to_array(expected)
        observed = self.p.array_to_oneliner(taxa, align)
        assert observed == expected

    def test_make_tree_name(self):
        d = {'chrm': 'chr1_1036'}
        observed = self.p.make_tree_name(d)
        expected = 'chrm=chr1_1036'
        assert observed == expected

    def test_split_oneliner(self):
        # send single oneliner
        expected = {'chrm': 'chr1_1036'}
        observed = self.p.split_oneliner(self.one, {})
        # dict
        assert observed[0] == expected
        # locus
        assert observed[1] == self.one.split(':')[1]

    def test_duplicate_oneliner(self):
        locus = self.prep_oneliner_array()
        d = {'chrm': 'chr1_1036'}
        tree_name = self.p.make_tree_name(d)
        oneliner = "%s:%s" % (tree_name, locus)
        expected = [(r, oneliner) for r in reversed(xrange(1, 6))]
        dupes = self.p.duplicate_oneliner(1, oneliner, 5)
        observed = [d for d in dupes]
        assert observed == expected

    def test_process_stats_file(self):
        pass

    def test_lines_to_oneliner(self):
        # I'm not entirely sure what this does
        pass

    def test_phyml(self):
        pass

    def test_get_genetrees_and_models_for_genetrees(self):
        obs_gen = self.p.get_genetrees_and_models(1, self.one, bin='../binaries', genetrees=True)
        exp_key, exp_tree = cPickle.load(open('pickles/expected_genetrees.pickle'))
        observed = [o for o in obs_gen]
        assert len(observed) == 1
        observed = observed[0]
        obs_key, obs_tree = observed
        assert obs_key == exp_key
        assert obs_tree == exp_tree

    def test_get_genetrees_and_models_for_oneliner(self):
        obs_gen = self.p.get_genetrees_and_models(1, self.one, bin='../binaries', genetrees=False)
        exp_key, exp_oneliner = cPickle.load(open('pickles/expected_genetree_oneliners.pickle'))
        observed = [o for o in obs_gen]
        assert len(observed) == 1
        obs_key, exp_oneliner = observed[0]
        assert obs_key == exp_key
        assert exp_oneliner == exp_oneliner


if __name__ == '__main__':
    unittest.main()
