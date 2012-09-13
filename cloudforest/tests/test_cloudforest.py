#!/usr/bin/env python
# encoding: utf-8
"""
File: test_cloudforest.py
Author: Brant Faircloth

Created by Brant Faircloth on 26 April 2012 16:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

"""

import copy
import cPickle
import unittest
import dendropy
from context import cloudforest as cl

# import pdb


class TestProcess(unittest.TestCase):

    def setUp(self):
        self.p = cl.Process()
        self.one = open('alignments/3.oneliners', 'rU').readline()

    def test_oneliner_to_phylip(self):
        """[Process] Oneliner to phylip"""
        expected = cPickle.load(open('pickles/expected_phylip.pickle'))
        observed = cl.oneliner_to_phylip(self.one)
        assert observed == expected

    def test_bootstrap(self):
        """[Process] Bootstrap"""
        # we are assuming sampling with replacement works as advertised
        # and just checking to make sure we re-order according to sample
        # here.
        expected = self.prep_oneliner_array()
        taxa, align = cl.oneliner_to_array(expected)
        bs, choices = cl.get_bootstraps(align, return_choices=True)
        for k, aln in enumerate(bs):
            assert (bs[k] == align[choices[k]]).all()

    def prep_oneliner_array(self):
        """[Process] Prep oneliner array"""
        locus = self.one.strip().split(';')
        locus = locus[0].split(':')[1]
        return locus

    def test_oneliner_to_array(self):
        """[Process] Oneliner to array"""
        exp_taxa = ['MusMuscu', 'GorGoril', 'PanTrogl']
        exp_align = cPickle.load(open('pickles/expected_align_to_array.pickle'))
        locus = self.prep_oneliner_array()
        obs_taxa, obs_align = cl.oneliner_to_array(locus)
        assert obs_taxa == exp_taxa
        # comparing arrays
        assert (obs_align == exp_align).all()

    def test_array_to_oneliner(self):
        """[Process] Array to oneliner"""
        expected = self.prep_oneliner_array()
        taxa, align = cl.oneliner_to_array(expected)
        observed = cl.array_to_oneliner(taxa, align)
        assert observed == expected

    def test_make_tree_name(self):
        """[Process] Tree name from args_dict"""
        d = {'chrm': 'chr1_1036'}
        observed = cl.make_tree_name(d)
        expected = 'chrm=chr1_1036'
        assert observed == expected

    def test_split_oneliner(self):
        """[Process] Split oneliner"""
        # send single oneliner
        expected = {'chrm': 'chr1_1036'}
        observed = cl.split_oneliner(self.one)
        # dict
        assert observed[0] == expected
        # locus
        assert observed[1] == self.one.split(':')[1]

    def test_split_oneliner_with_default_model(self):
        """[Process] Split oneliner and assign default model"""
        # send single oneliner
        expected = {'chrm': 'chr1_1036', 'model': 'GTR'}
        observed = cl.split_oneliner(self.one, default_model=True)
        # dict
        assert observed[0] == expected
        # locus
        assert observed[1] == self.one.split(':')[1]

    def test_duplicate_oneliner(self):
        """[Process] Duplicate oneliners"""
        locus = self.prep_oneliner_array()
        d = {'chrm': 'chr1_1036'}
        tree_name = cl.make_tree_name(d)
        oneliner = "%s:%s" % (tree_name, locus)
        expected = [(r, oneliner) for r in reversed(xrange(1, 6))]
        dupes = self.p.duplicate_oneliner(1, oneliner, 5)
        observed = [d for d in dupes]
        assert observed == expected

    def test_lines_to_oneliner(self):
        # I'm not entirely sure what this does
        pass

    def test_get_genetrees_gtr(self):
        """[Process] genetrees works with GTR"""
        obs_gen = self.p.get_genetrees(1, self.one, pth='../binaries', genetrees=True)
        exp_key, exp_tree = cPickle.load(open('pickles/expected_genetrees_gtr.pickle'))
        observed = [o for o in obs_gen]
        assert len(observed) == 1
        observed = observed[0]
        obs_key, obs_tree = observed
        assert obs_key == exp_key
        assert obs_tree == exp_tree

    def test_get_genetrees_hky(self):
        """[Process] genetrees works with models in onliner"""
        one_liner = open('alignments/3.oneliners', 'rU').readline()
        one_liner = copy.deepcopy(self.one).split(':')
        one_liner[0] += ",model=HKY"
        one_liner = ':'.join(one_liner)
        obs_gen = self.p.get_genetrees(1, one_liner, pth='../binaries', genetrees=True)
        exp_key, exp_tree = cPickle.load(open('pickles/expected_genetrees_hky.pickle'))
        observed = [o for o in obs_gen]
        assert len(observed) == 1
        observed = observed[0]
        obs_key, obs_tree = observed
        assert obs_key == exp_key
        assert obs_tree == exp_tree

    def test_get_genetrees_and_models_for_genetrees(self):
        """[Process] genetrees_and_models yields tree"""
        obs_gen = self.p.get_genetrees_and_models(1, self.one, pth='../binaries', genetrees=True)
        exp_key, exp_tree = cPickle.load(open('pickles/expected_genetrees_and_model.pickle'))
        observed = [o for o in obs_gen]

        assert len(observed) == 1
        observed = observed[0]
        obs_key, obs_tree = observed
        print obs_tree
        print exp_tree
        assert obs_key == exp_key
        assert obs_tree == exp_tree

    def test_get_genetrees_and_models_for_oneliner(self):
        """[Process] genetrees_and_models yields oneliner"""
        obs_gen = self.p.get_genetrees_and_models(1, self.one, pth='../binaries', genetrees=False)
        exp_key, exp_oneliner = cPickle.load(open('pickles/expected_genetree_oneliners.pickle'))
        observed = [o for o in obs_gen]
        assert len(observed) == 1
        obs_key, exp_oneliner = observed[0]
        assert obs_key == exp_key
        assert exp_oneliner == exp_oneliner


class TestPhymlMethods(unittest.TestCase):

    def setUp(self):
        self.phyml = cl.Phyml('alignments/phylip_primates/chr1_1036.phylip', pth='../binaries')

    def test_model_statements(self):
        """[Phyml] Model templates are correct"""
        expected = cPickle.load(open('pickles/phyml_models.pickle'))
        assert self.phyml.models == expected

    def test_model_numparams(self):
        """[Phyml] numparams is correct"""
        expected = cPickle.load(open('pickles/phyml_numparams.pickle'))
        assert self.phyml.numparams == expected

    def test_taxa_and_characters(self):
        """[Phyml] Parse taxa and nchar"""
        self.phyml._get_taxon_and_char_data()
        assert self.phyml.taxa == 10
        assert self.phyml.nchar == 430
        assert self.phyml.nbranch == (2 * 10) - 3

    def test_compute_aicc(self):
        """[Phyml] AICc formula"""
        # with loglik = 20 and params = 1 and nchar = 430
        expected = -2.335766423357664
        # set nchar and nbranch
        self.phyml.nchar = 430
        self.phyml.nbranch = 17
        # params of JC69I = 1
        aicc = self.phyml._compute_aicc('JC69I', 20)
        self.assertAlmostEqual(expected, aicc, 4)

    def get_tree_distances(self, newick1, newick2):
        tree1 = dendropy.Tree()
        tree2 = dendropy.Tree()
        tree1.read_from_string(newick1, 'newick')
        tree2.read_from_string(newick2, 'newick')
        return tree1.euclidean_distance(tree2)

    def test_slow_aicc_model(self):
        """[Phyml] AICc best models

        Run all these because it's slow and we don't want to recompute results
        to test every little thing.  The tests of trees are a little tricky
        because we need to ensure the trees are identical, but the ML optimization
        returns slightly different branch lengths all the time.
        So, compute tree distance, and assert that it's almost equal to 0.
        """
        # best model
        expected = 'GTR'
        observed = self.phyml.best_aicc_model()
        assert observed == expected
        # best tree
        expected = cPickle.load(open('pickles/best_aicc_tree.pickle'))
        observed = self.phyml.best_aicc_tree()
        distance = self.get_tree_distances(observed[1], expected[1])
        self.assertAlmostEqual(distance, 0.0, 2)
        # best model and tree
        expected = cPickle.load(open('pickles/best_aicc_model_and_tree.pickle'))
        observed = self.phyml.best_aicc_model_and_tree()
        assert observed[0] == expected[0]
        distance = self.get_tree_distances(observed[1], expected[1])
        self.assertAlmostEqual(distance, 0.0, 2)
        # all results
        expected = cPickle.load(open('pickles/best_aicc_model_results.pickle'))
        observed = self.phyml.aicc_model_results()
        observed_aicc = sorted(observed.keys())
        expected_aicc = sorted(expected.keys())
        for pos, aic in enumerate(observed_aicc):
            # aiccs are almost equal
            self.assertAlmostEqual(aic, expected_aicc[pos], 4)
            # model names are equal
            observed[aic][0] == expected[expected_aicc[pos]][0]
            # trees are almost equal
            distance = self.get_tree_distances(
                    observed[aic][1],
                    expected[expected_aicc[pos]][1]
                )
            self.assertAlmostEqual(distance, 0.0, 2)

    def test_run(self):
        """[Phyml] Phyml.run()"""
        expected = cPickle.load(open('pickles/gtr_lnl_and_model.pickle'))
        observed = self.phyml.run('GTR')
        distance = self.get_tree_distances(observed[1], expected[1])
        self.assertAlmostEqual(distance, 0.0, 2)
        self.assertAlmostEqual(float(observed[0]), float(expected[0]), 2)


class TestCoreFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def test_format_oneliner_from_dict(self):
        d = {
                'gor_gor': '---ATCATAGTTGAA',
                'oto_gar': 'GTAATCATAGTTGAA'
            }
        observed = cl.format_oneliner_from_dict(d, 'test')
        expected = 'chrm=test:gor_gor,---ATCATAGTTGAA,oto_gar,GTAATCATAGTTGAA;'
        assert observed == expected

    def test_phyml_to_oneliner(self):
        phylip = open('alignments/phylip_primates/chr1_1036.phylip', 'rU').read()
        observed = cl.phylip_to_oneliner(phylip, 'chr1_1036')
        expected = cPickle.load(open('pickles/expected_oneliner_from_phylip.pickle'))
        assert observed == expected


if __name__ == '__main__':
    unittest.main()
    #fast = unittest.TestLoader().loadTestsFromTestCase(TestCoreFunctions)
    #unittest.TextTestRunner().run(fast)
