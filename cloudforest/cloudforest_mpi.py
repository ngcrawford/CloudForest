#!/usr/bin/env python
# encoding: utf-8
"""
File: cloudforest_mpi.py
Author: Brant Faircloth

Created by Brant Faircloth on 12 May 2012 10:05 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import argparse
import numpy as np
from core import Phyml, is_dir, FullPaths, phylip_to_oneliner, split_oneliner, oneliner_to_array, get_bootstraps, array_to_oneliner, oneliner_to_phylip

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "input",
            type=is_dir,
            action=FullPaths,
            help="""Help text"""
        )
    parser.add_argument(
            "output",
            type=is_dir,
            action=FullPaths,
            help="""Help text""",
        )
    parser.add_argument(
            "type",
            choices=['genetrees', 'bootstraps', 'both'],
            default=None,
            help="""Help text""",
        )
    parser.add_argument(
            "phyml",
            action=FullPaths,
            help="""Help text""",
        )
    parser.add_argument(
            "--bootreps",
            type=int,
            default=5,
            help="""Help text""",
        )
    return parser.parse_args()


def make_tree_name(args_dict):
    """Converts dictionary of argument into a sorted string with the
    format:  """
    # join is faster and cleaner than concatenate
    # list comp is faster than for
    name = ','.join(['='.join([pair[0], pair[1]]) for pair in sorted(args_dict.items())])
    return name


def genetree_worker(params):
    locus, opts = params
    name, phylip = locus
    pth = opts
    args_dict = {}
    phyml = Phyml(phylip, pth)
    model, tree = phyml.best_aicc_model_and_tree()
    args_dict['chrm'] = name
    args_dict['model'] = model
    tree = "tree '%s' = [&U] %s" % (make_tree_name(args_dict), tree)
    return (name, model, tree)


def get_bootstrap_replicates(multilocus_bstrap):
    for locus in multilocus_bstrap:
        # split keys from alignments
        args_dict, locus = split_oneliner(locus)
        # convert alignments to arrays
        taxa, numpy_alignment = oneliner_to_array(locus.rstrip(';'))
        # transpose so we bootstrap by columns
        bases_by_col = np.column_stack(numpy_alignment)
        # bootstrap by columns
        shuffled = get_bootstraps(bases_by_col)
        # transpose back to rows of sequences
        shuffled = np.column_stack(shuffled)
        # convert back to oneliner
        oneliner = array_to_oneliner(taxa, shuffled)
        oneliner = "%s:%s" % (make_tree_name(args_dict), oneliner)
        yield oneliner


def bootstrap_worker(params):
    bootstrap_trees = []
    rep, pth, oneliners = params
    # first, bootstrap across loci
    multilocus_bstrap = get_bootstraps(oneliners)
    # second, bootstrap bases within loci
    bootstraps = get_bootstrap_replicates(multilocus_bstrap)
    for oneliner in bootstraps:
        args_dict, locus = split_oneliner(oneliner, default_model=True)
        phylip = oneliner_to_phylip(locus)
        phyml = Phyml(phylip, pth)
        # run phyml.  if no model, defaults to GTR
        # TOOD: Why do we need LnL?
        args_dict['lnL'], tree = phyml.run(args_dict['model'])
        #bootstrap_trees.append("tree '%s' = [&U] %s" % (make_tree_name(args_dict), tree))
        bootstrap_trees.append('''%s\t"%s"''' % (rep, tree))
    return bootstrap_trees


def generate_bootstrap_params(bootreps, phyml, oneliners):
    for i in xrange(bootreps):
        yield i, phyml, oneliners


def main():
    args = get_args()
    # read in phylip files
    alns = {}
    for f in glob.glob(os.path.join(args.input, '*.phy*')):
        alns[os.path.splitext(os.path.basename(f))[0]] = open(f, 'rU').read()
    opts = [args.phyml for i in range(len(alns))]
    if args.type == 'genetrees' or args.type == 'both':
        params = zip(alns.items(), opts)
        genetrees = map(genetree_worker, params)
        outf = open(os.path.join(args.output, 'genetrees.tre'), 'w')
        for tree in genetrees:
            outf.write("%s\n" % (tree[2]))
        outf.close()
    if args.type == 'both':
        models = dict([[tree[0], tree[1]] for tree in genetrees])
        oneliners = [phylip_to_oneliner(phylip, locus, models[locus]) \
                for locus, phylip in alns.iteritems()]
        # for every rep in boostraps, map loci onto worker that will
        # bootstrap, run phyml, and return bootstrap trees
        params = generate_bootstrap_params(args.bootreps, args.phyml, oneliners)
        bootreps = map(bootstrap_worker, params)
        # write
        outname = "%s-bootreps.tree" % (args.bootreps)
        outf = open(os.path.join(args.output, outname), 'w')
        for bootrep in bootreps:
            for tree in bootrep:
                outf.write("%s\n" % (tree))
    if args.type == 'bootstraps':
        # ensure we have bootstrap file
        pass


if __name__ == '__main__':
    main()
