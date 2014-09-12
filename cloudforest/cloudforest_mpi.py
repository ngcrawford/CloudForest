#!/usr/bin/env python
# encoding: utf-8
"""
File: cloudforest_mpi.py
Author: Brant Faircloth

Created by Brant Faircloth on 12 May 2012 10:05 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Compute genetrees and bootstrap replicates from a
directory of phylip files input on the CLI using MPI.

"""

import os
import sys
import glob
import argparse
import numpy as np
from time import localtime, strftime
from core import Phyml, is_dir, is_file, FullPaths, phylip_to_oneliner, split_oneliner, oneliner_to_array, get_bootstraps, array_to_oneliner, oneliner_to_phylip

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "input",
            type=is_dir,
            action=FullPaths,
            help="""The path to a directory of PHYLIP-formatted alignments."""
        )
    parser.add_argument(
            "output",
            type=is_dir,
            action=FullPaths,
            help="""The path to a directory in which to store the output.""",
        )
    parser.add_argument(
            "run",
            choices=['genetrees', 'bootstraps', 'both'],
            default=None,
            help="""The type of analysis to run.""",
        )
    parser.add_argument(
            "phyml",
            action=FullPaths,
            help="""The path (including executable name) to PhyML on your platform.""",
        )
    parser.add_argument(
            "--genetrees",
            action=FullPaths,
            type=is_file,
            default=None,
            help="""The path to the genetrees output by a previous run.  Needed only when bootstrapping a prior run.""",
        )
    parser.add_argument(
            "--bootreps",
            type=int,
            default=100,
            help="""The number of bootstrap replicates to run.""",
        )
    parser.add_argument(
            "--parallelism",
            choices=['mpi', 'multiprocessing', 'single'],
            default='mpi',
            help="""The type of parallelism to use.""",
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=7,
            help="""The number of compute cores to use.""",
        )

    args = parser.parse_args()

    if args.run == 'bootstraps' and args.genetrees is None:
        sys.exit("\nIf runnning only boostraps, you must pass a genetrees file")

    return args


def make_tree_name(args_dict):
    """Converts dictionary of argument into a sorted string with the
    format:  """
    # join is faster and cleaner than concatenate
    # list comp is faster than for
    name = ','.join(['='.join([pair[0], pair[1]]) for pair in sorted(args_dict.items())])
    return name


def generate_bootreps(bootreps, phyml, oneliners):
    """Replicate the data set bootrep numer of times, prior to bootstrapping"""
    # keep bootrep numbers indexed by one
    for i in xrange(1, bootreps + 1):
        yield i, phyml, oneliners


def get_models_from_genetrees(genetrees):
    """Iterate over a set of precomputed genetrees and return the models in metadata"""
    # iterate over genetrees to grab models
    models = {}
    for line, tree in enumerate(open(genetrees, 'rU')):
        ts = tree.split(' ')[1].strip("'")
        for item in ts.split(','):
            if item.startswith('chrm'):
                locus = item.split('=')[1]
            elif item.startswith('model'):
                model = item.split('=')[1]
        models[locus] = model
    assert len(models) == line + 1, "Not all loci have a model"
    return models


def get_bootstrap_replicates(multilocus_bstrap):
    """Boostrap the bases of all alignments in a resampled population of alignments"""
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


def genetree_worker(params):
    """Worker function to compute genetrees for individual loci"""
    locus, fullpth = params
    name, phylip = locus
    pth, exe = os.path.split(fullpth)
    args_dict = {}
    phyml = Phyml(phylip, pth=pth, exe=exe)
    model, tree = phyml.best_aicc_model_and_tree()
    args_dict['chrm'] = name
    args_dict['model'] = model
    tree = "tree '%s' = [&U] %s" % (make_tree_name(args_dict), tree)
    sys.stdout.write("[Info] {0} {1} genetree completed\n".format(
                strftime("%a, %d %b %Y %H:%M:%S", localtime()),
                name
            )
        )
    sys.stdout.flush()
    return (name, model, tree)


def bootstrap_worker(params):
    """Worker function to compute boostrap replicates of datasets and indiv. loci"""
    bootstrap_trees = []
    rep, fullpth, oneliners = params
    pth, exe = os.path.split(fullpth)
    # first, resample w/ replacement/bootstrap across loci
    multilocus_bstrap = get_bootstraps(oneliners)
    # Resample w/ replacement/boostrap bases within loci
    bootstraps = get_bootstrap_replicates(multilocus_bstrap)
    for oneliner in bootstraps:
        args_dict, locus = split_oneliner(oneliner, default_model=True)
        phylip = oneliner_to_phylip(locus)
        phyml = Phyml(phylip, pth=pth, exe=exe)
        # run phyml.  if no model, defaults to GTR
        # TOOD: Why do we need LnL?
        args_dict['lnL'], tree = phyml.run(args_dict['model'])
        #bootstrap_trees.append("tree '%s' = [&U] %s" % (make_tree_name(args_dict), tree))
        bootstrap_trees.append('''%s\t"%s"''' % (rep, tree))
    sys.stdout.write("[Info] {0} bootstrap completed\n".format(
            strftime("%a, %d %b %Y %H:%M:%S", localtime())
        )
    )
    sys.stdout.flush()
    return bootstrap_trees


def boostrap_all_loci(args, models, alns):
    """Compute trees from bootstrap replicates of a dataset"""
    oneliners = [phylip_to_oneliner(phylip, locus, models[locus]) \
                for locus, phylip in alns.iteritems()]
    # for every rep in boostraps, map loci onto worker that will
    # bootstrap, run phyml, and return bootstrap trees
    params = generate_bootreps(args.bootreps, args.phyml, oneliners)
    bootreps = mmap(bootstrap_worker, params)
    # write
    outname = "%s-bootreps.tree" % (args.bootreps)
    outf = open(os.path.join(args.output, outname), 'w')
    for bootrep in bootreps:
        for tree in bootrep:
            outf.write("%s\n" % (tree))


def main():
    args = get_args()
    # read in phylip files
    alns = {}
    for f in glob.glob(os.path.join(args.input, '*.phy*')):
        alns[os.path.splitext(os.path.basename(f))[0]] = open(f, 'rU').read()
    # replicate our options for passing to map()
    opts = [args.phyml for i in range(len(alns))]
    # compute genetrees
    if args.run == 'genetrees' or args.run == 'both':
        sys.stdout.write("Running genetrees...\n\n")
        params = zip(alns.items(), opts)
        genetrees = mmap(genetree_worker, params)
        # write genetrees to output file
        outf = open(os.path.join(args.output, 'genetrees.tre'), 'w')
        for tree in genetrees:
            outf.write("%s\n" % (tree[2]))
        outf.close()
    # compute bootreps on genetrees from above
    if args.run == 'both':
        sys.stdout.write("Running bootstraps of genetrees...\n")
        # get models for each locus based on genetrees in-memory
        models = dict([[tree[0], tree[1]] for tree in genetrees])
        boostrap_all_loci(args, models, alns)
    # compute bootreps on genetrees from a file
    if args.run == 'bootstraps':
        sys.stdout.write("Running boostraps...\n")
        # get models for each locus based on genetrees in genetree file
        models = get_models_from_genetrees(args.genetrees)
        boostrap_all_loci(args, models, alns)


if __name__ == '__main__':
    args = get_args()
    if args.parallelism == 'mpi':
        from deap.dtm import map as mmap
        from deap.dtm import start
        start(main)
    elif args.parallelism == 'multiprocessing':
        from multiprocessing import Pool
        pool = Pool(args.cores)
        mmap = pool.map
        main()
    elif args.parallelism == 'single':
        mmap = map
        main()
