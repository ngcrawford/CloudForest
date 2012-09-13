#!/usr/bin/env python
# encoding: utf-8

"""
File: cloudforest/cloudforest_emr.py

Created by Nicholas G. Crawford
Copyright (c) 2011-2012 Nicholas G. Crawford and Brant C. Faircloth. All rights reserved.

Description: Compute genetrees and bootstrap replicates from a
directory of phylip files input on the CLI using Amazon Elastic
Map Reduce (EMR) or a local Hadoop cluster.

"""

from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol

from core import Process


class ProcessPhyloData(MRJob, Process):
    """Main class for MrJob. Subclasses Process() from core.py"""
    def __init__(self, args):
        MRJob.__init__(self, args=args)
        # self.binaries = resource_filename(__name__, 'binaries')

    def configure_options(self):
        super(ProcessPhyloData, self).configure_options()
        self.add_passthrough_option(
                '--bootreps',
                dest='bootreps2run',
                default=False,
                type='int',
                help='Number of bootstrap replicates to generate'
            )

        self.add_passthrough_option(
                '--gene-trees',
                action="store_true",
                dest='gene_trees',
                default=False,
                help='Generate gene trees from alignments'
            )

        self.add_passthrough_option(
                '--full-analysis',
                action="store_true",
                default=False,
                help='Run full analysis'
            )

        self.add_passthrough_option(
                '--mraic',
                dest='mraic_opt',
                action='store_true',
                default=False,
                help='Use MrAIC to calculate models'
            )

        self.add_passthrough_option(
            '--constraint-tree', 
            type='str', 
            default='False',
            help="""Provide a constraint tree as quoted 
                    newick string. ex: --constraint-tree='((A,B),C);'""")

    def basic_reducer(self, key, line):
        """Do not reduce"""
        yield key, line

    def steps(self):
        """Job steps to run through hadoop/mrjob"""
        # Do full analysis
        if self.options.full_analysis == True:
            return [self.mr(self.get_genetrees_and_models, self.concatenate_oneliners),
                      self.mr(self.duplicate_oneliner, reducer=None),
                      self.mr(self.get_bootstrap_replicates, reducer=None),
                      self.mr(self.get_genetrees, reducer=None)]

        if self.options.gene_trees == True and self.options.mraic_opt == True:
            # TODO rewrite as a single fuction outside of steps.
            def output_protocol(self):
                return RawValueProtocol()
            return [self.mr(self.get_genetrees_and_models, reducer=None)]

        if self.options.gene_trees == True and self.options.mraic_opt == False:
            # TODO rewrite as a single fuction outside of steps.
            def output_protocol(self):
                return RawValueProtocol()
            return [self.mr(mapper=self.get_genetrees, reducer=None)]

if __name__ == '__main__':
    ProcessPhyloData.run()
