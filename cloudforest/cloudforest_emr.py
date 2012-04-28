
from mrjob.job import MRJob
from mrjob.protocol import RawValueProtocol

from core import Process


class ProcessPhyloData(MRJob, Process):
    """
    Main class for MrJob. Subclasses Process from cloudforest
    which is shared with cloudforest_mpi.
    """
    def __init__(self, args):
        MRJob.__init__(self, args=args)
        # self.binaries = resource_filename(__name__, 'binaries')

    def configure_options(self):
        super(ProcessPhyloData, self).configure_options()
        self.add_passthrough_option(
                '--bootreps',
                dest='bootreps2run',
                default=None,
                type='int',
                help='Number of bootstrap replicates to generate'
            )

        self.add_passthrough_option(
                '--gene-trees',
                action="store_true",
                dest='gene_trees',
                help='Generate gene trees from alignments'
            )

        self.add_passthrough_option(
                '--full-analysis',
                action="store_true",
                help='Run full analysis'
            )

        self.add_passthrough_option(
                '--mraic',
                dest='mraic_opt',
                action='store_true',
                default=False,
                help='Use MrAIC to calculate models'
            )

    def basic_reducer(self, key, line):
        """Do not reduce"""
        yield key, line

    def steps(self):
        """Job steps to run through hadoop/mrjob"""
        # Do full analysis
        if self.options.full_analysis == True:
            return [self.mr(self.mrAIC, self.lines2Oneliner),
                      self.mr(self.duplicateOneliners, reducer=None),
                      self.mr(self.bootstrapReplicates, reducer=None),
                      self.mr(self.phyml, reducer=None)]

        if self.options.gene_trees == True and self.options.mraic_opt == True:
            # TODO rewrite as a single fuction outside of steps.
            def output_protocol(self):
                return RawValueProtocol()
            return [self.mr(self.mrAIC, reducer=None)]

        if self.options.gene_trees == True and self.options.mraic_opt == None:
            # TODO rewrite as a single fuction outside of steps.
            def output_protocol(self):
                return RawValueProtocol()
            return [self.mr(mapper=self.phyml, reducer=None)]

if __name__ == '__main__':
    ProcessPhyloData.run()
