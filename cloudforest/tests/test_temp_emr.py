import unittest
from cloudforest.cloudforest_emr import ProcessPhyloData

import pdb


class TestCloudForestFunctions(unittest.TestCase):

    def get_models_from_mrjob_out(self, results):
        """split models from mrjob outstream"""
        return [r.split(" ")[1].strip("'").split('=')[-1] for r in results]

    def test_get_genetrees(self):
        """[MrJob] Generate genetrees"""

        test_data = open('alignments/3.oneliners', 'rU')
        mr_job = ProcessPhyloData([
                '-r', 'local',
                '--setup-cmd',
                'mkdir -p tmp',
                '--gene-trees',
                '--archive=../gzips/osx.phylo.tar.gz#bin',
                "-"]
            )
        mr_job.sandbox(stdin=test_data)
        results = []
        with mr_job.make_runner() as runner:
            runner.run()
            for line in runner.stream_output():
                # Use the job's specified protocol to read the output
                key, value = mr_job.parse_output_line(line)
                results.append(value)
        #TODO:better checking of tree structure here - do distance checking
        self.assertEqual(len(results), 3)

    def test_MrAIC(self):
        """[MrJob] Select model and generate genetrees"""
        test_data = open('alignments/3.oneliners', 'rU')
        mr_job = ProcessPhyloData([
                '-r', 'local',
                '--setup-cmd',
                'mkdir -p tmp',
                '--gene-trees',
                '--mraic',
                '--archive=../gzips/osx.phylo.tar.gz#bin',
                '-']
            )
        mr_job.sandbox(stdin=test_data)
        results = []
        with mr_job.make_runner() as runner:
            runner.run()
            for line in runner.stream_output():
                # Use the job's specified protocol to read the output
                key, value = mr_job.parse_output_line(line)
                results.append(value)
        # check models
        expected_models = ['HKY', 'GTR', 'HKY']
        observed_models = self.get_models_from_mrjob_out(results)
        assert observed_models == expected_models
        #TODO:better checking of tree structure here - do distance checking


if __name__ == '__main__':
    unittest.main()
