import unittest
from cloudforest.cloudforest_emr import ProcessPhyloData

class TestCloudForestFunctions(unittest.TestCase):

    def test_oneliner2phylip(self):
        pass

    def test_onelinerAlignment2Array(self):
        pass

    def test_array2OnelinerAlignmente(self):
        pass
         
    def test_duplicateOneliners(self):
        pass
    
    def test_GeneTrees(self):
        """Tests that creating genetrees with Cloudforest
            makes the appropriate number of trees.
        """

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
        self.assertEqual(len(results), 3)

    
    def test_MrAIC(self):
        """Tests that creating genetrees with the --mraic
            command set calculates the correct evolutionary models.
        """
        
        test_data = open('alignments/3.oneliners','rU')
        mr_job = ProcessPhyloData(['-r', 'local', '--setup-cmd', 
                                    'mkdir -p tmp','--gene-trees',
                                    '--mraic',
                                    '--archive=../gzips/osx.phylo.tar.gz#bin',
                                    "-"])
    
        mr_job.sandbox(stdin=test_data)
        results = []
        with mr_job.make_runner() as runner:
            runner.run()
            for line in runner.stream_output():
                # Use the job's specified protocol to read the output
                key, value = mr_job.parse_output_line(line)
                results.append(value) 
        
        # get models from output
        models = [ item.split(" ")[1].strip("'").split('=')[-1] for item in results]

        # check results.
        # result looks like: ['GTR', 'HKY', 'HKY'])
        self.assertIn('GTR', models)
        self.assertIn('HKY', models)
        self.assertEqual(len(results),3)


    def test_Bootstrapping(self):
        """Tests that bootstraps produces the correct number or replicates"""

        test_data = open('alignments/3.oneliners','rU')
        mr_job = ProcessPhyloData(['-r', 'local', 
                                    '--setup-cmd', 'mkdir -p tmp',
                                    '--full-analysis',
                                    '--bootreps=5',
                                    '--archive=../gzips/osx.phylo.tar.gz#bin',
                                    "-"])
    
        mr_job.sandbox(stdin=test_data)
        results = []
        with mr_job.make_runner() as runner:
            runner.run()
            for line in runner.stream_output():
                # Use the job's specified protocol to read the output
                key, value = mr_job.parse_output_line(line)
                results.append(value)

        self.assertEqual(len(results), 15)


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCloudForestFunctions)
    unittest.TextTestRunner(verbosity=3).run(suite)