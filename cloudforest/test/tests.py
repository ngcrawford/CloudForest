import unittest
from CloudForest.cloudforest import cloudforest

class TestCloudForestFunctions(unittest.TestCase):

    def test_oneliner2phylip(self):
        pass

    def test_onelinerAlignment2Array(self):
        pass

    def test_array2OnelinerAlignmente(self):
        pass
        
    def test_phyml(self):
        pass
        
    def test_mrAIC(self):
        pass
        
    def test_duplicateOneliners(self):
        pass

    # def basicReducer(self, key, line):
    #     """Yield lines from a generator."""
    #     for item in line:
    #         yield key, item  

    def test_basicReducer(self):

        test_data = open('/Users/nick/Desktop/Code/CloudForest/cloudforest/test/alignments/3.oneliners','rU')
        
        mr_job = cloudforest.ProcessPhyloData().sandbox(stdin=test_data)
        mr_job.options.gene_trees == True
        
        # mr_job.cloudforest.BootstrapAW().sandbox(fake_input)
        mr_job.run_mapper()
        stdin.close()
        pass
        
        
    def test_lines2Oneliner(self):
        pass
        
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestCloudForestFunctions)
    unittest.TextTestRunner(verbosity=2).run(suite)