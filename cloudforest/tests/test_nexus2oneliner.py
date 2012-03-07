"""
nexus2oneliners_test.py

Created by Nicholas Crawford and Brant C. Faircloth Copyright (c) 2010 Nicholas Crawford and
Brant C. Faircloth. All rights reserved.

Tests functions in nexus2oneliners.py
"""

import textwrap
import unittest
from cloudforest import nexus2oneliner

class TestNexus2OnelinereFunctions(unittest.TestCase):

	def setUp(self):
		"""Create necessary datasets for testing"""
		self.name = 'test_nexus'
		self.seq = "ATCT--ATCNN"
		self.no_data_seq = "?????----"
		self.nexus = """\
						#NEXUS
						begin data;
							dimensions ntax=10 nchar=430;
							format datatype=dna missing=? gap=-;
						matrix
						mus_musculus        ----TCCTAGCTGAACAGAGAAGGGTGATTAACGATAGCAATTTATT
						otolemur_garnetti   GTAATCATAGTTGAACCGAGAAAGGTGATTAACGATAGCAATTTATT
						callithrix_jacchus  GTAATCATAGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT
						gorilla_gorilla     ---ATCATAGTTGAACCAAGAAGAGTGATTAACGATAGCAATTTATT
						homo_sapiens        GTAATCATAGTTGAACCAAGAAGAGTGATTAACGATAGCAATTTATT
						pan_troglodytes     GTAATCACAGTTGAACCAAGAAGAGTGATTAACGATAGCAATTTATT
						pongo_abelii        GTAATCATAGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT
						nomascus_leucogenys GTAATCATAGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT
						papio_hamadryas     GTAATCATGGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT
						rhesus_macaque      GTAATCATGGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT
						;
						end;"""

		self.oneliner = """\
						chrm=test_nexus:\
						MusMuscu,----TCCTAGCTGAACAGAGAAGGGTGATTAACGATAGCAATTTATT,\
						OtoGarne,GTAATCATAGTTGAACCGAGAAAGGTGATTAACGATAGCAATTTATT,\
						CalJacch,GTAATCATAGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT,\
						GorGoril,---ATCATAGTTGAACCAAGAAGAGTGATTAACGATAGCAATTTATT,\
						HomSapie,GTAATCATAGTTGAACCAAGAAGAGTGATTAACGATAGCAATTTATT,\
						PanTrogl,GTAATCACAGTTGAACCAAGAAGAGTGATTAACGATAGCAATTTATT,\
						PonAbeli,GTAATCATAGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT,\
						NomLeuco,GTAATCATAGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT,\
						PapHamad,GTAATCATGGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT,\
						RheMacaq,GTAATCATGGTTGAACCAAGAAGGGTGATTAACGATAGCAATTTATT;\n"""

	def test_isSeqOK(self):
		ok_seq_result = nexus2oneliner.isSeqOK(self.seq)
		self.assertTrue(ok_seq_result)

		bad_seq_result = nexus2oneliner.isSeqOK(self.no_data_seq)
		self.assertFalse(bad_seq_result)

	def test_removeAmbiguousBases(self):
		cleaned_seq = nexus2oneliner.removeAmbiguousBases(self.seq)
		self.assertEqual(cleaned_seq, "ATCT--ATC--")

	def test_nexus2oneliner(self):
		# dedent and clean lines of intervening tabs
		nexus = textwrap.dedent(self.nexus)
		test_str = textwrap.dedent(self.oneliner).replace("\t", "")

		nexus_lines = nexus.split("\n")
		result = nexus2oneliner.nexus2oneliner(nexus_lines, self.name)
		self.assertEqual(result, test_str)


if __name__ == '__main__':
	suite = unittest.TestLoader().loadTestsFromTestCase(TestNexus2OnelinereFunctions)
	unittest.TextTestRunner(verbosity=3).run(suite)