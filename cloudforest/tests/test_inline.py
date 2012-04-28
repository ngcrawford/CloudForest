import os
import glob
import tempfile
import subprocess

import numpy as np

import pdb

DATA = ['chrm=chr1_1036,model=GTR:MusMuscu,----TCCTAGCTGAACAGAGAAGGGTGATTAACGATAGCAATTTATTGTATCTGAGGAGGAAAGAAAATAGGAGATGTTGTCGTTATCATAAATCATATCATAAATAACCAAAAAGAAGTCCTCGGAAAGCCTGTGTGGTGTCAGCAGCCAGCGCTTTCAGGCCCGTTAG-CCAGGATGTAGCCATAGAAACCCCAAAGTCATATTGTTATCTTGGTCTTGCTGCAGTTATGTGATACAACAGTATTAACATATAGATTAGGCCATTTATA-ATGCTAGGTGTCATTTCAGGTTGGCTTGC-TTATATCGC-ATTTTAGTTTCTTATTGTGCAGAATGTTAGTGGAATGGGGCTTTTTTACATGTCAGTGTGTGT-----AGGTATTTAT--AGCACATC---TTCACTAT--TGCTAATTTATCCTGTC,GorGoril,---ATCATAGTTGAACCAAGAAGAGTGATTAACGATAGCAATTTATTGTATCTGAGGAGGAGAGAAAATAGGAGATGTTGTCGTTATCATAAATCATATCATAAATAACCAAAAAGAAGCACTCAGAATGCTTGTGTGGTGTCAGCAGCCGGCCCTTTCGGGCACGTTAGCCCAGGATGTAGCCATAGAAACCCCAAAGTCATATTGTTATCTTAGTCTTGCTGCAGTTATGTGATACAACAGTATTAACATATAGATTAGGCCATTTATA-ATACTAGGTGTCATTTCAGGTTGGCTTGC-TTATATAGC-ATTTTAGTTTTTTATTGTGCAGAATGTTAGTGAAATAGTGCTTTTTTACAGGTCACTGTGTAT---CTAGGTATTTAT--AATGCATCATATTCACTCTGCTATTAATTTCTTGTG--,PanTrogl,GTAATCACAGTTGAACCAAGAAGAGTGATTAACGATAGCAATTTATTGTATCTGAGGAGGAGAGAAAATAGGAGATGTTGTCGTTATCATAAATCATATCATAAATAACCAAAAAGAAGCACTCAGAATGCTTGTGTGGTGTCAGCAGCCAGCCCTTTCAGGCACGTTAGCCCAGGATGTAGCCATAGAAACCCCAAAGTCATATTGTTATCTTAGTCTTGCTGCAGTTATGTGATACAACAGTATTAACATATAGATTAGGCCATTTATA-ATACTAGGTGTCATTTCAGGTTGGCTTGC-TTATATAGC-ATTTTAGTTTTTTATTGTGCAGAATGTTAGTGAAATAGTGCTTTTTTACAGGTCACTGTGTAT---CTAGGTATTTAT--AATGCATCATATTCACTCTACTATTAATTTCTTGTG--',
 'chrm=chr1_1039,model=HKY:MusMuscu,CATGAAAGGGTTCCAAGATAGAATTTAAATTTAAACTGATTCTTTAGAGATAAGCTTTTTAAAAATCCACCTCTTTGCTTCAAGTCAGGATTTTAAAAAT-AAATTTCCCCTTTAAATCTCTTGAGC-GATTTTCATTTTTTGCAAGG-CCGCACTGCTGGTGTCCTGGAAAAATGGAGAGTTAGGGCTTTATCAGTTGGTGCCCTGAGGTATGTGGCAAAGAAGCTAAATCCTTGAAGATTAGCAAAAAAAAATGCAGCACATGGCAATAAGGGGCTTATATGACTACTAGTGTTAGCTCAGGTGAAAGGAGTATTGCTTATACATAAAACTGGACAAATCCACTGACAATGTATTTTTAGTTAGCTACAAAGGAGTCTTATATTGCCGGACTTTCCCTC--TTTGGTGAATTGCGAAGCTTAGT,GorGoril,CATGAAAGAGTTCCAAGATAGAATTTAAATTTGAACTGATTCTTTATAGATAAGCTTTTTAAAAATCCACCTCTTTGCTTCAAGTCAGGATTTTAAAAAT-AAATTTCCCCTTTAAATCTCTTGAGC-GATTTTCATTTTTTGCAAGG-CCGCACTGCTGGTGTCCTGGAAAAATGGAGAGTTAGAGCTTTATCAGTTGGTGCCCTGAGGTATGTGGCGAAGAAGCTAAATCCTTGAAGATTAACAAAAAAAAATGCAGCACATGGCAATAAGGGGCTTATATGACTACTAGTATTAGTTCAGGTGAAAGAAGTATTGCTTCTACATAAAACTGGACAAATCCACTGACAATGTATTTTTAGTTAGCTACAAAGGAGTCTTATATTGCTGGACTTTACCTC--TTTGGTGAATTGGGAAGCTTAGT,PanTrogl,CATGAAAGAGTTCCAAGATAGAATTTAAATTTGAACTGATTCTTTATAGATAAGCTTTTTAAAAATCCACCTCTTTGCTTCAAGTCAGGATTTTAAAAAT-AAATTTCCCCTTTAAATCTCTTGAGC-GATTTTCATTTTTTGCAAGG-CCGCACTGCTGGTGTCCTGGAAAAATGGAGAGTTAGAGCTTTATCAGTTGGTGCCCTGAGGTATGTGGCGAAGAAGCTAAATCCTTGAAGATTAAC-AAAAAAAATGCAGCACATGGCAATAAGGGGCTTATATGACTACTAGTATTAGTTCAGGTGAAAGAAGTATTGCTTCTACATAAAACTGGACAAATCCACTGACAATGTATTTTTAGTTAGCTACAAAGGAGTCTTATATTGCTGGACTTTACCTC--TTTGGTGAATTGGGAAGCTTAGT',
 'chrm=chr1_1057,model=HKY:MusMuscu,--TACAAACTAGTTTCTGTGTTTGGGTAAAAAAGAGAGAATAAAAGTTTAATAAACGTGTGCTGTGTCGTGGATGACTTTGTGGCGACCTCTGTTGCAAATGGCCTATGCATGCGGAATAATGGCCTCCTTGCAGAGAGGGAAATAGCTTAGTGCTATCATCGTTGCCTCGGTAACCATCAGAG-TCCCCATCTAGCAAGAGA-GAAATGAGTTATG-AAAAGGCTTGAGTGAAGAAGGGA--TTAACACATGATCCCAGGGACAGTATGTCAATCAAAATCAACTGAG-AAATTACACT-CTCTATTTGAGAATTTCTCCGGTACCCAAAC-GGTAGTGAGGTT-AGCATTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTC-----------------TCTC-TCTCTCTC,GorGoril,-GTGCAAACTAGTTTCTGTG-TTGGGGGAAAAAAAGCGAATAAAAGTTTAATAAACGTGTGCTGTGTTGTGGATGACTTTGTGACGACCTCTGTTGCAAATGGCCTATGCATGCGGAATAATGGTCTCCTTGCAGAGAGGGAAATAGCTTAGTGCTATCATCGTTGCCTCGGTAACCATCAGAG-TCCCCATCTAGCAAGAGA-GAAATGAGTTATG-AAAAGGCTTGAGTGAAGAAGGGA--TTAACACATGATCCCAGGGACAGTGTGTCAATCAAAATCAACTGAG-AAATTTCACTGCTCTATTTGAGAATTTCTCCGGAACCCGCAT-GGTAGTGGAGTT-GGGTGCGTA--TTTTTTTTCCTCTCTCTCTTTTTTTCATTTGAGAAAGAAGGGGAAAGTATCATCGTGATT,PanTrogl,-GTGCAAACTAGTTTCTGTG-TTGGGGGAAAAAAAGCGAATAAAAGTTTAATAAACGTGTGCTGTGTTGTGGATGACTTTGTGACGACCTCTGTTGCAAATGGCCTATGCATGCGGAATAATGGTCTCCTTGCAGAGAGGGAAATAGCTTAGTGCTATCATCGTTGCCTCGGTAACCATCAGAG-TCCCCATCTAGCAAGAGA-GAAATGAGTTATG-AAAAGGCTTGAGTGAAGAAGGGA--TTAACACATGATCCCAGGGACAGTGTGTCAATCAAAATCAACTGAG-AAATTACACTGCTCTATTTGAGAATTTCTCCGGAACCCGCAT-GGTAGTGGAGTT-GGGTGCGTATTTTTTTTTTCCTCTCTCTCTTTTTTTCATTTGAGAAAGAAGGGGAAAGTATCATCGTGATT']

def mrAIC(self, key, line):
    """ Run mr-aic.pl on the each one-liner in the input file."""
    # weirdness parsing key, value
    oneliner = line.split("\t")[-1]
    if ":" in line:
        args_dict = self.split_oneliner(oneliner)
        # convert line to phylip
        phylip = self.oneliner_to_phylip(oneliner.split(":")[-1])
    else:
        # convert line to phylip
        phylip = self.oneliner_to_phylip(oneliner)
    phyml = Phyml(phylip)
    model, tree = phyml.best_model()
    args_dict = {'model': model}
    oneliner = "%s:%s" % (self.make_tree_name(args_dict), oneliner.split(":")[-1])
    if gt:
        for line in aic_fin:
            yield key, "tree '%s' = [&U] %s" % (self.make_tree_name(args_dict), tree)
    else:
        # give everything the same key so it can be reduced to a
        # 'oneliner' suitable for bootstrapping
        yield 1, oneliner

def get_bootstraps(sample, replicates=1):
    """Create boostrapped replicates of an a numpy array.

    Returns list

    """
    if not isinstance(sample, np.ndarray):
        try:
            sample = np.array(sample)
        except:
            raise TypeError("bootstrap() input must be list or numpy.ndarray")
    #replicates = int(replicates)
    size = len(sample)
    if replicates == 1:
        return sample[np.random.random_integers(0, size - 1, size)].tolist()
    else:
        return [sample[np.random.random_integers(0, size - 1, size)].tolist()
                        for i in xrange(replicates)]

def oneliner_to_array(line):
    """Convert oneliner to 2d numpy array.

    Returns tuple(list, array)

    """
    seqs = line.split(",")
    label_seqs = zip(seqs[:-1:2],seqs[1::2])
    taxa, bases = [], []
    for taxon, seq in label_seqs:
        bases.append(list(seq))
        taxa.append(taxon.strip())
    return taxa, np.array(bases)

def split_one_liner(line, args_dict={}, return_locus=False):
    """From a oneline, split locus/taxon info into dict and locus into string

    Returns dict or tuple(dict, locus)

    """
    args, locus = line.split(":")
    args = args.split(',')
    for item in args:
        dict_key, value = item.split("=")
        args_dict[dict_key] = value
    if not return_locus:
        return args_dict
    else:
        return args_dict, locus

def array_to_oneliner(taxa, bases):
    """Convert array of array of taxa and an array of bases to one-liner. 

    Returns string

    """
    # join is faster and cleaner than concatenate
    # list comp is faster than for
    # enumerate gets index of taxon from taxa
    oneliner = ','.join([','.join([taxa[count], seq.tostring()]) for count,seq in enumerate(bases)])
    return oneliner


def make_tree_name(args_dict):
    """Converts dictionary of arguement into a sorted string with the
    format:  """
    # join is faster and cleaner than concatenate
    # list comp is faster than for
    name = ','.join(['='.join([pair[0], pair[1]]) for pair in sorted(args_dict.items())])
    return name


def make_bootstrap_replicates(self, key, line):
    loci = line.strip().split(';')
    loci = loci[:-1]
    # first, bootstrap across loci
    multilocus_bstrap = get_bootstraps(loci)
    # second, bootstrap bases within loci
    for locus in multilocus_bstrap:
        # split keys from alignments
        args_dict, locus = split_one_liner(locus, return_locus=True)
        # convert alignments to arrays
        taxa, numpy_alignment = oneliner_to_array(locus)
        # transpose so we bootstrap by columns
        bases_by_col = np.column_stack(numpy_alignment)
        # bootstrap by columns
        shuffled = get_bootstraps(bases_by_col)
        # transpose back to rows of sequences
        shuffled = np.column_stack(shuffled)
        # convert back to oneliner
        oneliner = array_to_oneliner(taxa, shuffled)
        oneliner = "%s:%s" % (make_tree_name(args_dict), oneliner)
        yield key, oneliner

if __name__ == '__main__':
    make_bootstrap_replicates()