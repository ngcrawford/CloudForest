#!/usr/bin/env python
# encoding: utf-8

import numpy
import dendropy
import itertools


Tree1 = "(((((B,A2),C),A1),(D1,D2)),(F,E))"
Tree2 = "(((((B,C),A2),(D2,D1)),(A1,F)),E)"

dpTree1 = dendropy.Tree()
dpTree1.read_from_string(Tree1,'newick')

dpTree2 = dendropy.Tree()
dpTree2.read_from_string(Tree2,'newick')

    
def internode_distance(taxa_pair, tree):
    """Calculates the number of nodes that connect a pair of taxa."""
    
    counter_on = False
    counter = 0
    path = []
    start_node = None
    internode_distance = 0
    previous_node_level = None
    levels = []
    
    for nd in tree.preorder_node_iter(): 
        if nd.is_leaf():
            if nd.taxon.label in taxa_pair and nd.taxon.label != start_node and start_node != None:
                internode_distance = counter
                counter = False
        
            if nd.taxon.label in taxa_pair:
                start_node = nd.taxon.label
                counter_on = True
    
        if counter_on == True:
            levels.append(nd.level())
            counter += 1
    
    test = internode_distance(('B','F'),dpTree1)
    start = float(test[0])
    end = float(test[-1])
    mrca = float(min(test[1]))
    internode_distance = (start-mrca) + (end-mrca)
    
    return (internode_distance)

for pair in itertools.combinations(taxa,2):
    print pair, (internode_distance(pair,dpTree1) + internode_distance(pair,dpTree2))/2
    

    
