import dendropy

def topology_search(tree, query_tree, symmetric_distance=False):
	"""Given a tree and a query tree this function tests if the query tree is present or identical to the 
	tree. If 'symmetric_distance' is set to True the symmatric distance is returned. 

	Both trees should be tree objects and contain idential taxon labels.

	Example:
	--------

	>>> tree = dendropy.Tree() 
	>>> query_tree = dendropy.Tree()
	>>> tree.read_from_string('(((A,C),D),B,(E,F))', schema='newick')
	>>> query_tree.read_from_string('(E,F)', schema='newick')
	>>> print topology_search(tree, query_tree)
	True
	>>> query_tree.read_from_string('(D,(E,F))', schema='newick')
	>>> print contains_topology(tree, query_tree)
	False
	"""

	# there's probably a better/faster way to do this with bit_masks
	sub_tree_bitmask = tree.taxon_set.get_taxa_bitmask(labels=query_tree.taxon_set.labels())
	sub_tree_mrca_node = tree.mrca(split_bitmask=sub_tree_bitmask)
	sub_tree_newick = sub_tree_mrca_node.as_newick_string()
	sub_tree = dendropy.Tree()
	sub_tree.read_from_string(sub_tree_newick,schema='newick')
	sub_tree.ladderize(ascending=True) # may not be necessary
	sd = sub_tree.symmetric_difference(query_tree)
	
	if symmetric_distance == True:
		return sd

	elif sd == 0:
		return True
	
	else:
		return False

def subtree_percentage(trees, sub_tree):
    matches = 0
    for count, tree in enumerate(trees):
        if topology_search(tree, subtree, symmetric_distance=True) == 0:
            matches += 1
    percentage = float(matches)/float(len(trees))
    return (percentage, matches, trees)

def reroot_trees(trees, root):
    """Supply a treelist and a taxa label at which to root each tree. Returns a treelist with each tree 
    rerooted at the same tip/label."""
    new_tree_list = dendropy.TreeList()
    for tree in trees:
        node_root = tree.find_node_with_taxon_label(root)
        tree.reroot_at_edge(node_root.edge, update_splits=False)
        tree.ladderize(ascending=True)
        new_tree_list.append(tree)
    return new_tree_list

trees = dendropy.TreeList()
trees.read_from_path('/Users/ngcrawford/Downloads/reptiles-10-species.oneliners.AICc.genetrees.trees.nex','nexus')
trees = reroot_trees(trees,'HomSapie')


archosaurs_w_turtles = ['CroPoros','AllMissi','GalGallu','ZebFinch', 'ChrPicta', 'PelSubru']
birds_w_turtles = ['GalGallu','ZebFinch', 'ChrPicta', 'PelSubru']
crocs_w_turtles = ['CroPoros','AllMissi', 'ChrPicta', 'PelSubru']
lepidosaurs_w_turtles = ['AnoCarol', 'PanGutta', 'SphTuata','ChrPicta', 'PelSubru']
squamates_w_turtles = ['AnoCarol', 'PanGutta','ChrPicta', 'PelSubru']
mammals_w_turltes = ['HomSapie','ChrPicta', 'PelSubru' ]

lepidosaurs = ['AnoCarol', 'PanGutta', 'SphTuata']
archosaurs = ['CroPoros','AllMissi','GalGallu','ZebFinch']
squamates = ['AnoCarol', 'PanGutta']
birds = ['GalGallu','ZebFinch']
crocs = ['CroPoros','AllMissi']
turtles = ['ChrPicta', 'PelSubru']

print 'archosaurs_w_turtles', trees.frequency_of_split(labels=archosaurs_w_turtles)
print 'lepidosaurs_w_turtles', trees.frequency_of_split(labels=lepidosaurs_w_turtles)
print 'birds_w_turtles', trees.frequency_of_split(labels=birds_w_turtles)
print 'crocs_w_turtles', trees.frequency_of_split(labels=crocs_w_turtles)
print 'squamates_w_turtles', trees.frequency_of_split(labels=squamates_w_turtles)
print 'mammals_w_turltes', trees.frequency_of_split(labels=mammals_w_turltes)

print 'archosaurs', trees.frequency_of_split(labels=archosaurs)
print 'lepidosaurs', trees.frequency_of_split(labels=lepidosaurs)
print 'squamates', trees.frequency_of_split(labels=squamates)
print 'birds', trees.frequency_of_split(labels=birds)
print 'crocs', trees.frequency_of_split(labels=crocs)
print 'turtles', trees.frequency_of_split(labels=turtles)


# expected_tree_newick = '[&R] (HomSapie, ((((CroPoros,AllMissi),(GalGallu,ZebFinch)),(ChrPicta,PelSubru)),(SphTuata,(AnoCarol,PanGutta))))'
# expected_tree = dendropy.Tree()
# expected_tree.read_from_string(expected_tree_newick,'newick')
# expected_tree.as_ascii_plot()
# print expected_tree.as_ascii_plot()
# for node in expected_tree.preorder_internal_node_iter():
#     subtree = dendropy.Tree()
#     subtree_newick = node.as_newick_string()
#     subtree.read_from_string(subtree_newick,'newick')
#     print subtree_newick, subtree_percentage(trees, subtree)
#      
# 
# 
# #query_tree_newick = '(AnoCarol, PanGutta)'
# query_tree_newick = '((SphTuata,(AnoCarol,PanGutta)),(ChrPicta,PelSubru))'
# #query_tree_newick = '(SphTuata, (AnoCarol, PanGutta))'
# #query_tree_newick = '(((CroPoros,AllMissi),(GalGallu,ZebFinch)), (ChrPicta, PelSubru))'
# subtree = dendropy.Tree()
# subtree.read_from_string(query_tree_newick,'newick')
# print subtree_percentage(trees, subtree)








