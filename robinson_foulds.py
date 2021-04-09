from ete3 import Tree
import copy
from collections import Counter


#get leaf names of a tree (bipartition), for rooted trees,  add root as additional leaf to count clades
def get_leaves(nd, root=None):
    leaves = []
    for l in nd:
        leaves.append(l.name)
    if root:
        leaves.append(root)
    return leaves


#compare two unordered lists
def compare(s, t):
    return Counter(s) == Counter(t)


#get leafsets of all bipartitions of an unrooted tree
def get_leafset(t1):
    count2 = 0
    count3 = 0
    t1_leafsets = []

    for node in t1.traverse("levelorder"):
        count2 +=1

        #skip leaves
        if len(node.children) ==0:
            continue

        #rename children with no name (to keep track)
        for child in node.children:
            if child.name == "":
                count3 +=1
                child.name = str(count3)

        #append partitions of both children
        t1_leafsets.append(get_leaves(node.children[0]))
        t1_leafsets.append(get_leaves(node.children[1]))

        #skip root (only for unrooted trees)
        # if count2 ==1:
        #     continue

        tree_cpy1 = copy.deepcopy(t1)
        tree_cpy2 = copy.deepcopy(t1)
        #print("name", node.children[0].name)

        #break tree into two removing a child node
        ch1 = tree_cpy1.search_nodes(name=node.children[0].name)[0]
        ch2 = tree_cpy2.search_nodes(name=node.children[1].name)[0]
        ch1.detach()
        ch2.detach()
        #tree_cpy1 = tree_cpy1.remove_child(ch1)
        #tree_cpy2 = tree_cpy1.remove_child(ch2)
        t1_leafsets.append(get_leaves(tree_cpy1, "r"))
        t1_leafsets.append(get_leaves(tree_cpy2, "r"))

    #t1.show()
    return (t1_leafsets)
    # for leaf in t1:
    #     print(len(leaf.get_children()))


#get robinson foulds distance of one tree in comparison to the other
def robinson_foulds_one_way(lf1, lf2):
    c = 0
    for i in lf1:
        #print(len(i))
        found = False
        for j in lf2:
            if compare(i,j):
                found = True
                break
        if not found:
            c += 1
    return c


#get robinsond foulds distance between 2 trees
def robinson_foulds_dist(lf1, lf2):
    c1 = robinson_foulds_one_way(lf1, lf2)
    c2 = robinson_foulds_one_way(lf2, lf1)
    return ((c1 + c2)/2)


#compare a set of trees
def compare_trees(trees):
    for i,tree in enumerate(trees):
        for j,tree2 in enumerate(trees):
            if (i == j):
                continue
            lf1 = get_leafset(tree)
            lf2 = get_leafset(tree2)

            print("rooted robinson foulds distance between ", i, j, "is", robinson_foulds_dist(lf1, lf2))


# Load a tree structure from a newick file.
t0 = Tree("LysR family transcriptional regulator_sequences.nwk")
t1 = Tree("single-stranded DNA-binding protein.nwk")
t2 = Tree("efflux transporter outer membrane subunit_sequences.nwk")
t3 = Tree("helix-turn-helix domain-containing protein_sequences.nwk")


#remove protein name
for tree in [t0, t1, t2, t3]:
    for leaf in tree:
        leaf.name = leaf.name.split("_")[0] + "_" + leaf.name.split("_")[1]
        #print (leaf.name)

compare_trees([t0, t1, t2, t3])

