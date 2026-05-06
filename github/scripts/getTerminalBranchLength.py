#!/usr/bin/env python3

import sys
import os
import numpy as np
from scipy import stats
from ete3 import Tree
from Bio import SeqIO

def get_ancestors(tree):
    return list(set(leaf.up for leaf in tree.iter_leaves()))

def get_terminals(treefile):
    try:
        tetree = Tree(treefile)
    except:
        print(f'{treefile} could not be read')
        return
    ancestors = get_ancestors(tetree)
    erroneous = []
    for node in ancestors:
        pair = node.children
        dists = {pair[0].dist: pair[0], pair[1].dist: pair[1]}
        mind, maxd = min(dists.keys()), max(dists.keys())
        if (mind == 0 and maxd > 0.1):
            shorter = dists[mind]
            if shorter.name != '':
                erroneous.append(shorter.name)
        elif mind > 0 and maxd/mind > 100:
            shorter = dists[mind]
            if shorter.name != '':
                erroneous.append(shorter.name)
    prune_list = []
    for leaf in tetree.iter_leaves():
        if leaf.name not in erroneous:
            prune_list.append(leaf.name)
    tetree.prune(prune_list)
    return {(leaf.name, leaf.dist) for leaf in tetree.iter_leaves()}

def main():
    treefile = sys.argv[1]
    terminals = get_terminals(treefile)
    for i in terminals:
        print(i[0]+"\t"+str(i[1]))

if __name__ == '__main__':
    main()
