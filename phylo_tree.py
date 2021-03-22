#
# Phylogenetic tree generator
#
# Beth Ginsberg
# March 16, 2021
#
# This code generates a phylogenetic tree using one of the following algorithms:
# UPGMA (Unweighted Pair Group Method with Arithmetic Mean) by Sokal and Michener.
#   https://en.wikipedia.org/wiki/UPGMA
# WPGMA (Unweighted Pair Group Method with Arithmetic Mean) by Sokal and Michener.
#   https://en.wikipedia.org/wiki/WPGMA
# Neighbor Joining by Saitou and Nei
#   https://en.wikipedia.org/wiki/Neighbor_joining
#
# The tree is output to a file in GML format, which can then be visualized in Cytoscape.
#

import numpy as np
import sys
import networkx as nx
import argparse
import string


# Recursively calculate the total weight from the given node to a leaf.
# The method of tree construction ensures that it makes no difference which child gets traversed, so we keep left,
# like good South African drivers.
def distance_to_leaf(node, dist) -> float:
    if not tree.has_node(node):
        return dist
    edges = list(tree.out_edges(node, data=True))
    if len(edges) == 0:
        return dist
    edge = edges[0]
    wt = dist + edge[2]['weight']
    child = edge[1]
    return wt + distance_to_leaf(child, dist)


# Find the minimum non-zero entry of a matrix. Returns tuple (row, column).
def find_min_non_zero(inp_matrix: np.array) -> ():
    least = sys.maxsize
    position = (0, 0)
    for en, val in np.ndenumerate(inp_matrix):
        if (val != 0.0) and (val < least):
            least = val
            position = en
    return position


# Add nodes to the phylogeny tree.
def add_tree_nodes(inp_matrix: np.array, position: (), inp_labels: [], nj=False) -> {}:

    left = inp_labels[position[0]]
    right = inp_labels[position[1]]

    # The parent node will just be the concatenation of the its children.
    parent = left + right

    # Calculate the weight / distance between nodes.
    dist = inp_matrix[position[0], position[1]]
    if nj:
        n = inp_matrix.shape[0]
        if n == 2:
            tree.add_edge(left, right, weight=round(dist, 2))
            return {left: dist / 2.0, right: dist / 2.0}
        else:
            lw = (dist / 2.0) + (1.0 / (2 * (n - 2))) * (sum(inp_matrix[position[0]]) - sum(inp_matrix[position[1]]))
            rw = dist - lw
    else:
        half = dist / 2.0
        lw = half - distance_to_leaf(left, 0.0)
        rw = half - distance_to_leaf(right, 0.0)

    # Round to 2 decimal places for easy viewing.
    lw = round(lw, 2)
    rw = round(rw, 2)

    # Add nodes and new parent to the tree. Use labels to differentiate leaf and non-leaf nodes.
    if not tree.has_node(left):
        tree.add_node(left, leaf=True)
    if not tree.has_node(right):
        tree.add_node(right, leaf=True)
    tree.add_node(parent, leaf=False)

    # Add edges and store the weight.
    tree.add_edge(parent, left, weight=lw)
    tree.add_edge(parent, right, weight=rw)

    return {left: lw, right: rw}


# Combine rows and columns of closest sequences.
def collapse_rows(inp_matrix: np.array, position: (), inp_labels: [], weighted=True, nj=False, weights=None) -> (np.array, []):

    # Append a new row with the combined effect of the new node.
    new_row = np.empty((1, inp_matrix.shape[1]))
    columns = inp_matrix.shape[1]

    # Neighbor joining
    if nj:
        for i in range(columns):
            if i in position:
                new_row[0][i] = weights[inp_labels[i]]
            else:
                new_row[0][i] = (inp_matrix[position[0], i] + inp_matrix[position[1], i] - inp_matrix[position[0], position[1]]) / 2.0
    else:
        # Differentiate between weighted and unweighted algorithm
        if weighted:
            for i in range(columns):
                col = inp_matrix[:, i]
                new_row[0][i] = sum(col[list(position)]) / 2.0
        else:
            r1 = position[0]
            r2 = position[1]
            cluster1 = len(inp_labels[r1])
            cluster2 = len(inp_labels[r2])
            new_size = cluster1 + cluster2
            for i in range(columns):
                col = inp_matrix[:, i]
                new_row[0][i] = ((col[r1] * cluster1) + (col[r2] * cluster2)) / new_size

    inp_matrix = np.concatenate((inp_matrix, new_row), axis=0)

    # Append a new column that is the transpose of the new row
    new_row = np.append(new_row, [0.0])
    inp_matrix = np.concatenate((inp_matrix, np.transpose(new_row)[:, np.newaxis]), axis=1)

    # Delete the two old nodes.
    inp_matrix = np.delete(np.delete(inp_matrix, position, 0), position, 1)

    # Update the headers.
    new_label = ''.join(inp_labels[list(position)])
    inp_labels = np.delete(inp_labels, position)
    inp_labels = np.append(inp_labels, new_label)

    return inp_matrix, inp_labels


# Calculates the secondary matrix required by the neighbor-joining algorithm.
def calculate_joining_matrix(inp_matrix: np.array) -> np.array:

    ret_matr = np.zeros(inp_matrix.shape)
    n = inp_matrix.shape[0]
    for en, val in np.ndenumerate(inp_matrix):
        if en[0] == en[1]:
            continue
        if ret_matr[en] > 0:
            continue
        new_val = (n - 2) * val - sum(inp_matrix[:, en[0]]) - sum(inp_matrix[:, en[1]])
        ret_matr[en] = new_val
        ret_matr[en[1], en[0]] = new_val
    return ret_matr


def iterate_pgma(inp_matrix: np.array, inp_labels: [], weighted):

    while inp_matrix.shape[0] > 1:
        closest = find_min_non_zero(inp_matrix)
        add_tree_nodes(inp_matrix, closest, inp_labels)
        inp_matrix, inp_labels = collapse_rows(inp_matrix, closest, inp_labels, weighted)


# Perform UPGMA algorithm on the input matrix.
def upgma(inp_matrix: np.array, inp_labels: []):

    iterate_pgma(inp_matrix, inp_labels, False)


# Perform WPGMA algorithm on the input matrix.
def wpgma(inp_matrix: np.array, inp_labels: []):

    iterate_pgma(inp_matrix, inp_labels, True)


# Perform neighbor joining algorithm on the input matrix.
def neighbor_joining(inp_matrix: np.array, inp_labels: []):

    while inp_matrix.shape[0] > 1:
        matr = calculate_joining_matrix(inp_matrix)
        closest = find_min_non_zero(matr)
        wts = add_tree_nodes(inp_matrix, closest, inp_labels, True)
        inp_matrix, inp_labels = collapse_rows(inp_matrix, closest, inp_labels, nj=True, weights=wts)


if __name__ == '__main__':

    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Create a phylogenetic tree from a distance matrix using UPGMA or WPGMA.')
    parser.add_argument('inputfile', help='File to read. Must contain a valid distance matrix. This could be the output of distance_matrix.py.')
    parser.add_argument('outputfile', help='Filename for output. Output is in GML and can be rendered by Cytoscape.')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--upgma', dest='algorithm', const=upgma, default=upgma, nargs='?', help='Use UPGMA (default)')
    group.add_argument('--wpgma', dest='algorithm', const=wpgma, default=upgma, nargs='?', help='Use WPGMA. Default is UPGMA.')
    group.add_argument('--neighbor_joining', dest='algorithm', const=neighbor_joining, default=upgma, nargs='?', help='Use neighbor joining. Default is UPGMA.')
    args = parser.parse_args()

    # Read the input file into a matrix.
    with open(args.inputfile, 'r') as reader:
        orig_labs = reader.readline().split('||||')
        input_matrix = np.loadtxt(args.inputfile, skiprows=1)

    tree = nx.DiGraph()

    # Generate internal labels and map them back to the original labels.
    labs = np.array(list(string.ascii_uppercase)[0:np.shape(input_matrix)[0]])
    label_mapping = dict(zip(labs, orig_labs))

    # Call the specified algorithm on the input.
    args.algorithm(input_matrix, labs)

    # Replace original labels on the leaf nodes.
    nx.relabel_nodes(tree, label_mapping, copy=False)

    # Create the output file.
    nx.write_graphml(tree, args.outputfile)
