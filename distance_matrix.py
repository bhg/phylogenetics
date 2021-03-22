#
# Phylogenetic tree generator
#
# Beth Ginsberg
# March 16, 2021
#
# This code generates a distance matrix from an MSA using hamming distance.
#
# The distance matrix is output in a file, which can then be read by phylo_tree.py to create a phylogenetic tree.
#

import numpy as np
import argparse


# Calculates the hamming distance between sequences in a Multiple Sequence Alignment.
# Returns a square distance matrix, the size of which is equal to the number of sequences in the alignment.
# The distance matrix will always have zeros on the diagonal and be symmetric.
def hamming_distance_matrix(seq: []) -> np.array:
    ret_arr = np.zeros((len(seq), len(seq)))
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            if i == j:
                continue
            if ret_arr[i][j] > 0:
                continue
            seq1 = list(sequences[i])
            seq2 = list(sequences[j])
            if len(seq1) != len(seq2):
                raise Exception('Invalid MSA: sequence lengths differ.')
            ret_arr[i][j] = len([p for p in range(len(seq1)) if seq1[p] != seq2[p]])
            ret_arr[j][i] = ret_arr[i][j]
    return ret_arr


# Adjusts the distance matrix according to the Jukes Cantor model.
def jukes_cantor(inp_array: np.array) -> np.array:
    try:
        return -0.75 * sequence_length * np.log(1 - (4 * inp_array) / (3 * sequence_length))
    except ValueError:
        print("Insufficient conservation to apply Jukes Cantor model. Please try a better alignment.")
        exit(2)


sequence_length = 0

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Create a distance matrix from an MSA using hamming distance.')
    parser.add_argument('inputfile', help='File to read. Must contain a valid MSA.')
    parser.add_argument('outputfile', help='Filename for output of distance matrix')
    parser.add_argument('--jc', action='store_true', help='Adjust using the Jukes Cantor model')
    args = parser.parse_args()

    labels = []
    sequences = []

    with open(args.inputfile, 'r') as reader:
        line = reader.readline()
        while line != '':
            if line.startswith('>'):
                labels.append(line.rstrip('\n'))
                line = reader.readline()
            else:
                seq = ''
                while line != '' and not line.startswith('>'):
                    seq += line.rstrip('\n')
                    line = reader.readline()
                sequences.append(seq)
                sequence_length = len(seq)

    # Calculate the hamming distance.
    arr = hamming_distance_matrix(sequences)

    # If we're adjusting, apply Jukes Cantor model to update the matrix.
    if args.jc:
        arr = jukes_cantor(arr)

    # Write the output file
    with open(args.outputfile, 'w') as writer:
        writer.write('||||'.join(labels) + '\n')
        for line in arr:
            writer.write(' '.join(map(str, line)) + '\n')
