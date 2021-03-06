README.txt
========================
Phylogenetics
Beth Ginsberg
========================

This codebase consists of two Python files:
 - distance_matrix.py
 - phylo_tree.py

In combination, these two programs convert a DNA or RNA multiple sequence alignment into a phylogenetic tree.
They depend on the following packages:
 - argparse
 - numpy
 - networkx
 - sys
 - string
Please ensure these are installed before proceeding.

========================
distance_matrix.py
========================
Invoke python distance_matrix.py -h for usage.

Input file is an MSA in Pearson / FASTA format.
This can be generated by running ClustalO (https://www.ebi.ac.uk/Tools/services/web_clustalo/toolform.ebi) to perform the alignment.
Then convert the output format using MView (https://www.ebi.ac.uk/Tools/msa/mview/).
See alpha_globin.txt for an example input file.

The output will be a distance matrix, preceded by the names of the sequences in the alignment.
Distance is calculated using hamming distance, with the option to adjust using the Jukes Cantor model.

Example usage:
python distance_matrix.py --jc examples/real_life_alignment/alpha_globin.txt examples/real_life_alignment/alpha_globin_matrix.txt

========================
phylo_tree.py
========================
Invoke python phylo_tree.py -h for usage.

Input file is the output of distance_matrix.py (or follows the same format).

Example usage:
python phylo_tree.py examples/real_life_alignment/alpha_globin_matrix.txt examples/real_life_alignment/alpha_globin_tree.xml --neighbor_joining

========================
Cytoscape hints:
========================
The best way to view the GML output files of phylo_tree.py is with Cytoscape.
The recommended layout to use is YFiles Tree Layout, which can be installed inside the Cytoscape App Manager.
Use the included styles.xml to apply appropriate styling for phylogenetic trees.
