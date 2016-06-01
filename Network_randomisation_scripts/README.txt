Description and Basic usage for all scripts in Network_randomisation_scripts directory

########################################################################################################################################

correlation_randomisation_test.pl

DESCRIPTION:
Randomisation based test to determine if the mean correlation (across edges) of a gene network is >= than can be expected by chance,
controlling for feature (eg. intron) length.

This test works as follows:
1). Calculate the mean correlation value for the real gene network.
2). For each gene in the network, choose another at random, but matched for feature length class.
3). Calculate the mean correlation value for the network of re-sampled genes.
4). Repeat this N times.

This produces a p-value which represents the proportion of times the mean correlation of the randomised network is >= that of
the real network.

BASIC USAGE:
correlation_randomisation_test.pl -i Network.tsv -c Correlation_matrix.tsv -b Intron_positions.bed

INPUT FILE EXAMPLES:

Example of Network.tsv:
This file should be a list of edges, a tab delimited list where each row contains a gene-pair connected by an edge.

Example:
gene_A	gene_B
gene_A	gene_C
gene_C	gene_D

Example of Correlation_matrix.tsv
This should be a tab-delimited, square symmetrical matrix, or alternatively the lower triangle of the matrix.
The first column should contain gene IDs.

Example:
gene_A	0
gene_B	0.1	0
gene_C	-0.4	-0.03	0
gene_D	-0.2	-0.1	0.5	0

Example of Intron_positions.bed
A bed file containing genomic position of gene features to be used in this analysis (E.g. Intron positions).
This is used to calculate their length, which is used as a control in the randomisation test.

Requires BED fields are:
<Chromosome>    <Start_posiiton>  <End_position>   <Gene_ID>

Subsequent fields are ignored.

########################################################################################################################################