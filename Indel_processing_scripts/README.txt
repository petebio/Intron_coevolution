Description and Basic usage for all scripts in Indel_processing_scripts directory.

Example workflow for this analysis.

########################################################################################################################################

get_indels_from_emf.pl

DESCRIPTION:
Read an Ensembl Multi Format (EMF) file containing whole genome alignments, and find all insertion/deletion events that can
be mapped to a given phylogenetic tree.

Also provide a reference species for which genomic positions should be used.

This program was designed using the 15 Eutherian mammals EPO dataset from Ensembl Compara release 74.
Will probably work with other similar datasets, but is untested.

EMF files from Compare release 74 can be found here:
ftp://ftp.ensembl.org/pub/release-74/emf/ensembl-compara/epo_15_eutherian/

BASIC USAGE:
perl get_indels_from_emf.pl -i Alignment.emf -t Phylogeny.nwk -o Indels.bed

Alignment.emf is an EMF file from Compara (see link above).

Phylogeny.nwk is a text file containing the reference phylogenetic tree, in newick format. 
The 9 mammal tree used in this study is:

((((Fcat:0,Cfam:0):0,Ecab:0):0,(Oari:0,Btau:0):0):0,(((Ptro:0,Hsap:0):0,Ggor:0):0,Pabe:0):0):0;

Indels.bed is the output file (BED format).

########################################################################################################################################

map_indels_to_introns.pl

DESCRIPTION:
Map Insertion/Deletion events to intron positions.

Expects Indel file (BED formatted) produced by get_indels_from_emf.pl

Maps Indels in this file to Intron positions in another BED file. 
This file can be made using the convert_exon_to_intron_bed.pl script in the Misc_scripts directory of this repository.

BASIC USAGE:
perl map_indels_to_introns.pl -i Indels.bed -b Introns.bed -o Intronic_indels.bed

########################################################################################################################################

calculate_proportional_change.pl

DESCRIPTION:
Calculate the proportional change in intron content for all genes provided in the indel BED file.
The indel BED file can be created using map_indels_to_introns.pl

BASIC USAGE:
perl calculate_proportional_change.pl -i Intronic_indels.bed -b Introns.bed -t Phylogeny.nwk -o Proportional_change.tsv

Intronic_indels.bed is the BED file created by map_indels_to_introns.pl

Introns.bed is a BED file containing genomic positions for all introns. This is used to calculate intron lengths.
This file can be made using the convert_exon_to_intron_bed.pl script in the Misc_scripts directory of this repository.

Phylogeny.nwk is a text file containing the reference phylogenetic tree, in newick format. 
The 9 mammal tree used in this study is:

((((Fcat:0,Cfam:0):0,Ecab:0):0,(Oari:0,Btau:0):0):0,(((Ptro:0,Hsap:0):0,Ggor:0):0,Pabe:0):0):0;

Proportional_change.tsv is the output file.
This is a tab delimited matrix where rows are genes, and columns are branches of the phylogenetic tree.

########################################################################################################################################

partial_correlation_pearson.pl

DESCRIPTION:
Calculate gene-gene partial Pearson correlations for intron length change.
Controls for median change along each branch of the phylogenetic tree.

Expects matrix of proportional intron length changes from calculate_proportional_change.pl

BASIC USAGE:
perl partial_correlation_pearson.pl -i Proportional_change.tsv -o Pearson_correlation.tsv

Proportional_change.tsv is the tab delimited matrix of intron change produced by calculate_proportional_change.pl

Pearson_correlation.tsv is the output file.
The lower triangle of a symmetrical matrix containing gene-gene correlation values.
Retaining only the lower triangle saves on hard-disk space...

########################################################################################################################################

partial_correlation_spearman.pl

DESCRIPTION:
Calculate gene-gene partial Spearman correlations for intron length change.
Controls for median change along each branch of the phylogenetic tree.

Expects matrix of proportional intron length changes from calculate_proportional_change.pl

BASIC USAGE:
perl partial_correlation_spearman.pl -i Proportional_change.tsv -o Spearman_correlation.tsv

Proportional_change.tsv is the tab delimited matrix of intron change produced by calculate_proportional_change.pl

Spearman_correlation.tsv is the output file.
The lower triangle of a symmetrical matrix containing gene-gene correlation values.
Retaining only the lower triangle saves on hard-disk space...

########################################################################################################################################

correlation_to_simpson.pl

DESCRIPTION:
Calculate the Simpson coefficient from the gene-gene correlation matrix. 

BASIC USAGE:
perl correlation_to_simpson.pl -i Pearson_correlation.tsv -t 0.7 -o Simpson.tsv

Pearson_correlation.tsv is the gene-gene correlation matrix produced by partial_correlation_pearson.pl
Could also use spearman correlations from partial_correlation_spearman.pl

0.7 is the correlation threshold to use while calculating the Simpson coefficient. This can be specified by the user.

Simpson.tsv is the output file. An n x n square matrix of Simpson coefficient values.

########################################################################################################################################

EXAMPLE WORKFLOW:

=> Step 1: Infer all Insertion/Deletion events.
perl get_indels_from_emf.pl -i Alignment.emf -t Phylogeny.nwk -o Indels.bed

=> Step 2: Map Insertion/Deletion events to introns.
perl map_indels_to_introns.pl -i Indels.bed -b Introns.bed -o Intronic_indels.bed

=> Step 3: Calculate the proportional change in intron content along each branch of the phylogeny.
perl calculate_proportional_change.pl -i Intronic_indels.bed -b Introns.bed -t Phylogeny.nwk -o Proportional_change.tsv

=> Step 4: Calculate partial Pearson correlations for all gene pairs.
perl partial_correlation_pearson.pl -i Proportional_change.tsv -o Pearson_correlation.tsv

=> Step 5: Calculate Simpson coefficient for all gene pairs.
perl correlation_to_simpson.pl -i Pearson_correlation.tsv -t 0.7 -o Simpson.tsv

########################################################################################################################################