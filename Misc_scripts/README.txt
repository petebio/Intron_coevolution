Description and Basic usage for all scripts in Misc_scripts directory

########################################################################################################################################

convert_exon_to_intron_bed.pl

DESCRIPTION:
Convert a BED file that contains genomic positions for exons into a BED file containing intron positions.

BED file column order is:                                                                                                                      #
<Chromosome> <Exon start position> <Exon end position position> <Gene ID> <Transcript ID>

See script help (-h option) for file example. Columns after Transcript ID will be ignored. 

BASIC USAGE:
perl convert_exon_to_intron_bed.pl -i Exon_positions.bed -o Intron_positions.bed

########################################################################################################################################