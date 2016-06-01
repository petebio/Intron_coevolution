#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

########################################################################################################################################
#                                                                                                                                      #
# Convert a BED file that contains genomic positions for exons into a BED file containing intron positions.                            #
#                                                                                                                                      #
# BED file column order is:                                                                                                            #
# <Chromosome> <Exon start position> <Exon end position position> <Gene ID> <Transcript ID>                                            #
#                                                                                                                                      #
# See help for file example. Columns after Transcript ID will be ignored.                                                              #
#                                                                                                                                      #
# Author: Peter Keane.                                                                                                                 #
# Contact: peterakeane@gmail.com                                                                                                       #
#                                                                                                                                      #
########################################################################################################################################

my$exon_bed;
my$intron_bed;
my$help = 0;

GetOptions("infile=s" => \$exon_bed, "outfile=s" => \$intron_bed, "help|?" => \$help);

pod2usage(verbose => 99) if $help;

# Check required arguments.

my$args = 1;
if(not defined $exon_bed){
	print "\nMissing argument -i";
	$args = 0;
}

if(not defined $intron_bed){
	print "\nMissing argument -o";
	$args = 0;
}

# If required argument is missing, $args will be 0.
# Print usage and exit.
if($args == 0){
	print "\n\n"; # Some formatting.
	pod2usage(-verbose => 1);
}

########################################################################################################################################

# Read BED file of exon positions.
open(IN, $exon_bed) || die "Cannot open file $exon_bed: $!.\n";

my(%exon_tracking, %gene_transcript_map);
while(my$entry = <IN>){
	chomp $entry;
	my@fields = split("\t", $entry);
	
	my$chr = $fields[0];
	my$start = $fields[1];
	my$end = $fields[2];
	my$gene_id = $fields[3];
	my$transcript_id = $fields[4];
	
	$gene_transcript_map{$chr}{$gene_id}{$transcript_id} = 1;
	$exon_tracking{$transcript_id}{$start} = $end;
}

close IN;

########################################################################################################################################

open(OUT, ">", $intron_bed);

# Process each chromosome.
foreach my$chr (sort {$a cmp $b} keys %gene_transcript_map){
	# Process each gene.	
	foreach my$gene_id (keys %{$gene_transcript_map{$chr}}){
		# Process each transcript.
		foreach my$transcript_id (keys %{$gene_transcript_map{$chr}{$gene_id}}){
			my@exon_end_pos = sort {$a <=> $b} keys %{$exon_tracking{$transcript_id}};
			
			# Skip transcripts with only one exon (no introns).
			next if scalar @exon_end_pos == 1;
			
			# Start at the start position of the second exon.
			# This is the end of the first intron...
			for(my$i = 1; $i < scalar @exon_end_pos; $i++){
				# The start of the intron is the end of the upstream exon.
				my$upstream_exon_start = $exon_end_pos[$i - 1];
				my$upstream_exon_end = $exon_tracking{$transcript_id}{$upstream_exon_start};
				
				# Print intron position to output file.
				print OUT $chr,"\t",$upstream_exon_end,"\t",$exon_end_pos[$i],"\t",$gene_id,"\t",$transcript_id,"\n"; 
			}
		}
	}
}

close OUT;

########################################################################################################################################

exit;

########################################################################################################################################

__END__

=head1 DESCRIPTION

Script to convert a BED file containing exon positions into intron positions.

=head1 INPUT file example (Tab delimited)

=over 8

=item
1	935072	935552	ENSG00000188290	ENST00000428771

=item
1	934906	934993	ENSG00000188290	ENST00000428771

=item
1	934342	934812	ENSG00000188290	ENST00000428771

=item
1	935246	935491	ENSG00000188290	ENST00000304952

=item
1	935072	935167	ENSG00000188290	ENST00000304952

=item
1	934906	934993	ENSG00000188290	ENST00000304952

=item
1	934344	934812	ENSG00000188290	ENST00000304952

=item
1	934906	935476	ENSG00000188290	ENST00000481869

=item
1	934346	934812	ENSG00000188290	ENST00000481869

=item
1	935246	935361	ENSG00000188290	ENST00000484667

=item
1	934906	934993	ENSG00000188290	ENST00000484667

=item
1	934350	934812	ENSG00000188290	ENST00000484667

=back

=head1 AUTHOR

=over 8

=item
Peter Keane

=item
School of Mathematics, Statistics and Applied Mathematics

=item
National University of Ireland, Galway

=item
Email: peterakeane@gmail.com

=back

=head1 OPTIONS

=over 8

=item B<-i>
Exon bed file.

=item B<-o>
Output file.

=item B<-h>
Print help and exit.

=back
