#!/usr/bin/env perl
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

####################################################################################################################################
#                                                                                                                                  #
# Map Insertion/Deletion events to intron positions.                                                                               #
#                                                                                                                                  #
# Expects Indel file (BED formatted) produced by get_indels_from_emf.pl                                                            #
# Maps Indels in this file to Intron positions in another BED file.                                                                #
#                                                                                                                                  #
####################################################################################################################################


my$indel_file;
my$intron_file;
my$outfile;
my$help = 0;

GetOptions("infile=s" => \$intron_file, "bed=s" => \$intron_file, "outfile=s" => \$outfile, "help|?" => \$help);

pod2usage(-verbose => 99) if $help;

# Check required command line arguments are provided.

my$args = 1;
if(not defined $indel_file){
	print "\nMissing argument -i";
	$args = 0;
}

if(not defined $intron_file){
	print "\nMissing argument -b";
	$args = 0;
}

if(not defined $outfile){
	print "\nMissing argument -o";
	$args = 0;
}

# If a required argument is missing, $args will be equal to 0. In this case, print usage and exit.
if($args == 0){
	print "\n\n"; # Some formatting.
	pod2usage(-verbose => 1);
}

####################################################################################################################################

open(IND, $indel_file) || die "Cannot open file $indel_file: $!.\n";

my%indels;
while(my$line = <IND>){
	chomp $line;
	my($chr, $start, $end, $length, $branch, $type) = split("\t", $line);
	push(@{$indels{$chr}{$start}}, [($end, $length."\t".$branch."\t".$type)]);
}

close IN;

####################################################################################################################################

open(INT, $intron_file) || die "Cannot open file $intron_file: $!.\n";

my%introns;
while(my$line = <INT>){
	chomp $line;
	my@cols = split("\t", $line);

	my$chr = $cols[0];
	my$start = $cols[1];
	my$end = $cols[2];
	my$gene = $cols[3];

	$introns{$chr}{$start} = [($end, $gene)];
}

close IN;

####################################################################################################################################

open(OUT, ">", $outfile);

foreach my$chr (keys %indels){
	next if not exists $introns{$chr};

	my@indel_start_sort = sort {$a <=> $b} keys %{$indels{$chr}};
	my@intron_start_sort = sort {$a <=> $b} keys %{$introns{$chr}};

	while(scalar @indel_start_sort != 0 && scalar @intron_start_sort != 0){
		my$indel_start = $indel_start_sort[0];
		my$intron_start = $intron_start_sort[0];

		if($indel_start > $intron_start){
			my$intron_end = $introns{$chr}{$intron_start}[0];

			if($indel_start < $intron_end){
				foreach my$event (@{$indels{$chr}{$indel_start}}){
					my$indel_end = $event -> [0];

					if($indel_end < $intron_end){
						my$gene_id = $introns{$chr}{$intron_start}[1];
						print OUT $chr,"\t",$indel_start,"\t",$indel_end,"\t",$gene_id,"\t",$event -> [1],"\n";
					}
				}

				shift @indel_start_sort;
			}
			else{
				shift @intron_start_sort;
			}
		}
		else{
			shift @indel_start_sort;
		}
	}
}

close OUT;

####################################################################################################################################

exit;

####################################################################################################################################

__END__

=head1 DESCRIPTION

Map Insertion/Deletion events to intron positions. 

=head1 AUTHOR

=over 8

=item
Peter Keane

=item
School of Mathematics, Statistics and Applied Mathematics

=item 
National University of Ireland, Galway

=back

=head1 CONTACT

peterakeane@gmail.com

=head1 OPTIONS

=item B<-i>
Insertion/Deletion events in BED format.

=item B<-b>
Intron positions in BED format.

=item B<-o>
Output file (BED file)

=item B<-h>
Print usage and exit.

=back

