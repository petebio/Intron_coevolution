#!/usr/bin/env perl
use Statistics::Descriptive;
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

#############################################################################################################################
#                                                                                                                           #
# Calculate gene-gene partial Pearson correlations for intron length change.                                                #
# Controls for median change along each branch of the phylogenetic tree.                                                    #
#                                                                                                                           #
# Expects matrix of proportional intron length changes from calculate_proportional_change.pl                                #
#                                                                                                                           #
#############################################################################################################################

my$matrix_file;
my$outfile;
my$help = 0;

GetOptions("infile=s" => \$matrix_file, "outfile=s" => \$outfile, "help|?" => \$help);

pod2usage(-verbose => 99) if $help;

# Check required command line arguments are provided.

my$args = 1;
if(not defined $matrix_file){
	print "\nMissing argument -i";
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

#############################################################################################################################

# Read table of intron change values.
# Columns correspond to branches of the phylogenetic tree. Rows correspond to genes.

open(MAT, $matrix_file)||die "Cannot open file $matrix_file: $!.\n";

my(%gene_values, @header);
while(my$row = <MAT>){
	chomp $row;
	
	if($. == 1){
		@header = split("\t", $row);
		shift @header;
	}
	else{
		my@cols = split("\t", $row);
		my$gene_id = shift @cols;
		
		for(my$i = 0; $i < scalar @cols; $i++){
			my$col_id = $header[$i];
			$gene_values{$gene_id}{$col_id} = $cols[$i];
		}
	}
}

close MAT;

#############################################################################################################################

# Get the median value for each branch. This is the control vector for the partial correlation.

my%median_values;
foreach my$id (@header){
	my@values;
	foreach my$gene (keys %gene_values){
		push(@values, $gene_values{$gene}{$id});
	}
	
	my$stat = Statistics::Descriptive::Full -> new;
	$stat -> add_data(@values);
	my$median = $stat -> median;
	
	if($median != 0){
		$median_values{$id} = $median;
	}
}

#############################################################################################################################

open(OUT, ">", $outfile);

# Calculate partial correlations and write to file.
# Output is the lower triangle of a symmetrical matrix.
my@gene_ids = keys %gene_values;
for(my$i = 0; $i < scalar @gene_ids; $i++){
	my$id_x = $gene_ids[$i];
	print OUT $id_x;
	
	for(my$j = 0; $j < $i; $j++){
		my$id_y = $gene_ids[$j];
		
		my(@vector_x, @vector_y, @vector_z);
		foreach my$id (keys %median_values){
			push(@vector_x, $gene_values{$id_x}{$id});
			push(@vector_y, $gene_values{$id_y}{$id});
			push(@vector_z, $median_values{$id});
		}
		
		my$xy = cor(\@vector_x, \@vector_y);
		my$xz = cor(\@vector_x, \@vector_z);
		my$yz = cor(\@vector_y, \@vector_z);
		
		my$pcor = ($xy - $xz * $yz) / (sqrt(1 - $xz ** 2) * sqrt(1 - $yz ** 2));
		print OUT "\t",$pcor;
	}
	
	print OUT "\t0\n"; # Self correlations are set to zero.
}

#############################################################################################################################

exit;

#############################################################################################################################

##### SUBROUTINES #####

sub array_mean{
	my$vector = shift;
	my$stat = Statistics::Descriptive::Full -> new;
	$stat -> add_data($vector);
	return $stat -> mean;
}

# Calculate Pearson correlation for two numeric vectors.
sub cor{
	my($vector_x, $vector_y) = @_;
	
	my$x_bar = array_mean($vector_x);
	my$y_bar = array_mean($vector_y);
	
	my$x = 0;
	my$y = 0;
	my$xy = 0;
	
	for(my$i = 0; $i < scalar @{$vector_x}; $i++){
		my$x_calc = $vector_x -> [$i] - $x_bar;
		my$y_calc = $vector_y -> [$i] - $y_bar;
		
		$x += ($x_calc ** 2);
		$y += ($y_calc ** 2);
		$xy += ($x_calc * $y_calc);
	}
	
	if($x == 0 || $y == 0){
		return 0;
	}
	else{
		my$cor = $xy / (sqrt($x) * sqrt($y));
		return $cor;
	}
}

#############################################################################################################################

__END__

=head1 DESCRIPTION

Calculate gene-gene partial Pearson correlations for intron length change.
Controls for median change along each branch of the phylogenetic tree.

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

=head1 REQUIRED Perl modules

=over 8

=item
Statistics::Descriptive

=back

=head1 OPTIONS

=over 8

=item B<-i>
Matrix of intron length change values. Expects output from calculate_proportional_change.pl

=item B<-o>
Output file. The lower triangle of a symmetrical matrix.

=item B<-h>
Print usage and exit.