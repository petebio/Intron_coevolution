#!/usr/bin/env perl
use Statistics::Descriptive;
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

#####################################################################################################################################################
#                                                                                                                                                   #
# Randomisation based test to determine if the mean correlation (across edges) of a gene network is >= than can be expected by chance,              #
# controlling for feature (eg. intron) length.                                                                                                      #
#                                                                                                                                                   #
# This test works as follows:                                                                                                                       #
# 1). Calculate the mean correlation value for the real gene network.                                                                               #
# 2). For each gene in the network, choose another at random, but matched for feature length class.                                                 #
# 3). Calculate the mean correlation value for the network of re-sampled genes.                                                                     #
# 4). Repeat this N times.                                                                                                                          #
#                                                                                                                                                   #
# This produces a p-value which represents the proportion of times the mean correlation of the randomised network is >= that of                     #
# the real network.                                                                                                                                 #
#                                                                                                                                                   #
#####################################################################################################################################################

# Process command line arguments.

my$network_file;
my$cor_file;
my$bed_file;
my$n = 1000;
my$help = 0;

GetOptions("input=s" => \$network_file,
	"cor=s" => \$cor_file,
	"bed=s" => \$bed_file,
	"n=i" => \$n,
	"help|?" => \$help);

pod2usage(-verbose => 99) if $help;

# Check that all required arguments are provided.

my$args = 1;
if(not defined $network_file){
	print "\nMissing argument -i";
	$args = 0;
}

if(not defined $cor_file){
	print "\nMissing argument -c";
	$args = 0;
}

if(not defined $bed_file){
	print "\nMissing argument -b";
	$args = 0;
}

# If any required argument missing, $args is 0. Print usage and exit.
if($args == 0){
	print "\n\n"; # Some formatting...
	pod2usage(-verbose => 1);
}

#####################################################################################################################################################

# Read correlation values from $cor_file. 
# This should be a tab-delimited, square symmetrical matrix, or alternatively the lower triangle of the matrix.
# The first column should contain gene IDs.
#
# Example:
#
# gene_A	0
# gene_B	0.1	0
# gene_C	-0.4	-0.03	0
# gene_D	-0.2	-0.1	0.5	0

open(COR, $cor_file) || die "Cannot open file $cor_file: $!.\n";
print "Reading correlation matrix\n";

# %cor is a hash containing all pair-wise gene correlation values.
# @gene_order is used to keep track of rows/genes

my(%cor, @gene_order);
while(my$row = <COR>){
	chomp $row;
	my@cells = split("\t", $row);
	my$gene_id = shift @cells;
	
	for(my$i = 0; $i < scalar @gene_order; $i++){
		my$gene_id_b = $gene_order[$i];
		my$pair = join(":", sort {$a cmp $b} ($gene_id, $gene_id_b));
		
		next if $cells[$i] eq "NA"; # Avoid NA values.
		$cor{$pair} = $cells[$i];
	}
	
	push(@gene_order, $gene_id);
}

close COR;

# %cor_gene_lookup is used later to check if a gene is present in the correlation matrix.
my%cor_gene_lookup = map{$gene_order[$_] => 1}0..$#gene_order;

#####################################################################################################################################################

# Read the gene network from $network_file.
# This file should be a list of edges, a tab delimited list where each row contains a gene-pair connected by an edge.
#
# Example:
#
# gene_A	gene_B
# gene_A	gene_C
# gene_C	gene_D

open(NET, $network_file) || die "Cannot open file $network_file: $!";
print "Reading edge list\n";

# @edges is an array containing all gene-pairs contained in the network.
# %network_gene_lookup is used to keep track of which genes are present in the network.

my(@edges, %network_gene_lookup);
while(my$line = <NET>){
	chomp $line;
	my($node_a, $node_b) = split("\t", $line);
	
	# Check if both genes are contained in the correlation matrix. Keep only these edges.
	if(exists $cor_gene_lookup{$node_a} && exists $cor_gene_lookup{$node_b}){
		my$edge = join(":", sort {$a cmp $b} ($node_a, $node_b));
		push(@edges, $edge);
		
		$network_gene_lookup{$node_a} = 1;
		$network_gene_lookup{$node_b} = 1;
	}
}

close NET;

#####################################################################################################################################################

# Read genomic position of gene features to be used in this analysis (E.g. Intron positions).
# This is used to calculate their length, which is used as a control in the randomisation test.
#
# Requires BED fields are:
#
# Chromosome	Start_posiiton	End_position	Gene_ID
#
# Subsequent fields are ignored.

open(BED, $bed_file) || die "Cannot open file $bed_file: $!.\n";
print "Reading BED file\n";

# %feature_length is a hash containing all genes and their corresponding feature length.

my%feature_length;
while(my$entry = <BED>){
	chomp $entry;
	my@fields = split("\t", $entry);
	
	my$start = $fields[1];
	my$end = $fields[2];
	my$id = $fields[3];
	
	# Keep only genes that are found in the correlation matrix.
	if(exists $cor_gene_lookup{$id}){
		my$lenght = abs($end - $start);
		$feature_length{$id} += $lenght;
	}
}

close BED;

#####################################################################################################################################################

# Classify each gene into one of five classes based on their feature length.

# %gene_class is used to track which class a gene has been assigned to.
# %feature_class is used to track which genes have been assigned to each feature.

my(%gene_class, %feature_class);
foreach my$id (keys %feature_length){
	my$length = $feature_length{$id};
	
	my$class;
	if($length < 5000){
		$class = 1;
	}
	elsif($length >= 5000 && $length < 20000){
		$class = 2;
	}
	elsif($length >= 20000 && $length < 50000){
		$class = 3;
	}
	elsif($length >= 50000 && $length < 100000){
		$class = 4;
	}
	else{
		$class = 5;
	}
	
	$gene_class{$id} = $class;
	push(@{$feature_class{$class}}, $id)
}

#####################################################################################################################################################

# Calculate the mean correlation value across edges for the gene network.

my@network_cor_values;
foreach my$edge (@edges){
	if(exists $cor{$edge}){
		push(@network_cor_values, $cor{$edge});
	}
}

my$network_mean = array_mean(\@network_cor_values);

#####################################################################################################################################################

# Begin randomisation test.

print "Begin randomization test\n\n";

my$x = 0; # Keep track of the number of tests with a mean correlation >= that of the real network.

for(my$i = 1; $i <= $n; $i++){
	$|++;
	print "Running iteration $i of $n\r";

	# For each gene in the network, choose another from the same feature size class.

	my(%sampled_tracking, %sampled_gene);
	foreach my$gene (keys %network_gene_lookup){
		my$class = $gene_class{$gene};
		my$index = int(rand(scalar @{$feature_class{$class}}));
		my$sampled = $feature_class{$class}[$index];
		
		# Make sure each gene is only sampled once.
		if(not defined $sampled_gene{$sampled}){
			$sampled_tracking{$gene} = $sampled;
			$sampled_gene{$sampled} = 1;
		}
		else{
			redo;
		}
	}
	
	# Calculate the mean correlation value for the re-sampled network.
	
	my@sampled_cor_values;
	foreach my$edge (@edges){
		my($gene_a, $gene_b) = split(":", $edge);
		
		my$sampled_a = $sampled_tracking{$gene_a};
		my$sampled_b = $sampled_tracking{$gene_b};
		my$sampled_edge = join(":", sort {$a cmp $b} ($sampled_a, $sampled_b));
		
		if(exists $cor{$sampled_edge}){
			push(@sampled_cor_values, $cor{$sampled_edge});
		}
	}
	
	my$sampled_mean = array_mean(\@sampled_cor_values);
	
	# Compare to real network mean.
	if($sampled_mean >= $network_mean){
		$x++;
	}
}

print "\n\n";

#####################################################################################################################################################

# Calculate p-value and print results.

my$p = $x / $n;

print "Network mean = ",$network_mean,"\n";
print "P-value = ",$p,"\n";

#####################################################################################################################################################

exit;

#####################################################################################################################################################

#### SUBROUTINES ####

# Read an array of numbers and calculate their mean.
sub array_mean{
	my$vector = shift;
	my$stat = Statistics::Descriptive::Full -> new;
	$stat -> add_data($vector);
	my$mean = $stat -> mean;
	
	return $mean;
}

#####################################################################################################################################################

__END__

=head1 DESCRIPTION

Randomisation based test to determine if the mean correlation (across edges) of a gene network is >= than can be expected by chance.

=head1 CONTACT

=over 8

=item
Peter Keane

=item 
School of Mathematics, Statistics and Applied Mathematics

=item
National University of Ireland, Galway

=item
email: peterakeane@gmail.com

=back

=head1 REQUIRES Perl modules

=over 8

=item
Statistics::Descriptive

=back

=head1 OPTIONS

=over 8

=item B<-i>
Edge list for network to be analysed.

=item B<-c>
Gene-Gene correlation matrix.

=item B<-b>
BED file of feature positions.

=item B<-n>
Number of iterations for the randomisation test (default = 1000).

=item B<-h>
Print this helpful help message.

=back
