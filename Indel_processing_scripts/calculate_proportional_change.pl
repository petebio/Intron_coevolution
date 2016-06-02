#!/usr/bin/env perl
use Bio::TreeIO;
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

##################################################################################################################################################
#                                                                                                                                                #
# Calculate the proportional change in intron content for all genes provided in the indel BED file.                                              #
#                                                                                                                                                #
# Requires                                                                                                                                       #
# 1). a BED file containing all Insertion/Deletion events, as produced by map_indels_to_introns.pl                                               #
# 2). a BED file with intron positions, used to calculate intron length in the reference species.                                                #
# 3). The reference species to use. Default is Human.                                                                                            #
# 4). The reference phylogenetic tree.                                                                                                           #
# 5). Output file.                                                                                                                               #
#                                                                                                                                                #
##################################################################################################################################################

my$indel_file;
my$intron_file;
my$tree_file;
my$ref_species = "Hsap";
my$outfile;
my$help = 0;

GetOptions("infile=s" => \$indel_file,
	"bed=s" => \$intron_file,
	"tree=s" => \$tree_file,
	"species=s" => \$ref_species,
	"outfile=s" => \$outfile,
	"help|?" => \$help);

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

if(not defined $tree_file){
	print "\nMissing argument -t";
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

##################################################################################################################################################

# Read Insertion/Deletion events from $indel_file.

open(INDEL, $indel_file) || die "Cannot open file $indel_file: $!.\n";

my%indels;
while(my$line = <INDEL>){
	chomp $line;
	my@fields = split("\t", $line);
	
	my$gene_id = $fields[3];
	my$length = $fields[4];
	my$branch = $fields[5];
	my$type = $fields[6]; # $type = I for insertions, and = D for deletions.
	
	if($type eq "D"){
		$length = -$length;
	}
	
	push(@{$indels{$gene_id}{$branch}}, $length);
}

close INDEL;

##################################################################################################################################################

# Get intron length in reference species from $intron_file.

open(BED, $intron_file) || die "Cannot open file $intron_file: $!.\n";

my%intron_length;
while(my$line = <BED>){
	chomp $line;
	my@fields = split("\t", $line);
	
	my$start = $fields[1];
	my$end = $fields[2];
	my$length = $end - $start;
	
	my$gene_id = $fields[3];
	$intron_length{$gene_id} += $length;
}

close BED;

##################################################################################################################################################

# Read reference phylogenetic tree from $tree_file
my$tree_io = Bio::TreeIO -> new(-format => "newick", -file => $tree_file);
my$tree = $tree_io -> next_tree;

# Get parent/child relationships for all nodes.

my$root = $tree -> get_root_node;
my$root_id = get_node_id($root);

my%child_nodes;
foreach my$node ($root -> get_all_Descendents){
	my$parent = $node -> ancestor;

	my$node_id = get_node_id($node);
	my$parent_id = get_node_id($parent);

	push(@{$child_nodes{$parent_id}}, $node_id);
}

##################################################################################################################################################

# Calculate the proportional change in intron length for each branch of the phylogenetic tree.

# Define some variables.
my(%proportional_change, %node_tracking);

foreach my$gene_id (keys %indels){
	my$ref_species_length = $intron_length{$gene_id};

	# Get the intron length of the root taxon. 
	# Do this by subtracting all Insertions/Deletions from that of the reference species length.
	
	my$root_length = $ref_species_length;
	foreach my$node (keys %{$indels{$gene_id}}){
		if($node =~ m/$ref_species/){
			foreach my$indel (@{$indels{$gene_id}{$node}}){
				$root_length -= $indel;
			}
		}
	}
	
	# Now, infer the intron length for each node by adding all insertion/deletion events.
	# Start from the root, and work towards leaf nodes.
	
	my@nodes_to_process;
	push(@nodes_to_process, $root_id);
	
	my%inferred_length;
	$inferred_length{$root_id} = $root_length;
	
	while(scalar @nodes_to_process > 0){
		my$node = shift @nodes_to_process;
		my$node_length = $inferred_length{$node};
		
		foreach my$child (@{$child_nodes{$node}}){
			# If node is not a leaf node, add to list to process.
			if(exists $child_nodes{$child}){
				push(@nodes_to_process, $child);
			}
		
			# Initial length for intron at this node.
			$inferred_length{$child} = $node_length;
			
			# If no indels at this node, then the intron length is equal to that of its parent.
			next if not exists $indels{$gene_id}{$child};
			
			# Apply indels.
			foreach my$indel (@{$indels{$gene_id}{$child}}){
				$inferred_length{$child} += $indel;
			}
		}
	}
	
	# Now, calculate proportional change along each branch. 
	
	foreach my$node (keys %child_nodes){
		my$node_length = $inferred_length{$node};
		
		foreach my$child (@{$child_nodes{$node}}){
			my$child_length = $inferred_length{$child};
			if(not defined $child_length){print $gene_id,"\t",$child,"\n"; exit;};
			my$p = 2 * ($child_length - $node_length) / ($child_length + $node_length);
			$proportional_change{$gene_id}{$child} = $p;
			
			$node_tracking{$child} = 1;
		}
	}
}

##################################################################################################################################################

# Write output file.
my@column_order = sort {$a cmp $b} keys %node_tracking;
open(OUT, ">", $outfile);

print OUT "gene_id\t",join("\t", @column_order),"\n";

foreach my$gene_id (keys %proportional_change){
	my@output_vector;
	foreach my$node (@column_order){
		push(@output_vector, $proportional_change{$gene_id}{$node});
	}
	
	print OUT $gene_id,"\t",join("\t", @output_vector),"\n";
}

close OUT;

##################################################################################################################################################

exit;

##################################################################################################################################################

##### SUBROUTINES #####

sub get_node_id{
	my$node_obj = shift;
	
	my$node_id;
	if($node_obj -> is_Leaf){
		$node_id = $node_obj -> id;
	}
	else{
		my@leaf_nodes;
		foreach my$child ($node_obj -> get_all_Descendents){
			if($child -> is_Leaf){
				my$child_id = $child -> id;
				push(@leaf_nodes, $child_id);
			}
		}
		
		$node_id = join(",", sort {$a cmp $b} @leaf_nodes);
	}
	
	return $node_id;
}

##################################################################################################################################################

__END__

=head 1 DESCRIPTION

Calculate the proportional change in intron content for all genes provided in the indel BED file.

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

=head1 REQUIRED Perl/BioPerl modules

=over 8

=item
Bio::TreeIO

=back

=head1 VALID taxon codes - Ensembl Compara 74 (15 Eutherian mammals)

=over 8

=item
Hsap - Homo sapiens

=item
Ptro - Pan troglodytes

=item
Pabe - Pongo abelii

=item
Ggor - Gorilla gorilla

=item
Mmul - Macaca mulatta

=item
Cjac - Callithrix jacchus

=item
Fcat - Felis catus

=item
Cfam - Canis familiaris

=item
Ecab - Equus caballus

=item
Oari - Ovis aries

=item
Btau - Bos taurus

=item
Mmus - Mus musculus

=item
Rnor - Ratus norvegicus

=item
Sscr - Sus scrofa

=item
Ocun - Oryctolagus cuniculus

=back

=head1 OPTIONS

=item B<-i>
BED file of Indel positions.

=item B<-b>
BED file of intron positions.

=item B<-t>
Reference phylogenetic tree.

=item B<-s>
Reference taxon code (Default = Hsap)

=item B<-o>
Output file.

=item B<-h>
Print help message and exit.

=back