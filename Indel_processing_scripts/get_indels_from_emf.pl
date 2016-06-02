#!/usr/bin/env perl
use Array::Utils qw(intersect);
use Bio::TreeIO;
use IO::String;
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

##################################################################################################################################################
#                                                                                                                                                #
# Read an Ensembl Multi Format (EMF) file containing whole genome alignments, and find all insertion/deletion events that can                    #
# be mapped to a given phylogenetic tree.                                                                                                        #
#                                                                                                                                                #
# Also provide a reference species for which genomic positions should be used.                                                                   #
#                                                                                                                                                #
# This program was designed using the 15 Eutherian mammals EPO dataset from Ensembl Compara release 74.                                          #
# Will probably work with other similar datasets, but is untested.                                                                               #
#                                                                                                                                                #
##################################################################################################################################################

my$emf_file;
my$tree_file;
my$ref_taxon = "Hsap";
my$outfile;
my$help = 0;

GetOptions("infile=s" => \$emf_file, "tree=s" => \$tree_file, "species=s" => \$ref_taxon, "outfile=s" => \$outfile, "help|?" => \$help);

pod2usage(-verbose => 99) if $help;

# Check required command line arguments are provided.

my$args = 1;
if(not defined $emf_file){
	print "\nMissing argument -i";
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

# Read reference phylogenetic tree from $tree_file.
# Tree should be in newick format.

open(TR, $tree_file) || die "Cannot open file $tree_file: $!.\n";
chomp(my$ref_tree_raw = <TR>);
close TR;

my$ref_tree = read_newick($ref_tree_raw);

# Get a list of all taxa on the reference tree.
my@ref_taxa = get_taxa($ref_tree);

# Get parent/child relationships from phylogenetic tree.
my%parent_nodes = get_parent_nodes($ref_tree);

##################################################################################################################################################

# Process EMF file.

open(EMF, $emf_file) || die "Cannot open file $emf_file: $!.\n";

# Define some variables.
my(@taxa_order, $ref_chr, $ref_pos, %block_node_ids, %valid_nodes, @alignment);
my$valid = 0;

# Open output file.
open(OUT, ">", $outfile);

while(my$line = <EMF>){
	next if $line =~ m/^(:?#|DATA)/; # Skip these lines.
	chomp $line;

	if($line =~ m/^SEQ/){
		my($species_full, $chr, $pos) = $line =~ m/^SEQ\s+(\S+)\s+(\S+)\s+(\S+)\s+/;

		if($species_full !~ m/^anc/){
			my($first, $last) = $species_full =~ m/^(.)[^_]+_(.{3})/;
			my$species_id = uc($first).$last; # Reformat species name to match that of phylogenetic tree.

			push(@taxa_order, $species_id);

			if($species_id eq $ref_taxon){
				$ref_chr = $chr;
				$ref_pos = $pos;
			}
		}
		else{
			my$anc_id = "Aseq_".$chr;
			push(@taxa_order, $anc_id);
		}		
	}
	elsif($line =~ m/^TREE/){
		my($block_tree_raw) = $line =~ m/^TREE\s+(\S+)$/;
		$block_tree_raw =~ s/\[.\]//g; # re-format tree for compatability with Bio::TreeIO

		my$block_tree = read_newick($block_tree_raw);
		$valid = validate_tree($block_tree, \@ref_taxa, $ref_taxon);
		next if $valid == 0; # Skip invalid alignment block.

		# Get the proper node IDs for the phylogenetic tree from current alignment block.
		%block_node_ids = get_node_id_align($block_tree, \@ref_taxa);
	}
	elsif($line =~ m/^(:?A|G|T|C|-)/i){
		next if $valid == 0; # Skip invalid alignment block.
		push(@alignment, $line);
	}
	elsif($line eq "//"){
		# End of alignment block. Infer indels.
		# Skip invalid alignment block.
		if($valid == 1){
			infer_indels(\@alignment, \@taxa_order, \%parent_nodes, \%block_node_ids, $ref_taxon, $ref_chr, $ref_pos);
		}

		# Clear/reset variables for next block.
		undef @alignment;
		undef @taxa_order;
		undef %block_node_ids;
		undef %valid_nodes;
		undef $ref_chr;
		undef $ref_pos;
		$valid = 1;
	}
}

close OUT;

##################################################################################################################################################

exit;

##################################################################################################################################################

##### SUBROUTINES #####

# Read a string containing a phylogenetic tree in newick format.
# Return a tree object from Bio::TreeIO BioPerl module.
sub read_newick{
	my$raw_tree = shift;

	my$io = IO::String -> new($raw_tree);
	my$tree_io = Bio::TreeIO -> new(-format => "newick", -fh => $io);
	my$tree = $tree_io -> next_tree;
	
	return $tree;
}

# Get a list of all leaf nodes on the phylogenetic tree.
sub get_taxa{
	my$tree = shift;

	my@taxa;
	foreach my$taxon ($tree -> get_leaf_nodes){
		my$taxon_id = substr($taxon -> id, 0, 4);
		push(@taxa, $taxon_id);
	}

	return @taxa;
}

# Get the ID for a node in the phylogenetic tree.
# For tip nodes, this is the node label.
# For internal nodes, this is the list of all tip nodes descended from that node.
sub get_node_id{
	my$node = shift;

	my$node_id;
	if($node -> is_Leaf){
		$node_id = substr($node -> id, 0, 4);
	}
	else{
		my@leaf_nodes;
		foreach my$child ($node -> get_all_Descendents){
			if($child -> is_Leaf){
				my$id = substr($child -> id, 0, 4);
				push(@leaf_nodes, $id);
			}
		}

		$node_id = join(",", sort {$a cmp $b} @leaf_nodes);
	}

	return $node_id;
}

# Get child/parent relationships for each node of the phylogenetic tree.
sub get_parent_nodes{
	my$tree = shift;

	my$root_node = $tree -> get_root_node;

	my%parent_nodes;
	foreach my$child ($root_node -> get_all_Descendents){
		my$parent = $child -> ancestor;

		my$child_id = get_node_id($child);
		my$parent_id = get_node_id($parent);

		$parent_nodes{$child_id} = $parent_id;
	}

	return %parent_nodes;
}

# For each tree in the EMF file, we need to know which node in the tree corresponds to the relevant node in the reference tree.
# This function identifies such nodes, by finding the ancestral node that is closest to the last common ancestor of all taxa.
sub get_node_id_align{
	my($tree, $ref_taxa) = @_;
	
	my@taxa = @$ref_taxa;
	my%taxa_lookup = map{$taxa[$_] => 1}0..$#taxa;
	
	my$all_taxa = join(",", sort {$a cmp $b} @taxa);
	
	my$root = $tree -> get_root_node;
	push(my@nodes_to_visit, $root);

	my%new_nodes;
	my($root_id) = $root -> id =~ m/^(Aseq_Ancestor_\d+_\d+)_/;
	$new_nodes{$all_taxa} = $root_id;
	
	
	while(scalar @nodes_to_visit != 0){
		my$visit_node = shift @nodes_to_visit;
		
		foreach my$node ($visit_node -> each_Descendent){
			if($node -> is_Leaf){
				my$id = substr($node -> id, 0, 4);
				
				if(exists $taxa_lookup{$id}){
					$new_nodes{$id} = $id;
				}
			}
			else{
				my@leaf_nodes;
				foreach my$child ($node -> get_all_Descendents){
					if($child -> is_Leaf){
						my$id = substr($child -> id, 0, 4);
						
						if(exists $taxa_lookup{$id}){
							push(@leaf_nodes, $id);
						}
					}
				}
				
				next if scalar @leaf_nodes == 0;
				my$node_id = join(",", sort {$a cmp $b} @leaf_nodes);
				
				my($aseq_id) = $node -> id =~ m/^(Aseq_Ancestor_\d+_\d+)_/;
				$new_nodes{$node_id} = $aseq_id;
				
				push(@nodes_to_visit, $node);
			}
		}
	}
	
	return %new_nodes;	
}

# Perform checks to determine if the alignment block can be used. These tests are:
# 1). Check if the tree contains any duplicated species.
# 2). Check if the reference species is contained in the alignment block.
# 3). Check if all taxa present in the refernce tree provided are also in the alignment block.
# Any alignment block that fails any one of these tests returns 0 and is skipped.
# Those that pass, return 1.
sub validate_tree{
	my($tree, $taxa_list, $ref_taxon) = @_;

	my@ref_taxa = @$taxa_list;
	my@tree_taxa = get_taxa($tree);

	my%taxa_lookup = map{$tree_taxa[$_] => 1}0..$#tree_taxa;

	# Check if any taxa are duplicated.
	if(scalar keys %taxa_lookup != scalar @tree_taxa){
		return 0;
	}

	# Check if reference species is present in alignment block.
	if(not exists $taxa_lookup{$ref_taxon}){
		return 0;
	}

	# Check if all taxa from reference tree are also present in alignment block.
	my@common_taxa = intersect(@ref_taxa, @tree_taxa);
	if(scalar @common_taxa != scalar @ref_taxa){
		return 0;
	}

	# If all tests are passed, return 1.
	return 1;
}

# Find what alignment column corresponds to a given taxon.
sub get_col_num{
	my($order, $id) = @_;
	
	for(my$i = 0; $i < scalar @{$order}; $i++){
		if($order -> [$i] eq $id){
			return $i;
		}
	}
}

# Main function. Identify insertion/deletion events from the alignment.
sub infer_indels{
	my($align_ref, $order_ref, $parent_ref, $block_ref, $ref, $chr, $pos) = @_;
	
	my@alignment = @$align_ref;
	my%parent_nodes = %$parent_ref;
	my%block_node_ids = %$block_ref;
	
	my$ref_col = get_col_num($order_ref, $ref);
	
	foreach my$child (keys %parent_nodes){
		my$parent = $parent_nodes{$child};
		
		my$child_id = $block_node_ids{$child};
		my$parent_id = $block_node_ids{$parent};
		my$child_col = get_col_num($order_ref, $child_id);
		my$parent_col = get_col_num($order_ref, $parent_id);
		
		my$insertion_begin = 0;
		my$insertion_length = 0;
		my$deletion_begin = 0;
		my$deletion_length = 0;
		my$current_pos = $pos;
	
		foreach my$align (@alignment){			
			my$ref_base = substr($align, $ref_col, 1);
			
			if($ref_base ne "-"){
				$current_pos++;
			}
			
			my$child_base = substr($align, $child_col, 1);
			my$parent_base = substr($align, $parent_col, 1);
			
			if($child_base ne "-" && $parent_base eq "-"){
				# Is an insertion.
				if($insertion_length == 0){
					$insertion_begin = $current_pos;
				}
				
				$insertion_length++;
				
				# If previous line was a deletion, then report.
				if($deletion_length > 0){
					print OUT $chr,"\t",$deletion_begin,"\t",$current_pos + 1,"\t",$deletion_length,"\t",$child,"\tD\n";
					$deletion_length = 0;
				}
			}
			elsif($child_base eq "-" && $parent_base ne "-"){
				# Is a deletion.
				if($deletion_length == 0){
					$deletion_begin = $current_pos;
				}
				
				$deletion_length++;
				
				# If previous line was an insertion, then report.
				if($insertion_length > 0){
					print OUT $chr,"\t",$insertion_begin,"\t",$current_pos + 1,"\t",$insertion_length,"\t",$child,"\tI\n";
					$insertion_length = 0;
				}
			}
			elsif($child_base ne "-" && $parent_base ne "-"){
				# Is neither, positions are aligned. Report insertion/deletion.
				if($insertion_length > 0){
					print OUT $chr,"\t",$insertion_begin,"\t",$current_pos + 1,"\t",$insertion_length,"\t",$child,"\tI\n";
					$insertion_length = 0;
				}
				
				if($deletion_length > 0){
					print OUT $chr,"\t",$deletion_begin - 1,"\t",$current_pos + 1,"\t",$deletion_length,"\t",$child,"\tD\n";
					$deletion_length = 0;
				}
			}
		}
	}
}

##################################################################################################################################################

__END__

=head1 DESCRIPTION

Read an Ensembl Multi Format (EMF) file containing whole genome alignments, and find all insertion/deletion events that can
be mapped to a given phylogenetic tree.

Also provide a reference species for which genomic positions should be used.

This program was designed using the 15 Eutherian mammals EPO dataset from Ensembl Compara release 74. 
Will probably work with other similar datasets, but is untested.

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
Array::Utils

=item
IO::String

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

=over 8

=item B<-i>
EMF file.

=item B<-t>
Reference species tree.

=item B<-s>
Reference taxon code (Default = Hsap).

=item B<-o>
Outfile (BED file).

=item B<-h>
Print this helpful help message.

=back
