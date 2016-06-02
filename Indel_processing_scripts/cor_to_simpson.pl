#!/usr/bin/env perl
use Array::Utils (intersect);
use Getopt::Long;
use Pod::Usage;
use warnings;
use strict;

#######################################################################################################################################
#                                                                                                                                     #
# Read a matrix of correlation values, and calculate the Simpson coefficient for each gene pair.                                      #
#                                                                                                                                     #
#######################################################################################################################################

my$infile;
my$threshold;
my$outfile;
my$help = 0;

GetOptions("infile=s" => \$infile, "threshold=s" => \$threshold, "outfile=s" => \$outfile, "help|?" => \$help);

pod2usage(-verbose => 99) if $help;

# Check required command line arguments are provided.

my$args = 1;
if(not defined $infile){
	print "\nMissing argument -i";
	$args = 0;
}

if(not defined $threshold){
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

#######################################################################################################################################

open(IN, $infile) || die "Cannot open file $infile: $!.\n";
print "Reading correlation values from $infile\n";

my(@gene_order, %gene_lookup, %edges);
while(my$line = <IN>){
	chomp $line;
	my@cells = split("\t", $line);
	
	my$row_id = shift @cells;
	$gene_lookup{$row_id} = $.;
	
	for(my$i = 0; $i < scalar @gene_order; $i++){
		my$col_id = $gene_order[$i];
		
		if($cells[$i] > $threshold){
			push(@{$edges{$row_id}}, $gene_lookup{$col_id});
			push(@{$edges{$col_id}}, $gene_lookup{$row_id});
		}
	}
	
	push(@gene_order, $row_id);
}

close IN;

#######################################################################################################################################

print "Calculating Simpson coefficients\n";

my(%simpson, %valid_genes);
for(my$i = 0; $i < scalar @gene_order; $i++){
	my$id_a = $gene_order[$i];
	next if not exists $edges{$id_a};
	$valid_genes{$id_a} = 1;
	
	for(my$j = $i + 1; $j < scalar @gene_order; $j++){
		my$id_b = $gene_order[$j];
		next if not exists $edges{$id_b};
		$valid_genes{$id_b} = 1;
		
		my@sort_sizes = sort {$a <=> $b} (scalar @{$edges{$id_a}}, scalar @{$edges{$id_b}});
		my$min = $sort_sizes[0];
		
		my$shared = scalar intersect(@{$edges{$id_a}}, @{$edges{$id_b}});
		
		my$edge_id = join(":", sort {$a cmp $b} ($id_a, $id_b));
		
		$simpson{$edge_id} = $shared / $min;
	}
}

my@output_genes = keys %valid_genes;

#######################################################################################################################################

open(OUT, ">", $outfile);
print "Writing to file\n";

for(my$i = 0; $i < scalar @output_genes; $i++){
	my$id_a = $gene_order[$i];
	print OUT $id_a;
	for(my$j = 0; $j < scalar @output_genes; $j++){
		my$id_b = $gene_order[$j];
		
		if($id_a eq $id_b){
			print OUT "\t0";
		}
		else{
			my$edge_id = join(":", sort {$a cmp $b} ($id_a, $id_b));
			
			if(exists $simpson{$edge_id}){
				print OUT "\t",$simpson{$edge_id};
			}
			else{
				print OUT "\t0";
			}
		}
	}
	
	print OUT "\n";
}

close OUT;

#######################################################################################################################################

__END__

=head1 DESCRIPTION

Read a matrix of correlation values, and calculate the Simpson coefficient for each gene pair.

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
Array::Utils

=back

=head1 OPTIONS

=over 8

=item B<-i>
Correlation matrix file.

=item B<-t>
Correlation threshold.

=item B<-o>
Output file.

=item B<-h>
Print help message and exit.

=back