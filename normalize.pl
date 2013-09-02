#!/bin/env perl
$|=1;

use Getopt::Long;
use lib ".";
use Normalizer;


$verbose = 0;
$base = 2;
$infile = "sample.data.gz";
$format = "simple";
$code = "";
$method = "net";

%opts = (
	'verbose' => \$verbose,
	'base=i' => \$base,
	'infile=s' => \$infile,
	'format=s' => \$format,
	'help' => \$help,
	'code=s' => \$code,
	'outfile=s' => \$outfile,
	'method=s' => \$method,
	'guide=s' => \$guideGenes,
);
%defaults = (
	'verbose' => $verbose,
	'base' => $base,
	'infile' => $infile,
	'format' => $format,
	'help' => $help,
	'code' => $code,
	'outfile' => $outfile,
	'method' => $method,
);
GetOptions(%opts);
if ($help || !$base || !$infile) {
	usage();
	exit;
}
$logbase = log($base);
unless ($outfile) {
	($outfile) = $infile =~ /([^\.]+)/;
	$outfile .= ".normalized";
}





$norm = new Normalizer;
$norm->set_param("tmp", ".");
$norm->set_param("verbose", 1) if $verbose;
$norm->set_param("persistent_id", $code);
$norm->set_param("method", $method);
$norm->set_param("fraction", 1);

print "Reading $infile ...\n" if $verbose;
if (-e $infile) {
	if ($infile =~ /\.gz$/) {
		open DATA, "gunzip -c $infile |";
	} else {
		open DATA, $infile;
	}
	$_ = <DATA>;
	chomp;
	@samples = split /\t/; #first line names the samples
	if ($format eq "simple") {
		shift @samples; #...except for the first column
	} else {
		
		#...except for the first two columns
		foreach my $i (0..1) {
			shift @samples;
		}
	}
	
	print scalar @samples, " samples\n" if $verbose;
	while (<DATA>) {
		my(@values, $transcript, $gene, $class, $length);
		chomp;
		if ($format eq "simple") {
			($transcript, @values) = split /\t/;
		} else {
			($gene, $length, @values) = split /\t/;
			($gene, $transcript) = split /:/, $gene;
			$transcript = "$gene:$transcript";
			$gene{$transcript} = $gene;
			$class{$gene} = $class;
			$length{$transcript} = $length;
		}
		
		foreach my $i (0..$#samples) {
			if ((my $v = $values[$i]) > 0) {
				$counts{$transcript}{$samples[$i]} = $v;
				$nonzero++;
			}
		}
	}
	close DATA;
	print "nonzero: $nonzero\n" if $verbose;
} else {
	die "I don't know $infile\n";
}



$norm->set_data(\%counts, \%length);
$sampleCount = $norm->param("sample_count");
print "samples: $sampleCount, genes: ", scalar keys %counts, ", nonzero values: $nonzero\n" if $verbose;


if ($method ne 'quantile' && $method ne "none") {
	if ($method eq 'guide' || $method eq 'genorm') {
		open GUIDES, $guideGenes;
		while (<GUIDES>) {
			chomp;
			my($gene, $weight) = split /\t/;
			$weight{$gene} = ($weight || 1);
		}
		close GUIDES;
		die "Method '$method' requires a file of guide gene weights, specified via the 'guide' parameter.\n" unless %weight;
		$norm->set_gene_weights(%weight);
	} elsif ($method eq "total") {
		$norm->set_param("total_size_for_normalizing", 0);
	} elsif ($method eq "cpm") {
		$norm->set_param("total_size_for_normalizing", 1e6);
	} else {
		$ubiquitous = $norm->datarestore("ubiquitous_genes") || $norm->identify_ubiquitous_genes();
		foreach my $gene (@{$ubiquitous}) {
			$isUbiquitous{$gene} = 1;
		}
	}
	
	$fact = $norm->compute_scaling_factors();
	
	print "\nComputed scaling factors per sample:\n" if $verbose;
	foreach my $sample (@samples) {
		print join("\t", $sample, sprintf("%.5f", $fact->{$sample})), "\n";
	}
}
$normalized = $norm->normalize();
	

if ($outfile) {
	print "Saving $outfile ...\n" if $verbose;
	open OUTF, ">$outfile";
	print OUTF join("\t", "genes\\samples", @samples), "\n";
	while (my($gene) = each %counts) {
		$_ = $gene;
		s/\t/;/g;
		print OUTF $_;
		foreach my $sample (@samples) {
			print OUTF "\t";
			if (my $v = $normalized->{$gene}{$sample}) {
				print OUTF sprintf("%.3f", $v);
			}
		}
		print OUTF "\n";
	}
	close OUTF;
}


#######
sub usage {
	print "Usage: normalize.pl [options]\n";
	print "\tOption\tDefault\n";
	print "\t------\t-------\n";
	foreach my $opt (keys %defaults) {
		print join("\t", "", $opt, $defaults{$opt}), "\n";
	}
}

