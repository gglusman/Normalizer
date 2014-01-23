package Normalizer;
use strict;
use BaseClass;
our @ISA = ("BaseClass");
our $logtwo = log(2);

our $VERSION = '1.00';

our $lib = "/users/gglusman/scripts/Normalizer";
use lib "/users/gglusman/scripts/Normalizer";
use lib "/net/gestalt/cpan/Graph-0.91/lib";
use lib "/net/gestalt/cpan/Graph-Centrality-Pagerank-1.03/lib";

use Graph::Undirected;
use Graph::Centrality::Pagerank;
use Storable;

sub initialize {
	my($self) = @_;
	$self->SUPER::initialize(@_);
	
	$self->set_param('defaultScalingFactor', 1);
	$self->set_param('fraction', 1);
	$self->set_param('logmaxmincutoff', 2);
	$self->set_param('base', 2);
	$self->set_param("minimal_expression", 10);
	$self->set_param("limit_genes_to_compare_net", 2000);
	$self->set_param("limit_genes_to_compare_stability", 1000);
	$self->set_param("total_size_for_normalizing", 1e6);
	$self->set_param("method", "net");
	$self->set_param("CoV-cutoff", 0.25);
	$self->set_param("lower_trim_inclusive", 5);
	$self->set_param("upper_trim_inclusive", 95);
	$self->set_param("lower_trim", 30);
	$self->set_param("upper_trim", 85);
	my %info;
	open ML, "$lib/method.list";
	$_ = <ML>;
	chomp;
	my @fields = split /\t/;
	while (<ML>) {
		chomp;
		my(@values) = split /\t/;
		foreach my $i (0..$#fields) {
			$info{$fields[$i]} = $values[$i];
		}
		next;
		
		my $url = $info{'Method Link'};
		my @refs;
		my @pmid = split /;/, $info{'PMID'};
		my @fa = split /;/, $info{'First Author'};
		my @year = split /;/, $info{'Year'};
		my @cit = split /;/, $info{'Citation'};
		foreach my $i (0..$#pmid) {
			push @refs, "<a href=\"http://www.ncbi.nlm.nih.gov/pubmed?term=$pmid[$i]\">$fa[$i] ($year[$i]) $cit[$i]</a>";
		}
		
	}
	close ML;
	$self->set_param("method_info", \%info);
}

sub datastore {
	my($self, $ref, $part) = @_;
	
	my $name = $self->param("persistent_id");
	$name =~ s/\///g;
	if ($name) {
		store $ref, $self->param("outdir") . "/$name.$part";
	}
}

sub datarestore {
	my($self, $part) = @_;
	
	my $d = $self->{'index'}->{'data'};
	return $self->{$d}->{$part} if $self->{$d}->{$part};
	my $name = $self->param("persistent_id");
	$name =~ s/\///g;
	return unless $name;
	my $file = $self->param("outdir") . "/$name.$part";
	if (-s $file) {
		print "...(restoring $file)\n" if $self->param("verbose");
		return $self->{$d}->{$part} = retrieve $file;
	}
}

sub datareset {
	my($self, $part) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $name = $self->param("persistent_id");
	$name =~ s/\///g;
	if ($name) {
		unlink $self->param("outdir") . "/$name.$part";
	}
	delete $self->{$d}->{$part};
}

sub reset {
	my($self, @parts) = @_;
	
	foreach my $part (@parts) {
		$self->datareset($part);
	}
}


sub set_data {
	#expects a reference to a hash-of-hashes: counts{gene}{sample}
	my($self, $dataref, $lengthsref) = @_;
	
	my $logbase = log($self->param("base"));
	my %logv;
	my(%sampletotal, %genetotal, %sampleCount);
	print "...set data\n" if $self->param("verbose");
	my $d = $self->{'index'}->{'data'};
	$self->{$d}->{'values'} = $dataref;
	$self->{$d}->{'lengths'} = $lengthsref;
	foreach my $gene (keys %{$dataref}) {
		foreach my $sample (keys %{$dataref->{$gene}}) {
			my $v = $dataref->{$gene}{$sample};
			$logv{$gene}{$sample} = log($v)/$logbase;
			$v *= $lengthsref->{$gene}/1000 if defined $lengthsref->{$gene};
			$sampletotal{$sample} += $v;
			$genetotal{$gene} += $v;
		}
		$sampleCount{$gene} = scalar keys %{$dataref->{$gene}};
	}
	$self->{$d}->{'logvalues'} = \%logv;
	$self->set_param("samples", keys %sampletotal);
	$self->set_param("sample_count", scalar keys %sampletotal);
	$self->{$d}->{'sampletotal'} = \%sampletotal;
	print join("\t", %sampletotal), "\n" if $self->param("verbose");
	$self->{$d}->{'genetotal'} = \%genetotal;
	$self->{$d}->{'samplecount'} = \%sampleCount;
}

sub set_method {
	my($self, $method) = @_;
	
	$method = lc $method;
	my $alias = $self->{'alias'}{$method};
	if ($alias || $self->{'knownMethod'}{$method}) {
		$method = $alias || $method;
		$self->set_param("method", $method);
		return $method;
	} else {
		return;
	}
}


sub identify_values_at_percentiles {
	my($self) = @_;
	
	my(%rv, %zv, @nzq, @zq);
	my $d = $self->{'index'}->{'data'};
	my $v = $self->{$d}->{'values'};
	my @samples = $self->samples_to_use();
	while (my($gene) = each %{$v}) {
		foreach my $sample (@samples) {
			if ($v->{$gene}{$sample}) {
				$rv{$sample}{$gene} = $v->{$gene}{$sample};
			} else {
				$zv{$sample}++;
			}
		}
	}
	
	foreach my $sample (@samples) {
		my @v = values %{$rv{$sample}};
		my @sorted = sort {$a<=>$b} @v;
		my $n = scalar @v;
		$nzq[0]{$sample} = $sorted[0];
		foreach my $target (1..99) {
			my $uqn = $n*$target/100;
			my $f = $uqn-int($uqn);
			$nzq[$target]{$sample} = $sorted[$uqn]*(1-$f)+$f*$sorted[$uqn+1];
		}
		$nzq[100]{$sample} = $sorted[$#sorted];
		
		if ($zv{$sample}) {
			my $zn = $zv{$sample};
			$n += $zn;
			$zq[0]{$sample} = 0;
			foreach my $target (1..99) {
				my $uqn = $n*$target/100;
				if ($uqn<$zn) {
					$zq[$target]{$sample} = 0;
				} else {
					my $f = $uqn-int($uqn);
					$zq[$target]{$sample} = $sorted[$uqn-$zn]*(1-$f)+$f*$sorted[$uqn-$zn+1];
				}
			}
			$zq[100]{$sample} = $sorted[$#sorted];
		} else {
			foreach my $target (0..100) {
				$zq[$target]{$sample} = $nzq[$target]{$sample};
			}
		}
	}
	$self->datastore(\@nzq, "nonzeroquantiles");
	$self->datastore(\@zq, "zeroquantiles");
	print "...computed values at quantiles\n" if $self->param("verbose");
	$self->{$d}->{"nonzeroquantiles"} = \@nzq;
	$self->{$d}->{"zeroquantiles"} = \@zq;
	return \@nzq, \@zq;
}



sub use_samples {
	#specify which samples should be used
	#this is to make it easy to select subsets of a data set, without having to create separate data files for each subset
	#return value is the number of samples for which there's already some info
	my($self, @samples) = @_;
	
	my $d = $self->{'index'}->{'data'};
	$self->set_param('use_samples', @samples);
	my(%use, $already_seen);
	foreach my $sample (@samples) {
		$use{$sample} = 1;
		$already_seen++ if $self->{$d}->{'sampletotal'}{$sample};
	}
	$self->set_param('use_sample', \%use);
	print "...selected ", scalar @samples, " samples to use\n" if $self->param("verbose");
	
	my $dataref = $self->{$d}->{'values'};
	my $lengthsref = $self->{$d}->{'lengths'};
	my(%sampletotal, %genetotal, %sampleCount);
	foreach my $gene (keys %{$dataref}) {
		foreach my $sample (@samples) {
			my $v = $dataref->{$gene}{$sample};
			$v *= $lengthsref->{$gene}/1000 if defined $lengthsref->{$gene};
			$sampletotal{$sample} += $v;
			$genetotal{$gene} += $v;
			$sampleCount{$gene}++ if $v>0;
		}
	}
	$self->set_param("samples", keys %sampletotal);
	$self->set_param("sample_count", scalar keys %sampletotal);
	$self->{$d}->{'sampletotal'} = \%sampletotal;
	$self->{$d}->{'genetotal'} = \%genetotal;
	$self->{$d}->{'samplecount'} = \%sampleCount;
	
	#clear some cached computations...
	$self->reset(qw/minimal_nonmissing ubiquitous_genes gene_pairs gene_pairs_std graph ranked_genes gene_stability tmmnormweights tmmlogfold/);
	$self->set_param("minimal_nonmissing",0);
	return $already_seen;
}

sub samples_to_use {
	#get the list of samples that should be used, as set with use_samples(), or defaulting to all samples
	my($self) = @_;
	my $use_sample = $self->param('use_sample');
	return defined $use_sample ? keys %{$use_sample} : $self->param("samples");
}

sub compute_totalnorm_factors {
	#computes the normalization factors for the 'total expression' method (cpm)
	#if no target total is supplied, the average total sample expression is used, maintaining global counts
	my($self, $total) = @_;
	my(%norm, $count);
	
	my $d = $self->{'index'}->{'data'};
	my $tot = $self->{$d}->{'sampletotal'};
	my @samples = $self->samples_to_use();
	unless ($total) {
		foreach my $sample (@samples) {
			$total += $tot->{$sample};
			$count++;
		}
		$total /= $count;
	}
	foreach my $sample (@samples) {
		$norm{$sample} = $total/$tot->{$sample} if $tot->{$sample};
	}
	return %norm;
}

sub compute_total_in_ubiquitous {
	my($self, $total) = @_;
	my %norm;
	
	my $d = $self->{'index'}->{'data'};
	my @samples = $self->samples_to_use();
	my $v = $self->{$d}->{'values'};
	my $genes = $self->datarestore("ubiquitous_genes") || $self->identify_ubiquitous_genes();
	my %tot;
	
	foreach my $i (0..$#{$genes}-1) {
		my $gene = $genes->[$i];
		foreach my $sample (@samples) {
			$tot{$sample} += $v->{$gene}{$sample};
		}
	}
	
	return \%tot;
}

sub compute_prevnorm_factors {
	#computes the normalization factors for the 'total expression' method on ubiquitous genes
	#if no target total is supplied, the average total sample expression is used, maintaining global counts
	my($self, $total) = @_;
	my %norm;
	
	my $d = $self->{'index'}->{'data'};
	my @samples = $self->samples_to_use();
	my $v = $self->{$d}->{'values'};
	my $genes = $self->datarestore("ubiquitous_genes") || $self->identify_ubiquitous_genes();
	my %tot;
	
	foreach my $i (0..$#{$genes}-1) {
		my $gene = $genes->[$i];
		foreach my $sample (@samples) {
			$tot{$sample} += $v->{$gene}{$sample};
		}
	}
	
	$total ||= $self->avg(values %tot);
	foreach my $sample (@samples) {
		$norm{$sample} = $total/$tot{$sample} if $tot{$sample};
	}
	return %norm;
}


sub determine_minimal_nonmissing {
	#determine the minimal number of non-missing values to consider a gene
	my($self) = @_;
	
	my $fraction = $self->param("fraction");
	my @samples = $self->samples_to_use();
	my $sample_count = scalar @samples;
	my $minimal = int(0.5+$fraction*$sample_count);
	print "...determining minimal nonmissing: $minimal\n" if $self->param("verbose");
	return $self->set_param("minimal_nonmissing", $minimal);
}

sub identify_ubiquitous_genes {
	#identify genes observed in at least the desired fraction of samples, with a min/max expression percentile
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $verbose = $self->param("verbose");
	print "calling trim genes\n" if $verbose;
	my @genes = $self->trim_genes($self->param("lower_trim"), $self->param("upper_trim"));
	if (scalar @genes < 1000) {
		print "...identified too few ubiquitous genes (",scalar @genes,"), trying 5-95 range instead\n" if $verbose;
		@genes = $self->trim_genes(5, 95);
	}
	
	$self->datastore(\@genes, "ubiquitous_genes");
	print "...identified ", scalar @genes, " ubiquitous genes\n" if $verbose;
	return $self->{$d}->{"ubiquitous_genes"} = \@genes;
}
sub identify_ubiquitous_genes_inclusive {
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	print "calling trim genes\n" if $self->param("verbose");
	my @genes = $self->trim_genes($self->param("lower_trim_inclusive"), $self->param("upper_trim_inclusive"));
	$self->datastore(\@genes, "ubiquitous_genes_inclusive");
	print "...identified ", scalar @genes, " ubiquitous genes (inclusive cutoffs)\n" if $self->param("verbose");
	return $self->{$d}->{"ubiquitous_genes_inclusive"} = \@genes;
}

sub identify_trimmed_genes {
	#identify genes observed in at least the desired fraction of samples, with a min/max expression percentile
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my @genes = $self->trim_genes(5, 95);
	$self->datastore(\@genes, "trimmed_genes");
	print "...identified ", scalar @genes, " trimmed genes\n" if $self->param("verbose");
	return $self->{$d}->{"trimmed_genes"} = \@genes;
}

sub trim_genes {
	#identify genes observed in at least the desired fraction of samples, in the required percentile range
	my($self, $lower_trim, $upper_trim) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $dref = $self->{$d}->{'values'};
	my $minimal_nonmissing = $self->param("minimal_nonmissing") || $self->determine_minimal_nonmissing();
	
	my $q = $self->datarestore("nonzeroquantiles");
	($q) = $self->identify_values_at_percentiles() unless defined $q;
	my $qlower = $q->[$lower_trim];
	my $qupper = $q->[$upper_trim];
	
	print "...trimming genes ($lower_trim-$upper_trim)\n" if $self->param("verbose");
	my @samples = $self->samples_to_use();
	my %use;
	foreach my $sample (@samples) { $use{$sample}++; }
	my(%ubiquitous, %median);
	while (my($gene) = each %{$dref}) {
		my @counts;
		while (my($sample, $counts) = each %{$dref->{$gene}}) {
			next unless $use{$sample};
			next if $counts < $qlower->{$sample};
			next if $counts > $qupper->{$sample};
			push @counts, $counts;
		}
		my $seen = scalar @counts;
		if ($seen>=$minimal_nonmissing) {
			$ubiquitous{$gene} = $seen;
			@counts = sort {$a<=>$b} @counts;
			if ($minimal_nonmissing % 2) {
				$median{$gene} = $counts[$minimal_nonmissing/2];
			} else {
				$median{$gene} = ($counts[$minimal_nonmissing/2]+$counts[$minimal_nonmissing/2-1])/2;
			}
		}
	}
	
	return sort {$ubiquitous{$b}<=>$ubiquitous{$a} || $median{$b}<=>$median{$b}} keys %ubiquitous;
}

sub compare_genes_net {
	#compute logmaxmin for all pairs of ubiquitous genes
	my($self) = @_;
	my(@genes, %pair);
	
	my $d = $self->{'index'}->{'data'};
	my $v = $self->{$d}->{'values'};
	my $logv = $self->{$d}->{'logvalues'};
	my $verbose = $self->param("verbose");
	my $genes = $self->datarestore("ubiquitous_genes_inclusive") || $self->identify_ubiquitous_genes_inclusive();
	my $genes_to_use = $self->param("limit_genes_to_compare_net");
	$genes_to_use ||= scalar @{$genes};
	
	print "...comparing $genes_to_use genes for network\n" if $verbose;
	my $logmaxmincutoff = $self->param('logmaxmincutoff');
	my $minimal_nonmissing = $self->param("minimal_nonmissing") || $self->determine_minimal_nonmissing();
	my @samples = $self->samples_to_use();
	
	my $paircount;
	my %pairlmm;
	my $lowestlogmaxmin = 1e6;
	foreach my $i (0..$genes_to_use-2) {
		my $gene0 = $genes->[$i];
		foreach my $j ($i+1..$genes_to_use-1) {
			my $gene1 = $genes->[$j];
			#make sure the two genes share enough samples with non-zero values
			my %shared;
			foreach my $sample (@samples) {
				if ($v->{$gene0}{$sample} && $v->{$gene1}{$sample}) {
					$shared{$sample}[0] = $logv->{$gene0}{$sample};
					$shared{$sample}[1] = $logv->{$gene1}{$sample};
				}
			}
			my $shared = scalar keys %shared;
			next if $shared < $minimal_nonmissing;
			
			my @logratios;
			#compute, for this pair of genes, the observed expression ratios in all shared samples
			foreach my $sample (keys %shared) {
				push @logratios, $shared{$sample}[0] - $shared{$sample}[1];
			}
			
			#determine the highest and the lowest ratios, and compute the ratio between them (maxmin)
			@logratios = sort {$a<=>$b} @logratios;
			my $logmaxmin = $logratios[$#logratios]-$logratios[0];
			if ($logmaxmin<$lowestlogmaxmin) {
				print "...new lowest logmaxmin = $logmaxmin ($gene0 $gene1 $paircount)\n" if $verbose;
				$lowestlogmaxmin = $logmaxmin;
			}
			
			if ($logmaxmin>3) {
				#too dissimilar - discard
				next;
			}
			$pairlmm{"$gene0\t$gene1"} = $logmaxmin;
			$paircount++;
		}
	}
	print "...identified $paircount preliminary gene pairs\n" if $verbose;
	my @sorted = sort {$pairlmm{$a}<=>$pairlmm{$b}} keys %pairlmm;
	
	my $keep = 200;
	foreach my $i (0..$keep-1) {
		my($gene0, $gene1) = split /\t/, $sorted[$i];
		$pair{$gene0}{$gene1}{'similarity'} = 1/(1+$pairlmm{$sorted[$i]});
	}
	
	
	print "...keeping $keep gene pairs\n" if $verbose;
	$self->set_param("net_paircount", $paircount);
	$self->datastore(\%pair, "gene_pairs");
	return $self->{$d}->{'gene_pairs'} = \%pair;
}

sub compare_genes_stability {
	#compute std for all pairs of ubiquitous genes
	my($self, @genes) = @_;
	my($genes, %pair);
	
	
	
	my $d = $self->{'index'}->{'data'};
	my $v = $self->{$d}->{'values'};
	my $logv = $self->{$d}->{'logvalues'};
	my $genes_to_use;
	if (@genes) {
		$genes = \@genes;
		$genes_to_use = scalar @genes;
	} else {
		$genes = $self->datarestore("ubiquitous_genes_inclusive") || $self->identify_ubiquitous_genes_inclusive();
		$genes_to_use = $self->param("limit_genes_to_compare_stability");
		$genes_to_use ||= scalar @{$genes};
	}

	print "...comparing $genes_to_use genes for stability\n" if $self->param("verbose");
	my $minimal_nonmissing = $self->param("minimal_nonmissing") || $self->determine_minimal_nonmissing();
	my @samples = $self->samples_to_use();
	
	my $paircount;
	foreach my $i (0..$genes_to_use-1) {
		my $gene0 = $genes->[$i];
		foreach my $j (0..$genes_to_use-1) {
			next if $i==$j;
			my $gene1 = $genes->[$j];
			#make sure the two genes share enough samples with non-zero values
			my %shared;
			foreach my $sample (@samples) {
				if ($v->{$gene0}{$sample} && $v->{$gene1}{$sample}) {
					$shared{$sample}[0] = $logv->{$gene0}{$sample};
					$shared{$sample}[1] = $logv->{$gene1}{$sample};
				}
			}
			my $shared = scalar keys %shared;
			next if $shared < $minimal_nonmissing;
			
			my @logratios;
			#compute, for this pair of genes, the observed expression ratios in all shared samples
			foreach my $sample (keys %shared) {
				push @logratios, $shared{$sample}[0] - $shared{$sample}[1];
			}
			
			my($avg, $std) = $self->avgstd(@logratios);
			$pair{$gene0}{$gene1} = $std;
			$paircount++;
		}
	}
	print "...compared $paircount gene pairs\n" if $self->param("verbose");
	
	$self->datastore(\%pair, "gene_pairs_std") unless @genes;
	return $self->{$d}->{'gene_pairs_std'} = \%pair;
}

sub compute_gene_stability {
	my($self) = @_;
	my %stability;
	my $d = $self->{'index'}->{'data'};
	my $pairs = $self->datarestore('gene_pairs_std') || $self->compare_genes_stability();
	print "...computing gene stability\n" if $self->param("verbose");
	foreach my $gene0 (keys %{$pairs}) {
		my($sum, $n);
		while (my($gene1, $std) = each %{$pairs->{$gene0}}) {
			$sum += $std;
			$n++;
		}
		$stability{$gene0} = $sum/$n;
	}
	#$self->datastore(\%stability, "gene_stability");
	return $self->{$d}->{'gene_stability'} = \%stability;
}

sub create_graph {
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $pairs = $self->datarestore('gene_pairs') || $self->compare_genes_net();
	print "...creating graph\n" if $self->param("verbose");
	my $graph = Graph::Undirected->new();
	foreach my $gene0 (keys %{$pairs}) {
		foreach my $gene1 (keys %{$pairs->{$gene0}}) {
			$graph->add_weighted_edge($gene0, $gene1, $pairs->{$gene0}{$gene1}{'similarity'});
		}
	}
	$self->datastore($graph, "graph");
	return $self->{$d}->{'graph'} = $graph;
}

sub compute_gene_centrality {
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $graph = $self->datarestore('graph') || $self->create_graph();
	print "...computing gene centrality\n" if $self->param("verbose");
	my $ranker = Graph::Centrality::Pagerank->new();
	my $ranks = $ranker->getPagerankOfNodes('graph' => $graph, 'directed' => 0, 'useEdgeWeights' => 1);
	$self->datastore($ranks, "ranked_genes");
	my $ranked = scalar keys %{$ranks};
	print "...computed centrality for $ranked genes\n" if $self->param("verbose");
	$self->set_param("net_rankedgenes", $ranked);
	return $self->{$d}->{'ranked_genes'} = $ranks;
}




sub selectAtRandom {
	my($self, $target, @list) = @_;
	my(%use, $pick);
	my $n = scalar @list;
	while ($target > scalar keys %use) {
		$use{$list[rand $n]} = 1;
	}
	return \%use;
}



sub set_gene_weights {
	my($self, %weight) = @_;
	
	my $d = $self->{'index'}->{'data'};
	$self->datastore(\%weight, "gene_weights");
	print "...set weights for ", scalar keys %weight, " genes\n" if $self->param("verbose");
	return $self->{$d}->{'gene_weights'} = \%weight;
}

sub set_guide_genes {
	my($self, %guides) = @_;
	
	my $d = $self->{'index'}->{'data'};
	print "...set as guides: ", scalar keys %guides, " genes\n" if $self->param("verbose");
	return $self->{$d}->{'guide_genes'} = \%guides;
}

sub compute_gene_weights {
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $verbose = $self->param("verbose");
	my $method = $self->param("method");
	my %weight;
	if ($method eq "net") {
		my $ranks = $self->datarestore('ranked_genes') || $self->compute_gene_centrality();
		my $graph = $self->datarestore('graph') || $self->create_graph();
		my $genes = $self->set_param("genes_in_graph", scalar $graph->vertices());
		print "...computing gene weights\n" if $self->param("verbose");
		while (my($gene, $rank) = each %{$ranks}) {
			$rank *= $genes;
			$weight{$gene} = $rank-1 if $rank>1;
		}
	} elsif ($method eq "stability") {
		my $stability = $self->datarestore('gene_stability');
		unless (defined $stability) {
			$stability = $self->compute_gene_stability();
			$self->datastore($stability, "gene_stability");
		}
		#get stability parameter - max stability cutoff? number of genes to use?
		#filter genes by stability
		#translate to weights - equal weights for all genes that make the cut
		my @sorted = sort {$stability->{$a}<=>$stability->{$b}} keys %{$stability};
		my $use = $self->param("stabilityGenes") || 100;
		foreach my $i (0..$use-1) {
			my $gene = $sorted[$i];
			$weight{$gene} = 1;
		}
	} elsif ($method eq "genorm") {
		$self->compare_genes_stability(keys %{$self->{$d}->{'guide_genes'}});
		my $stability = $self->compute_gene_stability();
		
		my @sorted = sort {$stability->{$a}<=>$stability->{$b}} keys %{$stability};
		print "Using" if $verbose;
		my $use = $self->param("geNormGenes") || 3;
		foreach my $i (0..$use-1) {
			my $gene = $sorted[$i];
			$weight{$gene} = 1;
			print " ".$gene.":".$stability->{$gene} if $verbose;
		}
		print "\n" if $verbose;
	} else {
		#other methods shouldn't require a computation of gene weights...
		die "I don't know how to compute gene weights for method $method\n";
	}
	
	
	$self->datastore(\%weight, "gene_weights");
	my $gww = scalar keys %weight;
	print "...computed weights for $gww genes\n" if $self->param("verbose");
	$self->set_param("geneswithweight", $gww);
	return $self->{$d}->{'gene_weights'} = \%weight;
}


sub compute_logmedians {
	#median normalization
	my($self) = @_;
	
	my(%rv, %median);
	my $d = $self->{'index'}->{'data'};
	my $lv = $self->{$d}->{'logvalues'};
	my @samples = $self->samples_to_use();
	foreach my $gene (keys %{$lv}) {
		foreach my $sample (@samples) {
			$rv{$sample}{$gene} = $lv->{$gene}{$sample} if defined $lv->{$gene}{$sample};
		}
	}
	foreach my $sample (@samples) {
		my @sorted = sort {$a<=>$b} values %{$rv{$sample}};
		my $n = scalar @sorted;
		my $hn = int($n/2);
		if ($n % 2) {
			$median{$sample} = $sorted[$hn];
		} else {
			$median{$sample} = ($sorted[$hn+1]+$sorted[$hn])/2;
		}
	}
	
	return \%median;
}

sub compute_log_size_factors_ubiq {
	#DESeq using only ubiquitous genes
	my($self) = @_;
	
	my(%geomAvg, %sizeFactor);
	my $d = $self->{'index'}->{'data'};
	my $lv = $self->{$d}->{'logvalues'};
	my @samples = $self->samples_to_use();
	my $genes = $self->datarestore("ubiquitous_genes") || $self->identify_ubiquitous_genes();
	my $m = scalar @samples;
	foreach my $gene (@{$genes}) {
		my $sum;
		foreach my $sample (@samples) {
			$sum += $lv->{$gene}{$sample};
		}
		$geomAvg{$gene} = $sum/$m;
	}
	foreach my $sample (@samples) {
		my @v;
		foreach my $gene (@{$genes}) {
			push @v, $lv->{$gene}{$sample}-$geomAvg{$gene};
		}
		@v = sort {$a<=>$b} @v;
		my $n = scalar @v;
		my $hn = int($n/2);
		if ($n % 2) {
			$sizeFactor{$sample} = $v[$hn];
		} else {
			$sizeFactor{$sample} = ($v[$hn+1]+$v[$hn])/2;
		}
	}
	
	return \%sizeFactor;
}

sub compute_log_size_factors_nonzero {
	#DESeq using nonzero genes (in all samples)
	my($self) = @_;
	
	my(%geomAvg, %sizeFactor);
	my $d = $self->{'index'}->{'data'};
	my $lv = $self->{$d}->{'logvalues'};
	my @samples = $self->samples_to_use();
	my $m = scalar @samples;
	foreach my $gene (keys %{$lv}) {
		my @v;
		foreach my $sample (@samples) {
			push @v, $lv->{$gene}{$sample} if defined $lv->{$gene}{$sample};
		}
		if ($m==scalar @v) {
			my($avg) = $self->avgstd(@v);
			$geomAvg{$gene} = $avg;
		}
	}
	foreach my $sample (@samples) {
		my @v;
		foreach my $gene (keys %geomAvg) {
			push @v, $lv->{$gene}{$sample}-$geomAvg{$gene};
		}
		@v = sort {$a<=>$b} @v;
		my $n = scalar @v;
		my $hn = int($n/2);
		if ($n % 2) {
			$sizeFactor{$sample} = $v[$hn];
		} else {
			$sizeFactor{$sample} = ($v[$hn+1]+$v[$hn])/2;
		}
	}
	
	return \%sizeFactor;
}


sub compute_percentile {
	#Bullard's upper quartile normalization, extended to any percentile
	my($self, $target, $nonzero) = @_;
	
	$target ||= 0.75;
	
	my(%p, $q);
	my $logbase = log($self->param("base"));
	if ($nonzero) {
		$q = $self->datarestore("nonzeroquantiles");
		($q) = $self->identify_values_at_percentiles() unless defined $q;
	} else {
		$q = $self->datarestore("zeroquantiles");
		(undef, $q) = $self->identify_values_at_percentiles() unless defined $q;
	}
	foreach my $sample ($self->samples_to_use()) {
		if ($q->[100*$target]{$sample}) {
			$p{$sample} = log($q->[100*$target]{$sample})/$logbase;
		} else {
			print "No ", sprintf("%.0f", $target*100), " percentile value for sample $sample\n" if $self->param("verbose");
		}
	}
	return \%p;
}

sub compute_scaling_factors {
	#compute scaling factor for each sample
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $logv = $self->{$d}->{'logvalues'};
	my $base = $self->param("base");
	my $method = $self->param("method");
	my @samples = $self->samples_to_use();
	my(%corr);
	
	if ($method eq "none") {
		foreach my $sample (@samples) {
			$corr{$sample} = 1;
		}
	} elsif ($method eq "random") {
		foreach my $sample (@samples) {
			$corr{$sample} = 2**(rand(1)-.5);
		}
		%corr = $self->normalizeScalingFactors(%corr);
	} elsif ($method eq "cpm" || $method eq "total") {
		%corr = $self->compute_totalnorm_factors($self->param("total_size_for_normalizing"));
	} elsif ($method eq "ubiquitous") {
		%corr = $self->compute_prevnorm_factors();
	} elsif ($method eq "tmm") {
		%corr = $self->TMM();
	} elsif ($method eq "median") {
		my $medians = $self->compute_logmedians();
		my $avg = $self->avg(values %{$medians});
		foreach my $sample (@samples) {
			$corr{$sample} = $base**($avg-$medians->{$sample});
		}
	} elsif ($method eq "upperquartile") {
		my $upperquartiles = $self->compute_percentile(.75);
		my $avg = $self->avg(values %{$upperquartiles});
		foreach my $sample (@samples) {
			$corr{$sample} = $base**($avg-$upperquartiles->{$sample});
		}
	} elsif ($method eq "upperdecile") {
		my $upperdeciles = $self->compute_percentile(.9);
		my $avg = $self->avg(values %{$upperdeciles});
		foreach my $sample (@samples) {
			$corr{$sample} = $base**($avg-$upperdeciles->{$sample});
		}
	} elsif ($method eq "deseq") {
		my $sizefactors = $self->compute_log_size_factors_nonzero();
		my $avg = $self->avg(values %{$sizefactors});
		foreach my $sample (@samples) {
			$corr{$sample} = $base**($avg-$sizefactors->{$sample});
		}
	} elsif ($method eq "desequbiq") {
		my $sizefactors = $self->compute_log_size_factors_ubiq();
		my $avg = $self->avg(values %{$sizefactors});
		foreach my $sample (@samples) {
			$corr{$sample} = $base**($avg-$sizefactors->{$sample});
		}
	} else {
		my $weights;
		if ($method eq "genorm") {
			$weights = $self->compute_gene_weights();
		} elsif ($method eq 'guide') {
			$weights = $self->{$d}->{'guide_genes'};
		} elsif ($method eq '10ubiquitous') {
			my $genes = $self->datarestore("ubiquitous_genes") || $self->identify_ubiquitous_genes();
			$weights = $self->selectAtRandom(10, @{$genes});
		} elsif ($method eq '100ubiquitous') {
			my $genes = $self->datarestore("ubiquitous_genes") || $self->identify_ubiquitous_genes();
			$weights = $self->selectAtRandom(100, @{$genes});
		} elsif ($method eq 'allubiquitous') {
			my $genes = $self->datarestore("ubiquitous_genes") || $self->identify_ubiquitous_genes();
			my %w;
			foreach my $gene (@{$genes}) { $w{$gene} = 1; }
			$weights = \%w;
		} elsif ($method eq '100trimmed') {
			my $genes = $self->datarestore("trimmed_genes") || $self->identify_trimmed_genes();
			$weights = $self->selectAtRandom(100, @{$genes});
		} elsif ($method eq 'alltrimmed') {
			my $genes = $self->datarestore("trimmed_genes") || $self->identify_trimmed_genes();
			my %w;
			foreach my $gene (@{$genes}) { $w{$gene} = 1; }
			$weights = \%w;
		} else {
			$weights = $self->datarestore('gene_weights') || $self->compute_gene_weights();
		}
		
		#scale weights to a total of 1
		my $totalweight;
		while (my($gene, $weight) = each %{$weights}) {
			$totalweight += $weight;
		}
		my %w;
		while (my($gene, $weight) = each %{$weights}) {
			$w{$gene} = $weight/$totalweight;
		}
		
		#For each gene with some weight, compute the average of log-transformed expressions.
		#That average, minus the log-transformed expression in each sample, equals the
		# scaling factor to be used if normalization is based only on that gene.
		#Then combine, for each sample, the scaling factors suggested by all genes,
		# weighted by each gene's weight.
		my(%globalCorrection, %correctionWeight);
		while (my($gene, $weight) = each %w) {
			my %correction;
			foreach my $sample (@samples) {
				$correction{$sample} = $logv->{$gene}{$sample} if defined $logv->{$gene}{$sample};
			}
			my $avg = $self->avg(values %correction);
			foreach my $sample (@samples) {
				$globalCorrection{$sample} += ($avg - $correction{$sample})*$weight;
				$correctionWeight{$sample} += $weight;
			}
		}
		#...and scale so each sample's total correction weight is 1...
		my $totalCorrection;
		foreach my $sample (@samples) {
			next unless defined $correctionWeight{$sample};
			$globalCorrection{$sample} /= $correctionWeight{$sample};
			$totalCorrection += $globalCorrection{$sample};
		}
		#...and shift values such that the sum of all resulting scaling factors is zero.
		$totalCorrection /= scalar keys %globalCorrection;
		foreach my $sample (@samples) {
			next unless defined $correctionWeight{$sample};
			$globalCorrection{$sample} -= $totalCorrection;
			$corr{$sample} = $base**($globalCorrection{$sample});
		}
	}
	%corr = $self->normalizeScalingFactors(%corr) unless $method eq 'cpm';
	return $self->{$d}->{'scaling_factors'} = \%corr;
}

sub setScalingFactors {
	#set one or more scaling factors
	my($self, %data) = @_;
	
	my $d = $self->{'index'}->{'data'};
	while (my($sample, $factor) = each %data) {
		$self->{$d}->{'scaling_factors'}{$sample} = $factor if $factor>0;
	}
}

sub scalingFactor {
	#get one scaling factor
	my($self, $sample) = @_;
	
	my $d = $self->{'index'}->{'data'};
	return $self->{$d}->{'scaling_factors'}{$sample} || $self->param('defaultScalingFactor');
}

sub scalingFactors {
	#get all currently set scaling factors
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	return %{$self->{$d}->{'scaling_factors'}};
}

sub normalizeScalingFactors {
	my($self, %f) = @_;
	
	return unless %f;
	my $p;
	my @samples = keys %f;
	foreach my $sample (@samples) {
		$p += log($f{$sample});
	}
	$p /= scalar @samples;
	$p = exp($p);
	foreach my $sample (@samples) {
		$f{$sample} /= $p;
	}
	return %f;
}

sub normalize {
	#normalize all values
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $v = $self->{$d}->{'values'};
	my $method = $self->param("method");
	my @samples = $self->samples_to_use();
	
	if ($method eq "quantile") {
		return $self->quantile_normalization();
	}
	
	my $corr = $self->{$d}->{'scaling_factors'} || $self->compute_scaling_factors();
	#print "...scaling values\n" if $self->param("verbose");
	
	my %norm;
	foreach my $gene (keys %{$v}) {
		foreach my $sample (@samples) {
			$norm{$gene}{$sample} = $v->{$gene}{$sample} * $corr->{$sample} if defined $v->{$gene}{$sample};
		}
	}
	return $self->{$d}->{'normalized_values'} = \%norm;
}

sub quantile_normalization {
	#quantile normalization
	my($self) = @_;
	
	my(%rv, %sorted, @expr, @count, %norm);
	my $d = $self->{'index'}->{'data'};
	my $v = $self->{$d}->{'values'};
	my @samples = $self->samples_to_use();
	
	foreach my $gene (keys %{$v}) {
		foreach my $sample (@samples) {
			$rv{$sample}{$gene} = $v->{$gene}{$sample} if defined $v->{$gene}{$sample};
		}
	}
	foreach my $sample (@samples) {
		my @sorted = sort {$rv{$sample}{$b}<=>$rv{$sample}{$a}} keys %{$rv{$sample}};
		@{$sorted{$sample}} = @sorted;
		for (my $i=0;$i<scalar @sorted;$i++) {
			$expr[$i] += $rv{$sample}{$sorted[$i]};
			$count[$i]++;
		}
	}
	for (my $i=0;$i<scalar @expr;$i++) {
		$expr[$i] /= $count[$i] if $count[$i];
	}
	foreach my $sample (@samples) {
		for (my $i=0;$i<scalar @{$sorted{$sample}};$i++) {
			$norm{$sorted{$sample}[$i]}{$sample} = $expr[$i];
		}
	}
	return $self->{$d}->{'normalized_values'} = \%norm;
}

sub TMMnormWeights {
	my($self) = @_;
	
	my $d = $self->{'index'}->{'data'};
	my $v = $self->{$d}->{'values'};
	my $logv = $self->{$d}->{'logvalues'};
	my $logbase = log($self->param("base"));
	my $sampletotal = $self->{$d}->{'sampletotal'};
	my @samples = $self->samples_to_use();
	my %logst;
	foreach my $sample (@samples) {
		$logst{$sample} = log($sampletotal->{$sample})/$logbase if defined $sampletotal->{$sample};
	}
	my(%w, %m); #weights and log-fold changes
	my $genes = $self->datarestore("trimmed_genes") || $self->identify_trimmed_genes();
	
	foreach my $gene (@{$genes}) {
		foreach my $sample (@samples) {
			next unless my $gv = $v->{$gene}{$sample};
			my $st = $sampletotal->{$sample};
			$w{$gene}{$sample} = ($st-$gv)/$st/$gv;
			$m{$gene}{$sample} = $logv->{$gene}{$sample}-$logst{$sample};
		}
	}
	$self->datastore(\%m, 'tmmlogfold');
	$self->datastore(\%w, 'tmmnormweights');
	return \%m, \%w;
}

sub TMM {
	my($self) = @_;
	
	my $base = $self->param("base");
	my $w = $self->datarestore('tmmnormweights');
	my $m = $self->datarestore('tmmlogfold');
	($m, $w) = $self->TMMnormWeights() unless defined $w && defined $m;
	my @samples = $self->samples_to_use();
	my %factor;
	my $refsample = $self->param('tmmRefSample') || $samples[0];
	$factor{$refsample} = 1;
	
	foreach my $i (0..$#samples) {
		my $sample = $samples[$i];
		next if $sample eq $refsample;
		my %mgkr;
		foreach my $gene (keys %{$m}) {
			my $refm = $m->{$gene}{$refsample};
			my $sm = $m->{$gene}{$sample};
			next unless defined $refm && defined $sm;
			$mgkr{$gene} = $sm-$refm;
		}
		my @sorted = sort {$mgkr{$a}<=>$mgkr{$b}} keys %mgkr;
		my $n = scalar @sorted;
		@sorted = splice @sorted, $n*0.3, $n*0.4;
		my($tmm, $tmmw);
		foreach my $gene (@sorted) {
			my $wgkr = $w->{$gene}{$sample}+$w->{$gene}{$refsample};
			$tmm += $wgkr*$mgkr{$gene};
			$tmmw += $wgkr;
		}
		$factor{$sample} = $base**(-$tmm/$tmmw);
	}
	
	return %factor;
}

sub evaluate_normalization {
	my($self, $computeCutoffPercentile) = @_;
	my(@covHist, %uniform, %specific, $spec, @covs, %allCovs);
	my $d = $self->{'index'}->{'data'};
	my $normalized = $self->{$d}->{'normalized_values'};
	my @samples = $self->samples_to_use();
	my $halfTheSamples = (scalar @samples) / 2;
	my $covCutoff = $self->param("CoV-cutoff");
	while (my($gene) = each %{$normalized}) {
		my @vals = values %{$normalized->{$gene}};
		next if scalar @vals<scalar @samples;
		#next if $halfTheSamples>scalar @vals;
		$specific{$gene} = $spec if $spec = $self->jongeneel_specific(%{$normalized->{$gene}});
		my $cov = $self->coefficient_of_variation(@vals);
		$allCovs{$gene} = $cov;
		$uniform{$gene} = 1 if $cov<$covCutoff;
		if ($cov) {
			push @covs, $cov if $computeCutoffPercentile;
			$cov = log($cov)/$logtwo*10+50;
			$cov = 0 if $cov<0;
		}
		$covHist[$cov]++;
	}
	if ($computeCutoffPercentile) {
		@covs = sort {$a<=>$b} @covs;
		return $covs[$computeCutoffPercentile*scalar @covs];
	}
	
	return scalar keys %uniform, scalar keys %specific, \@covHist, \%uniform, \%specific, \%allCovs;
}

sub make_variant_on {
	my($self, $parent, $randomrange) = @_;
	my %new;
	$randomrange ||= .2;
	my $halfrandomrange = $randomrange/2;
	my @samples = $self->samples_to_use();
	my %proto;
	foreach my $sample (@samples) {
		$proto{$sample} = $parent->{'factors'}{$sample}*2**(rand($randomrange)-$halfrandomrange);
	}
	%proto = $self->normalizeScalingFactors(%proto);
	$self->setScalingFactors(%proto);
	$self->normalize();
	my($hom, $uniq, undef, $unif) = $self->evaluate_normalization();
	$new{'hom'} = $hom;
	$new{'uniq'} = $uniq;
	$new{'factors'} = \%proto;
	$new{'uniform'} = $unif;
	my($avg, $std) = $self->evaluate_gene_correlation();
	$new{'avg'} = $avg;
	$new{'std'} = $std;
	return \%new;
}

sub recombine_solutions {
	my($self, $p1, $p2) = @_;
	my(%new, @different_from);
	my @samples = $self->samples_to_use();
	foreach my $sample (@samples) {
		if (int(rand(2))) {
			$new{'factors'}{$sample} = $p1->{'factors'}{$sample};
			$different_from[2]++ if $new{'factors'}{$sample} != $p2->{'factors'}{$sample};
		} else {
			$new{'factors'}{$sample} = $p2->{'factors'}{$sample};
			$different_from[1]++ if $new{'factors'}{$sample} != $p1->{'factors'}{$sample};
		}
	}
	return unless $different_from[1] && $different_from[2];
	my %proto = %{$new{'factors'}};
	%proto = $self->normalizeScalingFactors(%proto);
	$self->setScalingFactors(%proto);
	$self->normalize();
	my($hom, $uniq, undef, $unif) = $self->evaluate_normalization();
	$new{'hom'} = $hom;
	$new{'uniq'} = $uniq;
	$new{'factors'} = \%proto;
	$new{'uniform'} = $unif;
	my($avg, $std) = $self->evaluate_gene_correlation();
	$new{'avg'} = $avg;
	$new{'std'} = $std;
	return 1, \%new;
}


sub evolution_strategy {
	my($self, $populationsize, $time_to_spend, $roundswithoutimprovement, @solution) = @_;
	my($previousbest, $roundsSinceImprovement, $time_started, $now, @best_series);
	my $time_to_end = ($time_started=time())+$time_to_spend;
	$roundswithoutimprovement ||= 10;
	$populationsize ||= 10;
	my @samples = $self->samples_to_use();
	my $verbose = $self->param("verbose");
	
	unless (scalar @solution) {
		foreach my $i (0..9) {
			foreach my $sample (@samples) {
				$solution[$i]{$sample} = 2**(rand(.2)-.1);
			}
		}
	}
	
	foreach my $i (0..$#solution) {
		my %proto = %{$solution[$i]};
		%proto = $self->normalizeScalingFactors(%proto);
		$solution[$i]{'factors'} = \%proto;
		$self->setScalingFactors(%proto);
		$self->normalize();
		my($hom, $uniq, undef, $unif) = $self->evaluate_normalization();
		$solution[$i]{'hom'} = $hom;
		$solution[$i]{'uniq'} = $uniq;
		$solution[$i]{'uniform'} = $unif;
		my($avg, $std) = $self->evaluate_gene_correlation();
		$solution[$i]{'avg'} = $avg;
		$solution[$i]{'std'} = $std;
	}
	
	my $iteration = 0;
	while (($now=time())<$time_to_end) {
		@solution = sort {
			$b->{'hom'} <=> $a->{'hom'} ||
			$a->{'uniq'} <=> $b->{'uniq'}
		} @solution;
		my $solutions = scalar @solution;
		if ($solutions>$populationsize) {
			@solution = @solution[0..$populationsize-1];
			$solutions = scalar @solution;
		}
		
		$iteration++;
		
		$best_series[($now-$time_started)/30] = join("\t", $solution[0]{'hom'}, $solution[0]{'uniq'});
		print join("\t", "(ES)", $now-$time_started, $solutions, $solution[0]{'hom'}."/".$solution[0]{'uniq'},
			$solution[$solutions-1]{'hom'}."/".$solution[$solutions-1]{'uniq'}
			), "\n" if $verbose;
		my $best = $solution[0]{'hom'};
		if ($best>$previousbest) {
			$roundsSinceImprovement = 0;
			$self->saveSolution("iter.".($iteration<10 ? "00$iteration" : ($iteration<100 ? "0$iteration" : $iteration)),
				{'uniform' => $solution[0]{'hom'}, 'specific' => $solution[0]{'uniq'},
					'avgcorr' => $solution[0]{'avg'}, 'stdcorr' => $solution[0]{'std'},
				},
				$solution[0]{'factors'}, undef, sort keys %{$solution[0]{'uniform'}});
		}
		$previousbest = $best;
		if ($roundsSinceImprovement>=$roundswithoutimprovement) {
			print "ES seems to have converged\n" if $verbose;
			last;
		}
		$roundsSinceImprovement++;
		
		foreach my $round (1..5) {
			push @solution, $self->make_variant_on($solution[0], .2);
			push @solution, $self->make_variant_on($solution[0], .1) if $roundsSinceImprovement>=5;
			push @solution, $self->make_variant_on($solution[0], .05) if $roundsSinceImprovement>=10;
			push @solution, $self->make_variant_on($solution[0], .01) if $roundsSinceImprovement>=20;
			push @solution, $self->make_variant_on($solution[1+int(rand(9))]);
			push @solution, $self->make_variant_on($solution[10+int(rand($solutions-10))]) if $solutions>10;
			
			my $pick1 = int(rand($solutions));
			my $pick2 = int(rand($solutions));
			while ($pick2 == $pick1) {
				$pick2 = int(rand($solutions));
			}
			my($success, $new) = $self->recombine_solutions($solution[$pick1], $solution[$pick2]);
			push @solution, $new if $success;
		}
	}
	
	$self->setScalingFactors(%{$solution[0]{'factors'}});
	$self->normalize();
	my($uniform, $specific);
	my($hom, $uniq, $covHist, $uniform, $specific, $allCov) = $self->evaluate_normalization();
	return $hom, $uniq, $covHist, $uniform, $specific, $allCov, \@best_series;
}

sub saveSolution {
	my($self, $name, $info, $factors, $allCov, @uniform) = @_;
	my $dir = $self->param("outdir") . "/solutions";
	mkdir $dir, 0755 unless -e $dir;
	open SOLDUMP, ">$dir/$name";
	foreach my $field (keys %{$info}) {
		print SOLDUMP join("\t", "#$field", $info->{$field}), "\n";
	}
	foreach my $sample (keys %{$factors}) {
		print SOLDUMP join("\t", "#factor", $sample, $factors->{$sample}), "\n";
	}
	print SOLDUMP join("\n", @uniform);
	close SOLDUMP;
}

sub evaluate_decorrelation {
	my($self) = @_;
	my $d = $self->{'index'}->{'data'};
	my $normalized = $self->{$d}->{'normalized_values'};
	my @samples = $self->samples_to_use();
	
	my $med = $self->compute_total_in_ubiquitous();
	my @sorted = sort {$med->{$a}<=>$med->{$b}} keys %{$med};
	my @sp;
	foreach my $gene (keys %{$normalized}) {
		my $sp = $self->spearman_presorted($normalized->{$gene}, \@sorted);
		push @sp, $sp;
	}
	my($avg, $std) = $self->mmad(@sp);
	return $avg, $std;
}

sub evaluate_gene_correlation {
	my($self, $rounds) = @_;
	$rounds ||= 10000;
	my $d = $self->{'index'}->{'data'};
	my $normalized = $self->{$d}->{'normalized_values'};
	my @samples = $self->samples_to_use();
	my %rank;
	my $genes;
	if (1) { #use all genes, or only ubiquitous?
		my @genes = keys %{$normalized};
		$genes = \@genes;
	} else {
		$genes = $self->datarestore("ubiquitous_genes") || $self->identify_ubiquitous_genes();
	}
	my $ngenes = scalar @{$genes};
	foreach my $gene (@{$genes}) {
		my @sorted = sort {$normalized->{$gene}{$a} <=> $normalized->{$gene}{$b}} @samples;
		foreach my $i (0..$#sorted) {
			$rank{$gene}{$sorted[$i]} = $i;
		}
	}
	my @sp;
	for (my $i=0;$i<$rounds;$i++) {
		my $g1 = $genes->[rand $ngenes];
		my $g2 = $genes->[rand $ngenes];
		while ($g1 eq $g2) {
			$g2 = $genes->[rand $ngenes];
		}
		my $sp = $self->spearman_double_presorted($rank{$g1}, $rank{$g2});
		push @sp, $sp;
	}
	
	my($avg, $std) = $self->avgstd(@sp);
	return $avg, $std, \@sp;
}


##########
sub avg {
	my($self, @values) = @_;
	my($value, $sum);
	
	return unless scalar @values;
	foreach $value (@values) {$sum += $value;}
	return $sum / scalar @values;
}
sub avgstd {
	my($self, @values) = @_;
	my($value, $sum, $avg, $devsqsum, $std);
	
	return unless 1<scalar @values;
	foreach $value (@values) {$sum += $value;}
	$avg = $sum / scalar @values;
	foreach $value (@values) {$devsqsum += ($value-$avg)**2;}
	$std = sqrt($devsqsum/(scalar @values-1));
	return $avg, $std;
}

sub mmad {
	my($self, @values) = @_;
	
	return unless my $n = scalar @values;
	@values = sort {$a<=>$b} @values;
	my $median = $values[$n/2]; #approx
	my @absdevs = map {abs($_-$median)} @values;
	@absdevs = sort {$a<=>$b} @absdevs;
	my $mad = $absdevs[$n/2]; #approx
	return $median, $mad;
}



sub jongeneel {
	my($self, @vals) = @_;
	
	my $total;
	@vals = sort {$b<=>$a} @vals;
	foreach my $i (1..$#vals) { $total += $vals[$i]; }
	return log(($vals[0]+1)/($total+1))/$logtwo;
}

sub jongeneel_specific {
	my($self, %vals) = @_;
	
	my($total, $max, $which);
	while (my($sample, $val) = each %vals) {
		$total += $val;
		if ($val>$max) {
			$max = $val;
			$which = $sample;
		}
	}
	return $which if $max>=$total/2;
}

sub coefficient_of_variation {
	my($self, @v) = @_;
	
	my($total, $s2);
	my $sampleCount = $self->param("sample_count");
	my $exp = 1/$sampleCount;
	foreach my $i (@v) {
		$total += $i;
	}
	return 1e6 unless $total>0;
	my $avg = $total/$sampleCount;
	foreach my $i (@v) {
		$s2 += ($i-$avg)**2;
	}
	if ($sampleCount > scalar @v) {
		#add deviations in unmentioned zero values
		my $avg2 = $avg**2;
		$s2 += ($sampleCount-scalar @v)*$avg2;
	}
	$s2 /= $sampleCount;
	return sqrt($s2)/$avg;
}

sub spearman_presorted {
	#$val is a reference to a hash with the normalized gene values, sample->value
	#$ref is a reference to an array of samples, sorted their characteristic value (e.g. median)
	my($self, $val, $ref) = @_;
	my $n = scalar @{$ref};
	my @sorted = sort {
		($val->{$ref->[$a]} <=> $val->{$ref->[$b]}) ||
		($a <=> $b)
	} (0..$n-1);
	my $c;
	foreach my $i (0..$n-1) { $c += ($sorted[$i]-$i)**2; }
	return 1 - 6*$c/$n/($n**2-1);
}

sub spearman_double_presorted {
	#$val and $ref are references to hashes of item->rank
	my($self, $val, $ref) = @_;
	my @items = keys %{$ref};
	my $n = scalar @items;
	my $c;
	foreach my $item (@items) { $c += ($val->{$item} - $ref->{$item})**2; }
	return 1 - 6*$c/$n/($n**2-1);
}



1;
__END__

=head1 NAME

Normalizer - Perl extension for normalizing gene expression tables

=head1 SYNOPSIS

	use Normalizer;
	$norm = new Normalizer();
	$norm->set_data(\%counts); #expected data structure: $counts{$gene_id}{$sample_id} = number
	
	$norm->set_param("method", "net"); #methods: none, cpm, total, ubiquitous, median, guide, net, stability, quantile
	$results = $norm->normalize();

=head1 ABSTRACT

  This should be the abstract for Normalizer.
  The abstract is used when making PPD (Perl Package Description) files.
  If you don't want an ABSTRACT you should also edit Makefile.PL to
  remove the ABSTRACT_FROM option.

=head1 DESCRIPTION

[documentation]

=head2 EXPORT

None by default.



=head1 SEE ALSO


=head1 AUTHOR

Gustavo Glusman, E<lt>Gustavo@SystemsBiology.orgE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright 2010 by Gustavo Glusman

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=cut

