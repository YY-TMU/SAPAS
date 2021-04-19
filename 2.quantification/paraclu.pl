#!/usr/bin/perl -w

use strict;
use List::Util qw(min max);

my %data;

warn "reading...\n";

while (<>) {
    s/#.*//;  # ignore comments
    next unless /\S/;  # skip blank lines

    my ($seq, $strand, $pos, $value) = split;
    my $key = "$seq $strand";
    push @{$data{$key}}, [ $pos, $value ];
}

warn "clustering...\n";

print "# sequence, strand, start, end, sites, sum of values, min d, max d\n";

for my $key (sort keys %data) {  # iterate over sequences / strands
    my ($seq, $strand) = split " ", $key;
    my $sites = $data{$key};

    @$sites = sort { $$a[0] <=> $$b[0] } @$sites;  # sort by position

    my $clusters = all_clusters($sites);

    for my $c (@$clusters) {
	my ($beg, $end, $tot, $sit, $min, $max) = @$c;
	my $beg_pos = $$sites[$beg][0];
	my $end_pos = $$sites[$end][0];
	printf "$seq\t$strand\t$beg_pos\t$end_pos\t$sit\t$tot\t%.3g\t%.3g\n",
	$min, $max;
    }
}

### Generic code to find clusters in a sparse sequence of values: ###

sub all_clusters {
    our $inf = 1e100;  # hopefully much bigger than any value in the input
    our $sites = shift;  # input: reference to array of site locations & values
    our $clusters = [];  # output: reference to array of clusters
    get_clusters(0, $#$sites, -$inf);
    return $clusters;
}

# get clusters of sites between beg and end with density > min_density
sub get_clusters {
    our ($clusters, $inf);
    my ($beg, $end, $min_density) = @_;

    my ($prefix, $pmin, $ptot, $psit) = weakest_prefix($beg, $end);
    my ($suffix, $smin, $stot, $ssit) = weakest_suffix($beg, $end);
    $ptot == $stot and $psit == $ssit or die "internal error!";
    my $max_density = min $pmin, $smin;

    unless ($max_density == $inf) {
	my $break = $pmin < $smin ? $prefix + 1 : $suffix;
	my $new_min = max $min_density, $max_density;
	get_clusters($beg, $break-1, $new_min);
	get_clusters($break, $end, $new_min);
    }

    push @$clusters, [ $beg, $end, $ptot, $psit, $min_density, $max_density ]
	if $max_density > $min_density;
}

# get least dense prefix (and total of values & sites)
sub weakest_prefix {
    our ($sites, $inf);
    my ($beg, $end) = @_;

    my $beg_pos = $$sites[$beg][0];
    my $min_density = $inf;
    my $min_prefix = $end;
    my $tot = 0;
    my $sit = 0;

    for (my $i = $beg; $i < $end; ++$i) {
	$tot += $$sites[$i][1];
	next if $$sites[$i][0] == $$sites[$i+1][0];  # idiot-proofing
	++$sit;
	my $dist = $$sites[$i+1][0] - $beg_pos;
	my $density = $tot / $dist;
	if ($density < $min_density) {
	    $min_prefix = $i;
	    $min_density = $density;
	}
    }

    $tot += $$sites[$end][1];
    ++$sit;
    return ($min_prefix, $min_density, $tot, $sit);
}

# get least dense suffix (and total of values & sites)
sub weakest_suffix {
    our ($sites, $inf);
    my ($beg, $end) = @_;

    my $end_pos = $$sites[$end][0];
    my $min_density = $inf;
    my $min_suffix = $beg;
    my $tot = 0;
    my $sit = 0;

    for (my $i = $end; $i > $beg; --$i) {
	$tot += $$sites[$i][1];
	next if $$sites[$i][0] == $$sites[$i-1][0];  # idiot-proofing
	++$sit;
	my $dist = $end_pos - $$sites[$i-1][0];
	my $density = $tot / $dist;
	if ($density < $min_density) {
	    $min_suffix = $i;
	    $min_density = $density;
	}
    }

    $tot += $$sites[$beg][1];
    ++$sit;
    return ($min_suffix, $min_density, $tot, $sit);
}
