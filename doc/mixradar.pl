#!/usr/bin/perl -w

die "Usage: $0 <string> <eofprob> radix1 radix2 radix3..." unless @ARGV > 2;
my ($text, $peof, @radix) = @ARGV;
$peof = eval($peof);

my %weight = ('0' => (1-$peof)/2,
	      '1' => (1-$peof)/2,
	      'x' => $peof);

my ($pmin, $pmax) = getRange ($text, %weight);
# warn "pmin=$pmin pmax=$pmax";
my @out = encode ($pmin, $pmax, @radix);
print "@out\n";

sub getRange {
    my ($seq, %weight) = @_;
    my @c = sort keys %weight;
    my $norm = 0;
    for my $c (@c) { $norm += $weight{$c} }
    my (%pmin, %pmax);
    my $p = 0;
    for my $c (@c) {
	$pmin{$c} = $p;
	$p += $weight{$c} / $norm;
	$pmax{$c} = $p;
#	warn "pmin{$c}=$pmin{$c} pmax{$c}=$pmax{$c}";
    }
    my ($pmin, $pmax) = (0, 1);
    for my $c (split(//,$seq)) {
	die "Allowed symbols: (@c)" unless exists $weight{$c};
	$new_pmin = $pmin + ($pmax - $pmin) * $pmin{$c};
	$new_pmax = $pmin + ($pmax - $pmin) * $pmax{$c};
#	warn "in=$c new_pmin=$new_pmin new_pmax=$new_pmax";
	($pmin, $pmax) = ($new_pmin, $new_pmax);
    }
    return ($pmin, $pmax);
}

sub encode {
    my ($pmin, $pmax, @radix) = @_;
    my $p = ($pmin + $pmax) / 2;
    my ($qmin, $qmax) = (0, 1);
    my @out;
    for my $radix (@radix) {
	my @q = map ($qmin + $_*($qmax-$qmin)/$radix, 0..$radix);
	my ($digit) = grep ($q[$_] <= $p && $q[$_+1] > $p, 0..$#q);
	push @out, $digit."_".$radix;
	($qmin, $qmax) = @q[$digit,$digit+1];
#	warn "out=${digit}_$radix new_qmin=$qmin new_qmax=$qmax q=(@q)";
	last if $qmin >= $pmin && $qmax <= $pmax;
    }
    die "Not enough radishes" unless $qmin >= $pmin && $qmax <= $pmax;
    return @out;
}
