#!/usr/bin/perl -w

my $mixrad = "../bin/mixradar.pl";

sub makeFrac {
    my ($n) = @_;
    if ($n =~ /^(\d+)\/(\d+)$/) {
# Leave as num/denom... looks better in table
#	$n = "\$\\frac{$1}{$2}\$";
    }
    return $n;
}

for my $blocklen (2..6) {
    my $eofprob = $blocklen <= 3 ? "0.01" : "0.001";
    my (%symsPerBit, $states, $transitions);
# Uncomment for verbose on stderr:
#    my @mixrad = `$mixrad -v -s $blocklen $eofprob`;
    my @mixrad = `$mixrad -s $blocklen $eofprob`;
    for my $line (@mixrad[-5..-1]) {
	if ($line =~ /^(\d+) (\d+)\/(\d+)/) {
	    my ($radix, $bits, $sym) = ($1, $2, $3);
	    $symsPerBit{$radix} = sprintf ("%.4g", $sym/$bits);
	} elsif ($line =~ /states: (\d+)/) {
	    $states = $1;
	} elsif ($line =~ /transitions: (\d+)/) {
	    $transitions = $1;
	}
    }
    print join (" & ", $blocklen, $eofprob, $states, $transitions, map ($symsPerBit{$_}, 2..4)), " \\\\\n";
}
