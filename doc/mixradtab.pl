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
    my $eofprob = $blocklen <= 3 ? "1/100" : "1/1000";
    my (%bitsPerSym, $states, $transitions);
    my @mixrad = `$mixrad -v -s $blocklen $eofprob`;
    for my $line (@mixrad[-5..-1]) {
	if ($line =~ /^(\d+) (\d+\/\d+)/) {
	    my ($radix, $bitsPerSym) = ($1, $2);
	    $bitsPerSym{$radix} = makeFrac ($bitsPerSym);
	} elsif ($line =~ /states: (\d+)/) {
	    $states = $1;
	} elsif ($line =~ /transitions: (\d+)/) {
	    $transitions = $1;
	}
    }
    print join (" & ", $blocklen, makeFrac($eofprob), $states, $transitions, map ($bitsPerSym{$_}, 2..4)), " \\\\\n";
}
