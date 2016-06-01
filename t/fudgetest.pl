#!/usr/bin/env perl

use warnings;
use File::Temp;

my $dryrun = @ARGV == 2 && $ARGV[0] eq "-n";
die "Usage: $0 [-n] <test>" unless @ARGV == 1 || $dryrun;
my $test = pop @ARGV;

open MAKE, "make -n $test |";
my %seen;
while (<MAKE>) {
    if (/testexpect.pl.* (\S+)$/) {
	my $out = $1;
	if ($seen{$out}++) {
	    print;
	    chomp;
	    system $_ unless $dryrun;
	} else {
	    s/ (\S+)$/ >$1/;
	    s/^t.testexpect.pl //;
	    print;
	    chomp;
	    system $_ unless $dryrun;
	}
    }
}
close MAKE;
