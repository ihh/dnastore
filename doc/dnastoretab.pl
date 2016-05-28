#!/usr/bin/perl -w

use strict;
use File::Temp;
use List::Util qw(sum);

my $dnastore = "../bin/dnastore";

my @params = ({ l=>2, c=>1 },
	      { l=>4, c=>2 },
	      { l=>6, c=>4 },
	      { l=>8, c=>4 },
	      { l=>10, c=>4 },
	      { l=>12, c=>4 });

for my $param (@params) {
    my $l = $param->{l};
    my $c = $param->{c};
    
    my $fh = File::Temp->new();
    my $filename = $fh->filename;

    my @line = `$dnastore -l $l -c $c -S $filename -R --print-controls`;

    my @controls;
    if ($line[0] =~ /Control words: (.*)/) {
	my $c = $1;
	@controls = split (/ /, $c);
    }

    my @tot;
    while ($line[1] =~ /\b[01]: ([0-9\.e\+\-]+)/g) {
	push @tot, $1;
    }
    my $rate = sum(@tot) / @tot;

    my $nTrans = 0;
    while (<$fh>) {
	while (/{([^{}]*\"to\":[^{}]*)}/g) {
	    my $trans = $1;
	    $trans =~ s/"//g;
	    if ($trans =~ /in:(0|1|EOF|FLUSH)/ || $trans !~ /in:/) {
		++$nTrans;
	    }
	}
    }

    my $nStates = 0;
    $fh->seek(0,SEEK_SET);
    while (<$fh>) {
	if (/"id":/) {
	    ++$nStates;
	}
    }

    print "\\\\hline\n";
    print join (" & ", $l, $c, $nStates, $nTrans, sprintf("%.4g",$rate), $controls[0]), " \\\\\n";
    for my $c (@controls[1..$#controls]) {
	print join (" & ", "", "", "", "", "", $c), " \\\\\n";
    }
}
