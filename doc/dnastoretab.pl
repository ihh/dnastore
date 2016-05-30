#!/usr/bin/perl -w

use strict;
use File::Temp;
use List::Util qw(sum);

my $dnastore = "../bin/dnastore";

my @noStart = ('-no-start' => '');
my @noEnd = ('-no-end' => '');
my @yesRepeats = ('repeats' => 1);
sub figref { my ($fig) = @_; return ('note' => '\subfigref{DNAStore}{' . $fig . '}') }
my @params = ({ l=>2, c=>0, @yesRepeats, @noStart, @noEnd, figref('a') },
	      { l=>2, c=>1, @yesRepeats, @noEnd, figref('b') },
	      { l=>2, c=>1, @yesRepeats, figref('c') },
	      { l=>2, c=>0, @noStart, @noEnd, figref('d') },
	      { l=>2, c=>1, @noEnd, figref('e') },
	      { l=>2, c=>1, figref('f') },
	      { l=>4, c=>0, @noStart, @noEnd },
	      { l=>4, c=>2 },
	      { l=>6, c=>0, @noStart, @noEnd },
	      { l=>6, c=>4 },
	      { l=>8, c=>0, @noStart, @noEnd },
	      { l=>8, c=>4 },
	      { l=>10, c=>0, @noStart, @noEnd },
	      { l=>10, c=>4 },
	      { l=>12, c=>0, @noStart, @noEnd },
	      { l=>12, c=>4 });

my %ignore = map (($_=>1), qw(note repeats l c));

for my $param (@params) {
    my $l = $param->{l};
    my $c = $param->{c};
    my $repeats = $param->{repeats} ? 'yes' : 'no';
    my $invreplen = $l >= 10 ? '4' : '';
    my $reparg = $param->{repeats} ? '-t 0' : '';
    my $note = $param->{note} || "";
    
    my $fh = File::Temp->new();
    my $filename = $fh->filename;

    my $args = join(" ",map("-$_ $param->{$_}",grep (!$ignore{$_}, keys %$param)));
    my $command = "$dnastore --length $l --controls $c --save-machine $filename --rate --print-controls $reparg $args";
    warn $command, "\n";
    my @line = `$command`;

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

    my @ttcontrols = map ("{\\tt $_}", @controls);
    if (@ttcontrols >= 2) {
	$ttcontrols[0] .= ' (start)';
	$ttcontrols[1] .= ' (end)';
    } elsif (@ttcontrols) {
	if (exists $param->{'-no-end'}) {
	    $ttcontrols[0] .= ' (start)';
	} else {
	    $ttcontrols[0] .= ' (start, end)';
	}
    }

    my @cols = ($l, $invreplen, $c, $repeats, $nStates, $nTrans, sprintf("%.4g",$rate));
    if (@controls <= 2) {
	print join (" & ", @cols, join(", ",@ttcontrols), $note), " \\\\\n";
    } else {
	for my $n (0..$#controls) {
	    print join (" & ",
			map ($n == 0 ? $_ : (" " x length($_)), @cols),
			$ttcontrols[$n] . ($n < $#controls ? ',' : ' '),
			$n == 0 ? $note : ''),
	    " \\\\\n";
	}
    }
}
