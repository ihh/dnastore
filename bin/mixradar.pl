#!/usr/bin/perl -w

use strict;
use Math::BigRat;
use Getopt::Long;
use List::Util qw(max);

my $usage = "Usage: $0 <msglen> <eofprob>\n";
my ($nomerge, $noprune, $norescale, $bigrat, $intervals, $getstats, $verbose);
GetOptions ("wide" => \$nomerge,
	    "deep" => \$noprune,
	    "pure" => \$norescale,
	    "rational" => \$bigrat,
	    "intervals" => \$intervals,
	    "stats" => \$getstats,
	    "verbose" => \$verbose)
    or die $usage;

$nomerge = $nomerge || $intervals;

sub newNumber {
    my ($val) = @_;
    return $bigrat ? Math::BigRat->new($val) : $val;
}

die $usage unless @ARGV == 2;
my ($msglen, $peofStr) = @ARGV;
my $peof = newNumber($peofStr);
my $pbit = (1 - $peof) / 2;
warn "pbit=$pbit peof=$peof" if $verbose;

# some constants
my ($epsilon, $div, $eof) = qw(e / $);
my %cprob = ('0' => $pbit,
	     '1' => $pbit,
	     $eof => $peof);
my @alph = (0, 1, $eof);
sub printable { my $word = shift; $word =~ s/\$/x/; return $word; }

my @radices = (2..4);

# generate prefix tree & input words
my @state = ( { word => '',
		dest => {},
		p => 1,
		start => 1 } );
my @prefixIndex = (0);
my @wordIndex;
while (@prefixIndex) {
    my $prefix = $state[shift @prefixIndex];
    for my $c (@alph) {
	my $childIndex = @state;
	my $child = { word => $prefix->{word} . $c,
		      dest => {},
		      p => $prefix->{p} * $cprob{$c} };
	$prefix->{dest}->{"$c$div$epsilon"} = $childIndex;
	push @state, $child;
	if ($c eq $eof || length($child->{word}) >= $msglen) {
	    $child->{input} = 1;
	    push @wordIndex, $childIndex;
	} else {
	    $child->{prefix} = 1;
	    push @prefixIndex, $childIndex;
	}
    }
}

# sort by probability & find intervals
my @sortedWordIndex = sort { $state[$b]->{p} <=> $state[$a]->{p}
			     || $state[$a]->{word} cmp $state[$b]->{word} } @wordIndex;
my $norm = 0;
for my $i (@sortedWordIndex) { $norm += $state[$i]->{p} }
for my $i (@sortedWordIndex) { $state[$i]->{p} /= $norm }
my $pmin = newNumber(0);
my $scale = newNumber(1);  # used to rescale probabilities after adjusting input intervals
my @allOutIndex = @sortedWordIndex;
my @finalIndex;
for my $i (@sortedWordIndex) {
    my $pmax = $pmin + $state[$i]->{p} * $scale;
    my $m = ($pmin + $pmax) / 2;
    # store
    $state[$i]->{A} = $pmin;
    $state[$i]->{B} = $pmax;
    $state[$i]->{m} = $m;
    $state[$i]->{D} = newNumber(0);
    $state[$i]->{E} = newNumber(1);
    $state[$i]->{outseq} = "";
    $pmin = $pmax;
    if ($verbose) {
	warn "P(", $state[$i]->{word}, ")=", $state[$i]->{p},
	" [A,B)=[", $state[$i]->{A}, ",", $state[$i]->{B}, ") m=", $m, "\n";
    }
    # generate tree
    my ($subtree, $final) = generateTree($i);
    push @allOutIndex, @$subtree;
    push @finalIndex, @$final;
    warn "Created ", @$final+0, " states to encode ", $state[$i]->{word}, "\n" if $verbose;
    # dynamically shrink input interval to just enclose all the output intervals actually used to encode it
    unless ($norescale) {
	my $new_pmax = max (map ($state[$_]->{E}, @$final));
	if ($new_pmax < $pmax) {
	    my $mul = (1 - $new_pmax) / (1 - $pmax);
	    warn "Shrinking B from $pmax to $new_pmax, increasing available space by factor of $mul\n" if $verbose;
	    $scale *= $mul;
	    $pmin = $new_pmax;
	}
    }
}

# subroutine to find a digit
sub findDigit {
    my ($m, $d, $e, $radix) = @_;
    my @d = map ($d + ($e - $d) * $_ / $radix, 0..$radix);
    my @digit = grep ($d[$_] <= $m && $d[$_+1] > $m, 0..$radix-1);
    if (@digit != 1) {
	warn "(D,E)=($d,$e) m=$m radix=$radix \@d=(@d)\n";
	die "Oops: couldn't find subinterval. Something's screwed";
    }
    my ($digit) = @digit;
    my ($new_d, $new_e) = @d[$digit,$digit+1];
    return ($digit, $new_d, $new_e);
}

# generate output trees
sub generateTree {
    my ($rootIndex) = @_;
    my @outputIndex = ($rootIndex);
    my (@subtree, @final);
    while (@outputIndex) {
	my $output = $state[shift @outputIndex];
	my ($a, $b, $m) = ($output->{A}, $output->{B}, $output->{m});
	for my $radix (@radices) {
	    my ($digit, $d, $e) = findDigit ($m, $output->{D}, $output->{E}, $radix);
	    my $outsym = $digit."_".$radix;
	    my $outseq = $output->{outseq} . (length($output->{outseq}) ? " " : "") . $outsym;
	    my $childIndex = @state;
	    my $child = { dest => {},
			  A => $a,
			  B => $b,
			  D => $d,
			  E => $e,
			  m => $m,
			  outseq => $outseq };
	    $output->{dest}->{"$epsilon$div$outsym"} = $childIndex;
	    push @state, $child;
	    if ($d >= $a && $e <= $b) {
		push @final, $childIndex;
	    } else {
		push @outputIndex, $childIndex;
	    }
	    push @subtree, $childIndex;
	}
    }
    return (\@subtree, \@final);
}

# if any nodes have a unique output sequence, remove all their descendants
my %nOutSeq;
for my $i (@allOutIndex) { ++$nOutSeq{$state[$i]->{outseq}} }
sub removeDescendants {
    my ($idx) = @_;
    while (my ($label, $destIdx) = each %{$state[$idx]->{dest}}) {
	warn "Removing #$destIdx (", $state[$destIdx]->{outseq}, ")\n" if $verbose;
	removeDescendants ($destIdx);
	$state[$destIdx]->{removed} = 1;
    }
    $state[$idx]->{dest} = {};
}

my @validOutIndex;
for my $i (@allOutIndex) {
    my $state = $state[$i];
    if (!$state->{removed}) {
	if (!$noprune && $nOutSeq{$state->{outseq}} == 1) {
	    warn "Pruning #$i: output sequence (", $state->{outseq}, ") is unique\n" if $verbose;
	    removeDescendants($i);
	}
	push @validOutIndex, $i;
    }
}

# do some analysis on the code
# find leaves for each codeword, collect statistics on number of digits per full codeword...
sub getLeaves {
    my ($idx) = @_;
    my $state = $state[$idx];
    my @kids = values %{$state->{dest}};
    return @kids ? map(getLeaves($_),@kids) : ($idx);
}
my @stats;
if ($getstats) {
    my %radixCount;
    for my $i (@sortedWordIndex) {
	my $word = $state[$i]->{word};
	my @outseqs = map ($state[$_]->{outseq}, getLeaves($i));
	grep (s/\d_//g, @outseqs);
	grep (s/ //g, @outseqs);
	@outseqs = sort { $a cmp $b } @outseqs;
	warn "Radices for ", $word, ": @outseqs\n" if $verbose;
	unless ($word =~ /\$/) {
	    for my $outseq (@outseqs) { ++$radixCount{$outseq} }
	}
    }
    my @radixSeqs = sort keys %radixCount;
    push @stats, "Frequencies of radix sequences for non-EOF codewords:\n";
    for my $radixSeq (@radixSeqs) {
	push @stats, $radixSeq, " ", $radixCount{$radixSeq}, "\n";
    }
    push @stats, "Mean bits/base for pure-radix sequences:\n";
    for my $radix (@radices) {
	my ($sum, $n) = (0, 0);
	for my $radixSeq (grep (/^$radix+$/, @radixSeqs)) {
	    $n += $radixCount{$radixSeq};
	    $sum += $radixCount{$radixSeq} * length($radixSeq);
	}
	push @stats, $radix, " ", $msglen/($sum/$n), "\n";
    }
}

# find output tree for each node, merge equivalence sets, assign IDs
my %equivIndex = ("()" => [0]);  # this takes care of the self-loop back to start
for my $outputIndex (reverse @validOutIndex) {
    my $output = $state[$outputIndex];
    my @destLabel = sort keys %{$output->{dest}};
    my @destIndex = map ($output->{dest}->{$_}, @destLabel);
    my @destSubtree = map ($state[$_]->{subtree}, @destIndex);
    my @destUndef = grep (!defined($destSubtree[$_]), 0..$#destIndex);
    if (@destUndef) { die "Oops. State $outputIndex child subtree(s) not defined (@destIndex[@destUndef]). Postorder?" }
    $output->{subtree} = '(' . join(',', map($destSubtree[$_].$destLabel[$_], 0..$#destLabel)) . ')';
    if ($nomerge && $output->{subtree} ne '()') {
	$output->{subtree} .= '[' . $output->{outseq} . ']';   # this will ensure uniqueness of every state
    }
    push @{$equivIndex{$output->{subtree}}}, $outputIndex;
}

for my $subtree (keys %equivIndex) {
    if (@{$equivIndex{$subtree}} > 1) {
	$equivIndex{$subtree} = [ sort { $a <=> $b } @{$equivIndex{$subtree}} ];
	warn "Merging (@{$equivIndex{$subtree}}) with subtree $subtree\n" if $verbose;
    }
}
for my $e (@{$equivIndex{"()"}}) {
    die $e if $e > 0 && keys(%{$state[$e]->{dest}});
}

my @equivIndex = map (exists($state[$_]->{subtree}) ? $equivIndex{$state[$_]->{subtree}}->[0] : $_, 0..$#state);
for my $state (@state) {
    for my $label (keys %{$state->{dest}}) {
	$state->{dest}->{$label} = $equivIndex[$state->{dest}->{$label}];
    }
}

my (%id, @uniqueState);
my ($nCodeStates, $nStates, $nTransitions) = (0, 0, 0);
for my $i (@equivIndex) {
    my $s = $state[$i];
    if (!$s->{removed}) {
	if (!exists $s->{id}) {
	    ++$nStates;
	    $nTransitions += keys(%{$s->{dest}});
	    if ($s->{start}) {
		$s->{id} = "S";
	    } elsif ($s->{prefix}) {
		$s->{id} = "P" . $s->{word};
	    } elsif ($s->{input}) {
		$s->{id} = "W" . printable($s->{word});
	    } else {
		$s->{id} = 'C' . (++$nCodeStates);
	    }
	    push @uniqueState, $s;
	}
    }
}

# print stats
if ($getstats) {
    print @stats;
    print "Number of states: ", $nStates, "\n";
    print "Number of transitions: ", $nTransitions, "\n";
    exit;
}

# print in dot format
print "digraph G {\n";
for my $state (@uniqueState) {
    my ($label, $shape, $style);
    if ($state->{start}) {
	$label = "START";
	$shape = "box";
	$style = "solid";
    } elsif ($state->{prefix}) {
	$label = $state->{word};
	$shape = "circle";
	$style = "solid";
    } elsif ($state->{input}) {
	$label = $state->{word};
	$shape = "doublecircle";
	$style = "solid";
    } else {
	$label = $nomerge ? $state->{outseq} : "";
	$shape = $nomerge ? "box" : "square";
	$style = "solid";
    }
    if ($intervals && exists $state->{A}) {
	$label .= "\n[A,B) = [" . ($state->{A}+0) . "," . ($state->{B}+0) . ")";
	$label .= "\n[D,E) = [" . ($state->{D}+0) . "," . ($state->{E}+0) . ")";
    }
    print " ", $state->{id}, " [style=$style;shape=$shape;label=\"$label\"];\n";
}
for my $src (@uniqueState) {
    while (my ($label, $destIndex) = each %{$src->{dest}}) {
	my $dest = $state[$destIndex];
	my $style = $label =~ /$epsilon$/
	    ? "dotted"
	    : ($label =~ /4$/
	       ? "bold"
	       : ($label =~ /3$/
		  ? "solid"
		  : ($label =~ /2$/
		     ? "dashed"
		     : "none")));
	my $srcStyle = $label =~ /^$epsilon/ ? "" : "dir=both;arrowtail=odot;";
	print " ", $src->{id}, " -> ", $dest->{id}, " [style=", $style, ";", $srcStyle, "arrowhead=empty;];\n";
    }
}
print "}\n";
