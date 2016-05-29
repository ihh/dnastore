#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use Math::Round qw(:all);
use File::Temp;

my $dnastore = "../bin/dnastore";

my ($bitseqlen, $codelen) = (8192, 8);
my ($duprates, $maxdupsize, $allowdupoverlaps) = (.01, 4, 0);
my ($delrates, $maxdelsize) = (0, 10);
my ($subrates, $ivratio) = (0, 10);
my $band = 100;
my $reps = 1;
my $rndseed = 123456789;
my ($keeptmp, $help);
my ($verbose, $dnastore_verbose) = (2, 2);

my $usage = "Usage: $0 [options]\n"
    . " -bits,-b <n>         number of random bits to encode (default $bitseqlen)\n"
    . " -codelen,-c <n>      codeword length in nucleotides (default $codelen)\n"
    . " -duprate,-d <n,n...> comma-separated list of duplication rates (default $duprates)\n"
    . " -maxdupsize,-m <n>   maximum length of duplications (default $maxdupsize)\n"
    . " -overlaps,-o         allow overlapping duplications\n"
    . " -delrate,-e <n,n...> comma-separated list of deletion rates (default $delrates)\n"
    . " -maxdelsize,-f <n>   maximum length of deletions (default $maxdelsize)\n"
    . " -subrate,-s <n,n...> comma-separated list of substitution rates (default $subrates)\n"
    . " -ivratio,-i <n>      transition/transversion ratio (default $ivratio)\n"
    . " -width,-w <n>        width of DP band for edit distance calculations (default $band)\n"
    . " -reps,-r <n>         number of repeat simulations at each (subrate,duprate,delrate) setting\n"
    . " -seed <n>            seed random number generator (default $rndseed)\n"
    . " -keeptmp,-k          keep temporary files\n"
    . " -verbose,-v          print log messages on stderr\n"
    . " -dnastore-verbose,-u set verbosity level for dnastore\n"
    . " -help,-h             print this help text\n"
    ;

GetOptions ("bits=i" => \$bitseqlen,
	    "codelen=i" => \$codelen,
	    "duprate|d=f" => \$duprates,
	    "maxdupsize|m=f" => \$maxdupsize,
	    "overlaps" => \$allowdupoverlaps,
	    "delrate|e=f" => \$delrates,
	    "maxdelsize|n=f" => \$maxdelsize,
	    "subrate|s=f" => \$subrates,
	    "ivratio=f" => \$ivratio,
	    "width=i" => \$band,
	    "reps=i" => \$reps,
	    "seed=i" => \$rndseed,
	    "keeptmp" => \$keeptmp,
	    "verbose=i" => \$verbose,
	    "dnastore-verbose|u=i" => \$dnastore_verbose,
	    "help" => \$help)
    or die $usage;

die $usage if $help;

my %transition = (A=>'G',C=>'T',G=>'A',T=>'C');
my %transversion = (A=>[qw(C T)],C=>[qw(A G)],G=>[qw(C T)],T=>[qw(A G)]);
my $nDupOverlap;

my $delext = 2 / $maxdelsize;  # hack

srand ($rndseed);

for my $subrate (split /,/, $subrates) {
    for my $duprate (split /,/, $duprates) {
	for my $delrate (split /,/, $delrates) {
	    warn "subrate=$subrate duprate=$duprate delrate=$delrate\n" if $verbose;
	    my @dist;
	    for my $rep (1..$reps) {
		my $bitseq = join ("", map (rand() < .5 ? '1' : '0', 1..$bitseqlen));
		warn $bitseq, "\n" if $verbose >= 5;
		my $cmd = "$dnastore -l $codelen -v $dnastore_verbose";
		my $encodeCmd = "$cmd -r -b $bitseq";
		warn "$cmd -r -b ...", "\n" if $verbose >= 2;
		my $origdna = `$encodeCmd`;
		chomp $origdna;

		my @origpos = (0..length($origdna)-1);
		my $dna = lc($origdna);

		$nDupOverlap = 0;
		$dna = evolve (\@origpos, $dna, $duprate, 1, $maxdupsize, \&dup, "duplication", $allowdupoverlaps);
		warn $nDupOverlap, " overlapping duplications\n" if $verbose >= 2;
		
		$dna = evolve (\@origpos, $dna, $subrate, 1, 1, \&subst, "substitution", 1);
		$dna = evolve (\@origpos, $dna, $delrate, 1, $maxdelsize, \&del, "deletion", 1);

		warn stockholm($origdna,$dna,\@origpos) if $verbose >= 3;
		
		warn $dna, "\n" if $verbose >= 5;
		$dna = uc($dna);

		my $fh = File::Temp->new (UNLINK => $keeptmp ? 0 : 1, DIR => "/tmp");
		my $filename = $fh->filename;

		print $fh ">seq\n$dna\n";
		my $decodeCmd = "$cmd -r -V $filename --error-global --error-sub-prob $subrate --error-dup-prob $duprate --error-del-open $delrate --error-del-ext $delext";
		warn $decodeCmd, "\n" if $verbose >= 2;
		my $decoded = `$decodeCmd`;
		chomp $decoded;
		$decoded =~ s/[\^\$]//g;

		warn $decoded, "\n" if $verbose >= 5;

		warn "Length(bitseq)=", length($bitseq), " Length(decoded)=", length($decoded), "\n" if $verbose;
		
		my $dist = editDistance ($bitseq, $decoded, $band);
		warn "Edit distance: $dist\n" if $verbose;
		push @dist, $dist;
	    }
	    my $distmean = sum(@dist) / @dist;
	    my $distsd = sqrt (sum(map($_*$_,@dist)) / @dist - $distmean**2);
	    print "$subrate $duprate $delrate $distmean $distsd\n";
	}
    }
}

sub evolve {
    my ($origpos, $seq, $rate, $minsize, $maxsize, $func, $desc, $allowoverlaps) = @_;
    my $nChanges = round ($rate * length($seq));
    warn "Introducing ", $nChanges, " ", $desc, ($nChanges == 1 ? "" : "s"), "\n" if $verbose >= 2;
    for (my $n = 0; $n < $nChanges; ++$n) {
	my ($pos, $size) = randcoords ($seq, $minsize, $maxsize, $allowoverlaps);
	warn "Mutating at (", $pos, "..", $pos+$size-1, ")\n" if $verbose >= 3;
	my ($newsub, $newpos) = &$func (substr ($seq, $pos, $size));
	die "Oops" unless @$newpos == length($newsub);
	substr ($seq, $pos, $size) = $newsub;
	splice (@{$origpos}, $pos, $size, map (defined($_) ? $origpos->[$pos+$_] : undef, @$newpos));
    }
    return $seq;
}

sub randcoords {
    my ($seq, $minsize, $maxsize, $allowoverlaps) = @_;
    my $len = length ($seq);
    my ($size, $pos);
    do {
	$size = int (rand() * ($maxsize + 1 - $minsize)) + $minsize;
	$pos = int (rand() * ($len + 1 - $size));
    } while (substr ($seq, $pos, $size) =~ /[A-Z]/ && !$allowoverlaps);
    return ($pos, $size);
}

sub dup {
    my ($seq) = @_;
    warn "Duplicating $seq\n" if $verbose >= 4;
    ++$nDupOverlap if $seq =~ /[ACGT]/;
    return ($seq . uc($seq), [0..length($seq)-1,map(undef,1..length($seq))]);
}

sub del {
    my ($seq) = @_;
    warn "Deleting $seq\n" if $verbose >= 4;
    return ("",[]);
}

sub subst {
    my ($base) = @_;
    my $newbase;
    if (rand() < 1 / (1 + $ivratio)) {
	$newbase = $transversion{uc $base}->[rand() * 2];
    } else {
	$newbase = $transition{uc $base};
    }
    warn "Substituting $base -> $newbase\n" if $verbose >= 4;
    return ($newbase,[0]);
}

sub stockholm {
    my ($orig, $new, $new2orig) = @_;
    my $stock = "# STOCKHOLM 1.0\n";
    my ($origrow, $newrow) = ("", "");
    my $origpos = -1;
    for (my $newpos = 0; $newpos < length($new); ++$newpos) {
	my $n2o = $new2orig->[$newpos];
	if (defined $n2o) {
	    while ($origpos + 1 < $n2o) {
		++$origpos;
		$origrow .= substr ($orig, $origpos, 1);
		$newrow .= "-";
	    }
	    $origrow .= substr ($orig, $origpos = $n2o, 1);
	} else {
	    $origrow .= "-";
	}
	$newrow .= substr ($new, $newpos, 1);
    }
    while ($origpos + 1 < length($orig)) {
	++$origpos;
	$origrow .= substr ($orig, $origpos, 1);
	$newrow .= "-";
    }
    my $width = 80;
    for (my $col = 0; $col < length($origrow); $col += $width) {
	$stock .= "\n" if $col > 0;
	$stock .= "old " . substr($origrow,$col,$width) . "\n";
	$stock .= "new " . substr($newrow,$col,$width) . "\n";
    }
    $stock .= "//\n";
    return $stock;
}

# Calculate Levenshtein edit distance
sub editDistance {
    my ($x, $y, $band) = @_;
    my $b2 = int($band/2);
    my @cell;
    my ($prev_jmin, $prev_jmax);
    my ($xlen, $ylen) = map (length($_), $x, $y);
    my $inf = $xlen + $ylen;  # max possible edit distance
    for (my $i = 0; $i <= $xlen; ++$i) {
	my $jmin = max (0, $i - $b2);
	my $jmax = min ($ylen, $i + $b2);
	if ($i >= 2) { shift @cell }
	push @cell, [map ($inf, $jmin..$jmax)];
	my $xi = $i > 0 ? substr($x,$i-1,1) : '?x';
	warn "cell[0] = (@{$cell[0]})\n" if @cell > 0 && $verbose >= 9;
	warn "cell[1] = (@{$cell[1]})\n" if @cell > 1 && $verbose >= 9;
	for (my $j = $jmin; $j <= $jmax; ++$j) {
	    my $yj = $j > 0 ? substr($y,$j-1,1) : '?y';
	    my $sc = $inf;
	    if ($i == 0 && $j == 0) {
		$sc = 0;
	    }
	    if ($i > 0 && $j <= $prev_jmax) {
		$sc = min ($sc, $cell[0]->[$j-$prev_jmin] + 1);
		if ($j > $prev_jmin) {
		    $sc = min ($sc, $cell[0]->[$j-$prev_jmin-1] + ($xi eq $yj ? 0 : 1));
		}
	    }
	    if ($j > $jmin) {
		$sc = min ($sc, $cell[-1]->[$j-$jmin-1] + 1);
	    }
	    $cell[-1]->[$j-$jmin] = $sc;
	    warn "i=$i xi=$xi j=$j yj=$yj sc=$sc\n" if $verbose >= 9;
	}
	($prev_jmin, $prev_jmax) = ($jmin, $jmax);
    }
    my $dist = $cell[1]->[-1];
    return defined($dist) ? $dist : $inf;
}
