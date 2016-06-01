#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use List::Util qw(min max sum);
use Math::Round qw(:all);
use File::Temp;
use File::Basename;
use Cwd qw(abs_path);

my $rootdir = abs_path (dirname($0) . "/..");
my $bindir = "$rootdir/bin";
my $datadir = "$rootdir/data";

my $dnastore = "$bindir/dnastore";
my $editdist = "$bindir/editdist";

my $h74path = "$datadir/hamming74.json";
my $m2path = "$datadir/mixradar2.json";
my $m6path = "$datadir/mixradar6.json";
my $flusherpath = "$datadir/flusher.json";
my $syncprefix = "$datadir/sync";
my $syncsuffix = ".json";

my ($bitseqlen, $codelen) = (8192, 8);
my ($duprates, $maxdupsize, $allowdupoverlaps) = (.01, 4, 0);
my ($delrates, $maxdelsize) = (0, 10);
my ($subrates, $ivratio) = (0, 10);
my $band;
my ($reps, $trainalign) = (1, 1);
my $rndseed = 123456789;
my ($hamming, $mixradar2, $mixradar6, $syncfreq, $exacterrs, $keeptmp, $help);
my ($verbose, $dnastore_verbose) = (2, 2);
my ($colwidth, $cmdwidth) = (80, 200);

my $usage = "Usage: $0 [options]\n"
    . " -bits,-b <n>         number of random bits to encode (default $bitseqlen)\n"
    . " -codelen,-c <n>      codeword length in nucleotides (default $codelen)\n"
    . " -mix2,-2             precompose with mixradar length-2 block code\n"
    . " -mix6,-6             precompose with mixradar length-6 block code\n"
    . " -hamming,-g          precompose with Hamming(7,4) error-correcting code\n"
    . " -syncfreq,-q <n>     precompose with synchronizer that inserts control word after every n bits\n"
    . " -duprate,-d <n,n...> comma-separated list of duplication rates (default $duprates)\n"
    . " -maxdupsize,-m <n>   maximum length of duplications (default $maxdupsize)\n"
    . " -overlaps,-o         allow overlapping duplications\n"
    . " -delrate,-e <n,n...> comma-separated list of deletion rates (default $delrates)\n"
    . " -maxdelsize,-f <n>   maximum length of deletions (default $maxdelsize)\n"
    . " -subrate,-s <n,n...> comma-separated list of substitution rates (default $subrates)\n"
    . " -ivratio,-i <n>      transition/transversion ratio (default $ivratio)\n"
    . " -width,-w <n>        limit width of DP band for edit distance calculations\n"
    . " -reps,-r <n>         number of repeat simulations at each (subrate,duprate,delrate) setting\n"
    . " -trainalign,-t <n>   number of pairwise alignments in training set for error model (default $trainalign)\n"
    . " -exacterrs,-x        don't train error model on data; cheat by giving it simulation arguments\n"
    . " -seed <n>            seed random number generator (default $rndseed)\n"
    . " -keeptmp,-k          keep temporary files\n"
    . " -verbose,-v <n>      set verbosity level (default $verbose)\n"
    . " -ds-verbose,-u <n>   set verbosity level for dnastore (default $dnastore_verbose)\n"
    . " -help,-h             print this help text\n"
    ;

GetOptions ("bits=i" => \$bitseqlen,
	    "codelen=i" => \$codelen,
	    "hamming|g" => \$hamming,
	    "mix2|2" => \$mixradar2,
	    "mix6|6" => \$mixradar6,
	    "syncfreq|q=i" => \$syncfreq,
	    "duprate|d=s" => \$duprates,
	    "maxdupsize|m=f" => \$maxdupsize,
	    "overlaps" => \$allowdupoverlaps,
	    "delrate|e=s" => \$delrates,
	    "maxdelsize|n=f" => \$maxdelsize,
	    "subrate|s=s" => \$subrates,
	    "ivratio=f" => \$ivratio,
	    "width=i" => \$band,
	    "reps=i" => \$reps,
	    "exacterrs|x" => \$exacterrs,
	    "trainalign=i" => \$trainalign,
	    "seed=i" => \$rndseed,
	    "keeptmp" => \$keeptmp,
	    "verbose=i" => \$verbose,
	    "ds-verbose|u=i" => \$dnastore_verbose,
	    "help" => \$help)
    or die $usage;

die $usage if $help;

my %transition = (A=>'G',C=>'T',G=>'A',T=>'C');
my %transversion = (A=>[qw(C T)],C=>[qw(A G)],G=>[qw(C T)],T=>[qw(A G)]);

my $delext = 2 / $maxdelsize;  # hack
srand ($rndseed);
sub tempfile { return File::Temp->new (UNLINK => $keeptmp ? 0 : 1, DIR => "/tmp") }

my $ncontrols = $codelen < 4 ? 1 : 4;
my $machine = tempfile();
my $cmdstub = "$dnastore --verbose $dnastore_verbose";
my $ctrlargs = " --controls $ncontrols";
my $hamargs = $hamming ? " --compose-machine $h74path" : "";
my $mixargs = $mixradar2 ? " --compose-machine $m2path" : ($mixradar6 ? " --compose-machine $m6path" : "");
die "You cannot use -syncfreq with a machine that does not allow 3 or more control words\n"
    unless $ncontrols >= 3;
die "You cannot use -syncfreq without a machine that disambiguates flushing (-hamming, -mix2 or -mix6)\n"
    if $syncfreq && !($hamming || $mixradar2 || $mixradar6);
die "You cannot use -syncfreq and -hamming with a sync period that is not a multiple of 4"
    if $syncfreq && $hamming && $syncfreq % 4 != 0;
my $syncpath;
if (defined $syncfreq) {
    $syncpath = $syncprefix . $syncfreq . $syncsuffix;
    die "You cannot use -syncfreq with a sync period whose transducer does not exist (try 16, 32, or 128)"
	unless -e $syncpath;
}
my $flushargs = $syncfreq ? " --compose-machine $syncpath --compose-machine $flusherpath" : "";
syswarn ("$cmdstub --length $codelen $ctrlargs $flushargs $hamargs $mixargs --save-machine $machine");
my $cmd = "$cmdstub --length $codelen --load-machine $machine";

print " ", join (" ", qw(SubProb DupProb DelProb MeanEditsPerBit StDevEditsPerBit)), "\n";
my $nRows = 0;
for my $subrate (split /,/, $subrates) {
    for my $duprate (split /,/, $duprates) {
	for my $delrate (split /,/, $delrates) {
	    warn "subrate=$subrate duprate=$duprate delrate=$delrate\n" if $verbose;

	    my $errmodfh = tempfile();
	    unless ($exacterrs) {
		my $trainfh = tempfile();
		my $trainlen = 8*$bitseqlen;  # crude way to match training seq len to simulated seqs
		for my $n (1..$trainalign) {
		    warn "Generating training sequence #$n (length $trainlen)\n" if $verbose;
		    my $orig = randseq ([qw(A C G T)], $trainlen);
		    my @origpos = (0..length($orig)-1);
		    my $seq = lc $orig;
		    $seq = evolve (\@origpos, $seq, $duprate, 1, $maxdupsize, \&dup, "duplication", $allowdupoverlaps);
		    $seq = evolve (\@origpos, $seq, $subrate, 1, 1, \&subst, "substitution", 1);
		    $seq = evolve (\@origpos, $seq, $delrate, 1, $maxdelsize, \&del, "deletion", 1);
		    my $trainstock = stockholm($orig,$seq,\@origpos);
		    print $trainfh $trainstock;
		    warn $trainstock if $verbose >= 3;
		}

		my $errmod = syswarn ("$cmd --fit-error $trainfh --error-global");
		print $errmodfh $errmod;
		warn "Estimated error model:\n", $errmod if $verbose >= 2;
	    }

	    my @dist;
	    for my $rep (1..$reps) {
		warn "Starting repetition $rep of $reps\n" if $verbose;
		my $bitseq = randseq ([0,1], $bitseqlen);

		warn $bitseq, "\n" if $verbose >= 5;

		my $origdna = syswarn ("$cmd --raw --encode-bits $bitseq");
		chomp $origdna;

		my @origpos = (0..length($origdna)-1);
		my $dna = lc($origdna);

		$dna = evolve (\@origpos, $dna, $duprate, 1, $maxdupsize, \&dup, "duplication", $allowdupoverlaps);
		$dna = evolve (\@origpos, $dna, $subrate, 1, 1, \&subst, "substitution", 1);
		$dna = evolve (\@origpos, $dna, $delrate, 1, $maxdelsize, \&del, "deletion", 1);

		warn stockholm($origdna,$dna,\@origpos) if $verbose >= 3;
		
		warn $dna, "\n" if $verbose >= 5;
		$dna = uc($dna);

		my $seqfh = tempfile();
		print $seqfh ">seq\n$dna\n";
		my $exactArgs =  "--error-global --error-sub-prob $subrate --error-dup-prob $duprate --error-del-open $delrate --error-del-ext $delext";
		my $errArgs = $exacterrs ? $exactArgs : "--error-file $errmodfh";
		my $decoded = syswarn ("$cmd --raw --decode-viterbi $seqfh $errArgs");
		chomp $decoded;
		$decoded =~ s/[\^\$]//g;

		warn $decoded, "\n" if $verbose >= 5;

		warn "Length(bitseq)=", length($bitseq), " Length(decoded)=", length($decoded), "\n" if $verbose;
		
		my $dist = editDistance ($bitseq, $decoded, $band);
		warn "Edit distance: $dist\n" if $verbose;
		push @dist, $dist/$bitseqlen;
	    }
	    my $distmean = sum(@dist) / @dist;
	    my $distsd = sqrt (sum(map($_*$_,@dist)) / @dist - $distmean**2);
	    print join (" ", ++$nRows, $subrate, $duprate, $delrate, $distmean, $distsd), "\n";
	}
    }
}

sub evolve {
    my ($origpos, $seq, $rate, $minsize, $maxsize, $func, $desc, $allowoverlaps) = @_;
    my $nChanges = round ($rate * length($seq));
    my ($nActual, $nOverlaps) = (0, 0);
    for (my $n = 0; $n < $nChanges; ++$n) {
	my ($pos, $size) = randcoords ($seq, $minsize, $maxsize);
	if (substr ($seq, $pos, $size) =~ /[A-Z]/) {
	    ++$nOverlaps;
	    next unless $allowoverlaps;
	}
	warn "Mutating at (", $pos, "..", $pos+$size-1, ")\n" if $verbose >= 5;
	my ($newsub, $newpos) = &$func (substr ($seq, $pos, $size));
	die "Oops" unless @$newpos == length($newsub);
	substr ($seq, $pos, $size) = $newsub;
	splice (@{$origpos}, $pos, $size, map (defined($_) ? $origpos->[$pos+$_] : undef, @$newpos));
	++$nActual;
    }
    warn
	"Introduced ", $nActual, " ", $desc, ($nActual == 1 ? "" : "s"),
	($allowoverlaps
	 ? " ($nOverlaps with overlap)"
	 : (" (rejected $nOverlaps due to overlap)")),
	"\n" if $verbose >= 2;
    return $seq;
}

sub randcoords {
    my ($seq, $minsize, $maxsize) = @_;
    my $len = length ($seq);
    my $size = int (rand() * (min($len,$maxsize) + 1 - $minsize)) + $minsize;
    my $pos = int (rand() * ($len + 1 - $size));
    return ($pos, $size);
}

sub dup {
    my ($seq) = @_;
    warn "Duplicating $seq\n" if $verbose >= 4;
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
    for (my $col = 0; $col < length($origrow); $col += $colwidth) {
	$stock .= "\n" if $col > 0;
	$stock .= "old " . substr($origrow,$col,$colwidth) . "\n";
	$stock .= "new " . substr($newrow,$col,$colwidth) . "\n";
    }
    $stock .= "//\n";
    return $stock;
}

sub randseq {
    my ($alph, $len) = @_;
    return join ("", map ($alph->[rand() * @$alph], 1..$len));
}

sub syswarn {
    my ($syscmd) = @_;
    if ($verbose) {
	my $cmdwarn = $syscmd;
	if (length $cmdwarn > $cmdwidth) {
	    $cmdwarn = substr ($cmdwarn, 0, $cmdwidth) . "...";
	}
	warn $cmdwarn, "\n" if $verbose;
    }
    return `$syscmd`;
}

# Calculate Levenshtein edit distance
# (reimplemented as separate C++ program for speed)
sub editDistance {
    my ($x, $y, $band) = @_;
    $band = "" unless defined $band;
    return `$editdist "$x" "$y" $band` + 0;
}
