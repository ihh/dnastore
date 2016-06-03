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
my $waterprefix = "$datadir/water";
my $codesuffix = ".json";

for my $dep ($dnastore, $editdist, $h74path, $m2path, $m6path, $flusherpath) {
    die "Dependency '$dep' not found" unless -e $dep;
}

my ($bitseqlen, $codelen) = (8192, 8);
my ($duprates, $maxdupsize, $allowdupoverlaps) = (.01, 4, 0);
my ($delrates, $maxdelsize) = (0, 10);
my ($subrates, $ivratio) = (0, 10);
my $mutrates = 1;
my ($reps, $trainalign) = (1, 1);
my $rndseed = 123456789;
my ($ldpcdir, $hamming, $mixradar2, $mixradar6, $syncfreq, $waterfreq, $exacterrs, $keeptmp, $help, @dnastore_opt);
my ($verbose, $dnastore_verbose) = (2, 2);
my ($colwidth, $cmdwidth) = (80, 200);

my $usage = "Usage: $0 [options]\n"
    . " -bits,-b <n>          number of random bits to encode (default $bitseqlen)\n"
    . " -length,-l <n>        codeword length in nucleotides (default $codelen)\n"
    . " -mix2,-2              precompose with mixradar length-2 block code\n"
    . " -mix6,-6              precompose with mixradar length-6 block code\n"
    . " -hamming,-g           precompose with Hamming(7,4) error-correcting code\n"
    . " -syncfreq,-q <n>      precompose with synchronizer that inserts control word after every n bits (n=16,32,128...)\n"
    . " -watermark,-w <n[.m]> precompose with watermark synchronizer with period n bits and subperiod m bits (n=128, m=16...)\n"
    . " -ldpc,-c <dir>        wrap with a Low Density Parity Check code using Radford Neal's LDPC package downloaded from https://github.com/radfordneal/LDPC-codes\n"
    . " -mutrate,-m <n>       comma-separated list of scale factors for all mutation rates (default $mutrates)\n"
    . " -duprate,-d <n>       comma-separated list of duplication rates (default $duprates)\n"
    . " -maxdupsize,-z <n>    maximum length of duplications (default $maxdupsize)\n"
    . " -overlaps,-o          allow overlapping duplications\n"
    . " -delrate,-e <n>       comma-separated list of deletion rates (default $delrates)\n"
    . " -maxdelsize,-f <n>    maximum length of deletions (default $maxdelsize)\n"
    . " -subrate,-s <n,n...>  comma-separated list of substitution rates (default $subrates)\n"
    . " -ivratio,-i <n>       transition/transversion ratio (default $ivratio)\n"
    . " -reps,-r <n>          number of repeat simulations at each (subrate,duprate,delrate) setting\n"
    . " -trainalign,-t <n>    number of pairwise alignments in training set for error model (default $trainalign)\n"
    . " -exacterrs,-x         don't train error model on data; cheat by giving it simulation arguments\n"
    . " -seed <n>             seed random number generator (default $rndseed)\n"
    . " -keeptmp,-k           keep temporary files\n"
    . " -verbose,-v <n>       set verbosity level (default $verbose)\n"
    . " -ds-verbose,-u <n>    set verbosity level for dnastore (default $dnastore_verbose)\n"
    . " -help,-h              print this help text\n"
    ;

GetOptions ("bits=i" => \$bitseqlen,
	    "length|l=i" => \$codelen,
	    "hamming|g" => \$hamming,
	    "mix2|2" => \$mixradar2,
	    "mix6|6" => \$mixradar6,
	    "syncfreq|q=i" => \$syncfreq,
	    "watermark|w=s" => \$waterfreq,
	    "ldpc|c=s" => \$ldpcdir,
	    "mutrate|m=s" => \$mutrates,
	    "duprate|d=s" => \$duprates,
	    "maxdupsize|z=f" => \$maxdupsize,
	    "overlaps" => \$allowdupoverlaps,
	    "delrate|e=s" => \$delrates,
	    "maxdelsize|n=f" => \$maxdelsize,
	    "subrate|s=s" => \$subrates,
	    "ivratio=f" => \$ivratio,
	    "reps=i" => \$reps,
	    "exacterrs|x" => \$exacterrs,
	    "trainalign=i" => \$trainalign,
	    "seed=i" => \$rndseed,
	    "keeptmp" => \$keeptmp,
	    "verbose=i" => \$verbose,
	    "ds-verbose|u=i" => \$dnastore_verbose,
	    "ds-opt|o=s" => \@dnastore_opt,
	    "help" => \$help)
    or die $usage;

die $usage if $help;

my %transition = (A=>'G',C=>'T',G=>'A',T=>'C');
my %transversion = (A=>[qw(C T)],C=>[qw(A G)],G=>[qw(C T)],T=>[qw(A G)]);

my $delext = 2 / $maxdelsize;  # hack
srand ($rndseed);
sub tempfile { return File::Temp->new (UNLINK => $keeptmp ? 0 : 1, DIR => "/tmp", @_) }

if ($maxdupsize > $codelen/2) {
    warn "\nWARNING: maximum duplication size (-maxdupsize) is greater than half of codeword length (-length).\nSome duplications will go undetected by the error decoder!\n\n";
}

# LDPC
my ($ldpc_pchk, $ldpc_gen);
if (defined $ldpcdir) {
    $ldpc_pchk = tempfile (SUFFIX => '.pchk');
    $ldpc_gen = tempfile (SUFFIX => '.gen');
    my ($ldpc_checks, $ldpc_bits) = ($bitseqlen, 2 * $bitseqlen);  # hardcoded
    my $ldpc_seed = 1;
    my $ldpc_checks_per_col = 3;
    syswarn ("$ldpcdir/make-ldpc $ldpc_pchk $ldpc_checks $ldpc_bits $ldpc_seed evenboth $ldpc_checks_per_col no4cycle");
    syswarn ("$ldpcdir/make-gen $ldpc_pchk $ldpc_gen dense");
}

sub ldpc_encode {
    my ($bitseq) = @_;
    my $ldpc_src = tempfile (SUFFIX => '.src');
    my $ldpc_enc = tempfile (SUFFIX => '.enc');
    print $ldpc_src $bitseq, "\n";
    syswarn ("$ldpcdir/encode $ldpc_pchk $ldpc_gen $ldpc_src $ldpc_enc");
    my @ldpc_enc = <$ldpc_enc>;
    my $encseq = join ("", @ldpc_enc);
    $encseq =~ s/[^01]//g;
    return $encseq;
}

sub ldpc_decode {
    my ($recseq, $subprob) = @_;
    my $ldpc_rec = tempfile (SUFFIX => '.rec');
    my $ldpc_dec = tempfile (SUFFIX => '.dec');
    my $ldpc_msg = tempfile (SUFFIX => '.msg');
    print $ldpc_rec $recseq, "\n";
    my $ldpc_max_iters = 250;
    syswarn ("$ldpcdir/decode $ldpc_pchk $ldpc_rec $ldpc_dec bsc $subprob prprp $ldpc_max_iters");
    syswarn ("$ldpcdir/extract $ldpc_gen $ldpc_dec $ldpc_msg");
    my @ldpc_msg = <$ldpc_msg>;
    my $msg = join ("", @ldpc_msg);
    $msg =~ s/[^01]//g;
    return $msg;
}

# DNASTORE
my $ncontrols = $codelen < 4 ? 1 : 4;
my $machine = tempfile (SUFFIX => '.code.json');
my $cmdstub = "$dnastore --verbose $dnastore_verbose @dnastore_opt";
my $ctrlargs = " --controls $ncontrols";
my $hamargs = $hamming ? " --compose-machine $h74path" : "";
my $mixargs = $mixradar2 ? " --compose-machine $m2path" : ($mixradar6 ? " --compose-machine $m6path" : "");
die "You cannot use -syncfreq with a machine that does not allow 3 or more control words\n"
    unless $ncontrols >= 3;
die "You cannot use -syncfreq without a machine that disambiguates flushing (-hamming, -mix2 or -mix6)\n"
    if $syncfreq && !($hamming || $mixradar2 || $mixradar6);
die "You cannot use -syncfreq and -hamming with a sync period that is not a multiple of 4"
    if $syncfreq && $hamming && $syncfreq % 4 != 0;
my $flushargs = "";
if (defined $syncfreq) {
    my $syncpath = $syncprefix . $syncfreq . $codesuffix;
    die "You cannot use -syncfreq with a sync period whose transducer does not exist (try 16, 32, or 128)"
	unless -e $syncpath;
    $flushargs = " --compose-machine $syncpath --compose-machine $flusherpath";
} elsif (defined $waterfreq) {
    my $waterpath = $waterprefix . $waterfreq . $codesuffix;
    die "You cannot use -watermark with a sync period whose transducer does not exist (try 128)"
	unless -e $waterpath;
    $flushargs = " --compose-machine $waterpath";
}
syswarn ("$cmdstub --length $codelen $ctrlargs $flushargs $hamargs $mixargs --save-machine $machine");
my $cmd = "$cmdstub --length $codelen --load-machine $machine";

print " ", join (" ", qw(MutProb SubProb DupProb DelProb MeanEditsPerBit StDevEditsPerBit MedianEditsPerBit UpperQuartileEditsPerBit LowerQuartileEditsPerBit)), "\n";
my $nRows = 0;
for my $mutrate (split /,/, $mutrates) {
    for my $subrate (map ($mutrate*$_, split (/,/, $subrates))) {
	for my $duprate (map ($mutrate*$_, (split /,/, $duprates))) {
	    for my $delrate (map ($mutrate*$_, (split /,/, $delrates))) {
		warn "mutrate=$mutrate subrate=$subrate duprate=$duprate delrate=$delrate\n" if $verbose;

		my $errmodfh = tempfile (SUFFIX => '.err.json');
		unless ($exacterrs) {
		    my $trainfh = tempfile (SUFFIX => '.stk');
		    my $trainlen = $bitseqlen;  # crude way to match training seq len to simulated seqs
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

		    my $encseq = defined($ldpcdir) ? ldpc_encode($bitseq) : $bitseq;
		    
		    my $origdna = syswarn ("$cmd --raw --encode-bits $encseq");
		    chomp $origdna;

		    my @origpos = (0..length($origdna)-1);
		    my $dna = lc($origdna);

		    $dna = evolve (\@origpos, $dna, $duprate, 1, $maxdupsize, \&dup, "duplication", $allowdupoverlaps);
		    $dna = evolve (\@origpos, $dna, $subrate, 1, 1, \&subst, "substitution", 1);
		    $dna = evolve (\@origpos, $dna, $delrate, 1, $maxdelsize, \&del, "deletion", 1);

		    warn stockholm($origdna,$dna,\@origpos) if $verbose >= 3;
		    
		    warn $dna, "\n" if $verbose >= 5;
		    $dna = uc($dna);

		    my $seqfh = tempfile (SUFFIX => '.fa');
		    print $seqfh ">seq\n$dna\n";
		    my $exactArgs =  "--error-global --error-sub-prob $subrate --error-dup-prob $duprate --error-del-open $delrate --error-del-ext $delext";
		    my $errArgs = $exacterrs ? $exactArgs : "--error-file $errmodfh";
		    my $recseq = syswarn ("$cmd --raw --decode-viterbi $seqfh $errArgs");
		    chomp $recseq;
		    $recseq =~ s/[\^\$]//g;

		    warn $recseq, "\n" if $verbose >= 5;

		    my $ldpc_err_prob = min (.999, max (1 - (1-$subrate)*(1-$duprate*($maxdupsize+1)/2)*(1-$delrate*($maxdelsize+1)/2), 1e-9));
		    my $decseq = defined($ldpcdir) ? ldpc_decode($recseq,$ldpc_err_prob) : $recseq;
		    
		    warn "Length(bitseq)=", length($bitseq), " Length(decseq)=", length($decseq), "\n" if $verbose;
		    
		    my $dist = editDistance ($bitseq, $decseq);
		    warn "Edit distance: $dist", defined($ldpcdir) ? (" (pre-LDPC: ", editDistance($encseq,$recseq), ")") : (), "\n" if $verbose;
		    push @dist, $dist/$bitseqlen;
		}
		@dist = sort { $a <=> $b } @dist;
		my $distmean = sum(@dist) / @dist;
		my $distsd = sqrt (sum(map($_*$_,@dist)) / @dist - $distmean**2);
		my ($distmedian, $distuq, $distlq) = @dist[@dist/4, @dist/2, 3*@dist/4];
		print join (" ", ++$nRows, $mutrate, $subrate, $duprate, $delrate, $distmean, $distsd, $distmedian, $distuq, $distlq), "\n";
	    }
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
    my ($x, $y) = @_;
    return `$editdist "$x" "$y"` + 0;
}
