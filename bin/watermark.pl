#!/usr/bin/perl -w

use Math::Random::MT;
use Math::Combinatorics;
use Digest::MD5 qw(md5_hex);
use Getopt::Long;

my ($subperiod, $subctrl, $nomixwater);
my $seedstr = "THEKID";
my $copies = 1;

my $usage = "Usage: $0 [options] <period>\n"
    .       "    -copies <n> number of copies of control word after each period\n"
    .       "       -sub <n> specify a sub-period, send a watermark bit after each sub-period\n"
    .       "       -subctrl send a control word after each sub-period, instead of a watermark bit\n"
    .       "    -nomixwater do not use mixed-radix encoding for watermark bits\n"
    .       " -seed <string> seed string for random-number generator (default is '$seedstr')\n";

GetOptions ("copies=i" => \$copies,
	    "sub=i" => \$subperiod,
	    "subctrl" => \$subctrl,
	    "nomixwater" => \$nomixwater,
	    "seed=s" => \$seedstr)
    or die $usage;

die $usage unless @ARGV == 1;
my ($period) = @ARGV;

die "-subctrl without -sub is useless\n" if defined($subctrl) && !defined($subperiod);
$subperiod = $period + 1 unless defined $subperiod;
$subperiod = int($subperiod);
die "Subperiod must be >1\n" unless $subperiod > 1;

my $seed = hex (substr (md5_hex ($seedstr), 0, 8));
my $gen = Math::Random::MT->new($seed);

my $states = $period + int($period/$subperiod);

print '{"state":[', "\n";
print ' {"n":0,"id":"B","trans":[{"in":"^","out":"^","to":1}]}', "\n";
for my $n (1..$states) {
    my $perm2 = $nomixwater ? 1 : int(1+$gen->rand(2));
    my $perm3 = $nomixwater ? 1 : int(1+$gen->rand(3*2));
    my $perm4 = $nomixwater ? 1 : int(1+$gen->rand(4*3*2));
    my $comb2 = Math::Combinatorics->new (count=>2, data=>[qw(i j)]);
    my $comb3 = Math::Combinatorics->new (count=>3, data=>[qw(x y z)]);
    my $comb4 = Math::Combinatorics->new (count=>4, data=>[qw(p q r s)]);
    my (@perm2, @perm3, @perm4);
    while ($perm2-- > 0) { @perm2 = $comb2->next_permutation() }
    while ($perm3-- > 0) { @perm3 = $comb3->next_permutation() }
    while ($perm4-- > 0) { @perm4 = $comb4->next_permutation() }
    print ' {"n":', $n, ',"id":"S', $n,
    '","trans":[',
    $n % $subperiod == 0
	? (defined($subctrl)
	   ? ('{"out":"B","to":', $n+1, '}')
	   : ('{"out":"',$perm2[0],'","to":', $n+1,
	      '},{"out":"',$perm3[0],'","to":', $n+1,
	      '},{"out":"',$perm4[0],'","to":', $n+1,
	      '}'))
	: ('{"in":"0","out":"',$perm2[0],'","to":', $n+1,
	   '},{"in":"1","out":"',$perm2[1],'","to":', $n+1,
	   '},{"in":"0","out":"',$perm3[0],'","to":', $n+1,
	   '},{"in":"1","out":"',$perm3[1],'","to":', $n+1,
	   '},{"in":"0","out":"',$perm4[0],'","to":', $n+1,
	   '},{"in":"1","out":"',$perm4[1],'","to":', $n+1,
	   '},{"in":"$","out":"$","to":', $states+$copies+1,
	   '}'),
	']}', "\n";
}
for my $c (1..$copies-1) {
    print ' {"n":', $states+$c, ',"id":"C', $c, '","trans":[{"out":"A","to":', $states+$c+1, '}]}', "\n";
}
print ' {"n":', $states+$copies, ',"id":"C', ($copies == 1 ? "" : $copies), '","trans":[{"out":"A","to":1}]}', "\n";
print ' {"n":', $states+$copies+1, ',"id":"E","trans":[]}', "\n";
print "]}\n";
