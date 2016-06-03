#!/usr/bin/perl -w

use Math::Random::MT;
use Math::Combinatorics;
use Digest::MD5 qw(md5_hex);

die "Usage: $0 <period> [<subperiod>]\n" unless @ARGV == 1 || @ARGV == 2;
my ($period, $subperiod) = @ARGV;

if (!defined($subperiod) && $period =~ /^(\d+)\.(\d+)$/) {
    # syntactic suuuuugar: "128.16" instead of "128 16"
    # to make the Makefile easier, basically
    ($period, $subperiod) = ($1, $2);
}
$subperiod = $period + 1 unless defined $subperiod;

my $seed = hex (substr (md5_hex ("THEKID"), 0, 8));
my $gen = Math::Random::MT->new($seed);

my $states = $period + int($period/$subperiod);

print '{"state":[', "\n";
print ' {"n":0,"id":"B","trans":[{"in":"^","out":"^","to":1}]}', "\n";
for my $n (1..$states) {
    my $perm2 = int(1+$gen->rand(2));
    my $perm3 = int(1+$gen->rand(3*2));
    my $perm4 = int(1+$gen->rand(4*3*2));
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
	? ('{"out":"',$perm2[0],'","to":', $n+1,
	   '},{"out":"',$perm3[0],'","to":', $n+1,
	   '},{"out":"',$perm4[0],'","to":', $n+1,
	   '}')
	: ('{"in":"0","out":"',$perm2[0],'","to":', $n+1,
	   '},{"in":"1","out":"',$perm2[1],'","to":', $n+1,
	   '},{"in":"0","out":"',$perm3[0],'","to":', $n+1,
	   '},{"in":"1","out":"',$perm3[1],'","to":', $n+1,
	   '},{"in":"0","out":"',$perm4[0],'","to":', $n+1,
	   '},{"in":"1","out":"',$perm4[1],'","to":', $n+1,
	   '},{"in":"$","out":"$","to":', $states+2,
	   '}'),
	']}', "\n";
}
print ' {"n":', $states+1, ',"id":"C","trans":[{"out":"A","to":1}]}', "\n";
print ' {"n":', $states+2, ',"id":"E","trans":[]}', "\n";
print "]}\n";
