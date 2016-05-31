#!/usr/bin/perl -w

die "Usage: $0 <period>\n" unless @ARGV == 1;
my ($period) = @ARGV;

print '{"state":[', "\n";
print ' {"n":0,"id":"B","trans":[{"in":"^","out":"^","to":1}]}', "\n";
for my $n (1..$period) {
    print ' {"n":', $n, ',"id":"S', $n,
    '","trans":[{"in":"0","out":"0","to":', $n+1,
    '},{"in":"1","out":"1","to":', $n+1,
    '},{"in":"$","out":"$","to":', $period+2,
    '}]}', "\n";
}
print ' {"n":', $period+1, ',"id":"C","trans":[{"out":"A","to":1}]}', "\n";
print ' {"n":', $period+2, ',"id":"E","trans":[]}', "\n";
print "]}\n";
