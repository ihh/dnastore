#!/usr/bin/perl -w

print "digraph G {\n";
print "rankdir=LR;\n";

for ($n = 0; $n < 2; ++$n) {
    trans ("S", bin($n,1), $n % 2, undef);
}

for ($n = 0; $n < 4; ++$n) {
    trans (bin($n>>1,1), bin($n,2), $n % 2, undef);
}

for ($n = 0; $n < 8; ++$n) {
    trans (bin($n>>1,2), bin($n,3), $n % 2, undef);
}

for ($n = 0; $n < 16; ++$n) {
    trans (bin($n>>1,3), bin($n,4), $n % 2, undef);
}

for ($n = 0; $n < 16; ++$n) {
    trans (bin($n,4), "p1_".bin($n,4), undef, parity($n,1,2,4));
}

for ($n = 0; $n < 16; ++$n) {
    trans ("p1_".bin($n,4), "p2_".bin($n,4), undef, parity($n,1,3,4));
}

for ($n = 0; $n < 16; ++$n) {
    trans ("p2_".bin($n,4), "d1_".bin($n,3), undef, bit($n,1));
}

for ($n = 0; $n < 8; ++$n) {
    trans ("d1_".bin($n,3), "p3_".bin($n,3), undef, parity($n,2,3,4));
}

for ($n = 0; $n < 8; ++$n) {
    trans ("p3_".bin($n,3), "d2_".bin($n,2), undef, bit($n,2));
}

for ($n = 0; $n < 4; ++$n) {
    trans ("d2_".bin($n,2), "d3_".bin($n,1), undef, bit($n,3));
}

for ($n = 0; $n < 2; ++$n) {
    trans ("d3_".bin($n,1), "S", undef, bit($n,4));
}


for ($n = 0; $n < 2; ++$n) {
    state (bin($n,1));
    state ("d3_".bin($n,1));
}

for ($n = 0; $n < 4; ++$n) {
    state (bin($n,2));
    state ("d2_".bin($n,2));
}

for ($n = 0; $n < 8; ++$n) {
    state (bin($n,3));
    state ("d1_".bin($n,3));
    state ("p3_".bin($n,3));
}

for ($n = 0; $n < 16; ++$n) {
    state (bin($n,4));
    state ("p1_".bin($n,4));
    state ("p2_".bin($n,4));
}

print "}\n";

sub trans {
    my ($src, $dest, $in, $out) = @_;
    # No latex, too cluttered!
    #    print "$src -> $dest [label=\$", defined($in) ? "${in}_2" : "\\epsilon", "/", defined($out) ? "${out}_2" : "\\epsilon", "\$];\n";
        print "$src -> $dest [label=\"", defined($in) ? "${in}" : "", "/", defined($out) ? "${out}" : "", "\"];\n";
}

sub state {
    my ($state) = @_;
    print "$state [label=\"\"];\n";
}

sub bin {
    my ($n, $d) = @_;
    my $s = "";
    for (my $k = 0; $k < $d; ++$k) {
	$s = (($n % 2) ? "1" : "0") . $s;
	$n = int ($n / 2);
    }
    return $s;
}

sub bit {
    my ($n, $pos) = @_;
    return ($n & (1 << (4-$pos))) ? 1 : 0;
}

sub parity {
    my ($n, @pos) = @_;
    my $p = 0;
    for my $pos (@pos) { $p += bit($n,$pos) }
    return $p % 2;
}
