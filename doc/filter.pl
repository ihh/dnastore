#!/usr/bin/perl -w

while (<>) {
    if (/\->/) {
	if (/\%/ || /EOF/ || /NULL/) {
	    s/"/\$/g;
	    s/\%/_/g;
	    s/([ACGT])/\\mbox\{$1\}/;
	    s/EOF/\\epsilon/g;
	    s/NULL/\\epsilon/g;
	    print;
	}
    } elsif (/\S/) {
	if (/"End/) {
	    s/"End.*"/"End"/;
	} elsif (/"Start/) {
	    s/"Start.*"/"Start"/;
	} else {
	    s/label="\S+ /label="/;
	    s/ +"/"/;
	    s/\*//g;
	}
	print;
    }
}
