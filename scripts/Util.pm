#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 
# -------------
# A collection of Utility methods
#
# -------------
use strict;
use warnings;
package Util;

# Sort method for sorting stuff like seqids that contain
# a combination of characters and numbers
sub sortseqids {
    return sort seqids @_;
}
sub seqids{
   (my $an = $a) =~ s[(\d+)][pack "N", $1]ge;
   (my $bn = $b) =~ s[(\d+)][pack "N", $1]ge;
   $an cmp $bn;
}

# Check for overlap between two sequences mapped to the same reference
sub overlaps{
   my($aStart, $aEnd, $bStart, $bEnd) = map @$_, @_;
   return 1 if $aEnd >= $bStart and $bEnd >= $aStart;
   return 0;
}
1;
