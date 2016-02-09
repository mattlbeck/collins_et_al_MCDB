#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 
# -------------
#
#
# -------------
use strict;
use warnings;
use Getopt::Std;
my %args;
getopts('l:', \%args);

my $length = $args{l} || 70;

while (<>){
    if(/^>/){
        print $_;
        my $seq = <>;
        chomp($seq);
        my @seqarr = split //, $seq;
        my $c = 0;
        my $newseq;
        foreach my $base (@seqarr){
            if($c == $length){
                print $newseq, "\n";    
                $newseq = "";
                $c = 0;
            }
            $c++;
            $newseq .= $base;
        }
        print $newseq,"\n";
    }
}
