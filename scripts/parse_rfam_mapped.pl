#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 14/03/2014
# -------------
#
#
# -------------
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Std;
my %args;
getopts('', \%args);

my %ucounts;
print(join(",",qw(seq hit sp mismatches)),"\n");
while(<>){
    chomp;
    my ($hit, $query, $start, $end, $strand, $mismatch) = split /\t/;
    print STDERR  "Description contains commas: $hit\n" if $hit =~ /,/;
    $hit =~ s/,/\./g;

    my ($id, $desc, $sp) = $hit =~ /(.+);(.+);.+\/\d+-\d+\s+\d+:(.+)/;
    die "Can't parse\n$hit" if !$id or !$desc or !$hit;

    print(join(",",$query,$desc,$sp,$mismatch),"\n")# if !$ucounts{$query};
    #$ucounts{$query} = 1;
}
