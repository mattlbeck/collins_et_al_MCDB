#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 
# -------------
# Combine seqs from multiple fasta files into one file
# Ignores all count stats etc
#
# -------------
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Std;
my %args;
getopts('c:h', \%args);

my $help = "Combine fastas into one file for efficient patman processing\n
USAGE combine_fasta.pl [fasta files] > output_fasta.fa\n";

die $help if $args{h};

my %unique;
while (<>){
    if (/^>/){
        my $seq = <>;
        chomp $seq;
        if ($args{c}){
            my ($read, $count) = $_ =~ /(\w+)\((\d+)\)/;
            $unique{$read} += $count;
        }
        else{
            $unique{$seq} = 1
        }
    }
}

foreach (keys %unique){
    if ($args{c}){
        print ">$_($unique{$_})\n$_\n";
    }
    else{
        print ">$_\n$_\n";
    }
}
