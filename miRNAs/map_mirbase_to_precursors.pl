#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 20/06/2015
# -------------
#
#
# -------------
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Std;
use Bio::DB::Fasta;
my %args;
getopts('p:m:g:', \%args);

my $prefile = $args{p};
my $maturefile = $args{m};
my $genome = $args{g};

my $genome = Bio::DB::Fasta->new($genome);

open(my $prefh, "<$prefile") or die "Can't open $prefile\n";
while(<$prefh>){
    my 
    my $chr = $genome->get_Seq_by_id($args{c});
    my $seq = $chr->subseq($args{s}, $args{e});
    if ($args{S} eq '-'){
        $seq =~ tr/ATCGU/TAGCA/;    
        $seq = reverse($seq);
    }
}
