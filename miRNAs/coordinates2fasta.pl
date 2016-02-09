#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 13/05/2015
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
getopts('g:', \%args);
my $genome = $args{g} || die "Supply a genome using -g\n";

my $fasta = Bio::DB::Fasta->new($genome);

<>;
while(<>){
    my ($uid,$mirname,$seqid,$start,$end,$strand) = split /,/;
    my $chrseq = $fasta->get_Seq_by_id($seqid);

    my $seq = $chrseq->subseq($start, $end);
    if ($strand eq '-'){
        $seq =~ tr/ATCG/TAGC/;
        $seq = reverse($seq);
    }
    
    print ">$uid\n$seq\n";
}

