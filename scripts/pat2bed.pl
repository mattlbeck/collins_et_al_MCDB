#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 
# -------------
# Converts patman output to bed format as long as
# query also contains count information in brackets
# also makes certain that read order is by chr, start, then end
#
# -------------
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Std;
use Sort::Naturally qw(ncmp);
my %args;
getopts('m:B', \%args);
my $mapping_threshold = $args{m};

# Are counts in the id column instead?
my $countcolreg = q/(\w+)\(([\d.]+)\)/;

my @reads_a;
my %map_count;
my %filtered;

while ( <> ){
    my %read_h;
    @read_h{qw(hit query start end strand mismatches)} = split /\t/;
    @read_h{qw(seq count)} = ($read_h{query}, 'NA');

    $map_count{$read_h{seq}}++; 
    
    $filtered{$read_h{seq}} = 1 if $mapping_threshold && $map_count{$read_h{seq}} > $mapping_threshold;

    next if $filtered{$read_h{seq}};


    if ($read_h{query} =~ /$countcolreg/){
        # counts are located in query column
        @read_h{qw(seq count)} = ($1, $2);
    }
    elsif ($read_h{hit} =~ /$countcolreg/){
        # counts are located in hit column
        @read_h{qw(hit count)} = ($1, $2);
    }

    my $seq = $read_h{seq};
    # If counts are nowhere, Reported as NA
    push @reads_a, \%read_h; 
}

# Print a sorted list of bed formatted reads
foreach (sort readcmp @reads_a){
    if (!$args{B}){
        print join("\t", @{$_}{qw(hit start end seq count strand)}),"\n" if !$filtered{$_->{seq}};
    }
    else{
        print join("\t", @{$_}{qw(hit query start end strand mismatches)}),"\n" if !$filtered{$_->{seq}};
    }
}

my $filtered = scalar(keys %filtered);
print STDERR "Filtered $filtered seqs that mapped more than $mapping_threshold times\n" if $mapping_threshold;

# sort a read by seqid (naturally) then end then start
sub readcmp{
    ncmp($a->{hit}, $b->{hit}) ||
    $a->{start} <=> $b->{start} ||
    $a->{end} <=> $b->{end};
}
