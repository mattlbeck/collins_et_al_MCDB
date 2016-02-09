#!/usr/bin/perl
# -------------
# file: merge_mirdeep.pl
# author: Matthew Beckers
# date: 16/04/2014
# -------------
# Takes as input all directories containing separate
# runs of miRDeep2 and merges all output files into
# one file
# -------------
use strict;
use warnings;
use FindBin qw($Bin);
use Getopt::Std;
my %args;
getopts('p:', \%args);

my $filepat = $args{p};
my $header=0;

foreach(@ARGV){
    if( -d $_){
        my $resultcsv = "${_}/${_}2.csv";
        my $samplename = $_;
        if ($filepat){
            ($samplename) = $_ =~ /$filepat/;
            die "Problem extracting sample name from directory name using $filepat on $_" if !$samplename;
        }
        open(my $thiscsv, $resultcsv) or die "Can't open result csv $resultcsv\n";
        my $thisheader = <$thiscsv>;
        if(!$header){
            print "sample,",$thisheader;
            $header = 1;
        }
        while(my $line = <$thiscsv>){
            print $samplename.',',$line;
        }
    }
}
