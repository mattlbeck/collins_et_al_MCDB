#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 
# -------------
# Accepts non-redundant fasta files and melts them
# to tables that can be easily read into R.
# Outputs a Read count table with read rows and sample columns
# and two length distribution tables
#
#
# -------------
use strict;
use warnings;
use FindBin qw($Bin);
use lib "./";
use Getopt::Std;
use Util;
my $help = '
-r supply a capturing regex to act on filenames OR a list of names you want to give the libraries (in the form name1,name2);
-t test the capturing regex only
-R custom header regex to retrieve counts using. The default is >(\w+)\((\d+)\). Currently only works for fasta files
-o write read table to specified output
-D write redudant length dist to specified output
-d write non-redundant length dist to specified output
-s the seperator to use when writing output. Default is ,
-F fastx style fasta formating (>id-count\nread)
';
my %args;
getopts("tr:FR:o:d:D:s:", \%args);
my $outfile = $args{o} || die "No output file specified using -o\n$help";

my @files = @ARGV;
die $help if !@files;

my @ids;
my $sampleregex;

my $sep = $args{s} || ',';
my $fasta_regex = $args{R} || '>(\w+)\((\d+)\)';

# determine whether r argument was regex or list
my $islist = 0;
if ($args{r} =~ /\(/ || !$args{r}){
    $sampleregex = $args{r} || "(.*)";
}
else{
    $islist = 1;
    @ids = split /,/, $args{r};
    die "Number of names does not match number of files$help" if @ids != @files;
}

my %reads;
my %rlengths;
my %nrlengths;
foreach (0..$#files){
    my $file = $files[$_];

    # extract name for file if regex argument was given
    my $id;
    if ($sampleregex){
        my @idbits = $file =~ /$sampleregex/;
        $id = join("",@idbits);
        push @ids, $id;
    }
    else{
       $id = $ids[$_];     
    }

    # Assess the type of file
    my ($type) = $file =~ m{\.([^/.]+)$};

    print STDERR "Found file $file for $id of type $type\n";

    next if $args{t};
    open (my $in, '<', $file) or die "Can't open $file\n";
    while (my $line = <$in>){
        my %unique;
        my ($seq, $count);
        if ($type =~ /fasta|fa/){
            if($args{F})
            {
                ($count) = $line =~ /\d+-(\d+)/;
                $seq = <$in>;
                chomp $seq;
            }
            else
            {
                if ($line =~ /$fasta_regex/){
                    ($seq, $count) = ($1, $2);
                }
            }
        }
        elsif ($type =~ /pat|patman/){
            my ($hit, $read) = split /\t/, $line;
            next if $unique{$read};
            $unique{$read} = 1;
            ($seq, $count) = $read =~ /(\w+)\((\d+)\)/;
        }
        else{ die "Bad input file $file\n" }
        if($seq){
            $reads{$seq}->{$id} = $count;
            $rlengths{length($seq)}->{$id} += $count;
            $nrlengths{length($seq)}->{$id}++;
        }
    }
}

# output read table 
write_table (\%reads, $outfile, \@ids, 'read'); 

# output redundant length dist
write_table (\%rlengths, $args{D}, \@ids, 'length', sub{$a <=> $b}) if $args{D};

# output non-redundant length dist
write_table (\%nrlengths, $args{d}, \@ids, 'length', sub{$a <=> $b}) if $args{d};

#---
# writes tables to specified file from the hash formatted as {rows} -> {columns}
# columns are written out by sorting in numerical order
#---
sub write_table{
    my ($table, $filename, $ids, $rowheader, $sort) = @_;
     
    my @sortedids;
    if(!$islist){
        @sortedids = Util::sortseqids( @$ids );
    }
    else{
        @sortedids = @ids;
    }
    open (my $out, '>', $filename) or die "Can't open $filename for writing table\n";
    print $out join($sep, $rowheader, @sortedids),"\n";
    my @rows = ($sort ? sort $sort keys %$table : keys %$table);
    foreach my $row (@rows){
        print $out join ($sep, $row, map { $_ || 0 } @{ $table->{$row} }{ @sortedids }),"\n"; 
    }
}
