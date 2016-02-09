#!/usr/bin/perl
# -------------
# file: preprocess2.pl
# author: Matthew Beckers
# date: 03/04/2013
# -------------
#
#
# -------------
use strict;
use warnings;
use Bio::SeqIO;
use FindBin qw($Bin);
use v5.10;
use Getopt::Std;

my $help = "
USAGE preprocess2.pl -a [3' adapter] -o [outputdir] [otheropts] fasta files
The script handles fasta files in batch
-c Use config file with key/values like key=value, the keys are shown in [] below
-C Print a template config file to working directory to fill out
-a [adapter] 3' adapter. Adding 'N' to the adapter will be used as a wildcard - useful for multiplexed adapters
-b [adapter] 5' adapter. Adding 'N' to the adapter will be used as a wildcard - useful for multiplexed adapters
-L [low]  min length
-H [high] max length
-d [dist] distributions to file (not specified)
-s [statsout] other stats to file specified
-o [outdir] output dir
";

my $commandstr = join(" ",@ARGV);
my %args;
getopts("Qc:a:b:L:H:dgs:o:C", \%args);

if($args{C}){
    my $fname = "preprocess.config";
    open(my $config, ">", $fname) or die "Couldn't create file $fname\n";
    print $config "adapter=
low=16
high=30
dist=
statsout=
outdir=
";
    print STDERR "Config file printed!\n";
    exit;
}

# If config file, use the parameters in that file
if($args{c}){
    open(my $infh, '<', $args{c}) or die "Can't open $args{c}\n";
    my %config;
    while(<$infh>){
        if (/^(\w+)=([\w\d\.\/]+)$/){
            $config{$1} = $2;
        }
        else{
            die "Bad config file line $_\n";
        }
    }
    @args{qw(Q a b L H d s g o)} = @config{qw(fastq adapter adapter5 low high dist statsout graph outdir)};
}

my $fastq = $args{Q};
my $adapter = $args{a} || die $help;
my $badapter = $args{b};
my $low = $args{L} || 16;
my $high = $args{H} || 35;
my $dist = $args{d}; #print out distributions to files
my $statsout = $args{s}; #print out other stats to one file specified
my $graph = $args{g};
my $outdir = $args{o} || "./";

# Change regex to fastq form if file is in fastq
my $inputformat = "fasta";
if ($fastq){
    $inputformat = "fastq";
    print STDERR "Input is fastq format";
}

# bail if there are no input files
die "No input files\n".$help if !@ARGV;

# append a / to the end of dir path if it doesn't have one yet
$outdir .= '/' if $outdir !~ /\/$/;

my $statsfh;
if($statsout){
    open ($statsfh, '>', $statsout) or die "Can't open stats out\n";
    print $statsfh join (",", qw/file total.reads no.adapter invalid low.complexity small large output unique/), "\n";
}
print STDERR join (",", qw/file total.reads no.adapter invalid low.complexity small large output/), "\n";

print STDERR "Running $commandstr\n";
# foreach dataset, remove adapters, filter for quality, and write to its own output file
foreach my $infile ( @ARGV ){
    my ($basefile) = $infile =~ m{([^/]+)\.\w+};
    my $outfile = "$outdir${basefile}_na.fasta";
    my $outdist = "$outdir${basefile}_lengthdist.tab";
    my $filein = Bio::SeqIO->new(-file => $infile, -format => $inputformat);

    print STDERR "Processing $infile to $outfile\n";

    # initialise stats
    my $totalreads = 0;
    my $trimmed;
    my $noadapter = 0;
    my $toohigh = 0;
    my $toolow = 0;
    my $invalid = 0; #contains Ns
    my $uncomplex = 0; #Only contains 1 or 2 different bases

    my %nrseqs;
    my %rseqlengths;
    my %nrseqlengths;

    while( my $seqObj = $filein->next_seq() ){
        my $read = $seqObj->seq(); #retrieves the line below the header.

        $totalreads++;

        # Trimming off the adapter
        my $seq = trimadapt($read, $adapter);
        $seq = trimadapt($seq, $badapter, 5) if $badapter;
        
        # Makes all the checks for an invalid sequence
        given($seq){
            when (!$seq)                      {$noadapter++}; # adapter wasn't found
            when (/N/)                        {$invalid++};   # error in sequence     
            when (/A/ + /T/ + /C/ + /G/ <= 2) {$uncomplex++}; # low complexity
            when (length $_ > $high)          {$toohigh++};   # too large 
            when (length $_ < $low)           {$toolow++};    # too small
            default {
                # default is a valid sequence
                $trimmed++;   
                 # add sequence to output seqs
                $nrseqs{$_}++;

                 # length distributions
                $rseqlengths{length($_)}++;
                $nrseqlengths{length($_)}++ if $nrseqs{$_} == 1;
            }
        }
    }
            
    # write the valid trimmed sequences to file in non-redudant form
    open(my $outfh, '>', $outfile) or die "Can't open $outfile\n";
    foreach my $seq ( keys %nrseqs ){
        print $outfh ">$seq($nrseqs{$seq})\n$seq\n"; 
    }

    #if required, write length distribution
    if($dist){
        open (my $distfh, '>', $outdist) or die "Can't open $outdist\n";
        print $distfh "length\tredundant\tunique\n";
        foreach my $len ( sort keys %rseqlengths ){
            print $distfh "$len\t$rseqlengths{$len}\t$nrseqlengths{$len}\n";
        }
    }

    if($graph){
        my $outfile = "$outdir${basefile}_lengthdist.pdf";
        printlengthdist ($outdist, $outfile, $basefile);
    }

    my $unique_seqs_out = keys %nrseqs;
    # if required, write to the stats file
    print $statsfh join(",", $infile, $totalreads, $noadapter, $invalid, $uncomplex, $toolow, $toohigh, $trimmed, $unique_seqs_out),"\n" if $statsfh;
    print STDERR join(",", $infile, $totalreads, atabp($totalreads, $noadapter, $invalid, $uncomplex, $toolow, $toohigh, $trimmed)),"\n";

}

# Given a sequence and adapter, remove the adapter and return
# the mature sequence
sub trimadapt{
    my ($sequence, $adapter, $prime) = @_;
    $prime = 3 if !$prime;
    $adapter =~ s/N/\\w/g; #Ns in the adapter count as wildcards
    my $mature;
    if($prime == 3){
        ($mature) = $sequence =~ /(\w+)$adapter/;
    }
    else{
         ($mature) = $sequence =~ /$adapter(\w+)/;
    }
    return $mature;
}

sub printlengthdist{
    my ($distfile, $outfile, $title) = @_;
    if($args{g}){
        my $cmd = "source(\"${Bin}/ggplot_lengthdist.R\");";
        $cmd .= "pdfLengthdist(\"$distfile\",\"$outfile\",\"$title\");";
        my $exit = `R -e '$cmd' 1>&2`;
        print "Error from R: $exit\n" if $?;
    }
}

# Return percentage arg1 of arg2
sub p{
    my ($v, $t) = @_;
    return sprintf("%.1f", $v/$t*100);
}

# return "value,percentage"
sub tabf{
    my $p = p(@_);
    return $_[0].','.$p;
}

sub atabf{
    my $total = shift;
    return map {tabf($_, $total)} @_;
}

sub atabp{
    my $total = shift;
    return map {p($_, $total)} @_;
}
