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
use FindBin qw($Bin);
use Getopt::Std;
use Bio::DB::Fasta;
my %args;
getopts('t:g:e:f:s:', \%args);

my $help = "Outputs a table and optional fasta file of mature tRNAs found by tRNA scan
-t trnascan output
-g genome file
-e number of base pairs to expand each tRNA sequence by. Default: 0
-f output to specified file a fasta file of mature sequences
-s structure file
";

my $trnas = $args{t} or die $help;
my $genomefile = $args{g} or die $help;
my $expand = $args{e} || 0;
my $structurefile = $args{s};

my $fasta;
if($args{f}){
    open($fasta, '>', $args{f}) or die "Can't open ".$args{f}."\n";
}

my $genome = Bio::DB::Fasta->new($genomefile);
open(my $tfh, '<', $trnas) or die "Can't open $trnas\n";
<$tfh> for(1..3);

my %trnas;
while(<$tfh>){
    my @cols = split /\s+/;
    my ($seqid, $num, $begin, $end, $type, $codon, $is, $ie, $cove) = @cols[0,1,2,3,4,5,6,7,8];
    $seqid =~ s/\s+$//;
    print $seqid,"|\n";
    
    my ($first, $last) = ($begin, $end)[$end < $begin, $end > $begin];
    my ($iposs, $ipose) = (0,0);
    
    my $chr = $genome->get_Seq_by_acc($seqid);
    my $seq = $chr->subseq($first, $last);
    my $strand = '+';
    if($end < $begin){
        $strand = '-';
        $seq =~ tr/ATCGatcg/TAGCtagc/;
        $seq = reverse($seq);
    }
    if ($is && $ie){
        my ($istart, $iend) = ($is - $first + 1, $ie-$first + 1);
        ($iposs, $ipose) = ($istart, $iend)[$iend < $istart, $iend > $istart];
        #print STDERR $seq,"\n";
        #print STDERR "$first $last $iposs $ipose $is $ie\n";
        
        #print STDERR "before $seq\n";
        
        substr($seq, $iposs-1, $ipose-$iposs+1) = "";
    }
    $seq .= 'CCA';

    # print join(",",$seqid.'_'.$num,$first,$last,$type,$codon,$seq),"\n";
    $trnas{$seqid.'_'.$num} = {seqid => $seqid, start => $first, end => $last, istart => $iposs, iend => $ipose, 
        type => $type, codon => $codon, seq => $seq, strand => $strand, cove => $cove}; 
    if($args{f}){
        print $fasta ">${seqid}_${num}\n$seq\n";
    }
}

if ($structurefile){
    open(my $sfh, '<', $structurefile) or die "Can't open $structurefile\n";
    my $currentid;
    my $current_trna;
    my ($is, $ie);
    while(<$sfh>){
        chomp;
        if(/(.+?)\s+\(\d+-\d+\)\s+Length/){
            # the id line            
            my $id = $1;
            $id =~ s/\.trna/_/;
            $currentid = $id;
            $current_trna = $trnas{$currentid};
        }
        elsif(/^Seq:\s+(\w+)/){
            # the sequence line
            my $seq = $1;

            my ($is, $ie) = @{$current_trna}{qw(istart iend)};
            if ($is || $ie){
                substr($seq, $is-1, $ie-$is+1) = "";
            }
            # currently do nothing with this sequence
            # print STDERR "$currentid $seq\n","$currentid ",$current_trna->{seq}."\n\n";
        }
        elsif(/^Possible intron:\s(\d+)-(\d+)/){
            ($is, $ie) = ($1, $2);
        }
        elsif(/^Str:\s+([><.]+)/){
            # the structure line
            my $struct = $1;
            #my ($is, $ie) = @{$current_trna}{qw(istart iend)}; 
            if ($is || $ie){
                #print STDERR $is,$ie,"\n$struct";
                substr($struct, $is-1, ($ie-$is)+1) = "";
            }
            $struct .= '...';
            $trnas{$currentid}->{struct} = $struct;
            ($is, $ie)=(0,0);
        }
    
    }
}

print join(",", qw(id seqid start end strand intron.start intron.end type codon cove.score mature.seq mature.structure)),"\n";
foreach(keys %trnas){
    my $trna = $trnas{$_}; 
    print join(",", $_,@{$trna}{qw(seqid start end strand istart iend type codon cove seq struct)}),"\n";
}
