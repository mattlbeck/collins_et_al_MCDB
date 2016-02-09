#!/usr/bin/perl
# -------------
# file:
# author: Matthew Beckers
# date: 25/11/2015
# -------------
# Input: predictions object from merge_miRNAs script, read counts for all libraries used by prediction tools
# Output: Specified fuzzyregions for probable five prime and three prime mature miRNAs in each hairpin. These can
# be used as alignment or overlapping targets to annotate reads
#
# -------------
use strict;
use warnings;
use FindBin qw($Bin);
use Data::Dumper;
use Getopt::Std;
use Storable;
use POSIX;
use Bio::DB::Fasta;
#$SIG{__WARN__} = sub { $DB::single = 1 }; # debugger breaks on warning
my %args;
getopts('p:c:g:f:', \%args);

my $predictions_file = $args{p} || die "Please run merge_miRNAs.pl script and then specify the predictions.obj file as -p here\n";
my $count_file = $args{c} || die "Please specify a csv count matrix with a 'read' column for all libs put through the predictions";
my $genomefile = $args{g} || die "Please specify genome using -g\n";
my $fastafile = $args{f};
my $genome = Bio::DB::Fasta->new($genomefile);

my $predictions = retrieve $predictions_file;


my %counts;
open (my $countfh, '<', $count_file) or die "Can't open $count_file\n";
my $header = <$countfh>;
chomp $header;
my @libnames = split /,/, $header;
shift @libnames;

my $fout;
if($fastafile){
    open($fout, ">", $fastafile) || die "can't open $fastafile";
}

while(<$countfh>)
{
    chomp;
    my @cols = split /,/;
    my $read = shift @cols;

    # add value to hash only if value is not 0
    for my $i (0..$#cols){
        if($cols[$i] > 0){
           $counts{$read}->{$libnames[$i]} = $cols[$i]; 
        }
    }
}

print join(",", qw(mirname seqid start end strand source arm region.start region.end region.seq size top.mature precursor.seq)),"\n";

# for each prediction
#   split hairpin in middle
#   find coverage pattern
#   identify largest steps to create regions
foreach my $prediction (values %$predictions){
    my ($start, $end) = @{$prediction}{qw(start end)};
    my $length = ($end - $start) + 1; # precursor length
    my $middle = ceil($length/2); # mid point - if even, base just downstream of mid. If odd, middle base
    my $mircat = $prediction->{mircat};
    my $mapmi = $prediction->{mapmi};
    my $mirdeep = $prediction->{mirdeep};

    my $seqid_seq = $genome->get_Seq_by_id($prediction->{seqid});
    my $strandarg = strandCode($prediction->{strand});
    my $pre = uc($genome->seq($prediction->{seqid}, $start, $end, $strandarg));

    my $mirname = $prediction->{mirname};

    my (@fivec, @threec);

    my $threelength = $length - $middle;

    # initialise with zeroes. One full length precursor for each strand
    @fivec[0..($length-1)] = ("0") x $length;
    my $isfive = 0;
    @threec[0..($length-1)] = ("0") x $length;
    my $isthree = 0;
    
    my %reads;
    my (%topfive, %topthree); # keep track of reads from mircat
    $topfive{c} = 0;
    $topthree{c} = 0;
    if($mircat){
        for my $p (@{$mircat->{miRNAs}}){
            my $fields = $p->{fields};
            my $read = $fields->{"mature.seq"};
            my $mcounts = $counts{$read}; # get counts for this mature

            # start site with reference to precursor 5' to 3'
            my $mstart =  $fields->{start} - $start;
            my $mend = $fields->{end} - $start;

            # if - strand
            if($fields->{strand} eq "-")
            {
                ($mstart, $mend) = ($length-1-$mend, $length-1-$mstart);
            }

            my $mid = $mstart + ceil((($mend - $mstart)+1)/2);

            my $cover = 0;
            for my $count (values %$mcounts){
                $cover += $count;
            }

            my $strand = $fields->{strand};

            if($cover){
                $reads{$read} = $cover; # mircat has used this read
                if($mid < $middle){
                    $_ += $cover for @fivec[($mstart)..($mend)];
                    $isfive=1;
                    @topfive{qw(s c)} = ($read,$cover) if !$topfive{c} || $topfive{c} < $cover;
                }
                else{
                    $_ += $cover for @threec[($mstart)..($mend)];
                    $isthree=1;
                    @topthree{qw(s c)} = ($read,$cover) if !$topthree{c} || $topthree{c} < $cover;
                }
            }
        }
    }
    if($mirdeep){
        for my $p (@{$mirdeep->{miRNAs}}){
            my $read = $p->{mature};

            # for some reason the mature mirdeep sequence is automatically flipped
            #  it is always the + strand. Flip it if minus
#            $read = flipstrand($read) if $p->{strand} eq "-";
            my $mstart = index(uc($pre), uc($read));
            my $mend = $mstart + (length($read)-1);

            my $mid = $mstart + ceil((($mend - $mstart)+1)/2);
            next if $reads{$read};
            
            my $mcounts = $counts{$read};
            
            my $cover = 0;
            for my $count (values %$mcounts){
                $cover += $count;
            }

            next if !$cover;
            $reads{$read} = $cover;

            my $strand = $p->{strand};
            if($mid < $middle){
                $_ += $cover for @fivec[($mstart)..($mend)];
                $isfive=1;
                @topfive{qw(s c)} = ($read,$cover) if !$topfive{c} || $topfive{c} < $cover;
            }
            else{
                $_ += $cover for @threec[($mstart)..($mend)];
                $isthree=1;
                @topthree{qw(s c)} = ($read,$cover) if !$topthree{c} || $topthree{c} < $cover;
            }
            
        }
    }
    if($mapmi){
        for my $p (@{$mapmi})
        {
            my $read = $p->{sequence};
            my $mstart = index(uc($pre), uc($read));
            my $mend = $mstart + (length($read) - 1);
            my $mid = $mstart + ceil((($mend - $mstart)+1)/2);

            next if $mstart == -1; # we are not interested if mapmi has found a precursor substantially longer than what we are estimating
            next if $reads{$read}; # we have already done this sequence

            my $mcounts = $counts{$read};
            next if !$mcounts; # not interested if there are no counts
            my $cover = 0;
            for my $count (values %$mcounts){
                $cover += $count;
            }
            next if !$cover;

            my $strand = $p->{strand};
            if($mid < $middle){
                $_ += $cover for @fivec[($mstart)..($mend)];
                $isfive=1;
                @topfive{qw(s c)} = ($read,$cover) if !$topfive{c} || $topfive{c} < $cover;
            }
            else{
                $_ += $cover for @threec[($mstart)..($mend)];
                $isthree=1;
                @topthree{qw(s c)} = ($read,$cover) if !$topthree{c} || $topthree{c} < $cover;
            }
            
            
        }
    }
    my ($fs, $ts, $fe, $te) = (-1,-1,-1,-1);
    my $line = join(",",@{$prediction}{qw(mirname seqid start end strand source)}).",";
    my $strand = $prediction->{strand};

    if($isfive){
        my $fiveregion = findRegion(\@fivec);
        ($fs, $fe) = @$fiveregion;
        ($fs, $fe) = ($length-1-$fe, $length-1-$fs) if $strand eq "-";
        ($fs, $fe) = ($start + $fs, $start + $fe);
        my $mseq = $genome->seq($prediction->{seqid}, $fs, $fe, strandCode($prediction->{strand}));
        print $line , join(",", "5_prime", $fs, $fe, uc($mseq), ($fe-$fs)+1, $topfive{s}, uc($pre)), "\n";
        print $fout ">",$mirname,"-5p", "\n", $mseq,"\n" if $fout;
    }
    if($isthree){
        my $threeregion = findRegion(\@threec);
        ($ts, $te) = @$threeregion;
        ($ts, $te) = ($length-1-$te, $length-1-$ts) if $strand eq "-";
        ($ts, $te) = ($start + $ts, $start + $te);
        my $mseq = $genome->seq($prediction->{seqid}, $ts, $te, strandCode($prediction->{strand}));
        print $line , join(",", "3_prime", $ts, $te, uc($mseq), ($te-$ts) + 1, $topthree{s}, uc($pre)), "\n";
        print $fout ">",$mirname,"-3p", "\n", $mseq,"\n" if $fout;
    }
    my $fline;
    my $tline;
    my $mfline;
    my $mtline;
    my $preline;
    for (0..$#fivec)
    {
        my ($flength, $tlength) = (length($fivec[$_]), length($threec[$_]));
        my $maxwidth = ($tlength , $flength)[$flength > $tlength];
        $fline .= $fivec[$_]. (" " x ($maxwidth-$flength+1));
        $tline .= $threec[$_]. (" " x ($maxwidth-$tlength+1));
        $preline .= substr($pre, $_, 1) . (" " x ($maxwidth));
        my $mf = ($start + $_ >= $fs && $start + $_ <= $fe) ? "v" : " ";
        my $mt = ($start + $_ >= $ts && $start + $_ <= $te) ? "^" : " ";
        $mfline .= (" " x ($maxwidth-1)) . $mf ." ";
        $mtline .= (" " x ($maxwidth-1)) . $mt ." ";
    }
    print STDERR $prediction->{mirname}, $prediction->{strand}, "\n", $mfline, "\n";
    print STDERR $fline, "\n";
    print STDERR $preline, "\n";
    print STDERR $tline, "\n";
    print STDERR $mtline, "\n\n";
    if(!$isfive && !$isthree && $mapmi){
        # use the mapmi information
        print $line, join(",", @{$mapmi->[0]}{qw(prime start end sequence)}, length($mapmi->[0]->{sequence}), @{$mapmi->[0]}{qw(sequence pre.seq)}), "\n";
    }
}

sub findRegion{
    my ($cover) = @_;
    # find largest step up and down
    my $biggestUp=$cover->[0];
    my $rstart = 0;
    my $biggestDown=-$cover->[-1];
    my $rend = $#{$cover};
    #print Dumper($cover);
    for my $b (1..$#{$cover})
    {
       my $diff = $cover->[$b] - $cover->[$b-1]; 
       if($diff > 0 && $diff > $biggestUp)
       {
            $biggestUp = $diff;
            $rstart = $b;
       }
       elsif($diff < 0 && $diff <= $biggestDown)
       {
           $biggestDown = $diff;
           $rend = ($b-1);
       }
    }
    return [$rstart, $rend];
}

# find middle of a seq
sub middle
{
    my ($start, $end) = @_;

}
sub strandCode
{
    return ((shift eq "-") ? -1 : 1);
}

sub flipstrand{
    my $seq = shift;
    $seq =~ tr/TCGA/AGCT/;
    return reverse($seq);
}
