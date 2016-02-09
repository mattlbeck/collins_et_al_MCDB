#!/usr/bin/perl
# -------------
# file: melt_mircat_csvs.pl
# author: Matthew Beckers
# date: 
# -------------
# Takes as input a list of directories representing mircat outputs from different runs.
# Reads in both the results csv and the hairpin file
# "melts" the results into two tables:
# * A table of unique predictions with all prediction information including hairpins
# * A table of mature reads and the runs that predicted them
# Both in csv form
#
# -------------
use strict;
use warnings;# FATAL => 'all';
use FindBin qw($Bin);
use Getopt::Std;
use Data::Dumper;
use Storable;
use Bio::DB::Fasta;
use Sort::Naturally;

my $genome = Bio::DB::Fasta->new("../genome/Bter20110317-genome-formatted.fa");
my %args;
my $help = " perl melt_mircat_csvs.pl [args] dirlist
-r filename pattern regex
-m modified mapmi file
-M mapmi output file
-D mirdeep merged csv file
-b bad hairpins file
";
getopts('r:m:M:b:D:', \%args);
my $mircat_csv_name = 'mircat_results.csv';
my $hairpin_file_name = 'miRNA_hairpins.txt';
my $filenamepat = $args{r} || '([^./]+)/?'; 
my $mapmifile = $args{m};
my $mapmioutfile = $args{M};
my $mirdeepfile = $args{D};


my $summary_miRNAs = "summary_predictions.csv";

my $mc_unique_out = "mircat_miRNAs.csv";
my $mc_predictions_out = "mircat_predictions.csv";
my @mc_colnames = qw(seqid start end strand abundance mature.seq mature.length genomic.hits hairpin.length hairpin.gc mfe amfe randfold star);
my @mc_unique_fields = qw(seqid start end strand mature.seq 
                       mature.length genomic.hits mfe amfe randfold star hairpin.length hairpin.gc);
my @mc_prediction_fields = qw(mature.seq abundance seqid start end strand);
my $stars_out = "mircat_stars.csv";
my $badhpfile = "mircat_badhps.csv";

my $md_unique_out = "mirdeep_miRNAs.csv";
my $md_predictions_out = "mirdeep_predictions.csv";
my @md_colnames = qw(run id score total.count mature.count loop.count star.count randfold.pval 
                     mature star precursor seqid start end strand);
my @md_unique_fields = qw(seqid start end strand mature star precursor);
my @md_prediction_fields = qw(run seqid start end strand score mature mature.count star star.count precursor loop.count randfold.pval);

my $mm_out = "mapmi_miRNAs.csv";
my $naming_conflicts_out = "miRNA_naming_conflicts.txt";

die $help if !@ARGV;


my $first = 0;
my @unique_fields = qw(mirid uid arm name seqid start end orientation sequence size genome.hits hp.GC MFE Adjusted.MFE Rfold.pval top.star hp.size hairpin hp.start hp.end structure);
my @prediction_fields = qw(uid name run mirid sequence abundance seqid start end orientation);
my @mapmi_header = qw(sequence mismatch seqid strand start end size pre.start pre.end score pre.seq id name );

my @bad_hp; # array fro hairpin predictions that do not fully contain the mature sequence.

my %star_predictions;
my @star_fields = qw(run mir.id star abundance);
my $id = 0;
my $uid_counter = 0;

# merge the mircat results into one hash
my $unique_mircat = merge_mircat_results(@ARGV);

# check for hairpin overlaps of unique sequences
$unique_mircat = overlap_mircat_hairpins($unique_mircat);

my $unique_mirdeep;
if($mirdeepfile){
    $unique_mirdeep = merge_mirdeep_results($mirdeepfile);
}

# join together mirdeep and mircat results
my $unique_mirs = merge_mirdeep_mircat($unique_mircat, $unique_mirdeep);

# If using mapmi file, corectly name overlapping mature sequences and any
# related sequence using mirbase names.
my %mapmi;
if($mapmifile){
    $unique_mirs = merge_mapmi($mapmifile, $unique_mirs);
}

print "Printing everything\n";

#open (my $unique_fh, '>', $unique_out) or die "Can't write to $unique_out";

# mircat predictions are individual mircat results per sample used
open (my $mcp_fh, '>', $mc_predictions_out) or die "Can't write to ";
my @mcp_out_cols = (qw(uid mirname arm run), @mc_prediction_fields);
print $mcp_fh join(",", @mcp_out_cols),"\n";

open (my $mdp_fh, '>', $md_predictions_out) or die "Can't write to ";
my @mdp_out_cols = (qw(uid mirname arm), @md_prediction_fields);
print $mdp_fh join(",", @mdp_out_cols),"\n";

open (my $mcm_fh, '>', $mc_unique_out) or die "Can't write to ";
my @mcm_out_cols = (qw(uid mirname arm), @mc_unique_fields, qw(precursor structure));
print $mcm_fh join(",", @mcm_out_cols),"\n";

open (my $mdm_fh, '>', $md_unique_out) or die "Can't write to ";
my @mdm_out_cols = (qw(uid mirname arm), @md_unique_fields );
print $mdm_fh join(",", @mdm_out_cols),"\n";

open (my $out_fh, '>', $summary_miRNAs) or die "Can't write to ";
my @out_cols = qw(uid mirname seqid start end strand source);
print $out_fh join(",", @out_cols),"\n";

open (my $mm_fh, '>', $mm_out) or die "Can't write to $mm_out";
print $mm_fh join(",", qw(uid mirname species), @mapmi_header),"\n";

my $predictions_obj_file = "predictions.obj";
store $unique_mirs, $predictions_obj_file;

foreach my $uid(keys %$unique_mirs){
    my $prediction = $unique_mirs->{$uid};
    my $mirname = $prediction->{mirname} || "NA";
    # top level miRNA summary
    print $out_fh join(",", $uid, $mirname, @{$prediction}{@out_cols[2..6]}),"\n";

    # mircat printing
    my $mircat = $prediction->{mircat};
    foreach my $miRNA (@{$mircat->{miRNAs}}){
        # miRNA is a hash of this mir (seqid,start,end) plus its predictions
        # print the unique mature miRNA prediction
        print $mcm_fh join(",", $uid, $mirname, $miRNA->{arm}, 
            @{$miRNA->{fields}}{@mc_unique_fields}, @{$miRNA}{qw(structure hairpin)}),"\n";

        # print the various predictions for each run
        foreach my $pred (@{$miRNA->{predictions}}){
            print $mcp_fh join(",", $uid, $mirname, $miRNA->{arm}, $pred->{run}, 
                @{$pred}{@mc_prediction_fields}),"\n";
        }
    }

    # mirdeep printing
    my $mirdeep = $prediction->{mirdeep};
    #$DB::single=1 if defined $mirdeep;
    foreach my $miRNA (@{$mirdeep->{miRNAs}}){
        # print the unique mature miRNA prediction
        print $mdm_fh join(",", $uid, $mirname, $miRNA->{arm}, @{$miRNA}{@md_unique_fields}),"\n"; 
        
        # print the various predictions for each run
        foreach my $pred (@{$miRNA->{predictions}}){
            print $mdp_fh join(",", $uid, $mirname, $miRNA->{arm}, 
                @{$pred}{@md_prediction_fields}),"\n";
        }
    }

    # mapmi printing
    my $mapmi = $prediction->{mapmi};
    foreach my $pred (@$mapmi){
        # print a line for each species the miRNA is conserved from
        foreach my $sp (@{$pred->{sp}}){
            print $mm_fh join(",", $uid, $mirname, $sp, @{$pred}{@mapmi_header}),"\n";
        }
    }
}

if (@bad_hp){
    print STDERR "There were predicted hairpins that did not properly contain the sequence\n";
    my $badfh = \*STDERR;

    if ($badhpfile){
        open($badfh, '>', $badhpfile) or die "Could not open $badhpfile\n";
    }

    foreach my $bhp (@bad_hp){
        print $badfh join(',', @$bhp),"\n";
    }
}

# Args: array of mircat result directories
# Returns: hash of unique mirs
#         {id}->{uid, hpstats, predictions}
sub merge_mircat_results{
    my %unique_mirs; 
    my $mir60013p = "Group13.5:4555800:4555820:-";
    foreach(@_){
        my ($run_name) = $_ =~ m{$filenamepat};
        die "Can't deduce dir name from $_" if !$run_name;

        print "Processing miRCat run $run_name\n";

        open(my $csv, '<', $_."/$mircat_csv_name") or die "Can't open $_\n";
        my $head = <$csv>; # do nothing with the head

        # example of a bad hp for debugging
        my $debugstr = "Group8.1-2737946-2737968";
        my $dbadded = 0;

        # compile unique mir hash and list of predicitons
        while (my $line = <$csv>){
            chomp $line;
            my @fields = split /,/, $line;

            my %lineh;

            # If a randfold was not carried out
            if (@fields == 13){
                my $stars = $fields[-1];
                $fields[-1] = "N/A";
                push @fields, $stars;
            }

            @lineh{@mc_colnames} = @fields;
            print STDERR $lineh{'mature.seq'}, ", ", $run_name, "\n" if $lineh{start} =~ /4555/;

            my $uid = join(":",@lineh{qw(seqid start end strand)});
            #print STDERR $run_name, " ",$uid,"\n";
            $dbadded = 1 if $uid eq $mir60013p;

            #$DB::single=1 if !$unique_mirs{$mir60013p} && $dbadded;
            
            # find top star prediction
            my %stars;
            unless ($lineh{star} eq 'NO'){
                %stars = split /[()\s]+/, $lineh{star};
                my @sortedstars = sort {$stars{$b} <=> $stars{$a}} keys %stars;
               $lineh{star} = $sortedstars[0]; 
            }

            my %ph;
            @ph{@mc_prediction_fields} = @lineh{@mc_prediction_fields};

            # mir already in unique hash? Use its id
            if($unique_mirs{$uid}){
                my $thisid = $unique_mirs{$uid}->{id};
                push @{$unique_mirs{$uid}->{predictions}}, {run => $run_name, 
                                                            id => $thisid, %ph};
            }
            else{
                # headers: Chromosome,Start,End,Orientation,Abundance,Sequence,sRNA length,# Genomic Hits,Hairpin Length,Hairpin % G/C content,Minimum Free Energy,Adjusted MFE,Randfold p-value,miRNA*
                #              0       1     2       3         4         5        6           6                 7                    8                 8                 9         10               11
                my %uh;
                @uh{@mc_unique_fields} = @lineh{@mc_unique_fields};
                $unique_mirs{$uid} = {id=>$id, 
                    fields => \%uh};

                push @{$unique_mirs{$uid}->{predictions}}, {run => $run_name, id => $id, %ph};
                $id++;
            }

            if (%stars){
                my $sid = "$run_name-$id"; 
                $star_predictions{$sid} = \%stars; 
            }
        }

        # Append hairpin structures to table
        open (my $hp, '<', $_."/$hairpin_file_name") or die "Can't open hairpin file $_";
        while (my $line = <$hp>){
            chomp $line;
            # header lines contain mature seq, seqid, start and end of mature seq
            if ($line =~ m!^>([ATCGU]+)_(.+)/(\d+)-(\d+)!){
                my ($mature, $seqid, $start, $end) = ($1,$2,$3,$4);
                my $hairpin = <$hp>; # next line is hairpin sequence
                chomp $hairpin;
                my $struct = <$hp>; # next line is hairpin structure
                chomp $struct;
                if(!$struct)
                {
                    print STDERR "No structure found for $line $run_name. Regenerating using RNAfold\n";
                    my @out = `echo $hairpin | RNAfold`; 
                    ($struct) = $out[1] =~ /([\(.\)]+)\s/;
                }
                chomp $hairpin;

                my $uid = join(":",$seqid,$start,$end);
                my $pos = $unique_mirs{$uid.':+'};
                my $neg = $unique_mirs{$uid.':-'};
                if((!$pos && !$neg) || ($pos && $neg)){
                    if($pos->{fields}->{"mature.seq"} eq $mature){
                        $uid .= ':+';
                    }
                    elsif($neg->{fields}->{"mature.seq"} eq $mature){
                        $uid .= ':-';
                    }
                    else{
                        die "Issue with finding prediction for hairpin $uid";
                    }
                }
                else{
                    $uid .= (($pos) ? ':+' : ':-');
                }
                #$DB::single = 1 if $uid eq 'GroupUn981:279375:279396:+';

                #$DB::single=1 if $uid eq $mir60013p;

                if($unique_mirs{$uid} && !$unique_mirs{$uid}->{hairpin}){
                    my $m_pos = index(uc($hairpin), uc($mature));
                    if ($m_pos < 0){
                        # attempt to find where at least part of the mature sequence is in the hairpin
                        my $mature2 = $mature;
                        my( $mstart, $mend);
                        my $numchomps = 0;
                        while(($mstart = index(uc($hairpin), uc($mature2))) < 0){
                            die "Mature sequence $mature has dropped to just length 2" if length($mature2) <=2;
                            $mature2 = substr($mature2, 1, length($mature2)-2);
                            $numchomps++;
                        }

                        # found a piece of mature sequence in the hairpin
                        # only allow it if its at the very end of the hairpin.
                        # If its in the middle somewhere then the match was totally erroneous.
                        my $left = 0;
                        if( $mstart + length($mature2)  == length($hairpin)){ # mstart is numchomps further from real start but the length is numchomps smaller so they cancel
                            # mend is at end of hairpin 
                        }
                        elsif( $mstart == 0){
                            # mstart is at beginning of hairpin
                            $left = 1;
                        }
                        else{ die "Found $mature whittled to $mature2 in hairpin $hairpin starting at $mstart. This is erroneous.";}

                        # use the start position of the mature sequence to find the hairpin coordinates
                        $mstart -= $numchomps;
                        my $hstart = $unique_mirs{$uid}->{fields}->{start} - $mstart;
                        my $hend = $hstart + length($hairpin) - 1;
                        my ($seqid, $strand) = @{$unique_mirs{$uid}->{fields}}{qw(seqid strand)};
                        $strand = int($strand.'1');

                        my (%hp1, %hp2);
                        @hp2{qw(start end)} = ($hstart - $numchomps, $hend + $numchomps);

                        # get two versions of the modified hairpin
                        if($left){
                            @hp1{qw(start end)} = ($hstart - $numchomps, $hend);
                        }
                        else
                        {
                            @hp1{qw(start end)} = ($hstart, $hend + $numchomps);
                        }

                        # Compute structures with RNAfold
                        %hp1 = seq_rnafold(%hp1, seqid => $seqid, strand => $strand);
                        %hp2 = seq_rnafold(%hp2, seqid => $seqid, strand => $strand);

                        my $newhp = (\%hp2, \%hp1)[$hp1{mfe} < $hp2{mfe}];
                        $hairpin = $newhp->{seq};
                        $struct = $newhp->{struct};

                        
                        # Mature sequence sometimes aligns outside of hairpin. Push these to a seperate table
                        push @bad_hp, [$uid, $seqid, $start, $end, $mature, $hairpin, $newhp->{seq}]; 
                        #delete $unique_mirs{$uid};
                        #next;
                        $m_pos = index(uc($hairpin), uc($mature));
                    }
                    # Need to grab orientation from hash
                    my ($h_startpos, $h_endpos);
                    my $orientation = $unique_mirs{$uid}->{fields}->{strand};
                    if ($orientation eq '-'){
                        my $m_endpos = $m_pos;
                        $h_endpos = $end + $m_endpos;
                        $h_startpos = ($h_endpos - length($hairpin)) + 1;
                    }
                    else{
                        my $m_startpos = $m_pos;
                        $h_startpos = $start - $m_startpos;
                        $h_endpos = ($h_startpos + length($hairpin)) - 1;
                    }

                    # Determine whether this is a 3-prime sequence or a 5-prime sequence with respect to hairpin structure
                    my $prime = findPrime($struct, $mature, $hairpin, $orientation);

                    @{$unique_mirs{$uid}}{qw(hairpin hp.start hp.end structure arm)} = ($hairpin, $h_startpos, $h_endpos, $struct, $prime) ;
                }
            }
        }
    }
    return \%unique_mirs;
}

sub seq_rnafold{
    my %obj = @_;
    $obj{seq} = $genome->seq(@obj{qw(seqid start end strand)}); 
    my @out = `echo $obj{seq} | RNAfold`;
    @obj{qw(struct mfe)} = $out[1] =~ /([\(.\)]+)\s\(([\s+\-\d\.]+)\)/;
    return %obj;
}

# Find the primeness of a mature sequence in a hairpin
sub findPrime{
    my ($struct, $mature, $hairpin, $strand) = @_;
    if($strand eq "-"){
        $struct = reverse $struct;
        $struct =~ tr/)(></()<>/;
        $mature = reverse $mature;
        $hairpin = reverse $hairpin;
    }
    my $m_pos = index(uc($hairpin), uc($mature));
    my $hairpinMiddle;
    my @loops = $struct =~ m!([(<{] [.\-=]+ [)>}])!xg;
    if(@loops > 1 || !$struct){
        # fall back to using center of sequence
        $hairpinMiddle = int(length($hairpin)/2)
    }
    else{
        $hairpinMiddle = int(($-[0] + $-[1]) / 2)
    }

    return "3_prime" if $m_pos > $hairpinMiddle;
    return "5_prime";
}


# args: hash of unique mirs
# return:
sub overlap_mircat_hairpins{
    my $unique_mirs = shift;
    my @ukeys = sort sort_uids keys %$unique_mirs;
    my @uid2mirids;
    my %mirids2uid;

    open (my $dbout, '>', "mircat_overlaps.txt") or die "Can't opent mircat_overlaps.txt";
    print "Assessing hairpin overlaps between miRCat files\n";
    foreach my $a (0..$#ukeys){
        my $ak = $ukeys[$a];
        my $aoverlaps = 0;
        foreach my $b ($a..$#ukeys){
            my $bk = $ukeys[$b];
            #$DB::single=1 if $unique_mirs->{$ak}->{fields}->{'mature.seq'} =~ /^GTAGGTAACGACTGATGGGAAC(A?)/
            #                 && $unique_mirs->{$bk}->{fields}->{'mature.seq'} =~ /^GTAGGTAACGACTGATGGGAAC(A?)/;
            #print Dumper($unique_mirs{$ak});

            # Don't compare self, check seqid
            if($a != $b && $unique_mirs->{$ak}->{fields}->{seqid} eq $unique_mirs->{$bk}->{fields}->{seqid}
                && $unique_mirs->{$ak}->{fields}->{strand} eq $unique_mirs->{$bk}->{fields}->{strand}){
                my($as, $ae) = @{$unique_mirs->{$ak}}{qw(hp.start hp.end)};
                my($bs, $be) = @{$unique_mirs->{$bk}}{qw(hp.start hp.end)};
                
                # check for coordinate overlap
                if($ae > $bs && $as < $be){
                    #$DB::single=1 if $unique_mirs->{$ak}->{fields}->{'mature.seq'} =~ /^GTAGGTAACGACTGATGGGAAC(A?)/;
                    
                    #$DB::single=1 if $unique_mirs->{$bk}->{fields}->{'mature.seq'} =~ /^GTAGGTAACGACTGATGGGAAC(A?)/;
                    my $overlap = ($ae, $be)[$ae > $be] - ($as, $bs)[$as < $bs];

                    # same miRNA if overlap is > maximum hairpin length
                    my ($al, $bl) = ($ae - $as, $be - $bs);
                    my $maxl = ($al,$bl)[$al < $bl];

                    if($overlap > $maxl/2){
                        $aoverlaps = 1;

                        # check $a for previous overlaps
                        my $thisuid = $mirids2uid{$ak};

                        # check $b for previous overlaps
                        my $thisuidb = $mirids2uid{$bk};
                        
                        if(defined $thisuid && defined $thisuidb && $thisuid != $thisuidb){
                            die "uid assignment mismatch between overlapping miRNAs\n", Dumper($unique_mirs->{$ak}), Dumper($unique_mirs->{$bk}),"\n";
                        }

                        # If both already have a defined uid, it must be the same due
                        # to passing the die condition above and we don't need to process
                        # this overlap again (this probably is never the case)
                        if(!defined $thisuidb){ 
                            if (defined $thisuid){
                                # previous overlap found, use the uid
                                push @{$uid2mirids[$thisuid]}, $bk;
                                $mirids2uid{$bk} = $thisuid; 
                                $unique_mirs->{$bk}->{uid} = $thisuid;
                            }
                            else{
                                # no previous overlap found, create new entry
                                push @uid2mirids, [$ak,$bk]; 
                                $mirids2uid{$ak} = $#uid2mirids;
                                $mirids2uid{$bk} = $#uid2mirids;
                                $unique_mirs->{$ak}->{uid} = $#uid2mirids;
                                $unique_mirs->{$bk}->{uid} = $#uid2mirids;
                            }
                        }
                    }
                    else
                    {
                        # Weakly overlapping entries. These seem to be mostly from 3 prime and 5 primes seqs that differ with the hairpin direction
                        # Attempt to reconcile at least these in to one prediction with a 5 prime and 3 prime seq

                        # Print this overlap out for debugging
                        #my $dbkey = "Group1.4:904361:904381:-";
                        #my $bdbkey = "Group1.4:904361:904383:-";
                        #$DB::single=1 if $ak eq $dbkey && $bk eq $bdbkey;
                        printOverlaps($dbout, $unique_mirs, $ak, $bk);
                        
                        # Do both mature sequences overlap one of the hairpins, just its own hairpin, or both overlap both?
                        my ($a, $b) = map {$unique_mirs->{$_}->{fields}} ($ak, $bk); 

                        my $aOverlapB = $a->{end} > $bs && $a->{start} < $be;
                        my $bOverlapA = $b->{end} > $as && $b->{start} < $ae;

                        if(!$aOverlapB && !$bOverlapA){
                            print STDERR "Weakly overlapping hairpins with no overlapping mature";
                            printOverlaps(\*STDERR, $unique_mirs, $ak, $bk);
                            next;
                        }

                        my $m;
                        my $prediction; # the hairpin that most likely can encorporate both predictions
                        my $predictionB; # secondary prediction where hairpin will be completely overwritten
                        if($aOverlapB && $bOverlapA) # bit of a tie, resolved by working out if one of the hairpins was made originally by miRCat
                        {
                            # does one of the structures contain some <>>><><><>
                            if($unique_mirs->{$bk}->{structure} =~ /\<|\>/){
                                $m = $a;
                                $prediction = $unique_mirs->{$bk};
                                $predictionB = $unique_mirs->{$ak};
                            }
                            else
                            {
                                $m = $b;
                                $prediction = $unique_mirs->{$ak};
                                $predictionB = $unique_mirs->{$bk};
                            }
                        }
                        elsif($aOverlapB) # mature seq of a overlaps hairpin B
                        {
                            # mature a overlaps hairpin B as well as mature b being in hairpin B
                            # Use hairpin B as a starting base
                            $m = $a;
                            $prediction = $unique_mirs->{$bk};
                            $predictionB = $unique_mirs->{$ak};
                        }
                        elsif($bOverlapA)
                        {
                            $m = $b;
                            $prediction = $unique_mirs->{$ak};
                            $predictionB = $unique_mirs->{$bk};
                        }

                        my %rnain = (seqid => $prediction->{fields}->{seqid},  strand => int($prediction->{fields}->{strand}."1"));
                        my $modded = 0;
                        if($m->{start} < $prediction->{'hp.start'}){ # mature overhang to the left 
                            @rnain{qw(start end)} = ($m->{start}, $prediction->{'hp.end'});
                            $modded++;
                        }
                        elsif($m->{end} > $prediction->{'hp.end'}){ # mature overhang to the right
                            @rnain{qw(start end)} = ($prediction->{'hp.start'}, $m->{end});
                            $modded++;
                        }
                        if($modded){ # the hairpin prediction needs to be extended to encorporate both predictions
                            my %rnafold = seq_rnafold(%rnain);
                            @{$prediction}{qw(hp.start hp.end hairpin structure)} = @rnafold{qw(start end seq struct)};
                        }
                        # update other prediction
                        @{$predictionB}{qw(hp.start hp.end hairpin structure)} = @{$prediction}{qw(hp.start hp.end hairpin structure)};
                        $predictionB->{arm} = findPrime($predictionB->{structure}, $predictionB->{fields}->{'mature.seq'}, $prediction->{hairpin}, $prediction->{fields}->{strand});    
                        my $newStruct;
                        my($mstart, $mend) = (( $predictionB->{fields}->{start} - $predictionB->{'hp.start'}), ($predictionB->{fields}->{end} - $predictionB->{'hp.start'}));
                        my $i=0;
                        my @struct = split("", $predictionB->{structure});

                        if ($predictionB->{fields}->{strand} eq "-"){
                            @struct = reverse @struct 
                        }
                        for (@struct)
                        {
                           if($i >= $mstart && $i <= $mend){
                               $_ =~ tr/()./<>-/;
                           }
                           else
                           {
                               $_ =~ tr/<>-/()./;
                           }
                           $newStruct .= $_;
                           $i++;
                        }
                        $newStruct = reverse $newStruct if $predictionB->{fields}->{strand} eq "-";
                        $predictionB->{structure} = $newStruct;

                        # check $a for previous overlaps
                        my $thisuid = $mirids2uid{$ak};

                        # check $b for previous overlaps
                        my $thisuidb = $mirids2uid{$bk};
                        
                        if(defined $thisuid && defined $thisuidb && $thisuid != $thisuidb){
                            die "uid assignment mismatch between overlapping miRNAs\n", Dumper($unique_mirs->{$ak}), Dumper($unique_mirs->{$bk}),"\n";
                        }

                        # If both already have a defined uid, it must be the same due
                        # to passing the die condition above and we don't need to process
                        # this overlap again (this probably is never the case)
                        if(!defined $thisuidb){ 
                            if (defined $thisuid){
                                # previous overlap found, use the uid
                                push @{$uid2mirids[$thisuid]}, $bk;
                                $mirids2uid{$bk} = $thisuid; 
                                $unique_mirs->{$bk}->{uid} = $thisuid;
                            }
                            else{
                                # no previous overlap found, create new entry
                                push @uid2mirids, [$ak,$bk]; 
                                $mirids2uid{$ak} = $#uid2mirids;
                                $mirids2uid{$bk} = $#uid2mirids;
                                $unique_mirs->{$ak}->{uid} = $#uid2mirids;
                                $unique_mirs->{$bk}->{uid} = $#uid2mirids;
                            }
                        }
                    }
                }
            }
        }

        # if no overlaps are found, store as a new uid
        if (!$aoverlaps && !defined $mirids2uid{$ak}){
            push @uid2mirids, [$ak];
            $mirids2uid{$ak} = $#uid2mirids;
            $unique_mirs->{$ak}->{uid} = $#uid2mirids;
        }
    }
    my %mircat_by_uids;
    # A uid groups all miRNAs into one precursor
    # Each prediction may have a variant of precursor
    # A summary is taken using min and max positions
    # ->{uid}->{minhp, maxhp, seqid, miRNAs->[{predictions}->[]]
    foreach my $uid (0..$#uid2mirids){
       my @mirids = @{$uid2mirids[$uid]};
       #      print $uid,": ",scalar(@mirids),"\n";
       my $min = (sort{$a <=> $b} map{$_->{"hp.start"}} @{$unique_mirs}{@mirids})[0];
       my $max = (sort{$b <=> $a} map{$_->{"hp.end"}} @{$unique_mirs}{@mirids})[0];
       my $seqid = $unique_mirs->{$mirids[0]}->{fields}->{seqid};
       my $strand = $unique_mirs->{$mirids[0]}->{fields}->{strand};
       $mircat_by_uids{getUid()} = {minhp => $min, maxhp => $max, seqid => $seqid, strand => $strand,
           miRNAs => [@{$unique_mirs}{@mirids}]}
    }
    return \%mircat_by_uids;
}

# Print the overlap of two hairpins from the unique_mirs hash
sub printOverlaps{
    my ($dbout, $unique_mirs, $ak, $bk) = @_;

    my ($aseq, $astruct, $as, $ae) = @{$unique_mirs->{$ak}}{qw(hairpin structure hp.start hp.end)};
    my ($bseq, $bstruct, $bs, $be) = @{$unique_mirs->{$bk}}{qw(hairpin structure hp.start hp.end)};
    my ($al, $bl) = (length($aseq), length($bseq));

    my $s = ($as, $bs)[$bs < $as];
    my $e = ($ae, $be)[$be > $ae];
    my $overlap = ($ae, $be)[$ae > $be] - ($as, $bs)[$as < $bs];
    my $l = ($al + $bl - $overlap);

    my $strand = $unique_mirs->{$ak}->{fields}->{strand};

    if($strand eq "-"){
        $bseq = reverse $bseq;
        $bstruct = reverse $bstruct;
        $aseq = reverse $aseq;
        $astruct = reverse $astruct;
    }

    print $dbout "--", $ak, " ", $bk, "--";
    print $dbout $unique_mirs->{$ak}->{fields}->{'mature.seq'}, " ", $unique_mirs->{$bk}->{fields}->{'mature.seq'}, "--\n";
    my @it = ($s..$e);
    my @aline;
    my @assline;
    my @bline;
    my @bssline;
    #$DB::single=1 if $aseq eq 'GAGTCACTGGTTGTATACAATGCATATGGGTGGCAACCATTTTGTGGCGGGATTTTGAATCGTGA' && $bseq eq 'GTGGCGGGATTTTGAATCGTGAGTATTGGATTAGATTTTGATTATGTTTCATACTCCGGCCAC';
    #print STDERR "\$aseq eq '",$aseq, "' && \$bseq eq '", $bseq, "'\n";
    #$DB::single=1 if length($aseq) != length($astruct);

    for my $i (@it)
    {
        my $n = (($i >= $as && $i <= $ae) ? substr($aseq,$i-$as, 1) : " ");
        my $ss = (($i >= $as && $i <= $ae) ? substr($astruct,$i-$as, 1) : " ");
        ($strand eq "-+") ? push @aline, $n : unshift @aline, $n;
        ($strand eq "-+") ? push @assline, $ss : unshift @assline, $ss;

        $n = (($i >= $bs && $i <= $be) ? substr($bseq,$i-$bs, 1) : " ");
        $ss = (($i >= $bs && $i <= $be) ? substr($bstruct,$i-$bs, 1) : " ");
        ($strand eq "-+") ? push @bline, $n : unshift @bline, $n;
        ($strand eq "-+") ? push @bssline, $ss : unshift @bssline, $ss;
    }
    print $dbout join("",@aline),"\n";
    print $dbout join("",@assline),"\n";
    print $dbout join("",@bline),"\n";
    print $dbout join("",@bssline),"\n";
    print $dbout  "\n";

    print $dbout "\n";
}

sub getUid{ my $thisuid = $uid_counter; $uid_counter++; return $thisuid}


sub merge_mirdeep_results{
    my $mirdeep_file = shift;

    my %unique_mds;
    open(my $mdfh, "<$mirdeepfile") or die "Can't open $mirdeepfile\n";
    my $colline = <$mdfh>;
    chomp $colline;
    my @colnames = split ",", $colline;
    while(<$mdfh>){
        chomp;
        my %md;
        @md{@md_colnames} = split ",";
        my $m_pos = index($md{"precursor"}, $md{"mature"});
        my ($start, $end) = @md{qw(start end)};
        die "Can't find mature in precursor\n", Dumper(\%md) if($m_pos < 0);

        my $strand = $md{"strand"};
        my($mstart, $mend);
        if($strand eq '+'){
            $mstart = $start + $m_pos;
            $mend = $mstart + length($md{"mature"});
        }else{
            $mstart = $end - ($m_pos + length($md{"mature"}));
            $mend = $mstart - length($md{"mature"});
        }

        # for mirdeep predictions we use the hairpin center to decide on the arm
        my $hairpinMiddle = int(length($md{precursor})/2);
        my $prime = "5_prime";
        $prime = "3_prime" if $m_pos > $hairpinMiddle;
        
        $md{"mstart"} = $mstart;
        $md{"mend"} = $mend;
       
        my $mdkey = join(":", $md{seqid}, $mstart, $mend, $strand);

        # value alterations
        # samples are codede "s#"
        $md{run} = 's'.$md{run};

        # sequence notation is capitalised and U->G
        $md{mature} =~ tr/atcguU/ATCGTT/;
        $md{precursor} =~ tr/atcguU/ATCGTT/;
        $md{star} =~ tr/atcguU/ATCGTT/;


        my %pmd;
        @pmd{@md_prediction_fields} = @md{@md_prediction_fields};
        my %umd;
        @umd{@md_unique_fields} = @md{@md_unique_fields};


        # assign a mirdeep uid and push to unique hash
        if($unique_mds{$mdkey}){
            my $uid = $unique_mds{$mdkey}->{id};
            push @{$unique_mds{$mdkey}->{predictions}}, \%pmd;
            $id++;
            # check and throw a warning if this hairpin is off from the stored hairpin
            if($md{start} != $unique_mds{$mdkey}->{'hp.start'} ||
               $md{end} != $unique_mds{$mdkey}->{'hp.end'}){
               print "Hairpins are slightly off:\n", 
                join("\n", map {join(",", @{$_}{@md_prediction_fields})} 
                            @{$unique_mds{$mdkey}->{predictions}}),
                "\n";
            }
        }
        else{
            # all common values are stored on this level
            
            $unique_mds{$mdkey} = {%umd, id => $id, 
                                   'hp.start' => $md{start},
                                   'hp.end' => $md{end},
                                   predictions =>[\%pmd],  
                                   arm => $prime};
        }

    }

    # group overlapping hairpins into the same mirid
    my @ukeys = sort sort_uids keys %unique_mds;
    my @uid2mirids;
    my %mirids2uid;

    foreach my $a (0..$#ukeys){
        my $ak = $ukeys[$a];
        my $aoverlaps = 0;
        foreach my $b ($a..$#ukeys){
            my $bk = $ukeys[$b];

            # if not same seq but on same seqid
            if($a != $b && $unique_mds{$ak}->{seqid} eq $unique_mds{$bk}->{seqid}
                        && $unique_mds{$ak}->{strand} eq $unique_mds{$bk}->{strand}){
                my($as, $ae) = @{$unique_mds{$ak}}{qw(hp.start hp.end)};
                my($bs, $be) = @{$unique_mds{$bk}}{qw(hp.start hp.end)};
                # check for coordinate overlap
                if($ae > $bs && $as < $be){
                    
                    my $overlap = ($ae, $be)[$ae > $be] - ($as, $bs)[$as < $bs];

                    # same precursor miRNA if overlap is > maximum hairpin length
                    my ($al, $bl) = ($ae - $as, $be - $bs);
                    my $maxl = ($al,$bl)[$al < $bl];
                    if($overlap > $maxl/2){
                        $aoverlaps = 1;

                        # check $a for previous overlaps
                        my $thisuid = $mirids2uid{$ak};

                        # check $b for previous overlaps
                        my $thisuidb = $mirids2uid{$bk};

                        if(defined $thisuid && defined $thisuidb && $thisuid != $thisuidb){
                            die "uid assignment mismatch between overlapping miRNAs\n", Dumper($unique_mds{$ak}), Dumper($unique_mds{$bk}),"";
                        }
                        
                        
                        if(!defined $thisuidb){ 
                            if (defined $thisuid){
                                # previous overlap found, use the uid
                                push @{$uid2mirids[$thisuid]}, $bk;
                                $mirids2uid{$bk} = $#uid2mirids;
                                $unique_mds{$bk}->{id} = $thisuid;
                            }
                            else{
                                # no overlap found, create new entry
                                push @uid2mirids, [$ak,$bk]; 
                                $mirids2uid{$ak} = $#uid2mirids;
                                $mirids2uid{$bk} = $#uid2mirids;
                                $unique_mds{$ak}->{id} = $#uid2mirids;
                                $unique_mds{$bk}->{id} = $#uid2mirids;
                            }
                        }
                    }
                }
            }
        }
    }

    # Checking this has worked by printing it back out right now
    open(my $mdout, ">mdout.csv") or die ":(";
    foreach my $mdkey (keys %unique_mds){
        my $umd = $unique_mds{$mdkey};
        my $uid = $umd->{id};
        foreach my $md (@{$umd->{predictions}}){
            print $mdout join(",", $uid, @{$md}{@md_prediction_fields}),"\n";
        }
    }

    my %mirdeep_by_uids;
    # A uid groups all miRNAs into one precursor
    # Each prediction may have a variant of precursor
    # A summary is taken using min and max positions
    # ->{uid}->{minhp, maxhp, seqid, miRNAs->[{predictions}->[]]
    foreach my $uid (0..$#uid2mirids){
       my @mirids = @{$uid2mirids[$uid]};
       #print $uid,": ",scalar(@mirids),"\n";
       my $min = (sort{$a <=> $b} map{$_->{"hp.start"}} @unique_mds{@mirids})[0];
       my $max = (sort{$b <=> $a} map{$_->{"hp.end"}} @unique_mds{@mirids})[0];
       my $seqid = $unique_mds{$mirids[0]}->{seqid};
       my $strand = $unique_mds{$mirids[0]}->{strand};
       $mirdeep_by_uids{getUid()} = {minhp => $min, maxhp => $max, seqid => $seqid, strand => $strand,
           miRNAs => [@unique_mds{@mirids}]};
    }

    return \%mirdeep_by_uids;
}

sub merge_mirdeep_mircat{
    my ($mircat, $mirdeep) = @_;

    my %predictions;
    my %overlapped_uids; # lists overlapped mapmi uids
    for my $mcuid (keys %$mircat){
        # this essentially a whole precursor with many miRNA variants
        my $mircatPrecursor = $mircat->{$mcuid};
        
        my $foundOverlap = 0;
        for my $mduid (keys %$mirdeep){
            my $mirdeepPrecursor = $mirdeep->{$mduid};
            my ($mcmin, $mcmax) = @{$mircatPrecursor}{qw(minhp maxhp)};
            my ($mdmin, $mdmax) = @{$mirdeepPrecursor}{qw(minhp maxhp)};
            
            # Assess overlap
            if($mircatPrecursor->{seqid} eq $mirdeepPrecursor->{seqid} &&
                $mircatPrecursor->{strand} eq $mirdeepPrecursor->{strand} &&
                overlaps_half($mcmin, $mcmax, $mdmin, $mdmax)){
                die "Too many overlaps between mircat/mirdeep precursors:\n",
                        Dumper($mirdeepPrecursor), Dumper($mircatPrecursor) if $foundOverlap;
                $overlapped_uids{$mduid} = 1;
                # Create an entry predicted by both tools
                # The Id's hairpin coordinates is the maximum of all found hairpins
                my $hpmin = ($mcmin, $mdmin)[$mcmin > $mdmin];
                my $hpmax = ($mcmax, $mdmax)[$mcmax < $mdmax];
                $predictions{$mcuid} = {start => $hpmin, end => $hpmax,
                    seqid => $mircatPrecursor->{seqid}, 
                    strand => $mircatPrecursor->{strand}, 
                    source => 1+2, # 1=mircat, 2=mirdeep, 4=mapmi
                    mircat => $mircatPrecursor,
                    mirdeep => $mirdeepPrecursor, mirname=>"MCD$mcuid"};
                $foundOverlap = 1; 
                print $mircatPrecursor->{seqid}, " ", $hpmin, " ", $hpmax, " ", $hpmax-$hpmin+1,"\n";
            }
            # keep going to check there is no funny chainy second overlap
        }
        if(!$foundOverlap){
            # no sequence found that overlaps with this mircat precursor. Add on its own

            $predictions{$mcuid} = {source => 1, mircat => $mircatPrecursor,
                                    seqid => $mircatPrecursor->{seqid}, strand => $mircatPrecursor->{strand},
                                    mirname => "MC$mcuid"};
            @{$predictions{$mcuid}}{qw(start end)} = 
                @{$mircatPrecursor}{qw(minhp maxhp)};
        }
    }

    # add all mirdeep seqs that were not found to overlap with mircat
    for my $mduid (keys %$mirdeep){
        if(!$overlapped_uids{$mduid}){
            # no sequence found that overlaps with this mirdeep precursor. Add on its own
            $predictions{$mduid} = {source => 2, mirdeep => $mirdeep->{$mduid},
                                    mirname => "MD$mduid"};
            @{$predictions{$mduid}}{qw(start end seqid strand)} = 
                @{$mirdeep->{$mduid}}{qw(minhp maxhp seqid strand)};
        }
    }
    return \%predictions;
}

sub merge_mapmi{
    my ($mapmifile, $predictions) = @_;
    my %mapmi;
    open(my $mdup, ">mapmidup.txt") or die "Can't open";
    open(my $mapmifh, '<', $mapmifile) or die "Couldn't open $mapmifile\n";
    print "Merging in Mapmi results\n";
    <$mapmifh>;

    # filter results down to the best unique matches
    while(<$mapmifh>){
        chomp;
        my @fields = split /,/;
        my %mm;
        @mm{@mapmi_header} = @fields[1..9, 12, 13, 15, 18];

        my $mmkey = join(":",@mm{qw(seqid start end strand)});

        my $hairpin = $mm{'pre.seq'};
        my $m_pos = index($hairpin, lc($mm{sequence}));
        my($start, $end) = @mm{qw(start end)};

        # find precursor coordinates
        my ($h_startpos, $h_endpos) = @mm{qw(pre.start pre.end)};
        my $orientation = $mm{strand};

        my ($sp) = $mm{id} =~ /^(\w{3})-/;
        die "Can't get species from $mm{id}" if !$sp;

        push @{$mm{sp}}, $sp;
#        if ($orientation eq '-'){
#            my $m_endpos = $m_pos;
#            $h_endpos = $end + $m_endpos;
#            $h_startpos = ($h_endpos - length($hairpin)) + 1;
#        }
#        else{
#            my $m_startpos = $m_pos;
#            $h_startpos = $start - $m_startpos;
#            $h_endpos = ($h_startpos + length($hairpin)) - 1;
#        }
#
        $mm{'pre.start'} = $h_startpos;
        $mm{'pre.end'} = $h_endpos;

        my $hairpinMiddle = int(length($mm{"pre.seq"})/2);
        my $prime = "5_prime"; 
        if($m_pos > $hairpinMiddle){
                $prime = "3_prime";
        }
        $mm{prime} = $prime;

        if(defined $mapmi{$mmkey}){
            # Assess the suitability of this entry rather than the current one
            next if $mapmi{$mmkey}->{id} eq $mm{id};
            my @splist;

            # add species ID to entry if the miRNA names are identical.
            #  Currently we "wipe" the species names if this is not the case
            #  This may cause issues with species names being "wiped" when they should stay
            if($mm{name} eq $mapmi{$mmkey}->{name})
            {
                @splist = (@{$mapmi{$mmkey}->{sp}}, $sp); 
            }

            if($mapmi{$mmkey}->{score} < $mm{score} || $mapmi{$mmkey}->{mismatch} > $mm{mismatch}){
                # This entry is better, replace current hash entry
                $mapmi{$mmkey} = \%mm;
            }
            elsif($mapmi{$mmkey}{score} == $mm{score} && $mapmi{$mmkey}->{mismatch} == $mm{mismatch}){
                # If there is a tie, we keep the current entry but print a warning
                print $mdup "Equally good mapmi entries:\n";
                print $mdup '1: ',join(",",@mm{@mapmi_header}),"\n";
                print $mdup '2: ',join(",",@{$mapmi{$mmkey}}{@mapmi_header}),"\n";
            }
            $mapmi{$mmkey}->{sp} = \@splist if @splist;
        }
        else{
            $mapmi{$mmkey} = \%mm;
        }
    }
    
   my $mir90uid; 
    # Check each mapmi entry against each prediction boundary (hpmin, hpmax)
    # If a mapmi entry overlaps a prediction, add to the prediction
    # If no matches, add it as a new mapmi-only prediction
    foreach my $mmkey ( keys %mapmi ){
        my %mm = %{$mapmi{$mmkey}};
        my $added = 0;
        # for each uid

        foreach my $uid (keys %$predictions){
            my $prediction = $predictions->{$uid};
            #next if $prediction->{source} == 4;

            my ($pseqid, $pstart, $pend) = @{$prediction}{qw(seqid start end)};
            my ($mseqid, $mstart, $mend) = @mm{qw(seqid start end)};
            # if this mapmi seq overlaps...
            if($mseqid eq $pseqid && 
                $pend > $mstart && $pstart < $mend){
                #overlaps_half($mm{start}, $mm{end}, $pstart, $pend)){

                # add to the prediction
                $prediction->{source} |= 4;
                push @{$prediction->{mapmi}} , \%mm;

                # flag as mapmi seq has been added
                $added++;
            }
        }
        # add as mapmi-only entry if not overlapping
        if(!$added){
            my $entry = {source => 4, mapmi => [\%mm]};
            @{$entry}{qw(seqid start end strand mirname)} = @mm{qw(seqid pre.start pre.end strand name)};
            my $thisUid = getUid();
            $mir90uid = $thisUid if $mm{name} eq "miR-190"; 
            $predictions->{$thisUid} = $entry;
        }
    }

    # resolve mapmi naming conflicts
    open(my $naming_conflicts_fh, ">$naming_conflicts_out") or die "Can't open $naming_conflicts_out";
    my %miRNA_names; # store used names to increment
    foreach my $uid (keys %$predictions){
        my $prediction = $predictions->{$uid};
        next if !$prediction->{mapmi};
        my %unique_names;
        foreach my $mm ( @{$prediction->{mapmi}}){
            $unique_names{$mm->{score}}->{$mm->{name}}++; 
        }
        my $thismismatch = (sort {$b <=> $a} keys %unique_names)[0];
        my $mmnames = $unique_names{$thismismatch};
        my @sortednames = sort {$mmnames->{$b} <=> $mmnames->{$a}} keys %$mmnames;
        my $namecount = $mmnames->{$sortednames[0]};
        my @topnames;
        foreach my $name (@sortednames){
            push @topnames, $name if $mmnames->{$name} == $namecount;
        }

        if(@topnames == 1){
            $prediction->{mirname} = $topnames[0]
        }
        else{
            print $naming_conflicts_fh "Naming conflict: \n", 
             join("\n", 
                 map{ join(",", @{$_}{qw(id sequence score)})} @{$prediction->{mapmi}}
             ),"";

            # Assign the topmost name anyway
            $prediction->{mirname} = $topnames[0];
            print $naming_conflicts_fh ", used ", $prediction->{mirname},"\n";
        }

        $miRNA_names{$prediction->{mirname}}++; # record number of times this name is used
    }

    # Finally, modify identical names by appending a -#
    my %mirnum;
    foreach my $k (sort {
                           return (ncmp($predictions->{$a}->{seqid}, $predictions->{$b}->{seqid}) 
                                        || $predictions->{$a}->{"start"} < $predictions->{$b}->{"start"} 
                                        || $predictions->{$a}->{"end"} < $predictions->{$b}->{"end"});
                        } keys %$predictions){
        my $p = $predictions->{$k};
        $mirnum{$p->{mirname}}++;
        my $suffix = $mirnum{$p->{mirname}};
        $p->{mirname} .= "-". $suffix if $miRNA_names{$p->{mirname}} && $miRNA_names{$p->{mirname}} > 1;
    }
    return $predictions;
}

sub sort_uids{
    my @sa = $a =~ /(.+):(\d+):(\d+):([+-])/;
    my @sb = $b =~ /(.+):(\d+):(\d+):([+-])/;
    die "Can't recover sid for sorting in comparison $a and $b" if !($a || $b);
    return $sa[0] cmp $sb[0] || $sa[1] <=> $sb[1] || $sa[2] <=> $sb[2] || $sa[3] cmp $sb[3];
}

sub sort_sids{
    my @sa = $a =~ /(.+)-(\d+)$/;
    my @sb = $b =~ /(.+)-(\d+)$/;
    die "Can't recover sid for sorting in comparison $a and $b" if !($a || $b);
    return $sa[0] cmp $sb[0] || $sa[1] <=> $sb[1];
}

# Checks if two coordinate systems overlap each other by more than
# half of the largest sequence
sub overlaps_half{
    my ($as, $ae, $bs, $be) = @_;
        
    if($ae > $bs && $as < $be){
        my $overlap = ($ae, $be)[$ae > $be] - ($as, $bs)[$as < $bs];

        # same miRNA if overlap is > maximum hairpin length
        my ($al, $bl) = ($ae - $as, $be - $bs);
        my $maxl = ($al,$bl)[$al < $bl];
        if($overlap > $maxl / 2){
            return 1;
        }
    }
    return 0;
}

