#!/usr/bin/perl
use strict; use warnings;
use Data::Dumper;
use Benchmark qw(:hireswallclock); # to check execution time 
use Getopt::Long;

#######################################################################################################
# This script takes the SAM format file as an input and output a new SAM file without PCR duplicates
# Input format can be both 'Single end' or 'Paired end' sequencing data
#                                                   written by Jungeui Hong 08/16/2014 Version v2.0
########################################################################################################

###########################################
# Take the starting time of this operation
###########################################
my $starttime = Benchmark->new;
my $finishtime;
my $timespent;

##################################################################
# Input arguments : input SAM format file / number of mismatches #
##################################################################
my $in;     # input NNN.sam
my $mm=1;   # number of mismatches in the random barcode  
GetOptions(
           'in=s'=>\$in,
           'mismatch=i'=>\$mm,
	   );

#########################################
# for FINAL SAM OUTPUT printing         #
#########################################
my @file_name=split(/\./, $in);
my $sam_out=$file_name[0].".rmdup.sam";
open OUT_SAM, ">$sam_out";

#########################################
# for Statistics OUTPUT                 #
#########################################
my $stats=$in.".stats";
open OUT, ">$stats";

###############################################################################################
# Read file and make an hash [chr number] of arrays [chr position] of arrays [aligned reads]  #
#  ,where each chr position has multiple reads according to its coordinate position           #
###############################################################################################
open(IN, "<$in") or die "No such a input duplicated SAM file!\n";
my %data=();
my $all_aln_counter=0;       # all alignments in the SAM file
my $good_aln_counter=0;      # good alignments in the SAM file
my $poor_aln_counter=0;      # poor alignments to the reference seq  : 5th col = '0' 

my $all_alignment_counter=0; # total number of alignments (including multiple alignment)
while(<IN>){
    chomp;
    if ($_ =~ /^HWI/ || $_=~/^M02455/){  # if this line is a new read : HWI - Hiseq, M - Miseq
	$all_aln_counter++;
	my @tmp_string=split(/\t/,$_);
	if($tmp_string[8] == 0 || $tmp_string[4] == 0){       # count multi-aligned reads   
	    $poor_aln_counter++;
	}else{  # only collect aligned reads
#	    my $tmp_key=$tmp_string[2].":".$tmp_string[3].":".$tmp_string[5].":".$tmp_string[6].":".$tmp_string[7].":".$tmp_string[8];
	    my $tmp_key=$tmp_string[2].":".$tmp_string[3].":".":".$tmp_string[6].":".$tmp_string[7].":".$tmp_string[8]; # make a key by concatenating chr_read1:position_read1:chr_read2:position_read2:input_size
	    push @{$data{$tmp_key}},$_;
	    $good_aln_counter++;
	}
    }else{ # print out header section
	print OUT_SAM "$_\n";
    }
}
#print Dumper\%data;die;
##################################
# PRINT OUT ALIGNMENT STATISITCS #
##################################
print OUT "This file is $in\n";
print OUT "Total num of all alignments is     $all_aln_counter\n";
print OUT "Total num of poor alignments is    $poor_aln_counter\n";
print OUT "Total num of good alignments is    $good_aln_counter\n";

print "This file is $in\n";
print "Total num of all alignments is     $all_aln_counter\n";
print "Total num of poor alignments is    $poor_aln_counter\n";
print "Total num of good alignments is    $good_aln_counter\n";

close IN;


########################################################################
# SELECT ONLY NON REDUNDANTLY ALIGNED READS FROM THE ORIGINAL SAM FILE #
#######################################################################

#my $dup_aln_counter=0;
my $uniq_aln_w_index=0;                         # all uniquely aligned reads when using index info 
my $uniq_aln_wo_index=scalar(keys %data);       # all redundantly aligned reads when not using index info

#print Dumper\%data;die;

foreach my $tmp_key (keys %data){ # for each unique key : based on alignment (coordinates) information : paired reads will be tested separately

    ###########################################################################################################
    # Take all reads aligned in the same coordinates                                                          #
    #   : sort them by reads name considering all reads sharing the same key are aligned in the same position #
    ###########################################################################################################
    my @all_reads= sort @{$data{$tmp_key}}; 

    my $dup_len=scalar @all_reads;
    if($dup_len > 1){  # if any read has the same alignment with any other reads

	foreach my $a (@all_reads){

	    my @tmp_reads=grep {$_ ne $a} @all_reads; # take all reads except $a

	    foreach my $b (@tmp_reads){  # compare $a and $b (every element other than $a)

		my @info1=split(/\t/,$a);my @info2=split(/\t/,$b); 
		my $len_info1=scalar(@info1);my $len_info2=scalar(@info2);
		my @tmp1=split(/\:/,$info1[$len_info1-1]);my @tmp2=split(/\:/,$info2[$len_info2-1]);
		my $len_tmp1=scalar(@tmp1);my $len_tmp2=scalar(@tmp2);

		my @ran_index1=split(//,$tmp1[$len_tmp1-1]); # 1st random 6mer
		my @ran_index2=split(//,$tmp2[$len_tmp2-1]); # 2nd random 6mer
		
		my $tmp_matches=0; # number of matches in the random 6mer
		for(my $count_matches=0;$count_matches<6;$count_matches++){
		    if ($ran_index1[$count_matches] eq $ran_index2[$count_matches]){
			$tmp_matches++;
		    }
		}
		if ($tmp_matches>=6-$mm){ # if two reads have the same random 6mers (- mismatch)
		    @all_reads=grep {$_ ne $b} @all_reads; # if two reads have the same random 6mer index, remove the later one from the @all_reads array
		}
	    }
	}
    }
    ###################################################
    # PRINT OUT NON-PCR DUPLICATE READS : FINAL OUTOUT
    ###################################################
    for(my $i=0;$i<scalar(@all_reads);$i++){
	print OUT_SAM "$all_reads[$i]\n";
	$uniq_aln_w_index++;
    }    
}

close OUT_SAM; 

my $per_pcr_dup_w_index=sprintf("%.1f",($good_aln_counter-$uniq_aln_w_index)/$good_aln_counter*100);
my $per_pcr_dup_wo_index=sprintf("%.1f",($good_aln_counter-$uniq_aln_wo_index)/$good_aln_counter*100);

print OUT "Total number of non-duplicated good alignments w.o index (pcr dup rate) is $uniq_aln_wo_index ($per_pcr_dup_wo_index%)\n";
print OUT "Total number of non-duplicated good alignments w index   (pcr dup rate) is $uniq_aln_w_index ($per_pcr_dup_w_index%)\n\n"; 
print "Total number of non-duplicated alignments w.o index (pcr dup rate) is $uniq_aln_wo_index ($per_pcr_dup_wo_index%)\n";
print "Total number of non-duplicated alignments w index   (pcr dup rate) is $uniq_aln_w_index ($per_pcr_dup_w_index%)\n\n"; 



###########################
# Take the finishing time
###########################
$finishtime = Benchmark->new;
$timespent = timediff($finishtime,$starttime);
print OUT "\nDone!\nSpent". timestr($timespent). "\n";
print "\nDone!\nSpent". timestr($timespent). "\n";
close OUT;


exit 0;
