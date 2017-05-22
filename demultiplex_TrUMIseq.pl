#!/usr/bin/perl

use strict;
use Data::Dumper;
use Getopt::Long;

my $fastq1;     # input fastq file1 (R1)
my $fastq2;     # input fastq file2 (R2)
my $library;    # library name and index
my $bq=10;      # base quality cutoff (default=10)
my $mm=1;       # number of allowable mismatches in the first 6 nucleotides

GetOptions(
	   'f1=s'=>\$fastq1,
	   'f2=s'=>\$fastq2,
	   'library=s'=>\$library,
	   'bq=i'=>\$bq,
	   'mismatch=i'=>\$mm,
 );


#########################################################
# Bulid a quality index based on Illumina 1.8 phred score 
#########################################################
my @qual_Illumina=split('',"\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/0123456789\:\;\<\=\>\?\@ABCDEFGHIJ");
my %bq_index=(); # base quality index

for(my $i=0;$i<scalar(@qual_Illumina);$i++){
    $bq_index{$qual_Illumina[$i]}=$i;
}

#print Dumper\%bq_index;die;
print "(1)\n";
########################################################
## separate each library from '$library' into an array
########################################################
my @lib=(); # library array
my $count_lib=0; # number of indice

open(IN, "<$library") or die("Can't open the library info file!: $!");
while(<IN>){
    chomp;
    my @tmp=split('\t',$_);
    $lib[$count_lib]->[0]=$tmp[0]; # library name
    $lib[$count_lib]->[1]=$tmp[1]; # index info
    $count_lib++;
}

close IN;
print "(2)\n";
##################################################################################################################
## Store original fastq file into an array of arrays : array1='the read tag', array2=sequence, array3=base quality
##################################################################################################################
open(IN, "<$fastq1") or die "No such a input fastq R1 file!\n";
my %fastq1;
while(<IN>){
    chomp;
    if($_=~/^\@HWI/ || $_=~/^\@M02455/){ # for finding read tags of HiSeq or MiSeq
	my @tmp_key=split(' ',$_);
	$fastq1{$tmp_key[0]}->[0]=$_;           # read tag (line 1) : unique id is the key
	chomp($fastq1{$tmp_key[0]}->[1]=<IN>);  # sequence (line 2)
	<IN>;                                   # skip the 3rd line : '+'
	chomp($fastq1{$tmp_key[0]}->[2]=<IN>);  # base quality (line 4) 
    }
}
close IN; 
my $total_count1=scalar(keys (%fastq1)); # total number of reads
print "There are in total $total_count1 reads in this $fastq1 fastq file!\n";

open(IN, "<$fastq2") or die "No such a input fastq R2 file!\n";
my %fastq2;
while(<IN>){
    chomp;
    if($_=~/^\@HWI/ || $_=~/^\@M02455/){ # for finding read tags of HiSeq of Miseq
	my @tmp_key=split(' ',$_);
	$fastq2{$tmp_key[0]}->[0]=$_;           # read tag (line 1) : unique id is the key
	chomp($fastq2{$tmp_key[0]}->[1]=<IN>);  # sequence (line 2)
	<IN>;                                   # skip the 3rd line : '+'
	chomp($fastq2{$tmp_key[0]}->[2]=<IN>);  # base quality (line 4) 
    }
}
close IN; 
my $total_count2=scalar(keys (%fastq2)); # total number of reads
print "There are in total $total_count2 reads in this $fastq2 fastq file!\n";



################################################################################
## read each read from the '$fastq1' into separate fastq files (demultiplexed)
################################################################################
my @out_file_names1; # output file (fastq) names
my @out_file_names2; # output file (fastq) names

for(my $i=0; $i<$count_lib;$i++){ # for output given the number of libraries
    $out_file_names1[$i]=$fastq1."_".$lib[$i][0]."_".$lib[$i][1];
    $out_file_names2[$i]=$fastq2."_".$lib[$i][0]."_".$lib[$i][1];
    my $num_reads=0; # number of reads for this sample index

    open OUT1, ">$out_file_names1[$i]";
    open OUT2, ">$out_file_names2[$i]";

    foreach my $key (keys %fastq1){

	my @seqs1=split('',$fastq1{$key}->[1]);  # sequences R1
	my @seqs2=split('',$fastq2{$key}->[1]);  # sequences R2
	my @quals1=split('',$fastq1{$key}->[2]); # qualities R1
	my @quals2=split('',$fastq2{$key}->[2]); # qualities R2
	    
	my @sam_index1=@seqs1[0..6]; # seq of the first 7 nu (R1)
	my @sam_index2=@seqs2[0..6]; # seq of the first 7 nu (R2)
	my @sam_index1_qual=@quals1[0..6]; #qual of the first 7 nu(R1)   
	my @sam_index2_qual=@quals2[0..6]; #qual of the first 7 nu(R2)
	    
	my @final_seqs1=@seqs1[7..scalar(@seqs1)-1];         # seq of real reads
	my @final_seqs2=@seqs2[7..scalar(@seqs2)-1];         # seq of real reads
	my @final_quals1=@quals1[7..scalar(@seqs1)-1];       # quality of real reads
	my @final_quals2=@quals2[7..scalar(@seqs2)-1];       # quality of real reads
	    
	my $num_matches1=0; # num of mismatch with > quality cutoff each index for R1
	my $num_matches2=0; # num of mismatch with > quality cutoff each index for R1
	    
	my @index_I_want=split('',$lib[$i][1]);
	
       
	for(my $num=0;$num<6;$num++){ # for each base, count num of matches
	    if(($sam_index1[$num] eq $index_I_want[$num]) && ($bq_index{$sam_index1_qual[$num]}>=$bq)){ # if base is matched and its qualit is greater than equal the bq cutoff
		$num_matches1++
		}
	    if(($sam_index2[$num] eq $index_I_want[$num]) && ($bq_index{$sam_index2_qual[$num]}>=$bq)){ # if base is matched and its qualit is greater than equal the bq cutoff
		$num_matches2++
		}
	}

	my $T_match=0; # the 7th position 'T' check
	if((($sam_index1[6] eq 'T') && ($bq_index{$sam_index1_qual[6]}>=$bq))||(($sam_index2[6] eq 'T') && ($bq_index{$sam_index2_qual[6]}>=$bq))){
	    $T_match=1;
	}

	if ((($num_matches1>=6-$mm)&& ($T_match==1)) || (($num_matches2>=6-$mm)&&($T_match==1))){ # matched in the first 6 bases and the T in the final 7th position
	    my $tmp_seq1=join('',@final_seqs1);
	    my $tmp_seq2=join('',@final_seqs2);
	    my $tmp_quals1=join('',@final_quals1);
	    my $tmp_quals2=join('',@final_quals2);
	    
	    print OUT1 "$fastq1{$key}->[0]\n$tmp_seq1\n\+\n$tmp_quals1\n";
	    print OUT2 "$fastq2{$key}->[0]\n$tmp_seq2\n\+\n$tmp_quals2\n";


	    #print OUT1 "$fastq1{$key}->[0]\n$fastq1{$key}->[1]\n\+\n$fastq1{$key}->[2]\n";
	    #print OUT2 "$fastq2{$key}->[0]\n$fastq2{$key}->[1]\n\+\n$fastq2{$key}->[2]\n";

	    delete($fastq1{$key});
	    delete($fastq2{$key});

	    $num_reads++;
	}	    
    }
    close OUT1;
    close OUT2;
    
    print "Total number of reads with index $lib[$i][1] ($lib[$i][0]) is $num_reads\n";
}

####################
# UNMATCHED READS
#####################

my $out_file_names1=$fastq1."_"."NNNNNN";
my $out_file_names2=$fastq2."_"."NNNNNN";
open OUT1, ">$out_file_names1"; #print out non matched reads (R1)
open OUT2, ">$out_file_names2"; #print out non matched reads (R1)
my $num_reads=0;
foreach my $key (keys %fastq1){
    my @seqs1=split('',$fastq1{$key}->[1]);  # sequences R1
    my @seqs2=split('',$fastq2{$key}->[1]);  # sequences R2
    my @quals1=split('',$fastq1{$key}->[2]); # qualities R1
    my @quals2=split('',$fastq2{$key}->[2]); # qualities R2
    
    my @final_seqs1=@seqs1[7..scalar(@seqs1)-1];         # seq of real reads
    my @final_seqs2=@seqs2[7..scalar(@seqs2)-1];         # seq of real reads
    my @final_quals1=@quals1[7..scalar(@seqs1)-1];       # quality of real reads
    my @final_quals2=@quals2[7..scalar(@seqs2)-1];       # quality of real reads

    $num_reads++;

    my $tmp_seq1=join('',@final_seqs1);
    my $tmp_seq2=join('',@final_seqs2);
    my $tmp_quals1=join('',@final_quals1);
    my $tmp_quals2=join('',@final_quals2);
	    
    print OUT1 "$fastq1{$key}->[0]\n$tmp_seq1\n\+\n$tmp_quals1\n";
    print OUT2 "$fastq2{$key}->[0]\n$tmp_seq2\n\+\n$tmp_quals2\n";
}
close OUT1;
close OUT2;
print "Total number of reads with no index match is $num_reads\n";

exit
    
