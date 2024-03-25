#!usr/bin/perl
use warnings;use strict;
my $input = shift @ARGV;

my (%rep_hash,@rep_array);

open REP_TOTAL, "$input" or die;
while (<REP_TOTAL>) {
chomp $_;
my @split = split "\t", $_;
my $kmer = $split[5]; 
$kmer=~ s/,.*//g;
my $kmerx2 = "$kmer$kmer";   
my @array = ();   
my $sort_counter = 0;     
my $kmer_length = length($kmer);  
while ($sort_counter < $kmer_length) {  
my $subkmer = (substr $kmerx2, $sort_counter, $kmer_length);
push @array, $subkmer;     
my $revcomp = reverse($subkmer);
$revcomp =~ tr/ACGTacgt/TGCAtgca/;
push @array, $revcomp;
$sort_counter ++;
}
@array = sort(@array); # alphabetizing all possible offsets of the repeat.
$kmer = $array[0];  
push @rep_array, "$split[0]\t$split[1]\t$split[2]\t$split[3]\t$split[4]\t$split[5]\t$kmer\n";
if ($rep_hash{$kmer}) {
next;  
} else {     
$rep_hash{$kmer} = 0;     
}       
}       
close REP_TOTAL;

#foreach my $rep (sort { $a cmp $b } keys %rep_hash) {
#        print "$rep $rep_hash{$rep}\n";
#}

foreach my $rep (@rep_array) {print "$rep";}
#
#
#
