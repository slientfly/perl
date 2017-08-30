#usr/bin/perl
use strict;
use warnings;
open F,"</disk1/DB/hg19/Sequence/WholeGenomeFasta/genome.fa";
my $all="";
while(<F>){
chomp;
$all.=$_;
}
my $count;
$count=$all=~tr/atgc/AGTC/;
print $count;
#print $all;
