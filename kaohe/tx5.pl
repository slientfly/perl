#usr/bin/perl
open F,'<chr1.fa';
use strict;
use warnings;
my $all="";
my $countG;
my $all_length;
my $GC;
my $countC;
while(<F>){
chomp;
$all.=$_;
}
$countG=$all=~tr/Gg/Gg/;
$countC=$all=~tr/Cc/Cc/;
$all_length=$all=~tr/AGTCatgc/ATGCatgc/;
$GC=($countG+$countC)/$all_length;
print "$GC\n";
