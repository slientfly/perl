#usr/bin/perl
use strict;
use warnings;
my $A='CTGCTCCCTGCATGAGGAGGACACCCAGAGACATGAGACCTACCACCAGCAGGGGCAGTGCCAGGTGCTGGTGCAGCGCTCGCCCTGGCTGATGATGCGGATGGGCATCCTCGGCCGTGGGCTGCAGGAGTACCAGCTGCCCTACCAGCT';
$A=~tr/ATGC/TACG/;
my $i;
my $w;
##for ($i=0;$i<=110;$i+=5){
##print  substr($A,$i,40);
##}
while($i<=110){
$w=substr($A,$i,40);

print " "x$i."$w\n";
$i+=5;
