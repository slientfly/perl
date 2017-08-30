#usr/bin/perl
use strict;
use warnings;
my $i;
open DMD,"DMD_exon";
open REGION,"/disk2/mygeno/chenyj/Python/test/region.txt";
my @reg=<REGION>;
close REGION;
while(<DMD>){
chomp;
my @a=split'\t';
foreach(@reg){
chomp;
my @b=split'\t';
if($b[0]eq$a[0]and $b[1]eq$a[1]){
$i=0;
last}
$i=1;
}
print "DMD:$a[2]: 被完全覆盖\n" if $i==0;
print "DMD:$a[2]:未被完全覆盖\n" if $i==1;
