#usr/bin/perl -w
use strict;
open F,"/disk2/mygeno/chenyj/Python/test/region.txt";
open D,"/disk1/DB/hg19/Sequence/WholeGenomeFasta/genome.fa";
my $chrx_seq;
while(<D>){
last if /^>chrX/}
while(<D>){
last if /^>chrY/;
$chrx_seq.=$_;
}
$chrx_seq=~s/Nn//g;
while(<F>){
next if$.==1;
chomp;
my @item=split'\t';
my $len=$item[1]-$item[0]+1;
print "chrX\t".$_."\t".substr($chrx_seq,$item[1],$len)."\n";
}
close F;
close D;
