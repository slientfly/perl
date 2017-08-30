#usr/bin/perl 
use strict;
use warnings;

open F,'</disk2/mygeno/chenyj/Python/test/SNP_INDEL.xls';
open D,'>new_file.xls';
my $line=<F>;

print D $line ;
while(<F>){
readline F;
chomp;
my @a=split'\t';
print D if $a[27]<0.3 and $a[29]<5;
}

close F;
close D;
