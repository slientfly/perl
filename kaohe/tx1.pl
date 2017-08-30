usr/bin/perl 
my $a='TTGCTAAAGGAAAGTGAAAGTGAAAGGAAGAGTCCTACGTCTGTCACTTTATGTCAA';
my $b=reverse $a;
$b=~tr/ATGC/TACG/;
print "$b";
