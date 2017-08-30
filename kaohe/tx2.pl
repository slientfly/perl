#usr/bin/perl
my @xulie=qw/A T G C/;
my $s=0;
my $length=50;
my @dna;
while($s<$length){
$x=$xulie[rand(@xulie)];
push(@dna,$x);
$s+=1;
}
print "@dna";

