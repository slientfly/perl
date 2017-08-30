#!usr/bin/perl -w
my $OMIM;
open (IN,"mim2gene.txt")||die ("failed");
$countchr1=0;
$countchr2=0;
$countchr3=0;
$countchr4=0;
$countchr5=0;
while (<IN>){
	chomp;
	@line=split/\t/;
	$seq=$line[2];
	if ($seq=~/NA/){
			next;
		}
	else{
		if($seq=~m/^1[a-z]/){
			$countchr1++;}
		elsif($seq=~m/^2[a-z]/){
			$countchr2++;}
		elsif($seq=~m/^3[a-z]/){
			$countchr3++}
		elsif($seq=~m/^13[a-z]/){
			$countchr4++;}
		elsif($seq=~m/mitochondria/){
			$countchr5++;}
		else{
			next;			
		}
		}
}	
	print"chr1:$countchr1\tchr2:$countchr2\tchr3:$countchr3\tchr13:$countchr4\tmitochondria:$countchr5\n";

close IN;
open (OF,"omim.txt")||die ("failed");
$omim1=0;
$omim2=0;
$omim3=0;
$omim4=0;
$omim5=0;
	while (<OF>){
		chomp;
		@list=split;
		$aa=$list[0];
		if($aa=~m/(^\*)/){
			$omim1++;}
		elsif($aa=~m/^\+/){
			$omim2++;}
		elsif($aa=~m/^#/){
			$omim3++;}
		elsif($aa=~m/^%/){
			$omim4++;}
		else{
			$omim5++;
		}
	}
	print "*:$omim1\t+:$omim2\t#:$omim3\t%:$omim4\tn:$omim5\n";