use strict;
use warnings;
use Data::Dumper;

my %pku;
open (IN,"/disk1/bed/LAMA2.bed");
while(<IN>){
	chomp;
	my @all=split/\t/;
	my $chr=$all[0];
	my @gene=split/:/,$all[4];
	my $gene=$gene[0];
	if(defined $pku{$gene}{'start'}){
		if($pku{$gene}{'end'}==$all[1]){
			$pku{$gene}{'start'}=$all[1];
			$pku{$gene}{'end'}=$all[2];
		}
		else{
		#	print Dumper(\%pku);
			print "$chr\t$pku{$gene}{'end'}\t$all[1]\n";
			$pku{$gene}{'start'}=$all[1];
			$pku{$gene}{'end'}=$all[2];
		}
	}
	else{	
		$pku{$gene}{'start'}=$all[1];
		$pku{$gene}{'end'}=$all[2];
	}
}
