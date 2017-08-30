#! usr/bin/perl 
open File,"</disk1/DB/hg19/Sequence/WholeGenomeFasta/genome.fa"||die("Could not open file");
while(<File>){
		chomp;
		if(/>(.*)/){
		$name=$1;}
		else{
		$hash{$name}.=$_;}
		#print $name;
		}
	#	print "chr1-----$hash{chr1}\n";
		foreach $key(sort keys %hash){
		#	print "$hash{$name}\n";
			$seq=$hash{$key};
		#	$seq=~s/[Nn]//g;
			while($seq=~m/([atgc]+)/g){
				my $list=$1;
				$len=length($1);
				$end=pos($seq);
				$start=$end-$len+1;
				print "$key\t$start\t$end\t$list\n";
			}
		}
close File;
