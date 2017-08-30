#! usr/bin/perl 
open File,"/disk2/mygeno/pengj/xutengfei/plum_0630.scafSeq.FG"||die("Could not open file");
while(<File>){
        chomp;
        if(/>(.*)/){
        $name=$1;}
        else{
        $hash{$name}.= $_;}
        #print "$name\n";
        }
       	foreach $key(sort keys %hash){
       		$value=$hash{$key};
       		$value=~s/N//g;
       		while($value=~m/([AGTC]+)/g){
       			$number=length($1);
       			#print "$key\t$number\n";
       			$i=0;
       			while($i<=$number){
       				$seq=substr($value,$i,250);
       				$i=$i+250;
       				$countG=$seq=~s/G/G/ig;
       				$countC=$seq=~s/C/C/ig;
       				$countGCrate=($countG+$countC)/250;
       				#$countGC=$seq=~tr/GC/GC/ig;
					   #$countGCrate=$countGC/250;
					print "$key\t$countGCrate\n";
       			}
       		}
      # 	print "$key\n";
       	}
