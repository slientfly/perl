#usr/bin/perl
#use strict;
#use warnings;
open F,"/disk1/bed/LAMA2.bed"|| die "not open the file";
my %hash;
#my $chr
my $list1;
my $list2;
#
while(<F>){
	chomp;
	my @all=split/\t/;
	my $chr=$all[0];
#	my $list1;
#	my $list2;
	if ($.==1){
			$list1=$all[1];
       	    $list2=$all[2];
		 	$hash{$list1}=$list2;
#print $hash{$list1};
		}

	else{
		if ($hash{$list1} ==  $all[1]){
	#	print "$hash{$list1}\n";			
			$list1=$all[1];
        	$list2=$all[2]; 
			$hash{$list1}=$list2;
			}
		else{
			print "$chr\t$hash{$list1}---$all[1]\t没有捕获\n";
			$list1=$all[1];
            $list2=$all[2];
			$hash{$list1}=$list2;
			}
		}
	}

	

    

