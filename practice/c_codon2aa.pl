#!usr/bin/perl
die "Usage: perl $0   <CDS.fasta> <codon.txt> " unless (@ARGV == 2);

open (CD,"$ARGV[1]") || die $!;
while(<CD>){
	chomp;
	@inf = split /'/;
	$codon2aa{$inf[1]}=$inf[3];
}
close CD;

#print $codon2aa{"ATG"};

open (FA,"$ARGV[0]") || die $!;
open (OUT,"> CDS2aa.fasta") || die $!;
while(<FA>){
        chomp;
        if (/^>/){
                $tag=$_;
                print OUT "$tag\n";
         }
        else{
                for($num=0;$num<length($_);$num=$num+3){
			if (length($_)-$num<3){last;}
                        $codon=substr($_,$num,3);
                        $aa=$codon2aa{$codon};
                        $pro.=$aa;

                }
                print OUT "$pro\n";
                undef $pro;
        }
}
close FA;
close OUT;                 

