#!usr/bin/perl -w
die "Usage: perl $0   <scaffod.gff> <chr.fasta> " unless (@ARGV == 2);

open (FA,"gzip -cd $ARGV[1] | ") || die $!;
while(<FA>){
        chomp;
        if (/^>/){
                next;
        }
        else{
                $fasta.=$_;
        }
}
close FA;

open (GFF,"gzip -cd $ARGV[0] | ") || die $!;
open (OUT,">CDS.fasta");
$count=1;
while(<GFF>){
        chomp;
        @inf=split;
        if ($inf[2]=~/CDS/){
                $seq=substr($fasta,$inf[3]-1,$inf[4]-$inf[3]+1);
                print OUT ">CDS_$count;$inf[0];$inf[8]\n$seq\n";
                $count++;
        }
}

close GFF;
close OUT;
