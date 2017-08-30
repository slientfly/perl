#!usr/bin/perl -w
die "Usage: perl $0   <scaffod.gff.gz> <chr.agp.gz> " unless (@ARGV == 2);

open (AGP,"gzip -cd $ARGV[1] | ") || die $!;
while(<AGP>){
	chomp;
	@inf = split;
	if ($inf[4] ne 'N'){
                $object{$inf[5]}{'s'}=$inf[1];
                $object{$inf[5]}{'e'}=$inf[2];
                $query{$inf[5]}{'s'}=$inf[6];
                $query{$inf[5]}{'e'}=$inf[7];
        }
}
close AGP;

open (IN, "gzip -cd $ARGV[0] | ") || die $!;
open (OUT,"| gzip > chromosome_gene.gff.gz");
while(<IN>){
        chomp;
        @inf=split;
        $changes=$object{$inf[0]}{'s'}+$inf[3]-$query{$inf[0]}{'s'};
        $inf[4]=$changes+$inf[4]-$inf[3];
        $inf[3]=$changes;
        $line=join("\t",@inf);
        print OUT "$line\n";
}
		
close OUT;
close IN;





