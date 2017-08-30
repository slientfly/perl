use strict;
use warnings;
use feature 'say';
use Data::Dumper;

my %h;
while(<DATA>){
    chomp;
    next unless $_;
    next if /^#/;
    my $depth;
    my ($sample,$lane,$chip) = split /\s+/;
    my $flg = undef;
    my ($r1,$r2) = (undef,undef);
    for my $path ("/disk1/data/MG-N111/17O003184_OC017_CapNGS","/Data/FQDATA4/","/Data/FQDATA3/","/Data/FQDATA2/","/Data/FQDATA/","/Data/FASTQ/","/Data/FQDATA5/tech_211","/Data_tmp/FASTQ/"){#/Data/FASTQ
        if (-d "$path/$lane"){
            if(-e "$path/$lane/$sample/${sample}_R1.fq.gz" && -e "$path/$lane/$sample/${sample}_R2.fq.gz" && -s "$path/$lane/$sample/${sample}_R1.fq.gz" > 0 && -s "$path/$lane/$sample/${sample}_R2.fq.gz"> 0){
                $r1 = "$path/$lane/$sample/${sample}_R1.fq.gz";
                $r2 = "$path/$lane/$sample/${sample}_R2.fq.gz";
                $flg = 1;last;
            }
            if(-d "$path/$lane/Sample_$sample"){
                chomp(my $tmp1 = `ls $path/$lane/Sample_$sample/*_R1_001.fastq.gz 2>/dev/null`);
                chomp(my $tmp2 = `ls $path/$lane/Sample_$sample/*_R2_001.fastq.gz 2>/dev/null`);
                if(-e $tmp1 && -e $tmp2 && -s $tmp1 > 0 && -s $tmp2 > 0){
                    $r1 = $tmp1;$r2 = $tmp2;
                    $flg = 1;last;
                }
            }
            chomp(my $tmp1 = `ls $path/$lane/$sample*_R1_001.fastq.gz 2>/dev/null`);
            chomp(my $tmp2 = `ls $path/$lane/$sample*_R2_001.fastq.gz 2>/dev/null`);
            if(-e $tmp1 && -e $tmp2 && -s $tmp1 > 0 && -s $tmp2 > 0){
                $r1 = $tmp1;$r2 = $tmp2;
                $flg = 1;last;
            }
            chomp($tmp1 = `ls $path/$lane/$sample\_R1.fq.gz 2>/dev/null`);
            chomp($tmp2 = `ls $path/$lane/$sample\_R2.fq.gz 2>/dev/null`);
            if(-e $tmp1 && -e $tmp2 && -s $tmp1 > 0 && -s $tmp2 > 0){
                $r1 = $tmp1;$r2 = $tmp2;
                $flg = 1;last;
            }
        }
    }
    unless ($flg){
        say STDERR "找不到 $sample 样本的fastq文件";
        next;
    }
    if($r1 && $r2){
        $h{$sample}{'R1'} = $r1;
        $h{$sample}{'R2'} = $r2;
        $h{$sample}{'lane'} = $lane;
    }
    open S,">$sample.sh";
    say S "cd /ssd1/mygeno/pengj/project;perl /disk2/mygeno/pengj/pipline/normal.pl -s $sample -r1 $r1 -r2 $r2 -c $chip\n\n";
    close S;
    print my $rst = `qsub -q all.q -l mem=15gb,walltime=100:00:00,nodes=1:ppn=1 $sample.sh`;
}
#print Dumper \%h;
#lane编号  样本编号   性别（F-女性，M-男性，不知道性别用 - 代替）
__DATA__
17C012763_H001_CapNGS MG-L626 MG_BIK_V1
#17C001237_MedE003_CapNGS    MG-L645  MG_One_V2
# zy-4-352    MG-N082 PKU
# zy-4-403    MG-N082 PKU
#17C016971_KY019_CapNGS  MG-L581 PKU
#17C016970_KY019_CapNGS  MG-L581 PKU
# 17C025548_KNGS003_CapNGS    MG-L615 PKU
# 17C025549_KNGS003_CapNGS    MG-L615 PKU
# 17C025547_KNGS003_CapNGS    MG-L615 PKU
# 17C025546_KNGS003_CapNGS    MG-L615 PKU
# 17C025545_KNGS003_CapNGS    MG-L615 PKU
# 17C025544_KNGS003_CapNGS    MG-L615 PKU
# 17C025543_KNGS003_CapNGS    MG-L615 PKU
# 17C025542_KNGS003_CapNGS    MG-L615 PKU
# 17C025541_KNGS003_CapNGS    MG-L615 PKU
# 17C025540_KNGS003_CapNGS    MG-L615 PKU
# 17C025539_KNGS003_CapNGS    MG-L615 PKU
# 17C025538_KNGS003_CapNGS    MG-L615 PKU
# 17C025537_KNGS003_CapNGS    MG-L615 PKU
# 17C025536_KNGS003_CapNGS    MG-L615 PKU
# 17C025535_KNGS003_CapNGS    MG-L615 PKU
# 17C025534_KNGS003_CapNGS    MG-L615 PKU
