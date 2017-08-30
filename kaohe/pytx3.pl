open FILE,"/disk1/DB/hg19/Sequence/WholeGenomeFasta/genome.fa";
$start=0;
while(<FILE>){
chomp;
if (/^>(chr\w+)/){
if ($seq){
$seq=~s/[Nn]//g;
if ($seq=~/(^[atcg]+)/){
$start=1;
$end=length $1;
print "$chr\t$start\t$end\t$1\n";
$start=$end+1;
}
while($seq=~/([ATCG]+)([atcg]+)/g){
$start+=length($1);
$end=$start+length($2)-1;
print "$chr\t$start\t$end\t$2\n";
$start=$end+1;
}
}
$chr=$1;
$seq='';
$start=0;
}
else{
$seq.=$_;
}
END{
$seq=~s/[Nn]//g;
if ($seq=~/(^[atcg]+)/){
$start=1;
$end=length $1;
print "$chr\t$start\t$end\t$1\n";
$start=$end+1;
}
while($seq=~/([ATCG]+)([atcg]+)/g){
$start+=length($1);
$end=$start+length($2)-1;
print "$chr\t$start\t$end\t$2\n";
$start=$end+1;
}
}
}
