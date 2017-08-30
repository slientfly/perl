open XLS,"/disk2/mygeno/chenyj/Python/test/SNP_INDEL.xls";
open Nxls,">/disk2/mygeno/pengj/perl/ceshi.xls";
$a=<XLS>;
@list=split '\t',$a;
foreach(@list){
chomp;
$index1=index($_,'MutRatio');
$index2=index($_,'MutCount');
}
print Nxls $a;
while(<XLS>){
@temp=split '\t';
print Nxls if $temp[$index1]<0.3 and $temp[$index2]<5;
}
