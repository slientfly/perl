use Excel::Writer::XLSX;
#use Smart::Comments;
use Encode;
$filename=shift;
$filename=~/(.*).xls/i or die "未指定xls文件!";
$filename_rmext=$1;
$filename_rmext=$1 if $filename_rmext=~/.*\/(.*)/;
$workbook=Excel::Writer::XLSX->new ("$filename_rmext"."_selectedgene.xlsx");
$worksheet1=$workbook->add_worksheet('RawData');
$worksheet2=$workbook->add_worksheet("Selected gene");
open(XLS1,$filename) || die $!;
$first=<XLS1>;
$first=decode("cp936",$first);
@firstcol=split '\t',$first;
$index=0;
$i=1;
$j=2;
@list=qw();
for (@firstcol){
$tag_index=$index if /Tag/i;
$publi_index=$index if /Publication/i;
$mutinn_index=$index if /MutinNormal/i;
$Eff_index=$index if  /Effect/i;
$mutratio_index=$index if /MutRatio/i;
$mutcount_index=$index if /MutCount/i;
$Generegion_index=$index if /Gene_Region/i;
$g2015apr_index=$index if /1000g2015aug_all/i;
$ESP_index=$index if /ESP6500/i;
$inhouse_index=$index if /Inhouse/i;
$Ex_ALL_index=$index if /ExAC_ALL/i;
$EX_EAS_index=$index if /ExAC_EAS/i;
$sift_index=$index if /SIFT_Predict/i;
$polyphen_index=$index if /PolyPhen_2_Predict/i;
$taster_index=$index if /MutationTaster_Predict/i;
$acid_index=$index if /Amino_Acid_Changes/i;
$grep_index=$index if /GERP\+\+_Predict/i;
$index++;
}

=head1
$mutinn_index
$Eff_index
$mutratio_index
$mutcount_index
$Generegion_index
$g2015apr_index
$ESP_index
$sift_index
$polyphen_index
$taster_index
$acid_index
$grep_index
=cut


$worksheet1->write_row(A1,\@firstcol);
push @firstcol,'risk';
$worksheet2->write_row(A1,\@firstcol);
while(<XLS1>){
$_=decode("cp936",$_);
@col=split '\t';
$worksheet1->write_row("A$j",\@col);
$j++;
next if not $col[$mutinn_index]=~m/N\/A/i;
next if not $col[$mutratio_index]>=0.3;
next if not $col[$mutcount_index]>=10;
next if $col[$Generegion_index]=~/intronic/i;
if ($col[$Generegion_index]=~/UTR/i){
next unless $col[$tag_index] eq 'DM';
};
push @list,$_;
$i++;
$worksheet2->write_row("A$i",\@col);
}
close XLS1;
open(XLS2,$filename) or die $!;
while(<XLS2>){
$_=decode("cp936",$_);
@col=split '\t';
next if $.==1;
next if $col[$publi_index]=~/#N\/A/i;

next if not $col[$mutratio_index]>=0.3;
next if not $col[$mutcount_index]>=10;
next if $col[$Generegion_index]=~/intronic/i;
if ($col[$Generegion_index]=~/UTR/i){
next unless $col[$tag_index] eq 'DM';
};
push @list,$_;

$i++;
$worksheet2->write_row("A$i",\@col);
}
close XLS2;
open (XLS3,$filename) or die $!;
while(<XLS3>){
$_=decode("cp936",$_);
@col=split '\t';
next if $.==1;
next unless ($col[$Eff_index]=~/^frameshift|^splicing|^stopgain|^stoploss/i );
next if not $col[$mutratio_index]>=0.3;
next if not $col[$mutcount_index]>=10;
next if $col[$Generegion_index]=~/intronic/i;
if ($col[$Generegion_index]=~/UTR/i){
next unless $col[$tag_index]=~/DM/i;
};
if ($col[$Generegion_index]=~/UTR/i){
next unless $col[$tag_index]=~/DM/i;
};
if ($col[$Generegion_index]=~/UTR/i){
next unless $col[$tag_index] eq 'DM';
};
push @list,$_;

$i++;
$worksheet2->write_row("A$i",\@col);
}
@uniq_list=grep {$hash{$_}++<1;} @list;
for (sort @uniq_list){
chomp;
@temp=split '\t';
$i++;
$badpre++ if $temp[$sift_index]=~/deleterious/i;
$badpre++ if $temp[$PolyPhen_2_Predict]=~/damaging/i;
$badpre++ if $temp[$taster_index]=~/deleterious/i;
$badpre++ if $temp[$grep_index]=~/^Conserved/i;
$norelt++ if $temp[$sift_index] eq '-';
$norelt++ if $temp[$PolyPhen_2_Predict] eq '-';
$norelt++ if $temp[$taster_index]  eq '-';
$norelt++ if $temp[$grep_index] eq '-';

if( $temp[$tag_index] eq 'DM'){
push @temp,'high_risk';
}
elsif($temp[$acid_index]=~/^p/ and $badpre>=2){
push @temp,'med_risk'}
elsif(($temp[$acid_index]=~/^p/ and $badpre<=1) or ($temp[$Eff_index]=~/unknown|nonsynonymous/i and $norelt==0)){
push @temp,'med_low_risk'}
else{
push @temp,'low_risk'}
$worksheet2->write_row("A$i",\@temp);
}
$worksheet1->freeze_panes(1,0);
$worksheet2->freeze_panes(1,0);
$workbook->close();
