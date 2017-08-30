use strict;
use Encode;
use Excel::Writer::XLSX;  
my $workbook = Excel::Writer::XLSX->new("xuesheng.xlsx");  #引号中为生成的excel的名称
#生成Excel表  
my $worksheet = $workbook->add_worksheet(decode('utf8',"学生信息"));  #引号中为excel工作簿中表的名称
#$xlsheet->freeze_panes(1, 0); #冻结首行
my $format1 = $workbook->add_format();#增加一种格式
$format1->set_bold();              #设置粗体
$format1->set_color( 'red' );      #设置颜色
$format1->set_align( 'center' );   #设置对齐方式（此处为居中）
$worksheet->set_column('A:A',14); #设置列的宽度 
$worksheet->set_column('B:B',10);  
$worksheet->set_column('C:C',10);  
$worksheet->set_column('D:D',20);
$worksheet->set_column('E:E',24);
#写表头（格式是使用上面添加的表头格式）  
$worksheet->write(0,0,decode('utf8',"学号"),$format1);  
$worksheet->write(0,1,decode('utf8',"姓名"),$format1);  
$worksheet->write(0,2,decode('utf8',"班级"),$format1);
$worksheet->write(0,3,decode('utf8',"身份证号"),$format1);
$worksheet->write(0,4,decode('utf8',"银卡号"),$format1);
my $format2 = $workbook->add_format();
$format2->set_align( 'center' );
open F,"</disk2/mygeno/ouyk/perl/8.28/xinxi.txt";
my $row=1;
while (<F>){
  chomp;
my @xinxi=split '\t';
for my $i (0..4){
$worksheet->write($row,$i,decode('utf8',$xinxi[$i]),$format2);
}
$row ++;
}
close F; 


