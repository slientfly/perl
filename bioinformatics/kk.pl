use strict;
use warnings;
use Data::Dumper;
while(<DATA>){
        chomp;
        next if /^$/;
        next if /^#/;
        my ($lane,$name,$sample,$special_sigle,$sex,$chip) = split /\s+/;
        unless(-d "/home/pengj/qsub/result1"){
                       `mkdir -p "/home/pengj/qsub/result1"`;
        }
       for($chip){
                for my $i(6..9){
                        if ($chip eq "MJNJK00$i"){
                                $chip="MJNJK002";
                                }
                        }
                for my $j(10..27,32..36){
                        if($chip eq "MJNJK0$j"){
                                $chip="MJNJK002";
                                }
                        }
       if($chip eq 'WJ003'||$chip eq 'MJNJK029'||$chip eq 'MJNJK030'||$chip eq 'MJNJK031'){
                        $chip="CanRisk";
                        } 
       if(-e "/disk2/mygeno/web/110SNP/health/project/$chip/$lane/$sample/$name-$special_sigle.genotype.xlsx"){
               `cp /disk2/mygeno/web/110SNP/health/project/$chip/$lane/$sample/$name-$special_sigle.genotype.xlsx /home/pengj/qsub/result1`;
       }else{
            print "不存在$name样本结果\n";
       }
    }
}
close DATA;

__DATA__

MG-N085         李仁秀  17H001221_RJK011_MPCR   17H001221       男      MJNJK032
MG-N044         2016-R-TJ009    2016-R-TJ009    2016-R-TJ009    男      MJNJK002
MG-N044         2016-R-TJ010    2016-R-TJ010    2016-R-TJ010    男      MJNJK002
MG-N044         2016-R-TJ011    2016-R-TJ011    2016-R-TJ011    男      MJNJK002
MG-N086         栾玉福  17H001142_RJK011_MPCR   17H001142       男      MJNJK032
MG-N109         孙全玉  17M001544_WJ004_MPCR    AA899543713     男      WJ004
MG-N109         曾令琼  17M001542_WJ004_MPCR    AA904918801     女      WJ004
MG-N108         张风雪  17M001539_WJ001_MPCR    AC192430676     女      WJ001
MG-N108         吕明鑫  17M001538_WJ001_MPCR    AC099979158     男      WJ001
MG-N108         石开鹏  17M001527_WJ001_MPCR    AC345378293     男      WJ001
MG-N108         李凤珍  17M001533_WJ001_MPCR    AC004226349     女      WJ001
MG-N106         邵艾    17H000350_MJNJK028_MPCR 17H000350       女      MJNJK028
MG-N106         卢晓平  17H000352_MJNJK028_MPCR 17H000352       男      MJNJK028
MG-N106         林婷    17H000351_MJNJK028_MPCR 17H000351       女      MJNJK028
MG-N106         陈俊英  17H001667       17H001667       女      MJNJK002
MG-N106         张立红  17H001665       17H001665       女      MJNJK002
MG-N106         许春富  17H001668       17H001668       男      MJNJK002
MG-N103         付颖颖  17H150036_MJNJK002_MPCR 17H150036       女      MJNJK002
MG-N103         温广洪  17M000947_MJNJK034_MPCR 17M000947       男      MJNJK034
MG-N103         张乃平  17M001428_WJ001_MPCR    AC655724953     女      WJ001
MG-N103         冷兴贵  17M004799_WJ001_MPCR    AC453505263     男      WJ001
MG-N103         虞治津  17M001423_WJ001_MPCR    AC695494695     男      WJ001
MG-N103         毕恩德  17M004797_WJ001_MPCR    AC416376235     男      WJ001
MG-N103         毕素菊  17M004796_WJ001_MPCR    AC539739734     女      WJ001
#MG-N103         申贞莉  17M001422_WJ001_MPCR    A
