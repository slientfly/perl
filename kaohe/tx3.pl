#usr/bin/perl 
#my $a='@ST-E00243:250:HKWLFALXX:4:1101:25469:1801';
##$a=~/@.*?(\d+):(\d+):.*?(\d+):(\d+):(\d+):(\d+)/;
##$new1="$1"+11111;
##$new2="$2"+111;
##$new3="$3"+1;
##$new4="$4"+1111;
##$new5="$5"+11101;
##$new6="$6"+1111;
my $c;
print"\@ST-E$new1:$new2:HKWLFALXX:$new3:$new4:$new5:$new6";
while($a=~/(.*?)(\d)/){
my $b=$2+1;
my $c=$c.($1.$b);
}
print "$c";

