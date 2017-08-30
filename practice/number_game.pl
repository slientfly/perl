#usr/bin/perl -w
#usr/bin/perl -w
use strict;
use utf8;
binmode(STDIN, ':encoding(utf8)');
binmode(STDOUT, ':encoding(utf8)');
binmode(STDERR, ':encoding(utf8)');
my $number=int(1+rand 100);
print "请输入一个1至100的数字：";
while (chomp(my $line=<STDIN>)){
    if ($line eq 'quit' ||$line eq 'exit'){
        last;
        }
    elsif($line<$number && $line>0){
        print "Too low\n";
        }
    elsif ($line>$number&& $line<100){
        print "Too hight\n";
        }
    elsif ($line==$number){
        print "binge\n";
        last;
        }
    else{
        print "请输入1至100以内的数！\n";
        }
}
