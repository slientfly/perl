#-*- coding:utf-8 -*-
#encode utf-8
def writr_file(txt):
    f=open('./new_file.xls','a')
    #txt=','.join(txt)
    f.write(txt)
    f.close()
with open('/disk2/mygeno/chenyj/Python/test/SNP_INDEL.xls','r') as f:
    #for line_1 in f.readlines()[0:1] :
    #   writr_file(line_1)
    #   print line_1
#   for line_last in f.readlines()[1:]:
    #writr_file(line_last)
#       print line_last
    lines=f.readlines()
    list_1=lines[0:1]
    list_1=list_1[0]
    writr_file(list_1)
    for line in lines[1:]:
        #writr_file(line)
        line=line.split('\t')
        if line[27]<0.3 and line[29]<5:
            writt_filr(line)

