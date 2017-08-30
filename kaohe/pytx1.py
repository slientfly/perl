import commands
import re
c=commands.getoutput('qnodes')
#pattern=re.compile(r'cu.*\n.*\n.*\n.*')
#w=pattern.findall(str(c))
pattern=re.findall(r'(cu.*\n.*?state.*\n.*?np.*\n.*?properties.*)\n.*?ntype.*',str(c),re.M)
for x in pattern:
    print x
