import re

inputfile = "count"
outputfile = "final"

out = open(outputfile,'w')
out.close()

with open(inputfile) as infile:
        for line in infile:
                count = re.split(':',line)[1]
                if int(count) > 2:
                        with open(outputfile,'a') as out:
                                out.write(re.split(':',line)[0]+'\t'+str(count))
