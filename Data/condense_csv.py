inFile = open('od_large.csv','r')
outFile = open('od_py_large.csv','w')

prev_line = ' '
count = 0

for line in inFile:
    if line == prev_line:
        count = count + 1
    else:
        output = prev_line.split('\n')[0] +","+str(count)+"\n"
        outFile.write(output)
        count = 1
        prev_line = line

output = prev_line.split('\n')[0] +","+str(count)+"\n"
outFile.write(output)

outFile.close()
inFile.close()
