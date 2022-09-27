import sys

values = []
with open(sys.argv[1], 'r') as f:
    for line in f :
        info = line.strip().split('\t')
        value = float(info[int(sys.argv[2])])
        values.append(value)
f.close()
print(round(sum(values)))