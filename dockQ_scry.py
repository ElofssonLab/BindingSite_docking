import sys
import pandas
pandas.set_option('display.max_rows', 10)

dockqs = [[line.split()[-2].split('/')[-1], line.split()[1]] for line in open(sys.argv[1])]

dic = {}
for code, score in dockqs:
    code = code.split('_')[0]
    if code not in dic: dic[code] = [float(score)]
    else: dic[code].append(float(score))

for code in dic: dic[code] = max(dic[code])
df = pandas.DataFrame(dic, index=[0]).T
pandas.set_option('display.max_rows', df.shape[0]+1)
print (df.sort_values(0, ascending=False))
