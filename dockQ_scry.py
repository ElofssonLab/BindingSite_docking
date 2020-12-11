import sys

dockqs = [[line.split[-2].split('/')[-1], line.split()[1]] for line in open(sys.argv[1])]

dic = {}
for code, score in dockqs:
    code = code.split('_')[0]
    if code not in dic: dic[code] = [float(score)]
    else: dic[code].append(float(score))

for code in dic: print('Max dockQ score for complex {}: {}'.format(code, max(dic[code])))
