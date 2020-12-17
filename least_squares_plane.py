import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def lsp(point_list):
    tmp_A = []
    tmp_b = []
    for i in range(len(point_list)):
        tmp_A.append([point_list[i][0], point_list[i][1], 1])
        tmp_b.append(point_list[i][2])
    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)
    return (A.T * A).I * A.T * b

def shortest_distance(x1, y1, z1, a, b, c, d):  
      
    d = abs((a * x1 + b * y1 + c * z1 + d))  
    e = (math.sqrt(a * a + b * b + c * c)) 
    print("Perpendicular distance is"), d/e 


plist = [[1,2,3], [3,3,3], [5,2,1], [9,2,11], [20,31,2]]

xs = [1,3,5,9,20]
ys = [2,3,2,2,31]
zs = [3,3,1,11,2]

fit = lsp(plist)

plt.figure()
ax = plt.subplot(111, projection='3d')
ax.scatter(xs, ys, zs, color='b')

# plot plane
xlim = ax.get_xlim()
ylim = ax.get_ylim()
X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
                  np.arange(ylim[0], ylim[1]))
Z = np.zeros(X.shape)
for r in range(X.shape[0]):
    for c in range(X.shape[1]):
        Z[r,c] = fit[0] * X[r,c] + fit[1] * Y[r,c] + fit[2]
ax.plot_wireframe(X,Y,Z, color='k')

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()

