import matplotlib
matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

fname = 'dt%d.csv'
x_i, y_i = np.loadtxt('initial.csv', delimiter=' ', usecols=(0,1), unpack=True)
plt.scatter(x_i, y_i, c='r', marker='o', s=50)

# for i in range(1,30):
    # x_f, y_f = np.loadtxt(fname%i, delimiter=' ', usecols=(0,1), unpack=True)
    # plt.scatter(x_f, y_f, c='b', marker='o')
# plt.show()

fig = plt.figure()
ax1 = plt.axes(xlim=(0, 100), ylim=(0,100))
line, = ax1.plot([], [], lw=1)
ax1.set_xlabel(r'Position $x(t)$')
ax1.set_ylabel(r'Position $y(t)$')

plotlays, plotcols = [3], ["black","red", "green"]
lines = []
for index in range(3):
    lobj = ax1.plot([],[],lw=1,color=plotcols[index])[0]
    lines.append(lobj)

def init():
    for line in lines:
        line.set_data([],[])
    return lines

x1,y1 = [],[]
x2,y2 = [],[]
x3,y3 = [],[]

def animate(i):
    if i is 0:
        x_f, y_f = np.loadtxt('initial.csv', delimiter=' ', usecols=(0,1), unpack=True)
    else:
        x_f, y_f = np.loadtxt(fname%(i), delimiter=' ', usecols=(0,1), unpack=True)

    x = x_f[0]
    y = y_f[0]
    x1.append(x)
    y1.append(y)

    x = x_f[1]
    y = y_f[1]
    x2.append(x)
    y2.append(y)

    x = x_f[2]
    y = y_f[2]
    x3.append(x)
    y3.append(y)

    xlist = [x1, x2, x3]
    ylist = [y1, y2, y3]

    for lnum,line in enumerate(lines):
        line.set_data(xlist[lnum], ylist[lnum]) 

    return lines

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=30, interval=100, blit=True)

anim.save('3body.mp4')
