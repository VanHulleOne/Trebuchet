
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

from fiffer import Simulator

def twoGraphs(S1):
    totalEnergy = S1.T + S1.V
    time = S1.time_array
    plt.figure(1)
    plt.plot(S1.projPos[X],S1.projPos[Y])
    plt.axis('equal')
    plt.xlabel('X position [m]')
    plt.ylabel('Y position [m]')
    plt.title('X-Y Position Projectile')
    
    plt.figure(2)
    plt.plot(time, S1.T, 'r', time, S1.V, 'g', time, totalEnergy, 'b')
    plt.xlabel('Time [sec]')
    plt.ylabel('Energy [Joules]')
    plt.title('System Energy vs Time')
    plt.legend(['Kinetic', 'Potential', 'Total'], loc='best')

#S1 = Simulator()

# Make data.
X = np.arange(1.5, 5, 0.05)
Y = np.arange(1.5, 7, 0.05)
Z = np.zeros([len(X)*len(Y)])
print('Total:', len(Z))

i = 0
for y in Y:
    for x in X:
        Z[i] = -S1.endRange([x,y])
        if i%100 == 0:
            print(i)
        i += 1
        

X, Y = np.meshgrid(X, Y)
Z = Z.reshape(X.shape)

xi, yi = np.unravel_index(Z.argmax(), Z.shape)
print('Arm:', X[xi,yi])
print('Sling:', Y[xi,yi])
print('Range:', Z[xi,yi])
#R = np.sqrt(X**2 + Y**2)
#Z = np.sin(R)
#Z = X*Y
fig = plt.figure()
ax = fig.gca(projection='3d')
# Plot the surface.
ax.plot_surface(X,Y,Z, cmap=cm.jet, rstride=2, cstride=2)

# Customize the z axis.
ax.set_zlim(100, 500)
ax.set_xlabel('Arm Length [m]')
ax.set_ylabel('Sling Length [m]')
ax.set_zlabel('Range [m]')
#ax.zaxis.set_major_locator(LinearLocator(50))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
