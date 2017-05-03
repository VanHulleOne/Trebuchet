import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import fiffer
from fiffer import Simulator

X,Y = 0,1

#S1 = Simulator()
fiffer.calcData(S1)

fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(7, 6.5)

ax = plt.axes(xlim=(-7, 5), ylim=(-5, 10))

arm = plt.Line2D([],'b-', lw=3)
sling = plt.Line2D([], 'r-', lw=2)
proj = plt.Circle([], 0.35, fc='r')

def init():
    arm.set_data([],[])
    sling.set_data([],[])
    proj.center = (0,0)
    plt.gca().add_line(arm)
    plt.gca().add_line(sling)
    plt.gca().add_patch(proj)
    return arm, sling, proj,

def animate(i):
    sling.set_data((S1.armTipPos[X][i], S1.projPos[X][i]),
                    (S1.armTipPos[Y][i], S1.projPos[Y][i]))
    arm.set_data((0,S1.armTipPos[X][i]), (0, S1.armTipPos[Y][i]))
    proj.center = (S1.projPos[X][i], S1.projPos[Y][i])
    return arm, sling, proj,

anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=len(S1.time_array)*2, 
                               interval=30,
                               blit=True)
#patch = plt.Circle((5, -5), 0.75, fc='y')
#
#def init():
#    patch.center = (5, 5)
#    ax.add_patch(patch)
#    return patch,
#
#def animate(i):
#    x, y = patch.center
#    x = 5 + 3 * np.sin(np.radians(i))
#    y = 5 + 3 * np.cos(np.radians(i))
#    patch.center = (x, y)
#    return patch,
#
#anim = animation.FuncAnimation(fig, animate, 
#                               init_func=init, 
#                               frames=360, 
#                               interval=20,
#                               blit=True)
#
#plt.show()