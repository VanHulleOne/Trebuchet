# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 07:56:54 2017

@author: Myself
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 08:34:51 2017

@author: Myself
"""

# Fiffer Simulation
import sympy as sym
from sympy.physics import mechanics as mech
import numpy as np
import matplotlib.pyplot as plt

sym.init_printing()

X, Y = 0, 1 # For matrix indexes

t, Mc, Mp, Ma, Ia, CWdrop, la, ls, r_cam, theta_Ai, g, space = sym.symbols('t M_c M_p M_a I_a CW_drop l_a l_s r_cam theta_Ai g space')
arm = 2.65 #3.178
sling = 4.23 #4.551
constants = {Mc:320, Mp:4.5, Ma:25, Ia:14, CWdrop:3.66, la:arm, ls:sling, r_cam:1,
             theta_Ai:-sym.pi/2, g:-9.8, space:0.5}

Yc = sym.Function('Yc')(t)
thS = sym.Function('theta_S')(t)
Ycd = sym.diff(Yc, t)
thSd = sym.diff(thS, t)

timeDependent = [Yc, Ycd, thS, thSd]

Yct, Ycdt, thSt, thSdt = sym.symbols('Yct, Ycdt, thSt, thSdt')
dummy_timedependent = (Yct, Ycdt, thSt, thSdt)

def lambdify_helper(variables, func, names=None, subConsts=True):
    if names is None:
        names = ['_dummy_' + str(i) for i in range(len(variables))]
    subs_dict = {var:name for var, name in zip(variables, names)}
    subs_func = func.subs(subs_dict)
    if subConsts:
        subs_func = subs_func.subs(constants)
    return sym.lambdify(names, subs_func, 'numpy', dummify=False)
    
print('Create Lagrangian')

Xcam = -CWdrop/2 - r_cam - space
Ycam = 0

Xpb = -CWdrop/2 + Yc**2/4
Ypb = Yc/2

theta_arm = theta_Ai - (sym.sqrt((Xcam - Xpb)**2 + (Ycam - Ypb)**2) - r_cam - space)/r_cam

armTipPos = la*sym.Matrix([sym.cos(theta_arm), sym.sin(theta_arm)])

Xproj = la*sym.cos(theta_arm) + ls*sym.cos(thS)
Yproj = la*sym.sin(theta_arm) + ls*sym.sin(thS)

projPos = sym.Matrix([Xproj, Yproj])

projVel = projPos.diff(t)
projAcc = projPos.diff(t, t)
projSpeedSq = projVel[0]**2 + projVel[1]**2

"""
L = T-V
where:
T = Kenetic Energy
V = Potential Energy
"""

V_c = Mc*-g*(Yc+CWdrop)
V_armcg = Ma*-g*armTipPos[Y]/2
V_p = Mp*-g*Yproj
V = V_c+V_armcg+V_p

T_c = Mc*Ycd**2/2
T_a = Ia*sym.diff(theta_arm,t)**2/2
T_p = Mp*projSpeedSq/2
T = T_c+T_a+T_p

L = T - V

LM = mech.LagrangesMethod(L, (Yc, thS))

eom = LM.form_lagranges_equations()
print('Solve for acc')
acc = LM.rhs()

print('Solved')

print('Create funcs')

Ycdd_func = lambdify_helper(timeDependent, acc[2])
thSdd_func = lambdify_helper(timeDependent, acc[3])
T_func = lambdify_helper(timeDependent, T, dummy_timedependent)
V_func = lambdify_helper((Yc, thS), V, (Yct, Ycdt))

projVel_func = lambdify_helper(timeDependent, projVel)

print('Start Loop')
dt = 0.01
end = constants[CWdrop]*(-0.96)
maxSteps = 2000

yc = [0]
ycd = [0]
ths = [0]
thsd = [0]
time = [0]



i = 0
while(i < maxSteps and yc[-1] > end):
    i += 1
    time.append(i*dt)
    y = yc[-1]
    yd = ycd[-1]
    th = ths[-1]
    thd = thsd[-1]
    ycdd = Ycdd_func(y, yd, th, thd)
    thsdd = thSdd_func(y, yd, th, thd)
    
    yc.append(y + yd*dt + 0.5*ycdd*dt**2)
    ycd.append(yd + ycdd*dt)
    ths.append(th + thd*dt + 0.5*thsdd*dt**2)
    thsd.append(thd + thsdd*dt)



    
Xp_func = sym.lambdify((Yc, thS),Xproj.subs(constants), 'numpy')
Yp_func = sym.lambdify((Yc, thS),Yproj.subs(constants), 'numpy')

projPos_func = sym.lambdify((Yc, thS),projPos.subs(constants), 'numpy')
projVel_func = lambdify_helper(timeDependent, projVel)

thA_func = sym.lambdify(Yc, theta_arm.subs(constants), 'numpy')

ycA = np.array(yc)
thsA = np.array(ths)
ycdA = np.array(ycd)
thsdA = np.array(thsd)

X = Xp_func(ycA, thsA)
Y = Yp_func(ycA, thsA)

T_array = T_func(ycA, ycdA, thsA, thsdA)
V_array = V_func(ycA, thsA)

totalEnergy = T_array + V_array

vel = projVel_func(yc[-1], ycd[-1], ths[-1], thsd[-1])
hangTime = (-vel[1,0] - (vel[1,0]**2 - 4*5*-9.8/2)**0.5)/(-9.8)
dist = vel[0,0]*hangTime
print('Vel:', vel)
print('Dist:', dist)
    
plt.figure(1)
plt.plot(X,Y)
plt.axis('equal')
plt.xlabel('X position [m]')
plt.ylabel('Y position [m]')
plt.title('X-Y Position Projectile')

plt.figure(2)
plt.plot(time, T_array, 'r', time, V_array, 'g', time, totalEnergy, 'b')
plt.xlabel('Time [sec]')
plt.ylabel('Energy [Joules]')
plt.title('System Energy vs Time')
plt.legend(['Kinetic', 'Potential', 'Total'], loc='best')

