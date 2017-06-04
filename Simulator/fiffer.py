# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 08:34:51 2017

@author: Myself
"""

# Fiffer Simulation
from collections import namedtuple as nt
import sympy as sym
from sympy.physics import mechanics as mech
import matplotlib.pyplot as plt
import numpy as np

from scipy import optimize

import fifferequations as feq

sym.init_printing()

X, Y = 0, 1
JIT, VECTORIZE = 0, 1

Data = nt('Data', 'Yc Ycd Ycdd thS thSd thSdd time l_arm l_sling')
data = Data._make([None]*9)

t, Mc, Mp, Ma, Ma_perMeter, Ia, CWdrop, la, ls, r_cam, theta_Ai, g, space = sym.symbols('t M_c M_p M_a Ma_perMeter I_a CW_drop l_arm l_sling r_cam theta_Ai g space')

la_init = 3.35
ls_init = 3.65

constants = {Mc:320, Mp:4.5, Ma_perMeter:12, CWdrop:3.66, r_cam:1.15,
             theta_Ai:-sym.pi/2, g:-9.8, space:0.5}

Yc = sym.Function('Yc')(t)
thS = sym.Function('theta_S')(t)
Ycd = sym.diff(Yc, t)
Ycdd = sym.diff(Yc,t,t)
thSd = sym.diff(thS, t)
thSdd = sym.diff(thS,t,t)

timeDependent = [Yc, Ycd, Ycdd, thS, thSd, thSdd]
timeIndepArgs = (la, ls)

Yct, Ycdt, Ycddt, thSt, thSdt, thSddt = sym.symbols('Yc, Ycd, Ycdd, thS, thSd, thSdd')
dummy_timedependent = (Yct, Ycdt, Ycddt, thSt, thSdt, thSddt)

timeSubsDict = {a:b for a,b in zip(timeDependent, dummy_timedependent)}

def numbify(func, name, _type):
    func, args = subFunc(func)
    signature = 'float64(' + 'float64,'*len(args) +')'
    if _type == JIT:
        lines = ['@jit(' + signature + ', ']
    elif _type == VECTORIZE:
        lines = ['@vectorize([' + signature + '], ']
    lines.append('nopython=True, cache=True)\n')
    lines.append('def ' + name + '(')
    lines.extend(str(s) + ', ' for s in args)
    lines.append('):\n')
    lines.append('\treturn ')
    lines.append(str(func)+'\n\n')
    return ''.join(lines)
        
    
def subFunc(func):
    func = func.subs(timeSubsDict)
    funcSet = set(sym.preorder_traversal(func))
    args = [arg for arg in dummy_timedependent + timeIndepArgs if arg in funcSet]
    func = func.subs(constants)
    return func, args    

class EquationBuilder:
    def __init__(self):
        self.symbolicEquations = {}
        self.accelerations = {}
        self.funcArgsMap = {}
        self.createLagrangian()
        self.formAccelerations()
        self.buildFuncArgMap()
        
    def __getattr__(self, attr):
        cls = type(self)
        if attr in self.funcArgsMap:
            return getattr(feq, attr)(*[getattr(data, arg) for arg in self.funcArgsMap[attr]])
        
        msg = '{.__name__!r} object has no attribute {!r}'
        raise AttributeError(msg.format(cls, attr))
        
    def plotter(self, *args, xAxis='time', ylabel=None, title=None):
        if ylabel is None:
            ylabel = ', '.join(args)
        if title is None:
            title = ylabel + ' vs ' + xAxis
        if xAxis == 'time':
            xData = data.time
        else:
            xData = getattr(self, xAxis)
        plt.figure(1)
        for arg in args:
            try:
                yData = getattr(self, arg)
            except Exception:
                yData = getattr(data, arg)
            plt.plot(xData, yData)
        plt.xlabel(xAxis)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.legend([arg for arg in args], loc='best')
    
    def buildFuncArgMap(self):
        for name, func in self.symbolicEquations.items():
            _, args = subFunc(func)
            self.funcArgsMap[name] = [str(arg) for arg in args]
        
    def formAccelerations(self):
        LM = mech.LagrangesMethod(self.symbolicEquations['L'],
                                  (Yc, thS),
                                  hol_coneqs=[self.symbolicEquations['projPosY']])
        
        LM.form_lagranges_equations()
        print('Solve for ground acc')
        acc = LM.rhs()
        self.accelerations['Ycdd_g'] = acc[2]
        self.accelerations['thSdd_g'] = acc[3]
        
        LM = mech.LagrangesMethod(self.symbolicEquations['L'], (Yc, thS))
        
        LM.form_lagranges_equations()
        print('Solve for general acc')
        acc = LM.rhs()
        self.accelerations['Ycdd'] = acc[2]
        self.accelerations['thSdd'] = acc[3]

    def createLagrangian(self):
        print('Create Lagrangian')
        
        Ma = Ma_perMeter*la
        Ia = Ma*la**2/3
        
        Xcam = -CWdrop/2 - r_cam - space
        Ycam = r_cam - r_cam # To make a symbolic zero
        
        Xpb = -CWdrop/2 + Yc**2/4
        Ypb = Yc/2
        
        theta_arm = theta_Ai - (sym.sqrt((Xcam - Xpb)**2 + (Ycam - Ypb)**2) - r_cam - space)/r_cam
        
        armTipPos = la*sym.Matrix([sym.cos(theta_arm), sym.sin(theta_arm)])
        
        projPos = sym.Matrix([la*sym.cos(theta_arm) + ls*sym.cos(thS),
                              la*sym.sin(theta_arm) + ls*sym.sin(thS)])
        
        projVel = projPos.diff(t)
        projAcc = projPos.diff(t, t)

        projSpeedSq = projVel[0]**2 + projVel[1]**2
        
        slingTension = sym.sqrt(projAcc[0]**2 + projAcc[1]**2)*Mp
        slingTensionY = -slingTension*sym.sin(thS)
        
        """
        L = T-V
        where:
        T = Kenetic Energy
        V = Potential Energy
        """
        
        V_c = Mc*-g*(Yc+CWdrop)
        V_armcg = Ma*-g*armTipPos[Y]/2
        V_p = Mp*-g*projPos[Y]
        V = V_c+V_armcg+V_p
        
        T_c = Mc*Ycd**2/2
        T_a = Ia*sym.diff(theta_arm,t)**2/2
        T_p = Mp*projSpeedSq/2
        T = T_c+T_a+T_p
        
        totalEnergy = T + V
        
        L = T - V
        
        for key, value in locals().items():
            if value is self:
                continue
            try:
                self.symbolicEquations[key+'X'] = value[X]
                self.symbolicEquations[key+'Y'] = value[Y]
            except Exception:
                self.symbolicEquations[key] = value
            
#@profile       
def endRange(la, ls=None, setData=False):
    if ls is None:
        la, ls = la
    dt = 0.001
    endHeight = constants[CWdrop]*(-0.96)
    maxSteps = 2000
    
    yc_array = np.zeros(maxSteps, dtype=np.float64)
    ycd_array = np.zeros(maxSteps, dtype=np.float64)
    ycdd_array = np.zeros(maxSteps, dtype=np.float64)
    ths_array = np.zeros(maxSteps, dtype=np.float64)
    thsd_array = np.zeros(maxSteps, dtype=np.float64)
    thsdd_array = np.zeros(maxSteps, dtype=np.float64)
    time_array = np.zeros(maxSteps, dtype=np.float64)
    
    arrays = (yc_array, ycd_array, ycdd_array,
              ths_array, thsd_array, thsdd_array,
              time_array)
    
    Ycdd_func = feq.Ycdd_g
    thSdd_func = feq.thSdd_g
    liftoff = False
    i = 0
    while(i < maxSteps and yc_array[i] > endHeight):
        time_array[i] = i*dt
        yc = yc_array[i]
        ycd = ycd_array[i]
        thS = ths_array[i]
        thSd = thsd_array[i]
        
        projVelX = feq.projVelX(yc, ycd, thS, thSd, la, ls)
        projVelY = feq.projVelY(yc, ycd, thS, thSd, la, ls)
        
        if projVelY > 0 and projVelX >= projVelY:
            break
        
        ycdd = Ycdd_func(yc, ycd, thS, thSd, la, ls)
        thSdd = thSdd_func(yc, ycd, thS, thSd, la, ls)
        
        if not liftoff:
            #projAccY = feq.projAccY(yc, ycd, ycdd, thS, thSd, thSdd, la, ls)   
            slingTensionY = feq.slingTensionY(yc, ycd, ycdd, thS, thSd, thSdd, la, ls)            
            armTipPosY = feq.armTipPosY(yc, la)
            
            if (slingTensionY > -constants[g]*constants[Mp] or la + armTipPosY >= ls):
                """
                If the upward acceleration is greater than gravity or
                the arm tip is higher than the sling length the projectile will
                leave the ground.
                """
                Ycdd_func = feq.Ycdd
                thSdd_func = feq.thSdd
                liftoff = True
        
        ycdd_array[i] = ycdd
        thsdd_array[i] = thSdd
        
        i += 1
        
        yc_array[i] = yc + ycd*dt + 0.5*ycdd*dt**2
        ycd_array[i] = ycd + ycdd*dt        
        
        ths_array[i] = thS + thSd*dt + 0.5*thSdd*dt**2
        thsd_array[i] = thSd + thSdd*dt
    else: # Excecuted if the break statement is not hit.
        '''
        If the sling is too long it will not come around fast enough and the
        while loop will exit because the CW has fallen all the way down. The
        real machine would continue running and possibly pull the CW back up a
        little before launching, but the simulation does not work that way so
        this is an attempt to handle that and smooth out the contour plot so the
        optimization can be more robust.
        '''
        cos45 = np.cos(np.pi/4)
        ideal = np.array([cos45, cos45])
        actual = np.array([projVelX, projVelY])
        velMagnitude = np.linalg.norm(actual)
        cosTheta = (np.dot(ideal, actual)/velMagnitude+1)/2
        projVelX, projVelY = ideal*velMagnitude*cosTheta
    
    
    # s = xi + vi*t + 1.2at^2
    # 0 = g/2*t^2 + vi*t + xi
    hangTime = (-projVelY - (projVelY**2 - 4*1*-9.8/2)**0.5)/(-9.8)
    dist = -projVelX*hangTime
    #projSpeed = np.sqrt(projVelY**2 + projVelX**2)
    if setData:
        global data
        data = Data._make([np.array(a[:i-1]) for a in arrays] + [la, ls])
    return dist


    
def opti(la = la_init, ls = ls_init):
    # 'Nelder-Mead' 'TNC' bounds=((1,7),(1,10)),'Nelder-Mead'
    res = optimize.minimize(endRange, [la, ls], method = 'Powell', tol=0.1)
    endRange(*res.x, setData=True)
    return res

def writeFunc(eqBuilder, fileName = 'fifferequations.py'):
    with open(fileName, 'w') as f:
        f.write('from math import cos, sin, sqrt, pi\n')
        f.write('from numba import jit, vectorize, float64\n')
        f.write('\n')
        for name, func in eqBuilder.symbolicEquations.items():
            print(name)
            f.write(numbify(func, name, VECTORIZE))
        
        for name, func in eqBuilder.accelerations.items():
            print(name)
            f.write(numbify(func, name, JIT))

if __name__ == '__main__':
    #print(endRange(2.7,3.1, setData=True))
    e1 = EquationBuilder()
    #opti()
