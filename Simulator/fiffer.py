# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 08:34:51 2017

@author: Myself
"""

# Fiffer Simulation
import sympy as sym
from sympy.physics import mechanics as mech
from sympy.utilities.lambdify import lambdastr
import numpy as np
import matplotlib.pyplot as plt

from scipy import optimize

sym.init_printing()

X, Y = 0, 1 # For matrix indexes

t, Mc, Mp, Ma, Ma_perMeter, Ia, CWdrop, la, ls, r_cam, theta_Ai, g, space = sym.symbols('t M_c M_p M_a Ma_perMeter I_a CW_drop l_a l_s r_cam theta_Ai g space')

la_init = 2.65
ls_init = 4.2

constants = {Mc:320, Mp:4.5, Ma_perMeter:12, CWdrop:3.66, r_cam:1,
             theta_Ai:-sym.pi/2, g:-9.8, space:0.5}

Yc = sym.Function('Yc')(t)
thS = sym.Function('theta_S')(t)
Ycd = sym.diff(Yc, t)
thSd = sym.diff(thS, t)

timeDependent = [Yc, Ycd, thS, thSd]

Yct, Ycdt, thSt, thSdt = sym.symbols('Yct, Ycdt, thSt, thSdt')
dummy_timedependent = (Yct, Ycdt, thSt, thSdt)
names = ('Yct', 'Ycdt', 'thSt', 'thSdt', 'l_a', 'l_s')


timeSubsDict = {a:b for a,b in zip(timeDependent, dummy_timedependent)}

def lambdify_helper(variables, func, names=None, subConsts=True):
    if names is None:
        names = ['_dummy_' + str(i) for i in range(len(variables))]
    subs_dict = {var:name for var, name in zip(variables, names)}
    subs_func = func.subs(subs_dict)
    if subConsts:
        subs_func = subs_func.subs(constants)
    return sym.lambdify(names, subs_func, 'numpy', dummify=False)

class Simulator:
    def __init__(self):
        self.symbolicEquations = {}
        self.functions = {}
        self.ground_symbolicEquations = {}
        self.ground_functions = {}
        self.distances = []
        #self.createLagrangian()
        #self.printFuncs()

    def createGroundLagrangian(self):
        print('Create Ground Lagrangian')
        
        Ma = Ma_perMeter*la
        Ia = Ma*la/12
        
        Xcam = -CWdrop/2 - r_cam - space
        Ycam = r_cam - r_cam # To make a symbolic zero
        
        Xpb = -CWdrop/2 + Yc**2/4
        Ypb = Yc/2
        
        theta_arm = theta_Ai - (sym.sqrt((Xcam - Xpb)**2 + (Ycam - Ypb)**2) - r_cam - space)/r_cam
        
        armTipPos = la*sym.Matrix([sym.cos(theta_arm), sym.sin(theta_arm)])
        armTipPosX = armTipPos[X]
        armTipPosY = armTipPos[Y]
        
        thS = -sym.asin(la*(1+sym.sin(theta_arm))/ls)
        thSd = sym.diff(thS,t)
        thSdd = sym.diff(thS,t,t)
        
        Xproj = la*sym.cos(theta_arm) + ls*sym.cos(thS)
        Yproj = la*sym.sin(theta_arm) + ls*sym.sin(thS)
        
        projPos = sym.Matrix([Xproj, Yproj])
        
        projVel = projPos.diff(t)
        projVelX = projVel[X]
        projVelY = projVel[Y]
        projAcc = projPos.diff(t, t)
        slingTension = sym.sqrt(projAcc[0]**2 + projAcc[1]**2)*Mp
        slingTensionY = -slingTension*sym.sin(thS)
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
        
        LM = mech.LagrangesMethod(L, (Yc,))
        
        LM.form_lagranges_equations()
        print('Solve for acc')
        vel = LM.rhs()
        Ycdd = sym.diff(vel[1],t)
        print('Solved')
        for key, value in locals().items():
            try:                
                self.ground_functions[key] = lambdify_helper(timeDependent, value, names)
                self.ground_symbolicEquations[key] = value
            except Exception as e:
                print('Key:', key)        

    def createLagrangian(self):
        print('Create Lagrangian')
        
        Ma = Ma_perMeter*la
        Ia = Ma*la/12
        
        Xcam = -CWdrop/2 - r_cam - space
        Ycam = r_cam - r_cam # To make a symbolic zero
        
        Xpb = -CWdrop/2 + Yc**2/4
        Ypb = Yc/2
        
        theta_arm = theta_Ai - (sym.sqrt((Xcam - Xpb)**2 + (Ycam - Ypb)**2) - r_cam - space)/r_cam
        
        armTipPos = la*sym.Matrix([sym.cos(theta_arm), sym.sin(theta_arm)])
        armTipPosX = armTipPos[X]
        armTipPosY = armTipPos[Y]
        
        Xproj = la*sym.cos(theta_arm) + ls*sym.cos(thS)
        Yproj = la*sym.sin(theta_arm) + ls*sym.sin(thS)
        
        projPos = sym.Matrix([Xproj, Yproj])
        
        projVel = projPos.diff(t)
        projVelX = projVel[X]
        projVelY = projVel[Y]
        projAcc = projPos.diff(t, t)
        slingTension = sym.sqrt(projAcc[0]**2 + projAcc[1]**2)*Mp
        slingTensionY = -slingTension*sym.sin(thS)
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
        
        LM.form_lagranges_equations()
        print('Solve for acc')
        acc = LM.rhs()
        Ycdd = acc[2]
        thSdd = acc[3]
#        self.symbolicEquations['Ycdd'] = acc[2]
#        self.symbolicEquations['thSdd'] = acc[3]
#        self.functions['Ycdd'] = lambdify_helper(timeDependent, acc[2], names)
#        self.functions['thSdd'] = lambdify_helper(timeDependent, acc[3], names)
        print('Solved')
        for key, value in locals().items():
            try:                
                self.functions[key] = lambdify_helper(timeDependent, value, names)
                self.symbolicEquations[key] = value
            except Exception as e:
                print('Key:', key)
#                print(e)
#                break

    def printFuncs(self):
        with open('trebfunctions.m', 'w') as f:
            f.write('classdef trebfunctions \n')
            f.write('methods(Static) \n')
            for key, value in self.symbolicEquations.items():
                if key == 'acc' or key == 'projVel' or key == 'projAcc':
                    continue
                f.write('function [' + key + '] = ' + key + '(Yct, Ycdt, thSt, thSdt, l_a, l_s) \n')
                func = value.subs(timeSubsDict)
                func = func.subs(constants)
                funcStr = lambdastr(names, func)
                funcStr = funcStr.replace('**', '^')
                f.write('\t' + key + ' = ' + funcStr[36:] + ';\n')
                f.write('end\n')
                f.write('\n\n')
            f.write('end\n')
            f.write('end\n')
#print('Start Loop')
#

#lengths_array = []
    def endRange(self, lengths):
    #    print('la:', la, ' ls:', ls)
        la, ls = lengths
        dt = 0.01
        end = constants[CWdrop]*(-0.96)
        maxSteps = 2000
        
        yc = [0]
        ycd = [0]
        ths = [0]
        thsd = [0]
        time = [0]
        
        functions = self.ground_functions
        switched = False
        i = 0
        while(i < maxSteps and yc[-1] > end):
            i += 1
            time.append(i*dt)
            y = yc[-1]
            yd = ycd[-1]
            th = ths[-1]
            thd = thsd[-1]
            ycdd = functions['Ycdd'](y, yd, th, thd, la, ls)
            thsdd = functions['thSdd'](y, yd, th, thd, la, ls)
            
            slingTensionY = functions['slingTensionY'](y, yd, th, thd, la, ls)
            armTipPosY = functions['armTipPosY'](y, yd, th, thd, la, ls)
            if (not switched
                and (slingTensionY > -constants['g']*constants['Mp']
                or la + armTipPosY >= ls)):
                functions = self.functions
            
            yc.append(y + yd*dt + 0.5*ycdd*dt**2)
            ycd.append(yd + ycdd*dt)
            ths.append(th + thd*dt + 0.5*thsdd*dt**2)
            thsd.append(thd + thsdd*dt)
            
            projVelX = functions['projVelX'](yc[-1], ycd[-1], ths[-1], thsd[-1], la, ls)
            projVelY = functions['projVelY'](yc[-1], ycd[-1], ths[-1], thsd[-1], la, ls)
            
            if projVelY > 0 and projVelX > projVelY:
                break
        
        # s = xi + vi*t + 1.2at^2
        # 0 = g/2*t^2 + vi*t + xi
        hangTime = (-projVelY - (projVelY**2 - 4*5*-9.8/2)**0.5)/(-9.8)
        dist = -projVelX*hangTime
#        print('dist:', dist, 'hangTime:', hangTime, 'la:', la, 'ls:', ls)
        self.distances.append(dist)
        return dist
    
    def opti(self, la = la_init, ls = ls_init):
        res = optimize.minimize(self.endRange, [la, ls], method = 'TNC', bounds=((1,7),(1,10)), tol=0.01)
        print(res)

#S1 = Simulator()


    
#opti()

def printFunc(func):
    func = func.subs({a:b for a,b in zip(timeDependent, dummy_timedependent)})
    func = func.subs(constants)
    funcStr = lambdastr(dummy_timedependent, func)
    funcStr = funcStr.replace('**', '^')
    print(funcStr)
#Xp_func = sym.lambdify((Yc, thS),Xproj.subs(constants), 'numpy')
#Yp_func = sym.lambdify((Yc, thS),Yproj.subs(constants), 'numpy')
#
#projPos_func = sym.lambdify((Yc, thS),projPos.subs(constants), 'numpy')
#
#thA_func = sym.lambdify(Yc, theta_arm.subs(constants), 'numpy')
#
#ycA = np.array(yc)
#thsA = np.array(ths)
#ycdA = np.array(ycd)
#thsdA = np.array(thsd)
#
#X = Xp_func(ycA, thsA)
#Y = Yp_func(ycA, thsA)
#
#T_array = T_func(ycA, ycdA, thsA, thsdA)
#V_array = V_func(ycA, thsA)
#
#totalEnergy = T_array + V_array
#
#plt.figure(1)
#plt.plot(X,Y)
#plt.axis('equal')
#plt.xlabel('X position [m]')
#plt.ylabel('Y position [m]')
#plt.title('X-Y Position Projectile')
#
#plt.figure(2)
#plt.plot(time, T_array, 'r', time, V_array, 'g', time, totalEnergy, 'b')
#plt.xlabel('Time [sec]')
#plt.ylabel('Energy [Joules]')
#plt.title('System Energy vs Time')
#plt.legend(['Kinetic', 'Potential', 'Total'], loc='best')
