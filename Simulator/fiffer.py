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
import testFuncs as tf

from scipy import optimize

sym.init_printing()

X, Y = 0, 1 # For matrix indexes

t, Mc, Mp, Ma, Ma_perMeter, Ia, CWdrop, la, ls, r_cam, theta_Ai, g, space = sym.symbols('t M_c M_p M_a Ma_perMeter I_a CW_drop l_a l_s r_cam theta_Ai g space')

la_init = 2.9
ls_init = 3.2

constants = {Mc:320, Mp:4.5, Ma_perMeter:12, CWdrop:3.66, r_cam:1.15,
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
        self.createLagrangian()
        self.createGroundLagrangian()
        self.la = None
        self.ls = None
        self.endRange((la_init, ls_init))
        #self.printFuncs()

    def createGroundLagrangian(self):
        print('Create Ground Lagrangian')
        
        Ma = Ma_perMeter*la
        Ia = Ma*la**2/3
        
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
        acc = LM.rhs()
        Ycdd = acc[1]
        print('Solved')
        for key, value in locals().items():
            if key == 'thSdd' or key == 'slingTensionY':
                ycdd = sym.diff(Yc,t,t)
                n = ('Yct', 'Ycdt', 'Ycddt', 'thSt', 'thSdt', 'l_a', 'l_s')
                tD = [Yc, Ycd, ycdd, thS, thSd]
                self.ground_functions[key] = lambdify_helper(tD, value, n)
                self.ground_symbolicEquations[key] = value
                continue
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
        
        LM = mech.LagrangesMethod(L, (Yc, thS))#, forcelist=[Ycd*0.1,])
        
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

    def printFuncs(self, eq_Dict=None):
        if eq_Dict is None:
            eq_Dict = self.symbolicEquations
        with open('trebfunctions.m', 'w') as f:
            f.write('classdef trebfunctions \n')
            f.write('methods(Static) \n')
            for key, value in eq_Dict.items():
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
            
    #@profile        
    def endRange(self, lengths):
        la, ls = lengths
        self.la, self.ls = lengths
        dt = 0.001
        end = constants[CWdrop]*(-0.96)
        maxSteps = 2000
        
        yc = [0]
        ycd = [0]
        ycdd_list = []
        ths = [0]
        thsd = [0]
        thsdd_list = []
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
            #ycdd = functions['Ycdd'](y, yd, th, thd, la, ls)
            if switched:
                ycdd = tf.ycdd(y, yd, th, thd, la, ls)
                thsdd = tf.thSdd(y, yd, th, thd, la, ls)
                #thsdd = functions['thSdd'](y, yd, th, thd, la, ls)
            else:
                ycdd = tf.g_ycdd(y, yd, th, thd, la, ls)
                thsdd = functions['thSdd'](y, yd, ycdd, th, thd, la, ls)
                
                slingTensionY = functions['slingTensionY'](y, yd, ycdd, th, thd, la, ls)
                armTipPosY = functions['armTipPosY'](y, yd, th, thd, la, ls)
                
                if (slingTensionY > -constants[g]*constants[Mp]
                    or la + armTipPosY >= ls):
                    functions = self.functions
                    switched = True
                    self.projLiftIndex = i
                    #print('Yc:', y)
            
            yc.append(y + yd*dt + 0.5*ycdd*dt**2)
            ycd.append(yd + ycdd*dt)
            ycdd_list.append(ycdd)
            
            ths.append(th + thd*dt + 0.5*thsdd*dt**2)
            thsd.append(thd + thsdd*dt)
            thsdd_list.append(thsdd)
            
            projVelX = functions['projVelX'](yc[-1], ycd[-1], ths[-1], thsd[-1], la, ls)
            projVelY = functions['projVelY'](yc[-1], ycd[-1], ths[-1], thsd[-1], la, ls)
            
            if projVelY > 0 and projVelX >= projVelY:
                break
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
            #print('VelMag:', velMagnitude, ' CosTheta:', cosTheta)
        
        
        # s = xi + vi*t + 1.2at^2
        # 0 = g/2*t^2 + vi*t + xi
        hangTime = (-projVelY - (projVelY**2 - 4*1*-9.8/2)**0.5)/(-9.8)
        dist = -projVelX*hangTime
        self.Yc_array = np.array(yc)
        self.Ycd_array = np.array(ycd)
        self.Ycdd_array = np.array(ycdd_list)
        self.thS_array = np.array(ths)
        self.thSd_array = np.array(thsd)
        self.thSdd_array = np.array(thsdd_list)
        self.time_array = np.array(time)
        return dist
    

    
def opti(self, la = la_init, ls = ls_init):
    # 'Nelder-Mead' 'TNC' bounds=((1,7),(1,10)),
    res = optimize.minimize(self.endRange, [la, ls], method = 'Nelder-Mead',  tol=1.0)
    print(res)

def calcData(self):
    length = len(self.time_array)
    self.Vp = np.zeros(length)
    self.Tp = np.zeros(length)
    self.projPos = np.zeros((2, length))
    self.armTipPos = np.zeros((2, length))
    self.T = np.zeros(length)
    self.V = np.zeros(length)
    firstSlice = slice(self.projLiftIndex)
    secondSlice = slice(self.projLiftIndex, length)
    
    for _slice, funcs in zip((firstSlice, secondSlice),
                             (self.ground_functions, self.functions)):
        data = (self.Yc_array[_slice],
            self.Ycd_array[_slice], 
            self.thS_array[_slice], 
            self.thSd_array[_slice], 
            self.la, self.ls)
        self.Vp[_slice] = funcs['V_p'](*data)
        self.Tp[_slice] = funcs['T_p'](*data)
        pos = funcs['projPos'](*data)
        self.projPos[X,_slice] = pos[X]
        self.projPos[Y,_slice] = pos[Y]
        armPos = funcs['armTipPos'](*data)
        self.armTipPos[X,_slice] = armPos[X]
        self.armTipPos[Y,_slice] = armPos[Y]
        self.T[_slice] = funcs['T'](*data)
        self.V[_slice] = funcs['V'](*data)

    
def printFunc(func):
    func = func.subs({a:b for a,b in zip(timeDependent, dummy_timedependent)})
    func = func.subs(constants)
    funcStr = lambdastr(dummy_timedependent, func)
    funcStr = funcStr.replace('**', '^')
    print(funcStr)
    
S1 = Simulator()
#S1.endRange((2.9,3.2))
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

