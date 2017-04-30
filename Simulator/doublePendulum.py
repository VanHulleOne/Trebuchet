# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 21:54:18 2017

@author: Myself
"""

import sympy as sym

m, l, th1, thd1, thdd1, th2, thd2, thdd2 = sym.symbols('m l theta_1 \\dot{\\theta_1} \\ddot{\\theta_1} theta_2 \\dot{\\theta_2} \\ddot{\\theta_2}')
g = sym.symbols('g')
e1 = 2*m*l*thdd1 + m*l*thdd2*sym.cos(th2-th1)
e2 = m*l*thd2**2*sym.sin(th2-th1)-2*m*g*sym.sin(th1)
e3 = l*thdd2 + l*thdd1*sym.cos(th2-th1)
e4 = -l*thd1**2*sym.sin(th2-th1)-g*sym.sin(th2)

e = [e1, e2, e3, e4]
f1 = sym.Eq(e2 - e1)
f2 = sym.Eq(e4 - e3)
res = sym.solve([f1, f2], [thdd1, thdd2])
