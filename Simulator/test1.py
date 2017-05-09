# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 08:58:40 2017

@author: Myself
"""

import numpy as np

aref = np.array([1,0])
for angle in range(0,314,10):
    angle /= 100
    atarg = np.array([np.cos(angle), np.sin(angle)])
    costheta = np.dot(aref, atarg)
    deg = np.arccos(costheta)/(2*np.pi)*360
    print('Atarg:', atarg, 'Cos_theta:', costheta, 'Angle:', deg)