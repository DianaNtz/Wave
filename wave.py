"""Runge Kutta method for the 2D wave equation"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import imageio.v2 as imageio
nfilenames = []
afilenames = []

#some initial values
nx=60
x0=0
xfinal=1
dx=(xfinal-x0)/(nx-1)
x=np.linspace(x0,xfinal,nx)

ny=60
y0=0
yfinal=1
dy=(yfinal-y0)/(ny-1)
y=np.linspace(y0,yfinal,ny)

nt=150
t0=0
tfinal=2/np.sqrt(2)
dt=(tfinal-t0)/(nt-1)
t=0
