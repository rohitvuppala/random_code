#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 22:53:26 2021

author: Rohit K S S Vuppala
         Graduate Student, 
         Mechanical and Aerospace Engineering,
         Oklahoma State University.

@email: rvuppal@okstate.edu

"""

import numpy as np
import matplotlib.pyplot as plt

xmin = 0.0
xmax = 500.0
nx   = 501
dx   = (xmax-xmin)/(nx-1) 
cfl  = 0.3
g    = 9.8

niter= 1500

#Set up initial stuff   
x = np.linspace(xmin,xmax,nx)
h = np.zeros(nx)
for i in range(nx):
    h[i] = 50.0 - 45.0*np.tanh((x[i]-300)/80)
    
h = h 

eta = np.zeros((nx))
m   = np.zeros((nx))

#initialize
for i in range(nx):
    #if(x[i]<1000):
    #    eta[i] = 1.0
    if(x[i] < 175 and x[i]>125):
        eta[i] = 5*np.exp(-(x[i]-150.0)**2/80)

m   = 100*eta

nplot = 100
#Time-stepping
for n in range(niter-1):
    hmax = np.amax(h)
    dt   = cfl*dx/np.sqrt(2*g*hmax)
    alp  = dt/dx
    #BCs
    eta[0] = 0.0
    eta[nx-1] = 0.0
    m[0]    = 0.0
    m[nx-1] = 0.0
    
    #Calculate
    for i in range(1,nx-2):
        eta_old = eta
        m_old = m
        
        eta[i] = eta_old[i] - 0.5*alp*(m_old[i+1]-m_old[i-1]) 
        m[i]   = m_old[i] -  0.5*alp*g*(eta_old[i] + h[i])*(eta_old[i+1] - eta_old[i-1])
    
    if(n%(nplot) == 0):
        print('iteration-'+str(n))
        plt.figure()
        plt.title("Tsunami Wave Travel")
        plt.ylim(-100,30)
        plt.xlabel("travel-direction")
        plt.ylabel("altitude")
        plt.plot(x,eta[:],'-',label=str(n))
        plt.plot(x,-h,'k--')
        plt.savefig("iteration"+str(n))
#%%
#plt.figure()
