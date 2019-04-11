#!/usr/bin/env python
#PyMultinest code for fitting absorption line with SPHOUT model
#Allows for fitting of shell thickness, maximum outflow velocity, velocity dispersion, optical depth,
# and a normalization constant for emission

import os
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import json
import sys
import scipy.stats, scipy
import pymultinest
import scipy
from scipy.stats import norm
import matplotlib.pyplot as plt
from sphout import makesphout
from uncertainties import unumpy as unp

#SOURCE PARAMETERS
filebase='B723N'       #Base file name      
zred=5.656             #Source redshift
TRANSITION='H2O3'      #Observed transition (see frtove below)
RMS=0.0003             #Observed RMS noise level per channel
cont=11.0E-3           #Observed continuum flux density
xlimits=[-1000,1950]   #Velocity limits for plotting

#FITTING PARAMETERS
fitwhich=[0,1,1,1,0]   #Flags for fitting: log(Shell thickness), vmax, FWHM, tau0, EFACT
dr=10.                 #Estimate for shell thickness
vmax=230.              #Estimate for outward speed of shell
FWHMkms=610.           #Estimate for FWHM of Gaussian to convolve with spectrum
tau0=0.00338           #Estimate for optical depth
EFACT=0.0              #Estimate for emission profile normalization

#ASSUMPTIONS
irrad=1.0          #Ratio of inner shell radius to IR source radius
gamma=-1           #Velocity & density power index: 0 for constant v, -2 for constant density
aper=1             #Fractional aperture of simulated observation
conangle=0.        #Opening angle of conical outflow (0 for sphere)

#Convert a frequency to an optical velocity
def frtove(freq,trans): 
	if trans=='HD': restfr=2674.986094/(1+zred)
	if trans=='CO21': restfr=230.538/(1+zred)
	if trans=='H2O3': restfr=2196.345756/(1+zred)
	if trans=='H2O4': restfr=2264.14965/(1+zred)
	temp=(3*10**5)*((restfr/freq)-1)
	return temp

#Uniform prior for PMN
def uniformprior(cube,howmany,a,b):
    return cube[howmany]*(b-a)+a

#log-uniform prior for PMN
def loguniformprior(cube,howmany,a,b):
    return 10**(cube[howmany]*(b-a)+a)

#Import data
N1=np.genfromtxt(filebase+'.txt',skip_header=9)
N1x=N1[:,0]           #Frequency [GHz]
N1y=1+N1[:,1]/cont 	  #1+Line/Cont
RMSy=np.zeros(len(N1x)) 

#Look at intial guess
N1v=np.zeros(len(N1x))
for i in range(len(N1v)): 
    N1v[i]=frtove(N1x[i],TRANSITION)
    RMSy[i]=RMS/cont
ANS_SPEC,ESP,ASP=makesphout(dr,irrad,vmax,FWHMkms,gamma,aper,conangle,N1v,tau0,EFACT,zred)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.plot(N1v,ANS_SPEC,color='b')
plt.step(N1v,ESP,color='g')
plt.step(N1v,ASP,color='c')
plt.step(N1v,N1y,color='r')
plt.fill_between(N1v,1+RMSy,1-1*RMSy,alpha=0.2,color='r')
plt.xlabel(r'v [km s$^{-1}$]')
plt.xlim(xlimits[0],xlimits[1])
plt.ylabel('Normalized Flux Density')
plt.show()

#--------------------------------
def prior(cube, ndim, nparams):
    howmany=0
    if fitwhich[0]:
        cube[howmany]=loguniformprior(cube,howmany,-2,1)
        howmany=howmany+1
    if fitwhich[1]:
        cube[howmany]=uniformprior(cube,howmany,50,1000)
        howmany=howmany+1
    if fitwhich[2]:
        cube[howmany]=uniformprior(cube,howmany,50,1000)
        howmany=howmany+1
    if fitwhich[3]:
        cube[howmany]=loguniformprior(cube,howmany,-7,1)
        howmany=howmany+1
    if fitwhich[4]:
        cube[howmany]=loguniformprior(cube,howmany,-8,0)
        howmany=howmany+1

    
def loglike(cube,ndim,nparams):
    totdif=0.0
    howmany=0
    if fitwhich[0]:
        global dr
        dr=cube[howmany]
        howmany=howmany+1
    if fitwhich[1]:
        global vmax
        vmax=cube[howmany]
        howmany=howmany+1
    if fitwhich[2]:
        global FWHMkms
        FWHMkms=cube[howmany]
        howmany=howmany+1
    if fitwhich[3]:
        global tau0
        tau0=cube[howmany]
        howmany=howmany+1
    if fitwhich[4]:
        global EFACT
        EFACT=cube[howmany]
        howmany=howmany+1
    modelspec=makesphout(dr,irrad,vmax,FWHMkms,gamma,aper,conangle,N1v,tau0,EFACT,zred)[0]
    for i in range(len(N1v)):
        totdif+=((N1y[i]-modelspec[i])/RMSy[i])**2+np.log10(2*np.pi*RMSy[i])**2
    return (-1/2)*totdif
    
#--------------------------------
parameters=[]
if fitwhich[0]:
    parameters.append('dr')
if fitwhich[1]:
    parameters.append('vmax')
if fitwhich[2]:
    parameters.append('FWHMkms')
if fitwhich[3]:
    parameters.append('tau0')
if fitwhich[4]:
    parameters.append('EFACT')

n_params = len(parameters)
datafile=whichone+'_sph'

# run MultiNest
pymultinest.run(loglike, prior, n_params, outputfiles_basename=datafile + '_1_', resume = False, verbose = True,max_iter=10000)

skippy=2*(n_params+3)-2 #n=5,14
ANS=np.genfromtxt(whichone+'_sph_1_stats.dat',skip_header=4,skip_footer=skippy,delimiter='   ')

howmany=0
if fitwhich[0]:
 dr=ANS[howmany,1]
 ddr=ANS[howmany,2]
 howmany=howmany+1
if fitwhich[1]:
 vmax=ANS[howmany,1]
 dvmax=ANS[howmany,2]
 howmany=howmany+1
if fitwhich[2]:
 FWHMkms=ANS[howmany,1]
 dFWHMkms=ANS[howmany,2]
 howmany=howmany+1
if fitwhich[3]:
 tau0=ANS[howmany,1]
 dtau0=ANS[howmany,2]
 howmany=howmany+1
if fitwhich[4]:
 EFACT=ANS[howmany,1]
 dEFACT=ANS[howmany,2]
 howmany=howmany+1

ANS_SPEC=makesphout(dr,irrad,vmax,FWHMkms,gamma,aper,conangle,N1v,tau0,EFACT,zred)[0]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.plot(N1v,ANS_SPEC,color='b')
plt.step(N1v,N1y,color='r')
plt.fill_between(N1v,1+RMSy,1-1*RMSy,alpha=0.2,color='r')
plt.xlabel(r'v [km s$^{-1}$]')
plt.xlim(xlimits[0],xlimits[1])
plt.ylabel('Normalized Flux Density')
plt.show()

if fitwhich[0]:
    dr=unp.uarray(dr,ddr)
if fitwhich[1]:
    vmax=unp.uarray(vmax,dvmax)
if fitwhich[2]:
    FWHMkms=unp.uarray(FWHMkms,dFWHMkms)
if fitwhich[3]:
    tau0=unp.uarray(tau0,dtau0)
if fitwhich[4]:
    EFACT=unp.uarray(EFACT,dEFACT)
    
#======================

a = pymultinest.Analyzer(n_params = n_params, outputfiles_basename=datafile + '_1_')
s = a.get_stats()

import json
with open('%sparams.json' % a.outputfiles_basename, 'w') as f:
	json.dump(parameters, f, indent=2)

with open('%sstats.json' % a.outputfiles_basename, mode='w') as f:
	json.dump(s, f, indent=2)
print()
print("-" * 30, 'ANALYSIS', "-" * 30)
print("Global Evidence:\n\t%.15e +- %.15e" % ( s['nested sampling global log-evidence'], s['nested sampling global log-evidence error'] ))

print '-----'
print 'dr = '+str(dr)
print 'vmax = '+str(vmax)+" km/s"
print "FWHMkms = "+str(FWHMkms)
print "tau0 = "+str(tau0)
print "EFACT = "+str(EFACT)
