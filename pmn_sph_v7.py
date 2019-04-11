#!/usr/bin/env python
#PMN code for fitting absorption line with SPHOUT model
#v5: with EFACT
#v6: updating continuum values
#v7: Adding in all three redshifts

import os
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
import json
import sys
import scipy.stats, scipy
sys.path.insert(0, '/Users/garethjones/PyMultiNest')
import pymultinest
import scipy
from scipy.stats import norm
import matplotlib.pyplot as plt
from sphout_vTEN import makesphout
#from sphout_vELEVEN import makesphout #Possible continuum normilzation issue
from uncertainties import unumpy as unp

#
whichone='B723N'       #B701N / B701SE / B723N / B723SE
z123='3'
fitwhich=[0,1,1,1,0]       #log(Shell thickness), vmax, FWHM, tau0, EFACT
dr=10.                #Shell thickness
vmax=230.              #Outward speed of shell
FWHMkms=610.           #FWHM of Gaussian to convolve with spectrum
tau0=0.00338
EFACT=0.0

#
irrad=1.0          #[1,0.5,0.1]
#n=2               #[2,3]
gamma=-1         #0 for cont v, -2 for const density ===================
aper=1           #[1,0.5]
conangle=0.       #pi/[2,3,6]

if z123=='1':
    zred=5.656
elif z123=='2':
    zred=5.652
elif z123=='3':
    zred=5.654
else: print "BAD REDSHIFT"

def frtove(freq,trans): 
	if trans=='HD': restfr=2674.986094/(1+zred)
	if trans=='CO21': restfr=230.538/(1+zred)
	if trans=='H2O3': restfr=2196.345756/(1+zred)
	if trans=='H2O4': restfr=2264.14965/(1+zred)
	temp=(3*10**5)*((restfr/freq)-1)
	return temp

def uniformprior(cube,howmany,a,b):
    return cube[howmany]*(b-a)+a

def loguniformprior(cube,howmany,a,b):
    return 10**(cube[howmany]*(b-a)+a)

if whichone=='B723N':
    N1=np.genfromtxt('/Volumes/Johto/2016/b7_old/calibrated_orig/B7_23_N_z.txt',skip_header=9)
    TRANSITION='H2O3'
    RMS=0.0003*np.sqrt(57/22)
    cont=11.0E-3
    xlimits=[-1000,1950]
else: print "BAD SOURCE"

N1x=N1[:,0]           #Frequency [GHz]
N1y=1+N1[:,1]/cont 	  #1+Line/Cont
RMSy=np.zeros(len(N1x)) 

#Look at intial guess
N1v=np.zeros(len(N1x))
for i in range(len(N1v)): 
    #N1v[i]=N1x[i]
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
#plt.axvline(vmax)
#plt.axvline(-1*vmax)
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


"""
#Write components
f=open("STUFF.txt",'w')
f.close()
f=open("STUFF.txt",'w')
for i in range(len(x2)):
    f.write(str(x2[i])+" "+str(np.log10(freefree(x[i])*sfrterm))+" "+str(np.log10(synch(x[i],lfnth,alpha)*sfrterm))+" "+str(np.log10(dust(x[i],td,beta)*sfrterm))+"\n")
f.close()
#plt.title(whichone+" "+str(fitwhich))
"""
    
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




