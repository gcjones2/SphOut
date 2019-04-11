#Model for spherical outflow
#Allows for aperture, opening angle, vacuous shell
#vTEN: includes EFACT

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import math
from scipy import signal

def makesphout(dr,irrad,vmax,FWHMkms,gamma,aper,conangle,vellist,tau0,EFACT,zred):

	#Assumptions
	innerrad=1.

	#Modeling constants
	aper=aper*(innerrad+dr)
	pi=np.pi
	numv=len(vellist)	#Number of velocity bins
	numr=100	#Number of radius bins
	#Velocity bin edges
	velbinedges=np.linspace(min(vellist)-np.abs(vellist[1]-vellist[0])/2.,max(vellist)+np.abs(vellist[numv-1]-vellist[numv-2])/2.,numv+1)	
	fullE=np.zeros((numr,numv)) #Spectrum for each angle
	spec=np.zeros(numv) #Integrated spectrum
	nphot=np.zeros(numr) #Total number of photons in a ring made by each angle
	
	dx=(innerrad+dr)/numr #Plane-based radius cellsize
	bigrad=np.linspace(dx/2.,innerrad+dr-dx/2.,numr)		 #Vertical distance array

	Espec=np.zeros(numv)
	Aspec=np.zeros(numv)

	#Returns unprojected velocity of cell at given radius from center - OK.
	def getvel(R):
		return -1*vmax*((R/innerrad)**gamma)

	#Returns fraction of flux blocked by cell at radius r - OK.
	def getblock(R,gamma):
		temp=1-np.exp(-tau0*((R/innerrad)**((-2*gamma)-1)))
		if temp>1: return 1
		elif temp<0: return 0
		else: return temp

	#Returns array of x-coordinates of a slice
	def sliceit(x0,x1):
		numdx=(x1-x0)/dx
		return np.linspace(x0,x1,numdx)

	#Make Gaussian for convolving
	def makegaus():
		sigmakms=FWHMkms/(2*np.sqrt(2*np.log(2)))
		x=np.arange(-3*FWHMkms, 3*FWHMkms, np.abs(vellist[1]-vellist[0]))
		gaussian=np.zeros(len(x))
		for i in range(len(x)): gaussian[i] = np.exp((-1/2.)*(x[i]/sigmakms)**2)
		return gaussian

	def getfrac(N1v,velspace,j,k):
	    try: 
	        dv=min(N1v[k+1]-N1v[k],N1v[k]-N1v[k-1])
	    except IndexError: 
	        try:
	            dv=min(N1v[k+1]-N1v[k],N1v[k+2]-N1v[k+1])
	        except IndexError: 
	            dv=min(N1v[k-2]-N1v[k-1],N1v[k-1]-N1v[k])
	    alpha=N1v[k]-dv/2.; beta=N1v[k]+dv/2.
	    a=velspace[j]; b=velspace[j+1]
	    if beta<a or alpha>b: return 0
	    else: return (min(b,beta)-max(a,alpha))/(beta-alpha)

#-------------

	contcounts=0.
	for i in range(numr):
		
		r=bigrad[i]

		#Check cone
		if r>(innerrad+dr)*np.cos(conangle): print 'CON_BREAK'; break

		#Check aperture
		if r>aper: print 'APE_BREAK'; break

		velbincounts=np.zeros(numv)

		#Inner section - all absorption
		if r<irrad:
			xshellmin=max(np.sqrt(innerrad**2-r**2),r*np.tan(conangle))
			xshellmax=np.sqrt((innerrad+dr)**2-r**2)	
			lilrad=sliceit(xshellmin,xshellmax)
			nphot[i]=pi*2*r/dx
			contcounts+=nphot[i]
			currentfrac=np.ones(len(velbincounts))
			for j in range(len(lilrad)):
				xx=lilrad[j]
				radx=np.sqrt(r**2+xx**2)
				velx=getvel(radx)*(xx/radx)
				for k in range(numv-1):	
					if velx>=velbinedges[k] and velx<velbinedges[k+1]:
						velbincounts[k]-=getblock(radx,gamma)*nphot[i]*currentfrac[k]
						Aspec[k]-=getblock(radx,gamma)*nphot[i]*currentfrac[k]
						currentfrac[k]*=(1-getblock(radx,gamma))
						if currentfrac[k]<0: currentfrac[k]=0

						
			for j in range(numv):	
				fullE[i,j]+=velbincounts[j] #Central absorption

	#-------------
	#Sum up the emission
	for i in range(numr): 
		for j in range(numv):
			spec[j]+=fullE[i,j]

	#Normalize by the continuum
	for j in range(numv):
		spec[j]=spec[j]/contcounts
		Espec[j]=Espec[j]/contcounts+1
		Aspec[j]=Aspec[j]/contcounts+1

	if contcounts==0: print dr,vmax,tau0,FWHMkms

	#Convolve
	spec2=np.zeros(len(spec))
	spec2=signal.convolve(spec, makegaus(), mode="same")
	
	#Add emission profile
	N1=np.genfromtxt('AVGPRO.txt')
	N1v=N1[:,0]
	N1F=N1[:,1]


	for i in range(len(N1v)):
		N1v[i]=(3E+5)*(((1+5.656)*(1+N1v[i]/(3E+5))/(1+zred))-1)
	
	for j in range(len(N1v)):
		for k in range(len(spec)-1):
			if getfrac(N1v,velbinedges,k,j)>0:
				spec2[k]+=EFACT*N1F[j]*getfrac(N1v,velbinedges,k,j)
				Espec[k]+=EFACT*N1F[j]*getfrac(N1v,velbinedges,k,j)
	
	#Add normalized continuum level
	for k in range(len(spec2)): 
		spec2[k]+=1.

	#Return spectrum
	return [list(reversed(spec2)),list(reversed(Espec)),list(reversed(Aspec))]

