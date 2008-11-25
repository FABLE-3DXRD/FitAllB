import numpy as n
from xfab import tools
import check_input


def	vars(inp):
	"""
	Calculated experimental errors of gexp for all peaks
	Jette Oddershede, July 2008
	"""
    
#	id = inp.id[0][0]
#	print inp.w[id],inp.dety[id],inp.detz[id],inp.Sww[id],inp.Syy[id],inp.Szz[id],inp.values['x0'],inp.values['y0'],inp.values['z0'],inp.values['L']
	check_input.set_globals(inp)
	inp.vars = []
	for i in range(inp.no_grains):
		inp.vars.append([])
		for j in range(inp.nrefl[i]):
			if i+1 in inp.fit['skip']:
				inp.vars[i].append([1,1,1])
			else:
				Sgg = error(inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],\
							inp.Sww[inp.id[i][j]],inp.Syy[inp.id[i][j]],inp.Szz[inp.id[i][j]],\
							inp.values['wx'],inp.values['wy'],inp.values['tx'],inp.values['ty'],inp.values['tz'],\
							inp.values['py'],inp.values['pz'],inp.values['cy'],inp.values['cz'],inp.values['L'],\
							inp.values['x%s' %i],inp.values['y%s' %i],inp.values['z%s' %i])
				inp.vars[i].append([Sgg[0,0],Sgg[1,0],Sgg[2,0]])


def error(w,dety,detz,Sww,Syy,Szz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z):
	"""
	Function for calculating the experimental errors of gexp for a single peak
	This is done by propagating the measured errors of the observables w, dety and detz
	Input:  observables: w,dety,detz
			corresponding experimental errors: Sww,Syy,Szz
			global parameters: wx,wy,tx,ty,tz,py,pz,cy,cz,L
			grain position: x,y,z
	Output:	vars, vector of propagated errors
	Jette Oddershede, July 2008
	"""
	
	hh = 1e-3
	dw = (gexp(w+hh/2,dety,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z)-gexp(w-hh/2,dety,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z))/hh
	dy = (gexp(w,dety+hh/2,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z)-gexp(w,dety-hh/2,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z))/hh
	dz = (gexp(w,dety,detz+hh/2,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z)-gexp(w,dety,detz-hh/2,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z))/hh
	vars = dw*dw*Sww + dy*dy*Syy + dz*dz*Szz
#	print 'Sww',sum(dw*dw),Sww,sum(dw*dw)*Sww
#	print 'Syy',sum(dy*dy),Syy,sum(dy*dy)*Syy
#	print 'Szz',sum(dz*dz),Szz,sum(dz*dz)*Szz, '\n'
	return vars 

	
def gexp(w,dety,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z):
	"""
	Experimental g-vector calculation
	Input:  observables: w,dety,detz
			global parameters: wx,wy,tx,ty,tz,py,pz,cy,cz,L
			grain position: x,y,z
	Output:	vars, vector of propagated errors
	Jette Oddershede, July 2008
	"""
	
	Omega = tools.quart2Omega(w,wx*n.pi/180,wy*n.pi/180)
	R = tools.detect_tilt(tx,ty,tz)
	d = n.dot(R,n.array([[0],[(dety+cy)*py],[(detz-cz)*pz]])) 
	d = d + n.array([[L],[0],[0]]) - n.dot(Omega,n.array([[x],[y],[z]]))
	gexp = n.dot(n.linalg.inv(Omega),(d/n.sqrt(n.sum(d**2)) - n.array([[1],[0],[0]])))
	return gexp 
	
	
