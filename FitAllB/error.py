

import numpy as n
from xfab import tools
from . import check_input
from . import reject
from copy import deepcopy

def vars_scale(inp):
    """
    Scale vars so that the expectation value for each refinement 
    becomes number of observations (3*nrefl) 
    Should be number of observations - number of parameters,
    but this is much more difficult and the scaling is only preliminary
    to ensure reasonable tolerances and up
     
    Jette Oddershede, December 2008
    """

    print('\n',inp.fit['goon'])
    print('present residuals')
    for i in range(inp.no_grains):
        print(i+1,n.sum(inp.residual[i]),len(inp.residual[i]),n.sum(inp.residual[i])/len(inp.residual[i]))
    # cakk if reject.residual is necessary to return to same reference point in as before first refinement
    reject.residual(inp,inp.fit['limit'][0],only=[])
    print('original reference residuals')
    for i in range(inp.no_grains):
        print(i+1,n.sum(inp.residual[i]),len(inp.residual[i]),n.sum(inp.residual[i])/len(inp.residual[i]))

    check_input.set_globals(inp)
    inp.vars = []
    for i in range(inp.no_grains):
        inp.vars.append([])
        fcn_current = n.sum(inp.residual[i]) 
        fcn_expected = 3.*len(inp.residual[i]) 
        for j in range(inp.nrefl[i]):
            if i+1 in inp.fit['skip']:
                inp.vars[i].append([1,1,1])
            else:
                Sgg = error(inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],\
                            inp.Sww[inp.id[i][j]]*fcn_current/fcn_expected,\
                            inp.Syy[inp.id[i][j]]*fcn_current/fcn_expected,\
                            inp.Szz[inp.id[i][j]]*fcn_current/fcn_expected,\
                            inp.values['wx'],inp.values['wy'],inp.values['tx'],inp.values['ty'],inp.values['tz'],\
                            inp.values['py'],inp.values['pz'],inp.values['cy'],inp.values['cz'],inp.values['L'],\
                            inp.values['x%s' %i],inp.values['y%s' %i],inp.values['z%s' %i])
#                inp.vars[i].append([Sgg[0,0],Sgg[1,0],Sgg[2,0]])
                inp.vars[i].append([fcn_current/fcn_expected,fcn_current/fcn_expected,fcn_current/fcn_expected])

                

def	vars(inp):
    """
    Calculated experimental errors of gexp for all peaks
    Jette Oddershede, July 2008
    """
    
    check_input.set_globals(inp)
    # start preparing for volume weighted errors
    volavg = []
    try:
        for i in range(inp.no_grains):
            if len(inp.volume[i]) > 0:
#                volavg.append(sum(inp.volume[i])/len(inp.volume[i])) #mean
                volavg.append(reject.median(inp.volume[i]))         #median
            else:
                volavg.append(0)
    except:
        volavg = [1.]*inp.no_grains
    volavg_significant = deepcopy(volavg)
    for i in range(inp.no_grains-1,-1,-1):
        if volavg_significant[i] == 0:
            volavg_significant.pop(i)
    if len(volavg_significant) > 0:
        volmedian = reject.median(volavg_significant)
    else:
        volmedian = 1.
    # end preparing for volume weighted errors
    
    inp.vars = []
    for i in range(inp.no_grains):
        inp.vars.append([])
        for j in range(inp.nrefl[i]):
#            if i+1 in inp.fit['skip']:
#                inp.vars[i].append([1,1,1])
#            else:
#                Sgg = error(inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],\
#                            inp.Sww[inp.id[i][j]],inp.Syy[inp.id[i][j]],inp.Szz[inp.id[i][j]],\
#                            inp.values['wx'],inp.values['wy'],inp.values['tx'],inp.values['ty'],inp.values['tz'],\
#                            inp.values['py'],inp.values['pz'],inp.values['cy'],inp.values['cz'],inp.values['L'],\
#                            inp.values['x%s' %i],inp.values['y%s' %i],inp.values['z%s' %i])
#                inp.vars[i].append([Sgg[0,0],Sgg[1,0],Sgg[2,0]]) # propagated errors
#                inp.vars[i].append([Sgg[0,0]+1e-9,Sgg[1,0]+1e-9,Sgg[2,0]+1e-9]) # propagated errors + constant contribution
            if volavg[i] != 0:
                inp.vars[i].append([4e-8*volmedian/volavg[i],4e-8*volmedian/volavg[i],1e-8*volmedian/volavg[i]]) # no error propagation
            else:
                inp.vars[i].append([4.,4.,1.]) # no error propagation


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
	
	hh = 1e-6
	dw = (gexp(w+hh/2,dety,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z)-gexp(w-hh/2,dety,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z))/hh
	dy = (gexp(w,dety+hh/2,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z)-gexp(w,dety-hh/2,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z))/hh
	dz = (gexp(w,dety,detz+hh/2,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z)-gexp(w,dety,detz-hh/2,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z))/hh
	vars = dw*dw*Sww + dy*dy*Syy + dz*dz*Szz
#	print 'Sww',sum(dw*dw),Sww,sum(dw*dw)*Sww
#	print 'Syy',sum(dy*dy),Syy,sum(dy*dy)*Syy
#	print 'Szz',sum(dz*dz),Szz,sum(dz*dz)*Szz, '\n'
	return vars 
#	return n.array([[1.],[1.],[1.]])
	
def gexp(w,dety,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z):
	"""
	Experimental g-vector calculation
	Input:  observables: w,dety,detz
			global parameters: wx,wy,tx,ty,tz,py,pz,cy,cz,L
			grain position: x,y,z
	Output:	vars, vector of propagated errors
	Jette Oddershede, July 2008
	"""
	
	Omega = tools.quart_to_omega(w,wx*n.pi/180,wy*n.pi/180)
	R = tools.detect_tilt(tx,ty,tz)
	d = n.dot(R,n.array([[0],[(dety+cy)*py],[(detz-cz)*pz]])) 
	d = d + n.array([[L],[0],[0]]) - n.dot(Omega,n.array([[x],[y],[z]]))
	gexp = n.dot(n.transpose(Omega),(d/n.sqrt(n.sum(d**2)) - n.array([[1],[0],[0]])))
	return gexp 
	
	
