import numpy as n
import check_input
import write_output
import reject
import fcn
import time
try:
    from iminuit import Minuit
except ImportError:
    from minuit import Minuit
import sys
import logging
from copy import deepcopy
logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')


class fit_minuit():
    def __init__(self,inp):
        self.inp = inp
			
    def refine(self):	
	"""
	Carry out one refinement cycle according to the order given by self.inp.fit['reforder']
	Reject reflections according to self.inp.fit['rej_resmean']
	Print and save refinement and rejection info and parameters
	
	Jette Oddershede, Risoe DTU, May 15 2008
	"""
		
        # initialise
        self.poor_value = []
        self.poor_nrefl = []

 		# create lists of parameters, global and for each grain
        self.globals = ["a","b","c","alpha","beta","gamma","wx","wy","tx","ty","tz","py","pz","cy","cz","L"]
        self.grains = []
        for i in range(self.inp.no_grains):
            self.grains.append(["x%s" %i,"y%s" %i,"z%s" %i,"rodx%s" %i,"rody%s" %i,"rodz%s" %i,
                                "epsaa%s" %i,"epsbb%s" %i,"epscc%s" %i,"epsbc%s" %i,"epsac%s" %i,"epsab%s" %i])

        #refinement update
        reload(fcn)
        self.m = Minuit(fcn.FCN)
        self.m.values = self.inp.values
        self.m.errors = self.inp.errors
        for entries in self.m.fixed:
            self.m.fixed[entries] = True

		# determine whether to refine
        self.ref = False
        if 'grain' in self.inp.fit['goon'] or 'final' in self.inp.fit['goon'] or 'rotpos' in self.inp.fit['goon']:
            self.ref = True
        elif 'start' in self.inp.fit['goon']:
            self.ref = False
        elif 'rod' in self.inp.fit['goon'] and self.inp.fit['rod'] != 0:
            self.ref = True		
        elif 'eps' in self.inp.fit['goon'] and self.inp.fit['eps'] != 0:
            self.ref = True		
        elif 'xyz' in self.inp.fit['goon'] and self.inp.fit['xyz'] != 0:
            self.ref = True		
		
    
		# carry out refinement
        if self.ref == True:
            self.mg = Minuit(fcn.FCNgrain)
            self.mg.values = self.m.values
            self.mg.errors = self.m.errors
            self.mg.fixed = self.m.fixed

            print '\n\n*****Now fitting %s*****' %self.inp.fit['goon']
            print 'newreject_grain', self.inp.fit['newreject_grain']
            # calculate starting values
            g = grain_values(self)
            self.fval = sum(g)
            print '\n%s starting value %e' %(self.inp.fit['goon'],self.fval)
            t1 = time.clock()
            self.mg = Minuit(fcn.FCNgrain)
            self.mg.values = self.m.values
            self.mg.errors = self.inp.errors
            for i in range(self.inp.no_grains):
                    if i+1 in self.inp.fit['skip']:
                        pass
                    elif 'final' in self.inp.fit['goon'] and i+1 not in self.inp.fit['newreject_grain']:
                        pass
                    elif 'rotpos' in self.inp.fit['goon'] and i+1 not in self.inp.fit['newreject_grain'] and abs(self.mg.errors['x%i' %i] - self.inp.param['y_size']/5.) > 1e-3:
                        pass
                    else:	
                        if 'grain' in self.inp.fit['goon']:
                            self.fitgrain(i)
                        elif 'final' in self.inp.fit['goon']:
                            self.fitgrain(i)
                        elif 'eps' in self.inp.fit['goon']:
                            self.fitepsgrain(i)
                        elif 'rod' in self.inp.fit['goon']:
                            self.fitrodgrain(i)
                        elif 'xyz' in self.inp.fit['goon']:
                            self.fitxyzgrain(i)
                        elif 'rotpos' in self.inp.fit['goon']:
                            self.fitrotposgrain(i)
                        if i == 0:
                            print 'Fit %s tolerance %e' %(self.inp.fit['goon'],self.mg.tol)
                        self.mg.values['i'] = i
                        print '\rRefining grain %i' %(i+1),
                        sys.stdout.flush()
                        self.mg.migrad()
                        scale_errors(self,i)
                        #print self.mg.edm, self.mg.ncalls
                        self.m.errors = self.mg.errors
                        write_output.write_cor(self,i)
                        write_output.write_cov(self,i)
                        write_output.write_errors(self,i)
                        self.m.values = self.mg.values
                        g[i] = self.mg.fval
				
            self.time = time.clock()-t1
            print '\nFit %s time %i s' %(self.inp.fit['goon'],self.time)
            self.fval = sum(g)
            print 'Fit %s value %e \n' %(self.inp.fit['goon'],self.fval)
			    
			
            # reject outliers and save cycle info	
            self.m.errors = self.inp.errors
            reject_outliers(self)
            write_output.write_values(self)
            write_output.write_rej(self.inp,message=self.inp.fit['goon'])
            write_output.write_log(self)


        if 'final' in self.inp.fit['goon'] and (self.inp.newreject > 0):
            self.inp.fit['goon'] = 'grain'+ self.inp.fit['goon'][5:]
        elif 'rotpos' in self.inp.fit['goon'] and (self.inp.newreject > 0):
            self.inp.fit['goon'] = 'start'+ self.inp.fit['goon'][6:]
        
		# move onto next refinement given by the reforder list	
        self.inp.fit['goon'] = self.inp.fit['reforder'][self.inp.fit['reforder'].index(self.inp.fit['goon'])+1]
	
        return
        
        

                		
    def fitstart(self):
	"""
	Set tolerance and fixed parameters for preliminary fit of the global parameters
	"""
        self.m.tol = self.inp.fit['tol_start']
        for entries in self.m.fixed:
            if entries[0]=='w' and self.inp.fit['w'] != 0:
                self.m.fixed[entries] = False
            elif entries[0]=='t' and self.inp.fit['tilt'] != 0:
                self.m.fixed[entries] = False
            elif 'p' in entries and len(entries) == 2 and self.inp.fit['pixel'] != 0:
                self.m.fixed[entries] = False
            elif entries[0]=='cy' and self.inp.fit['center'] != 0:
                self.m.fixed[entries] = False
            elif 'L' in entries and self.inp.fit['L'] != 0:
                self.m.fixed[entries] = False

		
    def fitrodgrain(self,i):
	"""
	Set tolerance and fixed parameters for fit of orientations for grain i
	"""
        self.mg.tol = self.inp.fit['tol_rod']
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for angles in self.grains[i]:
            if 'rod' in angles and self.inp.fit['rod'] != 0:
                self.mg.fixed[angles] = False
                self.mg.errors[angles] = self.mg.errors[angles] * 10.


    def fitxyzgrain(self,i):
	"""
	Set tolerance and fixed parameters for fit of positions for grain i
	"""
        self.mg.tol = self.inp.fit['tol_xyz']
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for pos in self.grains[i]:
            if (pos[0] == 'x' or pos[0] == 'y' or pos[0] == 'z') and self.inp.fit['xyz'] != 0:
                self.mg.fixed[pos] = False


    def fitepsgrain(self,i):
	"""
	Set tolerance and fixed parameters for fit of strains for grain i
	"""
        self.mg.tol = self.inp.fit['tol_eps']
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for strain in self.grains[i]:
            if 'eps' in strain and self.inp.fit['eps'] != 0:
                self.mg.fixed[strain] = False

				
    def fitgrain(self,i):
	"""
	Set tolerance and fixed parameters for fit of orientations, positions and strains for grain i
	"""
        
        if self.inp.fit['goon'] == 'grain':
            self.mg.tol = self.inp.fit['tol_grain']
        else:
            self.mg.tol = self.inp.fit['tol_grain']*0.1
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for entries in self.grains[i]:
            if entries[0]=='x' and self.inp.fit['xyz'] != 0 and self.inp.fit['fixx'] == 0:
                self.mg.fixed[entries] = False
            elif entries[0]=='y' and self.inp.fit['xyz'] != 0 and self.inp.fit['fixy'] == 0:
                self.mg.fixed[entries] = False
            elif entries[0]=='z' and self.inp.fit['xyz'] != 0 and self.inp.fit['fixz'] == 0:
                self.mg.fixed[entries] = False
            elif 'eps' in entries and self.inp.fit['eps'] != 0:
                self.mg.fixed[entries] = False
            elif 'rod' in entries and self.inp.fit['rod'] != 0:
                self.mg.fixed[entries] = False
                
                
    def fitrotposgrain(self,i):
        """
        Set tolerance and fixed parameters for fit of orientations and positions for grain i    
        """
        self.mg.tol = self.inp.fit['tol_rotpos']
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for entries in self.grains[i]:
            if entries[0]=='x' and self.inp.fit['xyz'] != 0 and self.inp.fit['constrx'] == 0 and self.inp.fit['fixx'] == 0:
                self.mg.fixed[entries] = False
            elif entries[0]=='y' and self.inp.fit['xyz'] != 0 and self.inp.fit['constry'] == 0 and self.inp.fit['fixy'] == 0:
                self.mg.fixed[entries] = False
            elif entries[0]=='z' and self.inp.fit['xyz'] != 0 and self.inp.fit['constrz'] == 0 and self.inp.fit['fixz'] == 0:
                self.mg.fixed[entries] = False
            elif 'rod' in entries and self.inp.fit['rod'] != 0:
                self.mg.fixed[entries] = False

            
def grain_values(lsqr):
        """
        Calculate the contributions from each grain
        For extreme contributions print a warning (*****)

		Jette Oddershede, Risoe DTU, May 15 2008
        """
        
        # rebuild function and load
        import build_fcn
        build_fcn.FCN(lsqr.inp)
        import fcn
        reload(fcn)
        # save values before making a new lsqr of minuit
#        temp1 = deepcopy(lsqr.m.values)		
#        temp2 = deepcopy(lsqr.m.errors)		
#        temp3 = deepcopy(lsqr.m.fixed)
#        temp4 = deepcopy(lsqr.mg.tol)
        g = n.zeros((lsqr.inp.no_grains))
        lsqr.inp.fit['poor'] = []
        lsqr.poor_value = []
        lsqr.poor_nrefl = []
        for i in range(lsqr.inp.no_grains):
            if i+1 in lsqr.inp.fit['skip']:
                pass
            else:		
                for j in range(lsqr.inp.nrefl[i]):
                    g[i] = g[i] + fcn.peak(lsqr.m.values['a'],lsqr.m.values['b'],lsqr.m.values['c'],lsqr.m.values['alpha'],lsqr.m.values['beta'],lsqr.m.values['gamma'],
                                     lsqr.inp.h[i][j],lsqr.inp.k[i][j],lsqr.inp.l[i][j],
                                     lsqr.inp.w[lsqr.inp.id[i][j]],lsqr.inp.dety[lsqr.inp.id[i][j]],lsqr.inp.detz[lsqr.inp.id[i][j]],
                                     lsqr.inp.vars[i][j], 
                                     lsqr.m.values['wx'],lsqr.m.values['wy'],
                                     lsqr.m.values['tx'],lsqr.m.values['ty'],lsqr.m.values['tz'],
                                     lsqr.m.values['py'],lsqr.m.values['pz'],
                                     lsqr.m.values['cy'],lsqr.m.values['cz'],
                                     lsqr.m.values['L'],
                                     lsqr.m.values['x%s' %i],lsqr.m.values['y%s' %i],lsqr.m.values['z%s' %i], 
                                     lsqr.inp.rod[i][0]+lsqr.m.values['rodx%s' %i],
                                     lsqr.inp.rod[i][1]+lsqr.m.values['rody%s' %i],
                                     lsqr.inp.rod[i][2]+lsqr.m.values['rodz%s' %i],
                                     lsqr.m.values['epsaa%s' %i],lsqr.m.values['epsab%s' %i],lsqr.m.values['epsac%s' %i], 
                                     lsqr.m.values['epsbb%s' %i],lsqr.m.values['epsbc%s' %i],lsqr.m.values['epscc%s' %i]) 
#        for i in range(lsqr.inp.no_grains):
#            if i+1 not in lsqr.inp.fit['skip']:
#                # make new lsqr of minuit        
#                lsqr.mg = Minuit(fcn.FCNgrain)
#                lsqr.mg.values = temp1
#                lsqr.mg.values['i'] = i
#                lsqr.mg.scan(("L",1,lsqr.mg.values['L']-1,lsqr.mg.values['L']+1)) # scan to set lsqr.m.fval, function starting value
#                g[i] = lsqr.mg.fval
        data = []
        poor = []
        for i in range(lsqr.inp.no_grains):
            if i+1 not in lsqr.inp.fit['skip']:
                data.append(g[i]/lsqr.inp.nrefl[i])
        reject.mad(data,poor,lsqr.inp.fit['rej_vol']**2)
        for i in range(lsqr.inp.no_grains):
            if i+1 not in lsqr.inp.fit['skip']:                
                print 'Grain %i %i: %e %f' %(i+1,lsqr.inp.nrefl[i],g[i],g[i]/lsqr.inp.nrefl[i])
		# give back old values	
#        lsqr.m.errors = temp2		
#        lsqr.m.fixed = temp3		
#        lsqr.mg.tol = temp4
            
        return g
			
			
def reject_outliers(lsqr):
        """
        Reject outliers peaks with a distance to the calculated peak position of
        more than lsqr.inp.fit['rej_resmean'] times the mean distance for the given grain	
		
		Jette Oddershede, Risoe DTU, May 15 2008
        """
		
        g = grain_values(lsqr)
        lsqr.inp.newreject = 0
        lsqr.inp.fit['newreject_grain'] = []
        #value = []
        new = 1
        while new == 1:
            new = 0
            for i in range(lsqr.inp.no_grains):
                #value.append([])
                if i+1 in lsqr.inp.fit['skip']:
                    pass
                else:		
                    for j in range(lsqr.inp.nrefl[i]-1,-1,-1): # loop backwards to make pop work
                        value = fcn.peak(lsqr.m.values['a'],lsqr.m.values['b'],lsqr.m.values['c'],lsqr.m.values['alpha'],lsqr.m.values['beta'],lsqr.m.values['gamma'],
                                        lsqr.inp.h[i][j],lsqr.inp.k[i][j],lsqr.inp.l[i][j],
                                        lsqr.inp.w[lsqr.inp.id[i][j]],lsqr.inp.dety[lsqr.inp.id[i][j]],lsqr.inp.detz[lsqr.inp.id[i][j]],
                                        lsqr.inp.vars[i][j], 
                                        lsqr.m.values['wx'],lsqr.m.values['wy'],
                                        lsqr.m.values['tx'],lsqr.m.values['ty'],lsqr.m.values['tz'],
                                        lsqr.m.values['py'],lsqr.m.values['pz'],
                                        lsqr.m.values['cy'],lsqr.m.values['cz'],
                                        lsqr.m.values['L'],
                                        lsqr.m.values['x%s' %i],lsqr.m.values['y%s' %i],lsqr.m.values['z%s' %i], 
                                        lsqr.inp.rod[i][0]+lsqr.m.values['rodx%s' %i],
                                        lsqr.inp.rod[i][1]+lsqr.m.values['rody%s' %i],
                                        lsqr.inp.rod[i][2]+lsqr.m.values['rodz%s' %i],
                                        lsqr.m.values['epsaa%s' %i],lsqr.m.values['epsab%s' %i],lsqr.m.values['epsac%s' %i], 
                                        lsqr.m.values['epsbb%s' %i],lsqr.m.values['epsbc%s' %i],lsqr.m.values['epscc%s' %i]) 
                        if value > lsqr.inp.fit['rej_resmean']*g[i]/lsqr.inp.nrefl[i]:
                            new = 1
                            print 'Rejected peak id %i from grain %i (hkl: %i %i %i, limit: %f): %f' %(lsqr.inp.id[i][j],i+1,lsqr.inp.h[i][j],lsqr.inp.k[i][j],lsqr.inp.l[i][j],lsqr.inp.fit['rej_resmean'],value*lsqr.inp.nrefl[i]/g[i])
                            reject.reject(lsqr.inp,i,j,value*lsqr.inp.nrefl[i]/g[i])
                        
        if 'globals' not in lsqr.inp.fit['goon']:
            lsqr.inp.mean_ia = []
            for i in range(lsqr.inp.no_grains):
                lsqr.inp.mean_ia.append([])
                for j in range(lsqr.inp.nrefl[i]):
                    lsqr.inp.mean_ia[i].append(1)
            reject.mean_ia(lsqr.inp,lsqr.inp.fit['rej_ia'])

            lsqr.inp.residual = []
            for i in range(lsqr.inp.no_grains):
                lsqr.inp.residual.append([])
                for j in range(lsqr.inp.nrefl[i]):
                    lsqr.inp.residual[i].append(1)
            reject.residual(lsqr.inp,lsqr.inp.fit['rej_resmedian'])

            lsqr.inp.volume = []
            for i in range(lsqr.inp.no_grains):
                lsqr.inp.volume.append([])
                for j in range(lsqr.inp.nrefl[i]):
                    lsqr.inp.volume[i].append(1)
            reject.intensity(lsqr.inp)

            reject.merge(lsqr.inp)
            reject.multi(lsqr.inp)
        
                        
        for i in range(lsqr.inp.no_grains):
            if lsqr.inp.nrefl[i] < lsqr.inp.fit['min_refl'] and i+1 not in lsqr.inp.fit['skip']:
                lsqr.inp.fit['skip'].append(i+1)
        lsqr.inp.fit['skip'].sort()
        
        
def scale_errors(lsqr,i=None):
        """
        Philosophy: Use const and near_const to tune final fval to approximately
                    3*sum(nrefl)-parameters, because:
                    1) Same const for a series facilitates evaluation of fit quality
                    2) fval is seen to decrease as the refinement proceeds
                    3) The tolerances depend on the scaling
        Scale the errors so that fval=3*sum(nrefl)-parameters
        This scale factor cannot be determined experimentally since it is detector
        specific and depends on for lsqr the gain.        
        """
        
        # remember only to apply correction to parameters refined in this particular cycle!!!!!!

        # parameters
        parameters = 0
        if i==None:
            for entries in lsqr.m.fixed:
                if lsqr.m.fixed[entries] == False:
                    parameters = parameters + 1
        else:
            for entries in lsqr.mg.fixed:
                if lsqr.mg.fixed[entries] == False:
                    parameters = parameters + 1

        #observations
        if i==None:
            observations = 0
            for j in range(lsqr.inp.no_grains):
                if j+1 in lsqr.inp.fit['skip']:
                    pass
                else:
                    observations = observations + lsqr.inp.nrefl[j]
        else:
            observations = lsqr.inp.nrefl[i]
              
        # expectation        
        expectation = 3*observations - parameters
        
        #correction
        if i==None:
            correction = lsqr.m.fval/expectation
            lsqr.m.up = correction
        else:
            correction = lsqr.mg.fval/expectation
            lsqr.mg.up = correction        
            
        # perform the  actual scaling task, NB must be done by calling hesse, with adjusted up, otherwise incorrect errors are estimated if the correct value of up is very far from 1
        if i==None:
            lsqr.m.hesse()
        else:
            lsqr.mg.hesse()
       
    


            
                
def refine(inp,killfile=None):
    check_input.interrupt(killfile)
    while inp.fit['goon'] != 'end':
        check_input.set_globals(inp)
        # build functions to minimise
        from FitAllB import build_fcn
        build_fcn.FCN(inp)
        import fcn
        reload(fcn)
        # minuit fitting
        from FitAllB import fit
        lsqr = fit.fit_minuit(inp)
        lsqr.refine()
        check_input.copy_globals(inp)
   

					
