

import numpy as n
from . import check_input
from . import write_output
from . import reject
import fcn
import time
import importlib
try:
    from iminuit import Minuit
except ImportError:
    from minuit import Minuit
import sys
import logging
from copy import deepcopy
from importlib import reload
logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')


class fit_minuit():
    def __init__(self,inp):
        self.inp = inp
        
			
    def refine(self):	
	"""
	Carry out one refinement cycle according to the order given by self.inp.fit['reforder']
	Reject reflection according to self.inp.fit['limit'][1]
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
        importlib.reload(fcn)

		# determine whether to refine
        self.ref = False
        if 'globals' in self.inp.fit['goon']:
            self.ref = True
		

		# carry out refinement
        if self.ref == True:
            print('\n\n*****Now fitting %s*****' %self.inp.fit['goon'])
#            print 'rerefine', self.inp.rerefine
            print('newreject_grain', self.inp.fit['newreject_grain'])
            # calculate starting values
            g = self.grain_values()
            self.fval = sum(g)
            print('\n%s starting value %e' %(self.inp.fit['goon'],self.fval))
            t1 = time.clock()
            try:
                self.mg = Minuit(fcn.FCNgrain,errordef=1)
            	self.mg.values = self.inp.values
            	self.mg.errors = self.inp.errors
            	self.mg.printMode = self.inp.fit['printmode']
            	self.mg.strategy = self.inp.fit['strategy']
            	for entries in self.mg.fixed:
                	self.mg.fixed[entries] = True
            except:
                self.mg = Minuit(fcn.FCNgrain,errordef=1,**self.inp.fitarg)
            global_parameters = []
            weight = []
            for i in range(self.inp.no_grains):
                    if i+1 in self.inp.fit['skip']:
                        pass
                    else:	
                        if i == 0:
                            print('Fit %s tolerance %e' %(self.inp.fit['goon'],self.mg.tol))
                        self.mg.values['i'] = i
                        try:
                            self.mg.fitarg['i'] = i
                        except:
                            pass
                        self.fitglobalgrain(i)
                        print('\rRefining grain %i' %(i+1), end=' ')
                        sys.stdout.flush()
                        self.mg.migrad()
#                        self.scale_errors(i)
                        g[i] = self.mg.fval
                        temp = []
                        for j in range(len(self.globals)):
                            temp.append(self.mg.values[self.globals[j]])
                        global_parameters.append(temp)
                        weight.append(len(self.inp.mean_ia[i])/n.sum(self.inp.mean_ia[i]))
                        self.m = self.mg
                        write_output.write_global(self)
				
            self.time = time.clock()-t1
            print('Fit %s time %i s' %(self.inp.fit['goon'],self.time))
            self.fval = sum(g)
            print('Fit %s value %e \n' %(self.inp.fit['goon'],self.fval))
			    
			# get global parameters and errors on these as average and spread of global_parameters
            global_parameters = n.array(global_parameters)
            weight = n.array(weight)
            weight = weight/n.sum(weight)
#            average_global_parameters = global_parameters.mean(0)
#            spread_global_parameters = global_parameters.std(0)
            average_global_parameters = n.zeros(len(self.globals))
            spread_global_parameters = n.zeros(len(self.globals))
            for i in range(len(global_parameters)):
                average_global_parameters = average_global_parameters + global_parameters[i] * weight[i] 
            for i in range(len(global_parameters)):
                spread_global_parameters = spread_global_parameters + (global_parameters[i] - average_global_parameters)**2 * weight[i] 
            spread_global_parameters = spread_global_parameters**.5
            
            print(global_parameters)
            print(weight)
            print(average_global_parameters)
            print(spread_global_parameters)
            
            for j in range(len(self.globals)):
                self.mg.values[self.globals[j]] = average_global_parameters[j]
                self.mg.errors[self.globals[j]] = spread_global_parameters[j]
#            print 'mg',self.mg.values['L']
#            print 'm',self.m.values['L']
#            print 'inp',self.inp.values['L']
            
            # reject outliers and save cycle info	
            self.inp.errors = self.mg.errors 
            self.inp.values = self.mg.values 
            self.m.errors = self.mg.errors 
            self.m.values = self.mg.values 
            self.reject_outliers()
            write_output.write_rej(self.inp,message=self.inp.fit['goon'])
            write_output.write_log(self)
            if 'globals' in self.inp.fit['goon']:
                write_output.write_par(self)

       
		# move onto next refinement given by the reforder list	
        self.inp.fit['goon'] = self.inp.fit['reforder'][self.inp.fit['reforder'].index(self.inp.fit['goon'])+1]
	
        return
        
        
    def scale_errors(self,i):
        """
        Philosophy: Use const and near_const to tune final fval to approximately
                    3*sum(nrefl)-parameters, because:
                    1) Same const for a series facilitates evaluation of fit quality
                    2) fval is seen to decrease as the refinement proceeds
                    3) The tolerances depend on the scaling
        Scale the errors so that fval=3*sum(nrefl)-parameters
        This scale factor cannot be determined experimentally since it is detector
        specific and depends on for instance the gain.        
        """
        
        # remember only to apply correction to parameters refined in this particular cycle!!!!!!

        example = 'L' 
        
        # parameters
        parameters = 0
        for entries in self.mg.fixed:
            if self.mg.fixed[entries] == False:
                parameters = parameters + 1
                example = entries
        #grains
        grains = 1
        #observations
        observations = self.inp.nrefl[i]
              
        # expectation        
        expectation = 3*observations - grains*parameters
        
        #correction
        correction = self.mg.fval/expectation
        self.mg.up = correction
            
        # perform the  actual scaling task, NB must be done by calling hesse, with adjusted up, otherwise incorrect errors are estimated if the correct value of up is very far from 1
        self.mg.hesse()
       
    
    def grain_values(self):
        """
        Calculate the contributions from each grain
        For extreme contributions print a warning (*****)

		Jette Oddershede, Risoe DTU, May 15 2008
        """
        
        g = n.zeros((self.inp.no_grains))
        for i in range(self.inp.no_grains):
            if i+1 not in self.inp.fit['skip']:                
                g[i] = n.sum(self.inp.residual[i])
                print('Grain %i %i: %e %f' %(i+1,self.inp.nrefl[i],g[i],g[i]/self.inp.nrefl[i]))
            
        return g
			
			
    def reject_outliers(self):
        """
        Reject outliers peaks with a distance to the calculated peak position of
        more than self.inp.fit['limit'][1] times the mean distance for the given grain	
		
		Jette Oddershede, Risoe DTU, May 15 2008
        """
		
        g = self.grain_values()
        self.inp.newreject = 0
        self.inp.fit['newreject_grain'] = []
#        self.inp.rerefine = []
        #value = []
        new = 1
        while new == 1:
            new = 0
            for i in range(self.inp.no_grains):
                #value.append([])
                if i+1 in self.inp.fit['skip']:
                    pass
                else:		
                    for j in range(self.inp.nrefl[i]-1,-1,-1): # loop backwards to make pop work
                        value = fcn.peak(lsqr.mg.values['a'],lsqr.mg.values['b'],lsqr.mg.values['c'],lsqr.mg.values['alpha'],lsqr.mg.values['beta'],lsqr.mg.values['gamma'],
                                        self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],
                                        self.inp.w[self.inp.id[i][j]],self.inp.dety[self.inp.id[i][j]],self.inp.detz[self.inp.id[i][j]],
                                        #n.array([self.inp.Syy[self.inp.id[i][j]],self.inp.Szz[self.inp.id[i][j]],self.inp.Sww[self.inp.id[i][j]]]),
                                        self.inp.vars[i][j], 
                                        self.mg.values['wx'],self.mg.values['wy'],
                                        self.mg.values['tx'],self.mg.values['ty'],self.mg.values['tz'],
                                        self.mg.values['py'],self.mg.values['pz'],
                                        self.mg.values['cy'],self.mg.values['cz'],
                                        self.mg.values['L'],
                                        self.mg.values['x%s' %i],self.mg.values['y%s' %i],self.mg.values['z%s' %i], 
                                        self.inp.rod[i][0]+self.mg.values['rodx%s' %i],
                                        self.inp.rod[i][1]+self.mg.values['rody%s' %i],
                                        self.inp.rod[i][2]+self.mg.values['rodz%s' %i],
                                        self.mg.values['epsaa%s' %i],self.mg.values['epsab%s' %i],self.mg.values['epsac%s' %i], 
                                        self.mg.values['epsbb%s' %i],self.mg.values['epsbc%s' %i],self.mg.values['epscc%s' %i]) 
                        if value > self.inp.fit['limit'][1]*g[i]/self.inp.nrefl[i]:
                            new = 1
                            print('Rejected peak id %i from grain %i (hkl: %i %i %i, limit: %f): %f' %(self.inp.id[i][j],i+1,self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],self.inp.fit['limit'][1],value*self.inp.nrefl[i]/g[i]))
                            reject.reject(self.inp,i,j,value*self.inp.nrefl[i]/g[i])
                        
        for i in range(self.inp.no_grains):
            if self.inp.nrefl[i] < self.inp.fit['min_refl'] and i+1 not in self.inp.fit['skip']:
                self.inp.fit['skip'].append(i+1)
        self.inp.fit['skip'].sort()

                		
    def fitglobalgrain(self,i):
	"""
	Set tolerance and fixed parameters for preliminary fit of the global parameters
	"""
        self.mg.tol = self.inp.fit['tol_global']
        
#        print self.mg.fixed
        for entries in self.mg.fixed:
            if entries=='wy' and self.inp.fit['w'] != 0:
                self.mg.fixed[entries] = False
#            elif entries=='wx' and self.inp.fit['w'] != 0:
#                self.mg.fixed[entries] = False
            elif entries[0]=='t' and self.inp.fit['tilt'] != 0:
                self.mg.fixed[entries] = False
            elif 'p' in entries and len(entries) == 2 and self.inp.fit['pixel'] != 0:
                self.mg.fixed[entries] = False
            elif entries=='cy' and self.inp.fit['center'] != 0:
                self.mg.fixed[entries] = False
            elif 'L' in entries and self.inp.fit['L'] != 0:
                self.mg.fixed[entries] = False
#            if self.mg.fixed[entries] == False:
#                print entries
             

                
def refine(inp):
    while inp.fit['goon'] != 'end':
        check_input.set_globals(inp)
        # calculate experimental errors using the present values 
        from FitAllB import error
        # build functions to minimise
        from FitAllB import build_fcn
        build_fcn.FCN(inp)
        import fcn
        importlib.reload(fcn)
        # minuit fitting
        from FitAllB import globals
        lsqr = globals.fit_minuit(inp)
        lsqr.refine()
        check_input.copy_globals(inp)
   

					
