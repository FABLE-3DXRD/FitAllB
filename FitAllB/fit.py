import numpy as n
import check_input
import write_output
import reject
import fcn
import time
import minuit
import sys
import logging
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
        self.globals = ["wx","wy","tx","ty","tz","py","pz","cy","cz","L"]
        self.grains = []
        for i in range(self.inp.no_grains):
            self.grains.append(["x%s" %i,"y%s" %i,"z%s" %i,"rodx%s" %i,"rody%s" %i,"rodz%s" %i,
                                "epsaa%s" %i,"epsbb%s" %i,"epscc%s" %i,"epsbc%s" %i,"epsac%s" %i,"epsab%s" %i])

        #refinement update
        reload(fcn)
        self.m = minuit.Minuit(fcn.FCN)
        self.m.values = self.inp.values
        self.m.errors = self.inp.errors
        self.m.printMode = self.inp.fit['printmode']
        self.m.strategy = self.inp.fit['strategy']
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
            self.mg = minuit.Minuit(fcn.FCNgrain)
            self.mg.values = self.m.values
            self.mg.errors = self.m.errors
            self.mg.fixed = self.m.fixed

            print '\n\n*****Now fitting %s*****' %self.inp.fit['goon']
#            print 'rerefine', self.inp.rerefine
            print 'newreject_grain', self.inp.fit['newreject_grain']
            # calculate starting values
            g = self.grain_values()
            self.fval = sum(g)
            print '\n%s starting value %e' %(self.inp.fit['goon'],self.fval)
            t1 = time.clock()
            self.mg = minuit.Minuit(fcn.FCNgrain)
            self.mg.values = self.m.values
            self.mg.errors = self.inp.errors
            for i in range(self.inp.no_grains):
                    if i+1 in self.inp.fit['skip']:
                        pass
#                    elif 'final' in self.inp.fit['goon'] and (i+1 not in self.inp.fit['newreject_grain'] or (self.inp.fit['newreject_grain'].count(i+1) == 1 and g[i]/self.inp.nrefl[i] < sum(g)/sum(self.inp.nrefl))):
                    elif 'final' in self.inp.fit['goon'] and i+1 not in self.inp.fit['newreject_grain']:# and i+1 not in self.inp.rerefine:
                        pass
                    elif 'xyz' in self.inp.fit['goon'] and i+1 not in self.inp.fit['newreject_grain'] and abs(self.mg.errors['x%i' %i] - self.inp.param['y_size']/5.) > 1e-3:# and i+1 not in self.inp.rerefine:
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
# if hesse != 0 covariance and errors from full hessian
#                        if self.inp.fit['hesse'] != 0:
#                            self.mg.hesse()
                        self.scale_errors(i)
                        #print self.mg.edm, self.mg.ncalls
                        self.m.errors = self.mg.errors
                        write_output.write_cor(self,i)
                        write_output.write_cov(self,i)
                        write_output.write_errors(self,i)
                        self.m.values = self.mg.values
                        g[i] = self.mg.fval
				
            self.time = time.clock()-t1
            print 'Fit %s time %i s' %(self.inp.fit['goon'],self.time)
            self.fval = sum(g)
            print 'Fit %s value %e \n' %(self.inp.fit['goon'],self.fval)
			    
			
            # reject outliers and save cycle info	
            self.m.errors = self.inp.errors
            self.reject_outliers()
            write_output.write_values(self)
            write_output.write_rej(self.inp,message=self.inp.fit['goon'])
            write_output.write_log(self)
#            write_output.write_gvectors(self)


        if 'final' in self.inp.fit['goon'] and (self.inp.newreject > 0):# or len(self.inp.rerefine) > 0):
            self.inp.fit['goon'] = 'grain'+ self.inp.fit['goon'][5:]
# delete rotpos rerefine
#        elif 'rotpos' in self.inp.fit['goon'] and (self.inp.newreject > 0):# or len(self.inp.rerefine) > 0):
#            self.inp.fit['goon'] = 'start'+ self.inp.fit['goon'][6:]
        elif 'xyz' in self.inp.fit['goon'] and (self.inp.newreject > 0):# or len(self.inp.rerefine) > 0):
            self.inp.fit['goon'] = 'start'+ self.inp.fit['goon'][3:]
        
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
        for entries in self.m.fixed:
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
        
        # rebuild function and load
        import build_fcn
        build_fcn.FCN(self.inp)
        import fcn
        reload(fcn)
        # save values before making a new instance of minuit
        temp1 = self.m.values		
        temp2 = self.m.errors		
        temp3 = self.m.fixed
        temp4 = self.mg.tol
        g = n.zeros((self.inp.no_grains))
        self.inp.fit['poor'] = []
        self.poor_value = []
        self.poor_nrefl = []
        for i in range(self.inp.no_grains):
            if i+1 not in self.inp.fit['skip']:
                # make new instance of minuit        
                self.mg = minuit.Minuit(fcn.FCNgrain)
                self.mg.values = temp1
                self.mg.values['i'] = i
                self.mg.scan(("L",1,self.mg.values['L']-1,self.mg.values['L']+1)) # scan to set self.m.fval, function starting value
                g[i] = self.mg.fval
        data = []
        poor = []
        for i in range(self.inp.no_grains):
            if i+1 not in self.inp.fit['skip']:
                data.append(g[i]/self.inp.nrefl[i])
        reject.mad(data,poor,self.inp.fit['mad'][1])
#        print max(data), data, poor
        for i in range(self.inp.no_grains):
            if i+1 not in self.inp.fit['skip']:                
                print 'Grain %i %i: %e %f' %(i+1,self.inp.nrefl[i],g[i],g[i]/self.inp.nrefl[i])
		# give back old values	
        self.m.errors = temp2		
        self.m.fixed = temp3		
        self.mg.tol = temp4
            
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
        self.inp.rerefine = []
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
                        value = fcn.peak(self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],
                                        self.inp.w[self.inp.id[i][j]],self.inp.dety[self.inp.id[i][j]],self.inp.detz[self.inp.id[i][j]],
                                        #n.array([self.inp.Syy[self.inp.id[i][j]],self.inp.Szz[self.inp.id[i][j]],self.inp.Sww[self.inp.id[i][j]]]),
                                        self.inp.vars[i][j], 
                                        self.m.values['wx'],self.m.values['wy'],
                                        self.m.values['tx'],self.m.values['ty'],self.m.values['tz'],
                                        self.m.values['py'],self.m.values['pz'],
                                        self.m.values['cy'],self.m.values['cz'],
                                        self.m.values['L'],
                                        self.m.values['x%s' %i],self.m.values['y%s' %i],self.m.values['z%s' %i], 
                                        self.inp.rod[i][0]+self.m.values['rodx%s' %i],
                                        self.inp.rod[i][1]+self.m.values['rody%s' %i],
                                        self.inp.rod[i][2]+self.m.values['rodz%s' %i],
                                        self.m.values['epsaa%s' %i],self.m.values['epsab%s' %i],self.m.values['epsac%s' %i], 
                                        self.m.values['epsbb%s' %i],self.m.values['epsbc%s' %i],self.m.values['epscc%s' %i]) 
                        if value > self.inp.fit['limit'][1]*g[i]/self.inp.nrefl[i]:
                            new = 1
                            print 'Rejected peak id %i from grain %i (hkl: %i %i %i, limit: %f): %f' %(self.inp.id[i][j],i+1,self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],self.inp.fit['limit'][1],value*self.inp.nrefl[i]/g[i])
                            reject.reject(self.inp,i,j,value*self.inp.nrefl[i]/g[i])
                        
        if 'final' in self.inp.fit['goon'] or 'grain' in self.inp.fit['goon']:# or 'rotpos' in self.inp.fit['goon']:
            self.inp.mean_ia = []
            for i in range(self.inp.no_grains):
                self.inp.mean_ia.append([])
                for j in range(self.inp.nrefl[i]):
                    self.inp.mean_ia[i].append(1)
            reject.mean_ia(self.inp,self.inp.fit['ia'])

        self.inp.residual = []
        for i in range(self.inp.no_grains):
            self.inp.residual.append([])
            for j in range(self.inp.nrefl[i]):
                self.inp.residual[i].append(1)
        reject.residual(self.inp,self.inp.fit['limit'][0])

        self.inp.volume = []
        for i in range(self.inp.no_grains):
            self.inp.volume.append([])
            for j in range(self.inp.nrefl[i]):
                self.inp.volume[i].append(1)
        reject.intensity(self.inp)

        reject.merge(self.inp)
        reject.multi(self.inp)
        
                        
        for i in range(self.inp.no_grains):
            if self.inp.nrefl[i] < self.inp.fit['min_refl'] and i+1 not in self.inp.fit['skip']:
                self.inp.fit['skip'].append(i+1)
        self.inp.fit['skip'].sort()

                		
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
        self.mg.printMode = 0
        
        if self.inp.fit['goon'] == 'grain':
            self.mg.tol = self.inp.fit['tol_grain']
        else:
            self.mg.tol = self.inp.fit['tol_grain']*0.1
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for entries in self.grains[i]:
            if (entries[0]=='x' or entries[0]=='y' or entries[0]=='z') and self.inp.fit['xyz'] != 0:
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
            if (entries[0]=='x' or entries[0]=='y' or entries[0]=='z') and self.inp.fit['xyz'] != 0:
                self.mg.fixed[entries] = False
            elif 'rod' in entries and self.inp.fit['rod'] != 0:
                self.mg.fixed[entries] = False

            
		                
                
def refine(inp):
    inp.rerefine = []
    while inp.fit['goon'] != 'end':
        check_input.set_globals(inp)
        # calculate experimental errors using the present values 
#        from FitAllB import error
#        error.vars_scale(inp)   # function to ensure correct relative scaling of variances between grains
#        error.vars(inp)
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
   

					
