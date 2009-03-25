import numpy as n
import check_input
import write_output
import reject
import fcn
import time
import minuit
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
	Reject reflection according to self.inp.fit['limit'][1]
	Print and save refinement and rejection info and parameters
	
	Jette Oddershede, Risoe DTU, May 15 2008
	"""
		
    # initialise
        self.poor_value = []
        self.poor_nrefl = []

 		# create lists of parameters, global and for each grain
        self.globals = ["wx","wy","tx","ty","tz","py","pz","cy","cz","L"]
#        for entry in self.globals:
#            print entry, self.inp.values[entry]
        self.grains = []
        for i in range(self.inp.no_grains):
            self.grains.append(["x%s" %i,"y%s" %i,"z%s" %i,"rodx%s" %i,"rody%s" %i,"rodz%s" %i,#"phia%s" %i,"PHI%s" %i,"phib%s" %i,
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
        if 'globals' in self.inp.fit['goon'] or 'final' in self.inp.fit['goon'] or 'rotpos' in self.inp.fit['goon']:
            self.ref = True
		

		# carry out refinement
        if self.ref == True:
            self.mg = minuit.Minuit(fcn.FCNgrain)
            self.mg.values = self.m.values
            self.mg.errors = self.m.errors
            self.mg.fixed = self.m.fixed

            print '\n\n*****Now fitting %s*****' %self.inp.fit['goon']
            print 'rerefine', self.inp.rerefine
            print 'newreject_grain', self.inp.fit['newreject_grain']
            # calculate starting values
            g = self.grain_values()
            self.fval = sum(g)
            print '\n%s starting value %e' %(self.inp.fit['goon'],self.fval)
            t1 = time.clock()
            self.fitglobals()
            print 'Fit %s tolerance %e' %(self.inp.fit['goon'],self.m.tol)
            self.m.errors = self.inp.errors
            self.m.migrad()
            if self.inp.fit['hesse'] != 0:
                self.mg.hesse()
            self.scale_errors(0)
            write_output.write_global(self)
				
            self.time = time.clock()-t1
            print 'Fit %s time %i s' %(self.inp.fit['goon'],self.time)
            print 'Fit %s value %e \n' %(self.inp.fit['goon'],self.m.fval)
			    
			
            # reject outliers and save cycle info	
            self.m.errors = self.inp.errors
            self.reject_outliers()
            self.mean_ia()
            write_output.write_values(self)
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
        
        example = 'L' 
        
        # parameters
        parameters = 0
        for entries in self.m.fixed:
                if self.m.fixed[entries] == False:
                    parameters = parameters + 1
                    example = entries
        
        grains = self.inp.no_grains - len(self.inp.fit['skip'])
        
        #observations
        observations = 0
        for j in range(self.inp.no_grains):
                if j+1 in self.inp.fit['skip']:
                    pass
                else:
                    observations = observations + self.inp.nrefl[j]
              
        # expectation        
        expectation = 3*observations - grains*parameters
        
        #correction
        correction = self.m.fval/expectation
        self.m.up = correction
            
        # perform the  actual scaling task
        self.m.hesse()
       
    
    def mean_ia(self):
        """
        Calculate mean internal angle for each grain and store in self.mena_ia array of length self.inp.no_grains
        Jette Oddershede Januar 2009
        """
        
        self.mean_ia = n.zeros(self.inp.no_grains)
        for i in range(self.inp.no_grains):
            if i+1 in self.inp.fit['skip']:
                pass
            else:
                rod = n.array([self.inp.rod[i][0]+self.inp.values['rodx%s' %i],self.inp.rod[i][1]+self.inp.values['rody%s' %i],self.inp.rod[i][2]+self.inp.values['rodz%s' %i]])
                rot = n.array([0,0,0])
#                print i, rod/n.linalg.norm(rod)
                for j in range(self.inp.nrefl[i]):
                    gexp = fcn.gexp(self.inp.w[self.inp.id[i][j]],self.inp.dety[self.inp.id[i][j]],self.inp.detz[self.inp.id[i][j]],
                                    self.inp.values['wx'],self.inp.values['wy'],self.inp.values['tx'],self.inp.values['ty'],self.inp.values['tz'],
                                    self.inp.values['py'],self.inp.values['pz'],self.inp.values['cy'],self.inp.values['cz'],self.inp.values['L'],
                                    self.inp.values['x%s' %i],self.inp.values['y%s' %i],self.inp.values['z%s' %i])
                    gcalc = fcn.gcalc(self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],
                                      self.inp.rod[i][0]+self.inp.values['rodx%s' %i],
                                      self.inp.rod[i][1]+self.inp.values['rody%s' %i],
                                      self.inp.rod[i][2]+self.inp.values['rodz%s' %i],
                                      self.inp.values['epsaa%s' %i],self.inp.values['epsab%s' %i],self.inp.values['epsac%s' %i],
                                      self.inp.values['epsbb%s' %i],self.inp.values['epsbc%s' %i],self.inp.values['epscc%s' %i])
                    (ia,rotj,ia2,rotj2) = self.IAforrod(n.transpose(gexp)[0],
                                                        n.transpose(gcalc)[0],
                                                        rod)

                    self.mean_ia[i] = self.mean_ia[i] + ia 
                self.mean_ia[i] = self.mean_ia[i] / self.inp.nrefl[i]  
    
    def grain_values(self):
        """
        Calculate the contributions from each grain
        For extreme contributions print a warning (*****)

		Jette Oddershede, Risoe DTU, May 15 2008
        """
        
        # save values before making a new instance of minuit
        temp1 = self.m.values		
        temp2 = self.m.errors		
        temp3 = self.m.fixed
        if 'start' in self.inp.fit['goon']:
            temp4 = self.m.tol
        else:
            temp4 = self.mg.tol
        # make new instance of minuit
        self.m = minuit.Minuit(fcn.FCN)
        self.m.values = temp1		
        self.m.scan(("L",1,self.m.values['L']-1,self.m.values['L']+1)) # scan to set self.m.fval, function starting value	
        g = n.zeros((self.inp.no_grains))
        self.inp.fit['poor'] = []
        self.poor_value = []
        self.poor_nrefl = []
        for i in range(self.inp.no_grains):
            if i+1 not in self.inp.fit['skip']:
                self.mg = minuit.Minuit(fcn.FCNgrain)
                self.mg.values = self.m.values
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
#                if g[i] > 2*self.m.fval/(self.inp.no_grains-len(self.inp.fit['skip'])):
#                if g[i]/self.inp.nrefl[i] > max(data):# and self.inp.nrefl[i] < 24:
#                    print 'Grain %i %i: %e %f *****' %(i+1,self.inp.nrefl[i],g[i],g[i]/self.inp.nrefl[i])
#                    self.inp.fit['poor'].append(i+1)
#                    self.poor_value.append(g[i]/self.inp.nrefl[i]*len(data)/sum(data))
#                    self.poor_nrefl.append(self.inp.nrefl[i])                    
#                else:	
#                    print 'Grain %i %i: %e %f' %(i+1,self.inp.nrefl[i],g[i],g[i]/self.inp.nrefl[i])
		# give back old values	
        self.m.errors = temp2		
        self.m.fixed = temp3		
        if 'start' in self.inp.fit['goon']:
            self.m.tol = temp4
        else:
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
                                        self.inp.w[self.inp.id[i][j]],self.inp.dety[self.inp.id[i][j]],self.inp.detz[self.inp.id[i][j]],self.inp.vars[i][j], 
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
#                            print 'Rejected peak id %i from grain %i (hkl: %i %i %i, limit: %f): %f' %(self.inp.id[i][j],i+1,self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],self.inp.fit['limit'][1],value*self.inp.nrefl[i]/g[i],g[i],self.inp.nrefl[i])
                            print 'Rejected peak id %i from grain %i (hkl: %i %i %i, limit: %f): %f' %(self.inp.id[i][j],i+1,self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],self.inp.fit['limit'][1],value*self.inp.nrefl[i]/g[i])
                            reject.reject(self.inp,i,j,value*self.inp.nrefl[i]/g[i])
                        
        for i in range(self.inp.no_grains):
            # rerefine if more than 10% change in fcn
            if self.inp.nrefl[i] < self.inp.fit['min_refl'] and i+1 not in self.inp.fit['skip']:
                self.inp.fit['skip'].append(i+1)
            if n.sum(self.inp.residual[i])/len(self.inp.residual[i]) < 2.9 and i+1 not in self.inp.fit['skip']: 
                self.inp.rerefine.append(i+1)
        self.inp.fit['skip'].sort()

                
    def volumes(self):
        """
        Determine the relative grain volumes based on the relative intensities.
        Volume set to 1 for the grain for which most reflections are assigned.
        Peaks with centres less than 10 pixels from the edge of the frame or less than 5 omega steps from the measured limits are disregarded
        Only use the middle one third of the volume fractions 
        """

        refmax = 0
        pixeltol = 10
        wsteptol = 5
        edge = []
        for i in range(self.inp.no_grains):
            edgegr = []
            if i+1 in self.inp.fit['skip']:
                pass
            elif self.inp.nrefl[i] > refmax:
                refmax = self.inp.nrefl[i]
                imax = i
            for j in range(self.inp.nrefl[i]):
                for w_limit in self.inp.fit['w_limit']:
                    if abs(w_limit-self.w[self.id[i][j]]) <= wsteptol*self.inp.fit['w_step']:
#                        print 'w', self.w[self.id[i][j]]
                        edgegr.append(j)
                        break
                if abs(self.dety[self.id[i][j]]) <= pixeltol or abs(self.dety[self.id[i][j]]) >= self.param['y_center']*2-pixeltol:
#                    print 'y', abs(self.dety[self.id[i][j]])
                    edgegr.append(j)
                elif abs(self.detz[self.id[i][j]]) <= pixeltol or abs(self.detz[self.id[i][j]]) >= self.param['z_center']*2-pixeltol:
#                    print 'z', abs(self.detz[self.id[i][j]])
                    edgegr.append(j)
                elif abs(n.sin(self.eta[self.id[i][j]]*n.pi/180.)) < 1e-4:
#                    print 'eta', n.sin(self.eta[self.id[i][j]]*n.pi/180.)
                    edgegr.append(j)
            edge.append(edgegr)
           
#        print 'edge', edge
        
        for i in range(self.inp.no_grains):
            self.volerr[i] = 0
            if i == imax or i+1 in self.inp.fit['skip']:
                self.volume[i] = 1
#                for j in range(self.inp.nrefl[i]):    
#                    print i,self.int[self.id[i][j]], self.h[i][j], self.k[i][j],self.l[i][j]
                pass
            else:
                self.volume[i] = 0
                data = []
                for j in range(self.inp.nrefl[i]):
                    if j in edge[i]:
                        pass
                    for jmax in range(self.inp.nrefl[imax]):
                        if jmax in edge[imax]:
                            pass
                        elif self.int[self.id[imax][jmax]] > 0 and self.h[i][j] == self.h[imax][jmax] and self.k[i][j] == self.k[imax][jmax] and self.l[i][j] == self.l[imax][jmax]:
                            vol = self.int[self.id[i][j]]*abs(n.sin(self.eta[self.id[i][j]]))/(self.int[self.id[imax][jmax]]*abs(n.sin(self.eta[self.id[imax][jmax]])))
                            data.append(vol)
                data.sort()
                ndata = 0
                for k in range(len(data)/3,len(data)*2/3,1):
                    self.volume[i] = self.volume[i] + data[k]
                    self.volerr[i] = self.volerr[i] + data[k]**2
                    ndata = ndata + 1
                if ndata != 0:
                    self.volume[i] = self.volume[i]/ndata
                    self.volerr[i] = n.sqrt(self.volerr[i]/ndata - self.volume[i]**2)
                else:
                    self.volume[i] = 1
                    self.volerr[i] = 1

 
		
    def fitglobals(self):
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
            elif entries=='cy' and self.inp.fit['center'] != 0:
                self.m.fixed[entries] = False
            elif 'L' in entries and self.inp.fit['L'] != 0:
                self.m.fixed[entries] = False

		

    def IAforrod(self,gv1,gv2,rod):
        """
        Calculates the internal angle ia between gvectors gv1 and gv2 relative to a 
        rotation axis given as a rodrigues vector rod.
        gv1,gv2,rod n.array(1x3)
        Returns ia in degrees
    
        Jette Oddershede, Jan 2009
        """
    
        # rotation axis must be projected onto the plane of gv1xgv2 and gv2+gv2
        gvcross = n.cross(gv1,gv2)
        gvcross = gvcross/n.linalg.norm(gvcross)
        gvadd = gv1+gv2
        gvadd = gvadd/n.linalg.norm(gvadd)
        # calculate normal of this plane
        gvnorm = n.cross(gvcross,gvadd)
        gvnorm = gvnorm/n.linalg.norm(gvnorm)
        
        # rotation vector in plane of gv1xgv2 and gv2+gv2 now given as
        rot = rod - n.dot(rod,gvnorm)*gvnorm
        rot = rot/n.linalg.norm(rot)
        
        #project gv1 and gv2 onto the plane perpendicular to rot
        gv1_proj = gv1 - n.dot(gv1,rot)*rot
        gv2_proj = gv2 - n.dot(gv2,rot)*rot
    
        # The desired angle is now the angle between gv1_proj and gv2_proj
        gv1_proj = gv1_proj/n.linalg.norm(gv1_proj)
        gv2_proj = gv2_proj/n.linalg.norm(gv2_proj)
        ia = n.arccos(n.dot(gv1_proj,gv2_proj))
        
        # ia2 angle between rod and rot around norm
        rod = rod/n.linalg.norm(rod)
        ia2 = n.arccos(n.dot(rod,rot))
        norm = n.cross(rod,rot)
        norm = norm/n.linalg.norm(norm)
    
        return (ia*180./n.pi,rot,ia2*180./n.pi,norm)
            
		                
                
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
        # minuit fitting
        from FitAllB import fitgg
        lsqr = fitgg.fit_minuit(inp)
        lsqr.refine()
        check_input.copy_globals(inp)
   
   

					
