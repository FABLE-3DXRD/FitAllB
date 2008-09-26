import numpy as n
from xfab import tools
from xfab import sg
from xfab import detector
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
            self.grains.append(["x%s" %i,"y%s" %i,"z%s" %i,"phia%s" %i,"PHI%s" %i,"phib%s" %i,
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
        if 'grain' in self.inp.fit['goon'] or 'final' in self.inp.fit['goon']:
            self.ref = True
        elif 'start' in self.inp.fit['goon'] and (self.inp.fit['w'] != 0 or self.inp.fit['tilt'] != 0 or self.inp.fit['pixel'] != 0 or self.inp.fit['center'] != 0 or self.inp.fit['L'] != 0):
            self.ref = True
        elif 'euler' in self.inp.fit['goon'] and self.inp.fit['euler'] != 0:
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
            # calculate starting values
            g = self.grain_values()
            fval = sum(g)
            print '\n%s starting value %e' %(self.inp.fit['goon'],fval)
            t1 = time.clock()
            if 'start' in self.inp.fit['goon']:
                self.fitstart()
                print 'Fit %s tolerance %e' %(self.inp.fit['goon'],self.m.tol)
                self.m.errors = self.inp.errors
                if self.inp.fit['hesse'] != 0:
                    self.mg.hesse()
                self.m.migrad()
                write_output.write_global(self)
            else:
                self.mg = minuit.Minuit(fcn.FCNgrain)
                self.mg.values = self.m.values
                self.mg.errors = self.inp.errors
                for i in range(self.inp.no_grains):
                    if i+1 in self.inp.fit['skip']:
                        pass
                    elif 'final' in self.inp.fit['goon'] and (i+1 not in self.inp.fit['newreject_grain'] or (self.inp.fit['newreject_grain'].count(i+1) == 1 and g[i]/self.inp.nrefl[i] < sum(g)/sum(self.inp.nrefl))):
                        pass
                    else:	
                        if 'grain' in self.inp.fit['goon'] or 'final' in self.inp.fit['goon']:
                            self.fitgrain(i)
                        elif 'eps' in self.inp.fit['goon']:
                            self.fitepsgrain(i)
                        elif 'xyz' in self.inp.fit['goon']:
                            self.fitxyzgrain(i)
                        elif 'euler' in self.inp.fit['goon']:
                            self.fiteulergrain(i)
                        if i == 0:
                            print 'Fit %s tolerance %e' %(self.inp.fit['goon'],self.mg.tol)
                        self.mg.values['i'] = i
                        print '\rRefining grain %i' %(i+1),
                        sys.stdout.flush()
                        # if hesse != 0 covariance and errors from full hessian
                        if self.inp.fit['hesse'] != 0:
                            self.mg.hesse()
                        self.mg.migrad()
                        #print self.mg.edm, self.mg.ncalls
                        self.m.errors = self.mg.errors
                        if 'grain' in self.inp.fit['goon'] or 'final' in self.inp.fit['goon']:
                            write_output.write_cor(self,i)
                            write_output.write_cov(self,i)
                            write_output.write_errors(self,i)
                        self.m.values = self.mg.values
#                        self.m.tol = self.mg.tol
                        g[i] = self.mg.fval
				
            self.time = time.clock()-t1
            print 'Fit %s time %i s' %(self.inp.fit['goon'],self.time)
            if 'start' in self.inp.fit['goon']:
                print 'Fit %s value %e \n' %(self.inp.fit['goon'],self.m.fval)
            else:
                fval = sum(g)
                print 'Fit %s value %e \n' %(self.inp.fit['goon'],fval)
			    
			
            # reject outliers and save cycle info	
            self.m.errors = self.inp.errors
            self.reject_outliers()
            write_output.write_values(self)
            write_output.write_rej(self.inp,message=self.inp.fit['goon'])
            write_output.write_log(self)

        if 'final' in self.inp.fit['goon'] and self.inp.newreject > 0:
            self.inp.fit['goon'] = 'grain'
        
		# move onto next refinement given by the reforder list	
        self.inp.fit['goon'] = self.inp.fit['reforder'][self.inp.fit['reforder'].index(self.inp.fit['goon'])+1]
	
        return
        
	
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
        reject.mad(data,poor,25)
#        print max(data), data, poor
        for i in range(self.inp.no_grains):
            if i+1 not in self.inp.fit['skip']:                
#                if g[i] > 2*self.m.fval/(self.inp.no_grains-len(self.inp.fit['skip'])):
                if g[i]/self.inp.nrefl[i] > max(data):# and self.inp.nrefl[i] < 24:
                    print 'Grain %i %i: %e %f *****' %(i+1,self.inp.nrefl[i],g[i],g[i]/self.inp.nrefl[i])
                    self.inp.fit['poor'].append(i+1)
                    self.poor_value.append(g[i]/self.inp.nrefl[i]*len(data)/sum(data))
                    self.poor_nrefl.append(self.inp.nrefl[i])                    
                else:	
                    print 'Grain %i %i: %e %f' %(i+1,self.inp.nrefl[i],g[i],g[i]/self.inp.nrefl[i])
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
        #value = []
        for i in range(self.inp.no_grains):
            #value.append([])
            if i+1 in self.inp.fit['skip']:
                pass
            else:				
                for j in range(self.inp.nrefl[i]-1,-1,-1): # loop backwards to make pop work
                    value = fcn.peak(self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],
                                     self.inp.w[self.inp.id[i][j]],self.inp.dety[self.inp.id[i][j]],self.inp.detz[self.inp.id[i][j]],fcn.vars[i][j], 
                                     self.m.values['wx'],self.m.values['wy'],
                                     self.m.values['tx'],self.m.values['ty'],self.m.values['tz'],
                                     self.m.values['py'],self.m.values['pz'],
                                     self.m.values['cy'],self.m.values['cz'],
                                     self.m.values['L'],
                                     self.m.values['x%s' %i],self.m.values['y%s' %i],self.m.values['z%s' %i], 
                                     self.m.values['phia%s' %i],self.m.values['PHI%s' %i],self.m.values['phib%s' %i], 
                                     self.m.values['epsaa%s' %i],self.m.values['epsab%s' %i],self.m.values['epsac%s' %i], 
                                     self.m.values['epsbb%s' %i],self.m.values['epsbc%s' %i],self.m.values['epscc%s' %i]) 
                    if value > self.inp.fit['limit'][1]*g[i]/self.inp.nrefl[i]:
                        print 'Rejected peak id %i from grain %i (hkl: %i %i %i, limit: %f): %f' %(self.inp.id[i][j],i+1,self.inp.h[i][j],self.inp.k[i][j],self.inp.l[i][j],self.inp.fit['limit'][1],value*self.inp.nrefl[i]/g[i])
                        reject.reject(self.inp,i,j,value*self.inp.nrefl[i]/g[i])
                        
        if self.inp.fit['goon'] == 'final':
            self.inp.residual = []
            self.inp.volume = []
            for i in range(self.inp.no_grains):
                self.inp.residual.append([])
                self.inp.volume.append([])
                for j in range(self.inp.nrefl[i]):
                    self.inp.residual[i].append(1)
                    self.inp.volume[i].append(1)
            reject.residual(self.inp,1,only=self.inp.fit['poor'])
            reject.intensity(self.inp)
            reject.merge(self.inp)
            reject.multi(self.inp)
        
                        
        for i in range(self.inp.no_grains):
            if self.inp.nrefl[i] < 12 and i+1 not in self.inp.fit['skip']:
                self.inp.fit['skip'].append(i+1)
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
            elif entries[0]=='c' and self.inp.fit['center'] != 0:
                self.m.fixed[entries] = False
            elif 'L' in entries and self.inp.fit['L'] != 0:
                self.m.fixed[entries] = False

		
    def fiteulergrain(self,i):
	"""
	Set tolerance and fixed parameters for fit of orientations and selected global parameters
	"""
        self.mg.tol = self.inp.fit['tol_euler']
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for angles in self.grains[i]:
            if ('phi' in angles or 'PHI' in angles) and self.inp.fit['euler'] != 0:
                self.mg.fixed[angles] = False


    def fitxyzgrain(self,i):
	"""
	Set tolerance and fixed parameters for fit of orientations and selected global parameters
	"""
        self.mg.tol = self.inp.fit['tol_xyz']
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for pos in self.grains[i]:
            if (pos[0] == 'x' or pos[0] == 'y' or pos[0] == 'z') and self.inp.fit['xyz'] != 0:
                self.mg.fixed[pos] = False


    def fitepsgrain(self,i):
	"""
	Set tolerance and fixed parameters for fit of orientations and selected global parameters
	"""
        self.mg.tol = self.inp.fit['tol_eps']
        for entries in self.mg.fixed:
            self.mg.fixed[entries] = True

        for strain in self.grains[i]:
            if 'eps' in strain and self.inp.fit['eps'] != 0:
                self.mg.fixed[strain] = False

				
    def fitgrain(self,i):
	"""
	Set tolerance and fixed parameters for fit of orientations and selected global parameters
	"""
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
            elif (entries[1]=='h' or entries[1]=='H') and self.inp.fit['euler'] != 0:
                self.mg.fixed[entries] = False
				

		
		
					
