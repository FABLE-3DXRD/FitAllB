

import numpy as n
from . import check_input
from . import write_output
from . import reject
from . import fit
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
#    """
#    Carry out one refinement cycle according to the order given by self.inp.fit['reforder']
#    Reject reflection according to self.inp.fit['limit'][1]
#    Print and save refinement and rejection info and parameters
#    
#    Jette Oddershede, Risoe DTU, May 15 2008
#    """
        
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
        
        try:
            for entries in self.inp.values:
                self.inp.fitarg["fix_%s" %entries] = True
        except:
            pass

        # carry out refinement
        if self.ref == True:
            print('\n\n*****Now fitting %s*****' %self.inp.fit['goon'])
            print('newreject_grain', self.inp.fit['newreject_grain'])
            # calculate starting values
            t1 = time.clock()
            self.m = Minuit(fcn.FCN,errordef=1,pedantic=False,print_level=-1,**self.inp.fitarg)
            try:
                self.m.values = self.inp.values
                self.m.errors = self.inp.errors
            except:
                self.m.fitarg = self.inp.fitarg
            self.mg = Minuit(fcn.FCNgrain,errordef=1,pedantic=False,print_level=-1,**self.inp.fitarg)
            try:
                self.mg.values = self.inp.values
                self.mg.errors = self.inp.errors
            except:
                self.mg.fitarg = self.inp.fitarg
            g = fit.grain_values(self)
            self.g_old = deepcopy(g)
            self.fval = sum(g)
            print('\n%s starting value %e' %(self.inp.fit['goon'],self.fval))
            try:
                for entries in self.mg.fixed:
                    self.mg.fixed[entries] = True
            except:
                pass
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
                            self.mg.fitarg['fix_i'] = True
                            self.mg.fitarg['error_i'] = 1.
                        except:
                            pass
                        self.fitglobalgrain(i)
                        print('\rRefining grain %i' %(i+1), end=' ')
                        sys.stdout.flush()
                        try:
                            errordef1 = self.g_old[i]/(3*self.inp.nrefl[i]-list(self.mg.fitarg.values()).count(False))   # Best alternative to scale_errors which is not possible by changing up with hesse in iminiut
                            self.mg = Minuit(fcn.FCNgrain,errordef=errordef1,pedantic=False,print_level=-1,**self.mg.fitarg)
                            self.mg.tol = self.mg.tol/errordef1
                        except:
                            pass
                        self.mg.migrad()
                        g[i] = self.mg.fval
                        temp = []
                        for j in range(len(self.globals)):
                            temp.append(self.mg.values[self.globals[j]])
                        global_parameters.append(temp)
                        #weight.append(len(self.inp.mean_ia[i])/n.sum(self.inp.mean_ia[i]))
#                        weight.append(float(sum(self.inp.volume[i])/self.inp.nrefl[i])) #mean
                        weight.append(float(reject.median(self.inp.volume[i])))         #median
                        self.m = self.mg
                        write_output.write_global(self)
                
            self.time = time.clock()-t1
            print('\nFit %s time %i s' %(self.inp.fit['goon'],self.time))
            self.fval = sum(g)
            print('Fit %s value %e \n' %(self.inp.fit['goon'],self.fval))
                
            # get global parameters and errors on these as average and spread of global_parameters
            global_parameters = n.array(global_parameters)
            weight = n.array(weight)
            weight = weight/n.sum(weight)
            average_global_parameters = n.zeros(len(self.globals))
            spread_global_parameters = n.zeros(len(self.globals))
            for i in range(len(global_parameters)):
                average_global_parameters = average_global_parameters + global_parameters[i] * weight[i] 
            for i in range(len(global_parameters)):
                spread_global_parameters = spread_global_parameters + (global_parameters[i] - average_global_parameters)**2 * weight[i] 
            spread_global_parameters = spread_global_parameters**.5
            
#            print global_parameters
#            print weight
#            print average_global_parameters
#            print spread_global_parameters
            
            for j in range(len(self.globals)):
                self.mg.values[self.globals[j]] = average_global_parameters[j]
                self.mg.errors[self.globals[j]] = spread_global_parameters[j]
            
            # reject outliers and save cycle info   
            self.inp.errors = self.mg.errors 
            self.inp.values = self.mg.values 
            try:
                self.m.errors = self.mg.errors
                self.m.values = self.mg.values
            except:
                for entries in list(self.mg.fitarg.keys()):
                    if "_i" in entries or entries == "i":
                        self.mg.fitarg.__delitem__(entries)
                self.m = Minuit(fcn.FCN,errordef=1,pedantic=False,print_level=-1,**self.mg.fitarg)
                self.inp.fitarg = self.m.fitarg
            fit.reject_outliers(self)
            write_output.write_rej(self.inp,message=self.inp.fit['goon'])
            write_output.write_log(self)
            write_output.write_par(self)

       
        # move onto next refinement given by the reforder list  
        self.inp.fit['goon'] = self.inp.fit['reforder'][self.inp.fit['reforder'].index(self.inp.fit['goon'])+1]
    
        return
        
        
    def fitglobalgrain(self,i):
#    """
#    Set tolerance and fixed parameters for preliminary fit of the global parameters
#    """
        self.mg.tol = self.inp.fit['tol_global']
        try:
            for entries in self.mg.fixed:
                if entries=='wy' and self.inp.fit['w'] != 0:
                    self.mg.fixed[entries] = False
                elif entries[0]=='t' and self.inp.fit['tilt'] != 0:
                    self.mg.fixed[entries] = False
                elif 'p' in entries and len(entries) == 2 and self.inp.fit['pixel'] != 0:
                    self.mg.fixed[entries] = False
                elif entries=='cy' and self.inp.fit['center'] != 0:
                    self.mg.fixed[entries] = False
                elif 'L' in entries and self.inp.fit['L'] != 0:
                    self.mg.fixed[entries] = False
        except:
            for entries in self.mg.fitarg:
                if "fix_" in entries:
                    if entries=='fix_wy' and self.inp.fit['w'] != 0:
                        self.mg.fitarg[entries] = False
                    elif entries[0:5]=='fix_t' and self.inp.fit['tilt'] != 0:
                        self.mg.fitarg[entries] = False
                    elif 'fix_p' in entries and self.inp.fit['pixel'] != 0:
                        self.mg.fitarg[entries] = False
                    elif entries=='fix_cy' and self.inp.fit['center'] != 0:
                        self.mg.fitarg[entries] = False
                    elif 'fix_L' in entries and self.inp.fit['L'] != 0:
                        self.mg.fitarg[entries] = False
                    else:                               
                        self.mg.fitarg[entries] = True

                
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
        from FitAllB import fitgg
        lsqr = fitgg.fit_minuit(inp)
        lsqr.refine()
        check_input.copy_globals(inp)
   

                    
