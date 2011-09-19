#!/usr/bin/env python

#
# Checking input  
#
import ImageD11.columnfile as ic
from string import split
import sys, os 
import write_output 
import conversion
from xfab import tools
from xfab import symmetry
from xfab import detector
import numpy as n
import logging
import minuit
from copy import deepcopy


logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')

def interrupt(killfile):
    if killfile is not None and os.path.exists(killfile):
        print 'keyboard interrupt'
        raise KeyboardInterrupt
            
class parse_input:
    def __init__(self,input_file = None):
        self.filename = input_file
        self.files = {}
        self.fit = {}
        self.param = {}
        self.param_near = {}
        
        self.needed_items = {
                    'w_step': 'Missing input: omega step size in deg',
                    'crystal_system': 'Missing input: crystal_system\n possible values:\n isotropic, cubic, hexagonal, trigonal_high, trigonal_low,\n tetragonal_high, tetragonal_low, orthorhombic, monoclinic, triclinic',
                    'log_file' : 'Missing input: grainspotter log file',
                    'flt_file' : 'Missing input: peaksearch filtered peaks file',
                    'par_file' : 'Missing input: ImageD11 detector.par file'
                    }
        self.optional_items = {
            #misc options
            'title': 'Title',
            'bg': 1000,
            'near_bg': 100,
            'const': 1,
            'near_const': 1,
            #fitgloball and fitglobalgrain options
            'cycle': 20,
            'near_cycle': 20,
            #optional files
            'structure_file': None,
            'near_flt_file': None,
            'near_par_file': None,
            # Optional space group
            'sgno' : 1,
            #experimental setup
            'dety_size': 2048,
            'detz_size': 2048,
            'near_dety_size': 1536,
            'near_detz_size': 1024,
            'w_limit': None,
            'beampol_factor':1.,
            'beampol_direct':0.0,
            #fit parameters
            'w': 0,
            'tilt': 0,
            'pixel': 0,
            'center': 0,
            'L': 0,
            'd0': 0,
            'rod': 1,
            'xyz': 1,
            'eps': 1,
            'constrx': 0,
            'constry': 0,
            'constrz': 0,
            'near_w': 0,
            'near_tilt': 0,
            'near_pixel': 0,
            'near_center': 0,
            'near_L': 0,
            'near_rod': 0,
            'near_xyz': 0,
            'near_eps': 0,
            #outlier rejection
            'skip': [],
            'rej_vol': 5,
            'rej_resmean': 10,
            'rej_resmedian': 5,
            'rej_ia': 0.2,
            'rej_multi': 1,
            'min_refl': 9,
            'overlap': 1.,
            'near_rej_vol': 5,
            'near_rej_resmean': 10,
            'near_rej_resmedian': 5,
            'near_rej_ia': 0.5,
            'near_rej_multi': 1,
            'near_min_refl': 3,
            'overlap': 1.,
            #tolerances
            'tol_global': 1e-2,
            'tol_rotpos': 1e-2,
            'tol_grain': 1e-3,
            'tol_fw_proj': 2,
            #resume refinement option
            'res_file': None,
            'rej_file': None,
            'near_rej_file': None,
            #strain to stress conversion
            'stress': 0,
            'c11': None,
            'c12': None,
            'c13': None,
            'c14': None,
            'c15': None,
            'c16': None,
            'c22': None,
            'c23': None,
            'c24': None,
            'c25': None,
            'c26': None,
            'c33': None,
            'c34': None,
            'c35': None,
            'c36': None,
            'c44': None,
            'c45': None,
            'c46': None,
            'c55': None,
            'c56': None,
            'c66': None
            }

        self.newreject = 0
        self.fit['outliers'] = 0
        self.fit['rejectgrain'] = []
        self.fit['rejectdet'] = []
        self.fit['rejectid'] = []
        self.fit['hh'] = []
        self.fit['kk'] = []
        self.fit['ll'] = []
        self.fit['rejectvalue'] = []
            
    def read(self):     
        try:
            f = open(self.filename,'r')
        except IOError:
            logging.error('No file named %s' %self.filename)
            raise IOError
        
        input = f.readlines()
        f.close()

        for lines in input:
            if lines.find('#') != 0:
                if lines.find('#') > 0:
                    lines = split(lines,'#')[0]
                line = split(lines)
                if len(line) != 0:
                    key = line[0]
                    val = line[1:]

                    valtmp = '['
                    if key == 'title':
                        valtmp = ''
                        for i in val:
                            valtmp = valtmp + i + ' '
                        val = valtmp
                    elif len(val) > 1 or key == 'skip':
                        for i in val:
                            valtmp = valtmp + i +','
                            
                        val = valtmp + ']'
                    else:
                        val = val[0]

                    # save input file names in self.files and fitting info in self.fit
                    if 'file' in key:
                        self.files[key] = val
                    else:
                        try:
                            self.fit[key] = eval(val)
                        except:
                            self.fit[key] = val
        
        stem = split(self.filename,'.')[0]      
        self.fit['stem'] = stem
        self.fit['direc'] = deepcopy(stem)
        try:
            print self.fit['title']
        except:
            pass
            
                
                        
    def check(self):
        # Needed items
        self.missing = False
        for item in self.needed_items:
            if item == 'log_file':
                if item not in self.files and 'res_file' not in self.files:
                    print 'Either log_file or res_file should be supplied.'
                    print 'If both are given, res_file takes priority.'
                    self.missing = True
            else:
                if item not in self.files:
                    if item not in self.fit:
                        print self.needed_items[item]
                        self.missing = True

            
            
    def initialize(self): 
        # Does output directory exist?
        if not os.path.exists(self.fit['stem']):
            os.mkdir(self.fit['stem'])
        sys.path.insert(0,self.fit['stem'])
        #print sys.path[0]
        #sys.exit()
        # Set default options
        for item in self.optional_items:
            if 'file' in item and item not in self.files:
                self.files[item] = self.optional_items[item] 
            elif 'file' not in item and item not in self.fit:
                self.fit[item] = self.optional_items[item] 
                
        # calculate stiffness tensor
        if self.fit['stress'] == 0:
            print 'No strain to stress conversion'
            self.C = n.zeros((6,6))
        else:
            print 'Both strain and stress given in output'
            self.C = conversion.formStiffnessMV(self.fit['crystal_system'],
                                                c11=self.fit['c11'],c12=self.fit['c12'],c13=self.fit['c13'],c14=self.fit['c14'],c15=self.fit['c15'],c16=self.fit['c16'],
                                                                    c22=self.fit['c22'],c23=self.fit['c23'],c24=self.fit['c24'],c25=self.fit['c25'],c26=self.fit['c26'],
                                                                                        c33=self.fit['c33'],c34=self.fit['c34'],c35=self.fit['c35'],c36=self.fit['c36'],
                                                                                                            c44=self.fit['c44'],c45=self.fit['c45'],c46=self.fit['c46'],
                                                                                                                                c55=self.fit['c55'],c56=self.fit['c56'],
                                                                                                                                                    c66=self.fit['c66'])
        #print self.C                       
           
           
    def read_par(self,par_file): # read detector.par
        try:
            f=open(par_file,'r')
        except IOError:
            logging.error('par_file: no file named %s' %par_file)
            raise IOError
        
        input = f.readlines()
        f.close()

        for lines in input:
            if lines.find('#') != 0:
                if lines.find('#') > 0:
                    lines = split(lines,'#')[0]
                line = split(lines)
                if len(line) != 0:
                    key = line[0]
                    val = line[1]

                    # evaluate and store parameters in self.param 
                    try:
                        self.param[key] = eval(val)
                    except:
                        self.param[key] = val
                        
        self.unit_cell = n.array([self.param['cell__a'],self.param['cell__b'],self.param['cell__c'],self.param['cell_alpha'],self.param['cell_beta'],self.param['cell_gamma']])
        self.param['unit_cell'] = self.unit_cell
        (dety_center, detz_center) = detector.xy_to_detyz([self.param['z_center'],self.param['y_center']],
                                                          self.param['o11'],self.param['o12'],self.param['o21'],self.param['o22'],
                                                          self.fit['dety_size'],self.fit['detz_size'])
        self.param['y_center'] = dety_center
        self.param['z_center'] = detz_center
        #read Jons wedge convention and convert to Sorens
        self.param['wedge'] = -1.*self.param['wedge']
                
                
    def read_flt(self,flt_file): # read peaks_t##.flt and calculate experimental variances Sww,Syy,Szz
        # create parameters, must be lists in order to append
        try:
            f=open(flt_file,'r')
        except IOError:
            logging.error('flt_file: no file named %s' %flt_file)
            raise IOError
        
        input = f.readlines()
        f.close()
        
        #read as columnfile to avoid problems if peaksearch output is changed
        flt = ic.columnfile(flt_file)
        self.int = flt.getcolumn('sum_intensity')
        intmax = flt.getcolumn('IMax_int')
        pixels = flt.getcolumn('Number_of_pixels')
        
        self.w = flt.getcolumn('omega')
        sigw = flt.getcolumn('sigo')
        spot = flt.getcolumn('spot3d_id')        
        sc = flt.getcolumn('sc')
        fc = flt.getcolumn('fc')
        sigs = flt.getcolumn('sigs')
        sigf = flt.getcolumn('sigf')
        covsf = flt.getcolumn('covsf')
        self.dety = []
        self.detz = []
        sigy = []
        sigz = []
        covyz = covsf*n.linalg.det([[self.param['o11'],self.param['o12']],[self.param['o12'],self.param['o22']]])
        
        smin = flt.getcolumn('Min_s')
        smax = flt.getcolumn('Max_s')
        sstep = smax - smin + 1
        fmin = flt.getcolumn('Min_f')
        fmax = flt.getcolumn('Max_f')
        fstep = fmax - fmin + 1
        wmin = flt.getcolumn('Min_o')
        wmax = flt.getcolumn('Max_o')
        wstep = wmax - wmin + self.fit['w_step']
        ystep = []
        zstep = []
        ymin = []
        ymax = []
        zmin = []
        zmax = []
        eta = n.zeros(flt.nrows)
        rr = n.zeros(flt.nrows)
        self.labels = None
        if self.files['res_file'] != None:
            try:
                self.labels = flt.getcolumn('labels')
                if 0 in self.labels:    
                    self.labels = self.labels+1
                h = flt.getcolumn('h')
                k = flt.getcolumn('k')
                l = flt.getcolumn('l')
                self.tth = flt.getcolumn('tth_per_grain')
                self.eta = flt.getcolumn('eta_per_grain')
                self.nrefl = []
                self.id = []
                self.h = []
                self.k = []
                self.l = []
                for no in self.grainno:
                    idgr = []
                    hgr = []
                    kgr = []
                    lgr = []
                    for i in range(flt.nrows):
                        if self.labels[i] == no:
                            idgr.append(int(spot[i]))
                            hgr.append(int(h[i]))
                            kgr.append(int(k[i]))
                            lgr.append(int(l[i]))
                    self.nrefl.append(len(idgr))
                    self.id.append(idgr)        
                    self.h.append(hgr)
                    self.k.append(kgr)
                    self.l.append(lgr)
                            
                    
            except:
                pass
     
        for i in range(flt.nrows):
            (dety,detz) = detector.xy_to_detyz([sc[i],fc[i]],
                                               self.param['o11'],
                                               self.param['o12'],
                                               self.param['o21'],
                                               self.param['o22'],
                                               self.fit['dety_size'],
                                               self.fit['detz_size'])
            self.dety.append(dety)              
            self.detz.append(detz)
            (eta[i],rr[i]) = detector.detyz_to_eta_and_radpix([dety*self.param['y_size'],detz*self.param['z_size']], 
                                                               self.param['y_center']*self.param['y_size'], 
                                                               self.param['z_center']*self.param['z_size'])

            (sz,sy) = n.dot(n.array([[abs(self.param['o11']),abs(self.param['o12'])],[abs(self.param['o21']),abs(self.param['o22'])]]),
                               n.array([sigs[i],sigf[i]]))
            sigy.append(sy)             
            sigz.append(sz)
            
            
            (zs,ys) = n.dot(n.array([[abs(self.param['o11']),abs(self.param['o12'])],[abs(self.param['o21']),abs(self.param['o22'])]]),
                               n.array([sstep[i],fstep[i]]))
            ystep.append(ys)                
            zstep.append(zs)
            (zmi,ymi) = n.dot(n.array([[abs(self.param['o11']),abs(self.param['o12'])],[abs(self.param['o21']),abs(self.param['o22'])]]),
                               n.array([smin[i],fmin[i]]))
            ymin.append(ymi)                
            zmin.append(zmi)
            (zma,yma) = n.dot(n.array([[abs(self.param['o11']),abs(self.param['o12'])],[abs(self.param['o21']),abs(self.param['o22'])]]),
                               n.array([smax[i],fmax[i]]))
            ymax.append(yma)                
            zmax.append(zma)

        # convert into arrays so sorting according to spotid is possible 
        Syy = n.array(sigy)**2-1
        Szz = n.array(sigz)**2-1
        Syz = n.array(covyz)*n.sqrt(Syy*Szz)
        Srr =  (self.param['z_size']**4*(self.detz-self.param['z_center'])**2*Szz+
                self.param['y_size']**4*(self.dety-self.param['y_center'])**2*Syy+
                2*self.param['z_size']**2*self.param['y_size']**2*(self.dety-self.param['y_center'])*(self.detz-self.param['z_center'])*Syz)/rr**2
        self.sig_eta = n.zeros(flt.nrows)
        for i in range(flt.nrows):
            if -2*(self.detz[i]-self.param['z_center'])*(self.dety[i]-self.param['y_center'])*Syz[i] >0:
                self.sig_eta[i] = 180/n.pi/rr[i]**2*self.param['y_size']*self.param['z_size']*n.sqrt((self.detz[i]-self.param['z_center'])**2*Syy[i]+(self.dety[i]-self.param['y_center'])**2*Szz[i]-2*(self.detz[i]-self.param['z_center'])*(self.dety[i]-self.param['y_center'])*Syz[i])
            else:
                self.sig_eta[i] = 180/n.pi/rr[i]**2*self.param['y_size']*self.param['z_size']*n.sqrt((self.detz[i]-self.param['z_center'])**2*Syy[i]+(self.dety[i]-self.param['y_center'])**2*Szz[i])
        self.sig_tth = n.zeros(flt.nrows)
        for i in range(flt.nrows):
            if Srr[i] > 0:
                self.sig_tth[i] = 180*self.param['distance']/n.pi/(self.param['distance']**2+rr[i]**2)*n.sqrt(Srr[i])
            else:
                self.sig_tth[i] = 0

        self.dety = n.array(self.dety)
        self.detz = n.array(self.detz)  
        sigy = n.array(sigy)
        sigz = n.array(sigz)
        ymin = n.array(ymin)
        ymax = n.array(ymax)
        zmin = n.array(zmin)
        zmax = n.array(zmax)
        # do the sorting
        index = n.argsort(spot)
        self.w = self.w[index]
        self.dety = self.dety[index]
        self.detz = self.detz[index]
        sigw = sigw[index]
        sigy = sigy[index]
        sigz = sigz[index]
        self.int = self.int[index]
        intmax = intmax[index]
        wmin = wmin[index]
        wmax = wmax[index]
        ymin = ymin[index]
        ymax = ymax[index]
        zmin = zmin[index]
        zmax = zmax[index]
        self.sig_eta = self.sig_eta[index]
        self.sig_tth = self.sig_tth[index]
        try:
            self.tth = self.tth[index]
            self.eta = self.eta[index]
            self.labels = self.labels[index]
        except:
            pass

        # we now have arrays on length len(self.spot), but in the end we want lists of length maxspotno+1
        # therefore create temporary lists with zero values of length maxspotno+1
        self.param['total_refl'] = int(max(spot))+1 #necessary if spots are not consequtively numbered
        tw = [0.]*self.param['total_refl']
        twmin = [0.]*self.param['total_refl']
        twmax = [0.]*self.param['total_refl']
        tymin = [0.]*self.param['total_refl']
        tymax = [0.]*self.param['total_refl']
        tzmin = [0.]*self.param['total_refl']
        tzmax = [0.]*self.param['total_refl']
        tdety = [0.]*self.param['total_refl']
        tdetz = [0.]*self.param['total_refl']
        tsigw = [-1.]*self.param['total_refl']
        tsigy = [-1.]*self.param['total_refl']
        tsigz = [-1.]*self.param['total_refl']
        tint = [0.]*self.param['total_refl']
        tintmax = [1.]*self.param['total_refl']
        tsigeta = [1.]*self.param['total_refl']
        tsigtth = [1.]*self.param['total_refl']
        ttth = [-1.]*self.param['total_refl']
        teta = [-1.]*self.param['total_refl']
        tlabels = [-1.]*self.param['total_refl']
        missing = 0
        # update temporary lists for all read reflections
        for i in range(self.param['total_refl']):
            if i in spot:
                tw[i] = self.w[i-missing]
                twmin[i] = wmin[i-missing]
                twmax[i] = wmax[i-missing]
                tymin[i] = ymin[i-missing]
                tymax[i] = ymax[i-missing]
                tzmin[i] = zmin[i-missing]
                tzmax[i] = zmax[i-missing]
                tdety[i] = self.dety[i-missing]
                tdetz[i] = self.detz[i-missing]
                tsigw[i] = sigw[i-missing]
                tsigy[i] = sigy[i-missing]
                tsigz[i] = sigz[i-missing]
                tint[i] = self.int[i-missing]
                tintmax[i] = intmax[i-missing]
                tsigeta[i] = self.sig_eta[i-missing]
                tsigtth[i] = self.sig_tth[i-missing]
                try:
                    ttth[i] = self.tth[i-missing]
                    teta[i] = self.eta[i-missing]
                    tlabels[i] = self.labels[i-missing]
                except:
                    pass
            else:
                missing = missing+1

        #copy temporary lists to variables        
        self.w = tw
        self.dety = tdety
        self.detz = tdetz
        sigw = tsigw
        sigy = tsigy
        sigz = tsigz
        self.int = tint
        intmax = tintmax
        wmin = twmin
        wmax = twmax
        ymin = tymin
        ymax = tymax
        zmin = tzmin
        zmax = tzmax
        self.sig_eta = tsigeta
        self.sig_tth = tsigtth
        try:
            self.tth = ttth
            self.eta = teta
            self.labels = tlabels
        except:
            pass
        
        if self.fit['w_limit'] == None:
            self.fit['w_limit'] = [min(self.w),max(self.w)]
        else:
            assert len(self.fit['w_limit']) % 2 == 0, 'An even number of omega-limits must be given'
            self.fit['w_limit'].sort()

        # set default variances
        self.Sww = [self.fit['w_step']**2/12.]*self.param['total_refl']
        self.Syy = [1.]*self.param['total_refl']
        self.Szz = [1.]*self.param['total_refl']
        for j in range(len(self.int)):

            if self.int[j] > 0:
                if intmax[j] > 2**16-2*self.fit['bg']:# and intmax[j]*pixels[j]/self.int[j] < 8.:
                    self.Sww[j] = -10
                else:
                    if sigw[j] > 1:
#                        self.Sww[j] = (.01+1./self.int[i])*(sigw[j]**2-1) + (self.fit['w_step']**2/12.)*(self.fit['w_step']/wstep[i])**2
                        self.Sww[j] = (1./self.int[i])*(sigw[j]**2-1) 
                    if sigy[j] > 1:
#                        self.Syy[j] = (.001+1./self.int[i])*(sigy[j]**2-1) + (1./12.)*(1./ystep[i])**2
                        self.Syy[j] = (1./self.int[i])*(sigy[j]**2-1) 
                    if sigz[j] > 1:
#                        self.Szz[j] = (.001+1./self.int[i])*(sigz[j]**2-1) + (1./12.)*(1./zstep[i])**2
                        self.Szz[j] = (1./self.int[i])*(sigz[j]**2-1) 
                for limit in self.fit['w_limit']:
                    if abs(limit-wmin[j]) < self.fit['w_step'] or abs(limit-wmax[j]) < self.fit['w_step']:
                        self.Sww[j] = -100
                if ymin[j] < 2 or ymax[j] > self.fit['dety_size']-2:
                        self.Syy[j] = -100
                if zmin[j] < 2 or zmax[j] > self.fit['detz_size']-2:
                        self.Szz[j] = -100
                        
        if self.files['res_file'] != None:
            try:
        # calculate self.F2vol which is the intensity divided the Lorentz factor, thus the squared structure factor times the volume
                self.F2vol = [0]*self.param['total_refl']
                for i in range(self.param['total_refl']):
                    rho = n.pi/2.0 + self.eta[i]*n.pi/180. + self.fit['beampol_direct']*n.pi/180.0 
                    P = 0.5 * (1. + n.cos(self.tth[i]*n.pi/180.)**2 + self.fit['beampol_factor']*n.cos(2*rho)*n.sin(self.tth[i]*n.pi/180.)**2)
                    Linv = (n.sin(self.tth[i]*n.pi/180.)*abs(n.sin(self.eta[i]*n.pi/180.)))
                    self.F2vol[i] = self.int[i]*Linv/P
                self.param['theta_min'] = min(self.tth)/2.
                self.param['theta_max'] = max(self.tth)/2.
            except:
                pass
                        

        
 
    def read_log(self): # read grainspotter.log
        self.nrefl = []
        self.rod = []
        self.h = []
        self.k = []
        self.l = []
        self.id = []
        self.x = []
        self.y = []
        self.z = []
        self.eta = [0]*self.param['total_refl']
        self.tth = [0]*self.param['total_refl']
        ia = []
        B = tools.form_b_mat(self.unit_cell)
        Binv = n.linalg.inv(B)
        
        try:
            f=open(self.files['log_file'],'r')
        except IOError:
            logging.error('log_file: no file named %s' %self.files['log_file'])
            raise IOError
        
        input = f.readlines()
        f.close()

        self.no_grains = int(split(input[0])[1])
        self.grainno = range(1,self.no_grains+1)
        nn = 23 # jumping to first grain

        for gr in range(self.no_grains):
            nn = nn + 1
            self.nrefl.append(int(split(input[nn])[1]))
            nn = nn + 1
            # read grain positions from new grainspotter output
            if len(split(input[nn])) >= 4: 
                ia.append(eval(split(input[nn])[0]))
                self.x.append(eval(split(input[nn])[1]))
                self.y.append(eval(split(input[nn])[2]))
                self.z.append(eval(split(input[nn])[3]))
            nn = nn + 9
            self.rod.append([eval(split(input[nn])[0]),eval(split(input[nn])[1]),eval(split(input[nn])[2])])
            nn = nn + 5
#            nn = nn + 2
#            self.euler.append([eval(split(input[nn])[0]),eval(split(input[nn])[1]),eval(split(input[nn])[2])])
#            nn = nn + 3
            idgr = []
            h = []
            k = []
            l = []
            for refl in range(self.nrefl[gr]):
                nn = nn + 1
                idgr.append(int(split(input[nn])[2]))
                h.append(int(split(input[nn])[3]))
                k.append(int(split(input[nn])[4]))
                l.append(int(split(input[nn])[5]))
                self.tth[int(split(input[nn])[2])]=float(split(input[nn])[12])
                self.eta[int(split(input[nn])[2])]=float(split(input[nn])[18])

            self.id.append(idgr)                
            self.h.append(h)                
            self.k.append(k)                
            self.l.append(l)
            nn = nn + 2
        
        # calculate self.F2vol which is the intensity divided the Lorentz factor, thus the squared structure factor times the volume
        self.F2vol = [0]*self.param['total_refl']
        for i in range(self.param['total_refl']):
            rho = n.pi/2.0 + self.eta[i]*n.pi/180. + self.fit['beampol_direct']*n.pi/180.0 
            P = 0.5 * (1. + n.cos(self.tth[i]*n.pi/180.)**2 + self.fit['beampol_factor']*n.cos(2*rho)*n.sin(self.tth[i]*n.pi/180.)**2)
            Linv = (n.sin(self.tth[i]*n.pi/180.)*abs(n.sin(self.eta[i]*n.pi/180.)))
            self.F2vol[i] = self.int[i]*Linv/P
        self.param['theta_min'] = min(self.tth)/2.
        self.param['theta_max'] = max(self.tth)/2.
        
        # delete grains with an internal angle above the set threshold
#        for i in range(self.no_grains):
#            if ia[i] > self.fit['ia']:
#                self.fit['skip'].append(i+1)
                
        # Check for equal grains and convert orientations to the fundamental zone if crystal system is given
        # Fundamental zone definition: trace(U) is maximal
        # NB! Only possible for U1*B*hkl1 = U2*B*hkl2, where U2 = U1*p[k] and hkl2 = pinv[k]*hkl1 if either B or p is diagonal, hence the checks
        cs = 1
        if 'isotropic' in self.fit['crystal_system'] or 'cubic' in self.fit['crystal_system']:
            cs = 7
        elif 'hexagonal' in self.fit['crystal_system']:
            cs = 6
        elif 'trigonal' in self.fit['crystal_system']:
            cs = 5
        elif 'tetragonal' in self.fit['crystal_system']:
            cs = 4
        elif 'orthorhombic' in self.fit['crystal_system']:
            cs = 3
        elif 'monoclinic' in self.fit['crystal_system']:
            cs = 2
        for i in range(self.no_grains):
            Ui = tools.rod_to_u([self.rod[i][0],self.rod[i][1],self.rod[i][2]])
            p = symmetry.permutations(cs)
            t = Ui.trace()
            Ut = Ui.copy()
            pt = n.eye(3,3) 
            for k in range(len(p)):
                Urot = n.dot(Ui,n.dot(B,n.dot(p[k],Binv)))
                trace = Urot.trace()
                if trace > t:
                    t = trace
                    Ut = Urot
                    pt = n.linalg.inv(p[k])
            for j in range(self.nrefl[i]):
                [self.h[i][j],self.k[i][j],self.l[i][j]] = n.dot(pt,n.array([self.h[i][j],self.k[i][j],self.l[i][j]]))
            [self.rod[i][0],self.rod[i][1],self.rod[i][2]] = tools.u_to_rod(Ut)
            
            for j in range(i):
                Uj = tools.rod_to_u([self.rod[j][0],self.rod[j][1],self.rod[j][2]])
                Umis = symmetry.Umis(Ui,Uj,1)
                mis = 180.
                for k in range(len(Umis)):
                    if Umis[k][1] < mis:
                        mis = Umis[k][1]
                if mis < 5:
                    dist = n.sqrt((self.x[i]-self.x[j])**2+(self.y[i]-self.y[j])**2)
                    if  dist < 0.1:
                        print i+1,j+1,mis,self.x[i],self.y[i],self.z[i],self.x[j],self.y[j],self.z[j],dist
            
            
        
            

    def read_res(self): # read file of positions, orientations and strain tensor components to resume refinement
        if self.files['res_file'] != None:
            try:
                f=open(self.files['res_file'],'r')
                f.close()
                print 'Resume refinement'
                res = ic.columnfile(self.files['res_file'])
                self.grainno = res.getcolumn('grainno')
                self.grainno = self.grainno.astype(n.int)
                if 0 in self.grainno:
                    self.grainno = self.grainno + 1
                self.grainno = self.grainno.tolist()
                self.x = res.getcolumn('x')
                self.y = res.getcolumn('y')
                self.z = res.getcolumn('z')
                U11 = res.getcolumn('U11')            
                U12 = res.getcolumn('U12')            
                U13 = res.getcolumn('U13')            
                U21 = res.getcolumn('U21')            
                U22 = res.getcolumn('U22')            
                U23 = res.getcolumn('U23')            
                U31 = res.getcolumn('U31')            
                U32 = res.getcolumn('U32')            
                U33 = res.getcolumn('U33')
                U = n.zeros((len(U11),3,3))
                for i in range(len(U11)):
                    U[i][0][0] = U11[i]
                    U[i][0][1] = U12[i]
                    U[i][0][2] = U13[i]
                    U[i][1][0] = U21[i]
                    U[i][1][1] = U22[i]
                    U[i][1][2] = U23[i]
                    U[i][2][0] = U31[i]
                    U[i][2][1] = U32[i]
                    U[i][2][2] = U33[i]
                self.rodx = n.zeros((len(U11)))
                self.rody = n.zeros((len(U11)))
                self.rodz = n.zeros((len(U11)))
                for i in range(len(U11)):
                    [self.rodx[i],self.rody[i],self.rodz[i]] = tools.u_to_rod(U[i])
                try:
                    self.eps11 = res.getcolumn('eps11')
                    self.eps22 = res.getcolumn('eps22')
                    self.eps33 = res.getcolumn('eps33')
                    self.eps23 = res.getcolumn('eps23')
                    self.eps13 = res.getcolumn('eps13')
                    self.eps12 = res.getcolumn('eps12')
                except:
                    self.eps11 = [0.0]*len(self.grainno)
                    self.eps22 = [0.0]*len(self.grainno)
                    self.eps33 = [0.0]*len(self.grainno)
                    self.eps23 = [0.0]*len(self.grainno)
                    self.eps13 = [0.0]*len(self.grainno)
                    self.eps12 = [0.0]*len(self.grainno)
                    
                try:
                    self.ia = res.getcolumn('mean_IA')
                    self.grainvolume = res.getcolumn('grainvolume')
                    self.eps11_s = res.getcolumn('eps11_s')
                    self.eps22_s = res.getcolumn('eps22_s')
                    self.eps33_s = res.getcolumn('eps33_s')
                    self.eps23_s = res.getcolumn('eps23_s')
                    self.eps13_s = res.getcolumn('eps13_s')
                    self.eps12_s = res.getcolumn('eps12_s')
                    self.sig11_s = res.getcolumn('sig11_s')
                    self.sig22_s = res.getcolumn('sig22_s')
                    self.sig33_s = res.getcolumn('sig33_s')
                    self.sig23_s = res.getcolumn('sig23_s')
                    self.sig13_s = res.getcolumn('sig13_s')
                    self.sig12_s = res.getcolumn('sig12_s')  
                except:
                    pass
            except IOError:
                logging.error('res_file: no file named %s' %self.files['res_file'])
                raise IOError
        else:
            print 'Start refinement from scratch' 
            return


                
    def read_rej(self,rej_file): # read file containing rejected peaks to resume refinement
        try:
            f=open(rej_file,'r')
        except:
            print 'Start refinement without apriori information about peak rejection' 
            return
        print 'Use apriori information about peak rejections' 
        
        input = f.readlines()
        f.close()

        # build rejection list in rigth format
        rejectid = []
        for i in range(self.no_grains):
            rejectid.append([])
                    
        # read parameters by appending
        for line in input:
            if 'Rejected peak id' in line:
                try:
                    rejectid[int(split(line)[7])-1].append(int(split(line)[4]))
                    self.fit['rejectid'].append(int(split(line)[4]))
                    self.fit['rejectgrain'].append(int(split(line)[7]))
                    self.fit['hh'].append(int(split(line)[9]))
                    self.fit['kk'].append(int(split(line)[10]))
                    self.fit['ll'].append(int(split(line)[11]))
                    try:
                        self.fit['rejectvalue'].append(eval(split(line)[13]))
                    except:
                        self.fit['rejectvalue'].append(split(line)[13])
                    self.fit['outliers'] = self.fit['outliers'] + 1
                except:
                    pass
            if 'Skip grains' in line:
                string = ''
                for i in range(2,len(split(line))):
                    string = string+split(line)[i]
                self.fit['skip'].extend(eval(string))
        for i in range(len(self.fit['skip'])-1,-1,-1):
            if self.fit['skip'][i] > self.no_grains:
                self.fit['skip'].pop(i)
        try:
            for i in range(self.no_grains):
                for j in range(self.nrefl[i]-1,-1,-1): # loop backwards to make pop work
                    if self.id[i][j] in rejectid[i]:
                        self.id[i].pop(j)
                        self.h[i].pop(j)
                        self.k[i].pop(j)
                        self.l[i].pop(j)
                        self.nrefl[i] = self.nrefl[i] - 1
        except:
            pass # exception for analyse_reject.py


    def set_start(self): # build fcn, initiate minuit and set starting values and errors
        self.values = {}
        # grain values
        if self.files['res_file'] != None:
            self.no_grains = max(len(self.grainno),max(self.grainno))
            self.rod = []
            for i in range(self.no_grains):
                self.rod.append([0.0,0.0,0.0])
            self.param['theta_min'] = 0.0
            self.param['theta_max'] = 7.5
        for i in range(self.no_grains):
            self.values['x%s' %i] = 0.0
            self.values['y%s' %i] = 0.0
            self.values['z%s' %i] = 0.0
            self.values['epsaa%s' %i] = 0.0 
            self.values['epsab%s' %i] = 0.0  
            self.values['epsac%s' %i] = 0.0
            self.values['epsbb%s' %i] = 0.0 
            self.values['epsbc%s' %i] = 0.0 
            self.values['epscc%s' %i] = 0.0
            self.values['rodx%s' %i] = 0.0
            self.values['rody%s' %i] = 0.0
            self.values['rodz%s' %i] = 0.0
        # grain values for resuming refinement
        if self.files['res_file'] != None:
            for i in range(self.no_grains):
                if i+1 in self.grainno:
                    self.values['x%s' %i] = 1000.*self.x[self.grainno.index(i+1)]
                    self.values['y%s' %i] = 1000.*self.y[self.grainno.index(i+1)]
                    self.values['z%s' %i] = 1000.*self.z[self.grainno.index(i+1)]
                    self.values['epsaa%s' %i] = self.eps11[self.grainno.index(i+1)] 
                    self.values['epsab%s' %i] = self.eps12[self.grainno.index(i+1)] 
                    self.values['epsac%s' %i] = self.eps13[self.grainno.index(i+1)]
                    self.values['epsbb%s' %i] = self.eps22[self.grainno.index(i+1)]
                    self.values['epsbc%s' %i] = self.eps23[self.grainno.index(i+1)]
                    self.values['epscc%s' %i] = self.eps33[self.grainno.index(i+1)]
                    self.rod[i][0] = self.rodx[self.grainno.index(i+1)]
                    self.rod[i][1] = self.rody[self.grainno.index(i+1)]
                    self.rod[i][2] = self.rodz[self.grainno.index(i+1)]
                else:
                    self.fit['skip'].append(i+1)
        # else if start from scratch with new grainspotter log file use positions from this
        elif len(self.x) == self.no_grains:
            for i in range(self.no_grains):
                self.values['x%s' %i] = 1000.*self.x[i]
                self.values['y%s' %i] = 1000.*self.y[i]
                self.values['z%s' %i] = 1000.*self.z[i]
        
        for i in range(len(self.fit['skip'])-1,-1,-1):
            if self.fit['skip'][i] > self.no_grains:
                self.fit['skip'].pop(i)

        self.errors = {}
        # global errors
        self.param['a_error'] = 0.001
        self.param['b_error'] = 0.001
        self.param['c_error'] = 0.001
        self.param['alpha_error'] = 0.001
        self.param['beta_error'] = 0.001
        self.param['gamma_error'] = 0.001
        self.param['chi_error'] = 0.001
        self.param['wedge_error'] = 0.001
        self.param['tilt_x_error'] = 0.001
        self.param['tilt_y_error'] = 0.001
        self.param['tilt_z_error'] = 0.001
        self.param['y_size_error'] = 0.1
        self.param['z_size_error'] = 0.1
        self.param['y_center_error'] = 0.1
        self.param['z_center_error'] = 0.1
        self.param['distance_error']  = 10.
        self.param['i_error']  = 1
        self.param['j_error']  = 1
        # grain errors
        for i in range(self.no_grains):
            self.errors['x%s' %i] = self.param['y_size']/5.
            self.errors['y%s' %i] = self.param['y_size']/5.
            self.errors['z%s' %i] = self.param['z_size']/20.
            self.errors['epsaa%s' %i] = 0.001 
            self.errors['epsab%s' %i] = 0.001
            self.errors['epsac%s' %i] = 0.001
            self.errors['epsbb%s' %i] = 0.001
            self.errors['epsbc%s' %i] = 0.001
            self.errors['epscc%s' %i] = 0.001
            lenrod = n.linalg.norm(self.rod[i])
            ia = 0.1*n.pi/180.
            errorscale = ia*(1+lenrod*lenrod)/(lenrod*(1-0.25*ia*ia*lenrod*lenrod))
            if errorscale < 0:
                errorscale = 1
            self.errors['rodx%s' %i] = errorscale*abs(self.rod[i][0])
            self.errors['rody%s' %i] = errorscale*abs(self.rod[i][1])
            self.errors['rodz%s' %i] = errorscale*abs(self.rod[i][2])
    

        self.fit['newreject_grain'] = []
        
            
    def reject(self): # carry out initial rejections

        import reject
        print '\n\nNumber of read grains', self.no_grains
        print 'Number of assigned reflections: ', n.sum(self.nrefl)
        # set starting values
        self.newreject = 0
        self.fit['newreject_grain'] = []
        self.residual = []
        self.volume = []
        self.mean_ia = []
        self.spr_eta = []
        self.spr_tth = []
        for i in range(self.no_grains):
            self.residual.append([])
            self.volume.append([])
            self.mean_ia.append([])
            self.spr_eta.append([])
            self.spr_tth.append([])
            for j in range(self.nrefl[i]):
                self.residual[i].append(1)
                self.volume[i].append(1)
                self.mean_ia[i].append(1)
                self.spr_eta[i].append(1)
                self.spr_tth[i].append(1)
        # do the actual rejections
        reject.overflow(self)
        reject.edge(self)
        reject.intensity(self)
        reject.peak_spread(self)
        reject.mean_ia(self,2*self.fit['rej_ia'])
        reject.residual(self,self.fit['rej_resmedian'])
        reject.multi(self)
        reject.merge(self)
        #reject.friedel(self)

        
    def write_rej(self): # write the rejected peaks to rejection file
    
        import reject
        reject.unique_list(self.fit['skip'])
        print 'Skip the following grains:', self.fit['skip']
        print 'Actual number of grains in fit', self.no_grains - len(self.fit['skip'])
        observations = 0
        for i in range(self.no_grains):
            if i+1 in self.fit['skip']:
                pass
            else:
                observations = observations + self.nrefl[i]
        print 'Total number of reflections in these', observations,'\n' 
       
        write_output.write_rej(self,message=('%s\n\ncheck_input' %self.fit['title']))
        
def set_globals(inp):
        # global values
        inp.values['a'] = deepcopy(inp.param['cell__a'])
        inp.values['b'] = deepcopy(inp.param['cell__b'])
        inp.values['c'] = deepcopy(inp.param['cell__c'])
        inp.values['alpha'] = deepcopy(inp.param['cell_alpha'])
        inp.values['beta'] = deepcopy(inp.param['cell_beta'])
        inp.values['gamma'] = deepcopy(inp.param['cell_gamma'])
        inp.values['wx'] = deepcopy(inp.param['chi'])
        inp.values['wy'] = deepcopy(inp.param['wedge'])
        inp.values['tx'] = deepcopy(inp.param['tilt_x'])
        inp.values['ty'] = deepcopy(inp.param['tilt_y'])
        inp.values['tz'] = deepcopy(inp.param['tilt_z'])
        inp.values['py'] = deepcopy(inp.param['y_size'])
        inp.values['pz'] = deepcopy(inp.param['z_size'])
        inp.values['cy'] = deepcopy(inp.param['y_center'])
        inp.values['cz'] = deepcopy(inp.param['z_center'])
        inp.values['L']  = deepcopy(inp.param['distance'])
        # global errors
        inp.errors['a'] = deepcopy(inp.param['a_error'])
        inp.errors['b'] = deepcopy(inp.param['b_error'])
        inp.errors['c'] = deepcopy(inp.param['c_error'])
        inp.errors['alpha'] = deepcopy(inp.param['alpha_error'])
        inp.errors['beta'] = deepcopy(inp.param['beta_error'])
        inp.errors['gamma'] = deepcopy(inp.param['gamma_error'])
        inp.errors['wx'] = deepcopy(inp.param['chi_error'])
        inp.errors['wy'] = deepcopy(inp.param['wedge_error'])
        inp.errors['tx'] = deepcopy(inp.param['tilt_x_error'])
        inp.errors['ty'] = deepcopy(inp.param['tilt_y_error'])
        inp.errors['tz'] = deepcopy(inp.param['tilt_z_error'])
        inp.errors['py'] = deepcopy(inp.param['y_size_error'])
        inp.errors['pz'] = deepcopy(inp.param['z_size_error'])
        inp.errors['cy'] = deepcopy(inp.param['y_center_error'])
        inp.errors['cz'] = deepcopy(inp.param['z_center_error'])
        inp.errors['L']  = deepcopy(inp.param['distance_error'])
        inp.errors['i']  = deepcopy(inp.param['i_error'])
        inp.errors['j']  = deepcopy(inp.param['j_error'])
    
    
def copy_globals(inp):
        # Necessary to save copies of global parameters in param when switching between near and farfiel detectors
        # global values
        inp.param['cell__a'] = deepcopy(inp.values['a']) 
        inp.param['cell__b'] = deepcopy(inp.values['b']) 
        inp.param['cell__c'] = deepcopy(inp.values['c']) 
        inp.param['cell_alpha'] = deepcopy(inp.values['alpha']) 
        inp.param['cell_beta'] = deepcopy(inp.values['beta']) 
        inp.param['cell_gamma'] = deepcopy(inp.values['gamma']) 
        inp.param['chi'] = deepcopy(inp.values['wx']) 
        inp.param['wedge'] = deepcopy(inp.values['wy']) 
        inp.param['tilt_x'] = deepcopy(inp.values['tx']) 
        inp.param['tilt_y'] = deepcopy(inp.values['ty']) 
        inp.param['tilt_z'] = deepcopy(inp.values['tz']) 
        inp.param['y_size'] = deepcopy(inp.values['py']) 
        inp.param['z_size'] = deepcopy(inp.values['pz']) 
        inp.param['y_center'] = deepcopy(inp.values['cy']) 
        inp.param['z_center'] = deepcopy(inp.values['cz']) 
        inp.param['distance'] = deepcopy(inp.values['L'])  
        # global errors
        inp.param['a_error'] = deepcopy(inp.errors['a']) 
        inp.param['b_error'] = deepcopy(inp.errors['b']) 
        inp.param['c_error'] = deepcopy(inp.errors['c']) 
        inp.param['alpha_error'] = deepcopy(inp.errors['alpha']) 
        inp.param['beta_error'] = deepcopy(inp.errors['beta']) 
        inp.param['gamma_error'] = deepcopy(inp.errors['gamma']) 
        inp.param['chi_error'] = deepcopy(inp.errors['wx']) 
        inp.param['wedge_error'] = deepcopy(inp.errors['wy']) 
        inp.param['tilt_x_error'] = deepcopy(inp.errors['tx']) 
        inp.param['tilt_y_error'] = deepcopy(inp.errors['ty']) 
        inp.param['tilt_z_error'] = deepcopy(inp.errors['tz']) 
        inp.param['y_size_error'] = deepcopy(inp.errors['py']) 
        inp.param['z_size_error'] = deepcopy(inp.errors['pz']) 
        inp.param['y_center_error'] = deepcopy(inp.errors['cy']) 
        inp.param['z_center_error'] = deepcopy(inp.errors['cz']) 
        inp.param['distance_error'] = deepcopy(inp.errors['L'])  
        inp.param['i_error'] = deepcopy(inp.errors['i'])  
        inp.param['j_error'] = deepcopy(inp.errors['j'])  
                   
        
        
 
