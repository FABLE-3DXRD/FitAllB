#!/usr/bin/env python

# Modules to import 
import sys
from FitAllB import check_input
import logging
from optparse import OptionParser
from copy import deepcopy
import numpy as n
import os

logging.basicConfig(level=logging.INFO,format='\n%(levelname)s: %(message)s')
parser = OptionParser()
parser.add_option("-i", "--input", action="store",
                  dest="filename", type="string",
                  help="Name of the file containing the input parameters")
options , args = parser.parse_args()
if options.filename == None:
    parser.print_help()
    print "\nNo input file supplied [-i filename]\n"
    sys.exit()

    
#Read and check input
far = check_input.parse_input(input_file=options.filename)  # Make instance of parse_input class
far.read()              # read input file
far.check()             # check validity of input
if far.missing == True: # if problem exit
    logging.info('MISSING ITEMS')
    sys.exit()
far.initialize()        # if ok initialize
# mandatory farfield info 
far.read_par(far.files['par_file'])          # read detector.par file
far.read_flt(far.files['flt_file'])          # read peaks_t##.flt file
far.read_res()          # read paramters file to resume refinement                              NB! optional
if far.files['res_file'] == None:
    far.read_log()          # read grainspotter.log file
far.read_rej()          # read file containing rejected peaks to resume refinement   NB! optional
far.set_start()         # set values and errors for refinement start
check_input.set_globals(far)


# optional nearfield info
if far.files['near_flt_file'] != None:
    assert far.files['near_par_file'] != None, 'near_par_file parameter file for near-field detector is missing'
    near = deepcopy(far)
    # take special care of near-field keywords (eg copy near_dety_size to dety_size)
    for key in near.fit.keys():
        if key[0:5] == 'near_': 
            near.fit[key[5:len(key)]] = near.fit[key]
    near.fit['stem'] = far.fit['stem'] + '_near'
    print 'far.fit', far.fit
    print 'near.fit', near.fit
    near.read_par(near.files['near_par_file'])
    near.read_flt(near.files['near_flt_file'])
    keys = ['cell__a','cell__b','cell__c','cell_alpha','cell_beta','cell_gamma','cell_lattice_[P,A,B,C,I,F,R]','chi','omegasign','wavelength']
    for key in keys:
        assert near.param[key] == far.param[key], '%s is different for near- and far-field detectors' %key
    # in case of different wedge use farfield value and refine
    if near.param['wedge'] != far.param['wedge']:
        near.param['wedge'] = far.param['wedge']
#        near.fit['w'] = 1
#        far.fit['w'] = 1
    near.values = far.values
    near.errors = far.errors
    near.fit['skip'] = far.fit['skip']
    check_input.set_globals(near)
    # match peaks on nearfield detector and reject outliers
    from FitAllB import near_field
    near_field.find_refl(near)
    near_field.match(near)
    from FitAllB import error
    error.vars(near)
    from FitAllB import build_fcn
    build_fcn.FCN(near)
    near.reject()
    near.write_rej()
    near.fit['reforder'] = ['starta','xyza','end']
#    near.fit['reforder'] = ['starta','rotposa','end']
    if near.fit['near_resume'] != None:
        near.fit['goon'] = near.fit['near_resume']
    else:
        near.fit['goon'] = near.fit['reforder'][0]
    # nearfield refinement
    from FitAllB import fit
    fit.refine(near)


#  Farfield outlier rejection
if far.fit['resume'] == None: # do outlier rejection
    if far.files['res_file'] != None:
        from FitAllB import near_field
        near_field.find_refl(far)
        near_field.match(far)
#    from FitAllB import near_field
#    near_field.find_refl(far)
#    near_field.match(far)
    from FitAllB import error
    error.vars(far)
    from FitAllB import build_fcn
    build_fcn.FCN(far)
    far.reject()            
else:                          # if refinement is resumed build residual and volume arrays 
    far.residual = []
    far.volume = []
    far.mean_ia = []
    for i in range(far.no_grains):
        far.residual.append([])
        far.volume.append([])
        far.mean_ia.append([])
        for j in range(far.nrefl[i]):
            far.residual[i].append(1)
            far.volume[i].append(1)
            far.mean_ia[i].append(1)
    from FitAllB import reject
    reject.intensity(far)       # necessary to get correct volumes in output file, very few peaks actually rejected
far.write_rej()         # write the files rejected during friedel and merge to rejection file
 
# farfield refinement
#far.fit['reforder'] = ['start','eps','grain','start1','final','end']
#far.fit['reforder'] = ['start','grain','start1','final','end']
far.fit['reforder'] = ['start','grain','final','end'] # scaling of errors does not work if start1 is present and no global parameters are refined
if far.fit['resume'] != None:
    far.fit['goon'] = far.fit['near_resume']
else:
    far.fit['goon'] = far.fit['reforder'][0]
from FitAllB import fit
fit.refine(far)


    
    
    
os.remove('%s/fcn.py' %far.fit['direc'])
os.remove('%s/fcn.pyc' %far.fit['direc'])
sys.exit()
# program ends here after deleting fcn.py and fcn.pyc




#nearfield refinement
try: 
#    near.fit['reforder'] = ['startb','rotposb','end']
    near.fit['reforder'] = ['startb','xyzb','end']
    if near.fit['near_resume'] != None:
        near.fit['goon'] = near.fit['near_resume']
    else:
        near.fit['goon'] = near.fit['reforder'][0]
    # set errors back to standard values to ensure refinement of all grains
    for i in range(near.no_grains):
        near.errors['x%s' %i] = near.param['y_size']/5.
        near.errors['y%s' %i] = near.param['y_size']/5.
        near.errors['z%s' %i] = near.param['z_size']/10.
    # nearfield refinement
    from FitAllB import fit
    fit.refine(near)

# farfield refinement
    far.fit['reforder'] = ['grainb','start1','finalb','end']
    far.fit['goon'] = 'grainb'
    from FitAllB import fit
    fit.refine(far)
    
except:
    pass

    

