import sys
from FitAllB import check_input
import logging
from optparse import OptionParser
from copy import deepcopy
import numpy as n
import os

def get_options(parser):
    parser.add_option("-i", "--input", action="store",
                      dest="filename", type="string",
                      help="Name of the file containing the input parameters")
    options, args = parser.parse_args()
    options.killfile = None

    return options

def run(options):
    if options.filename == None:
        raise ValueError, "\nNo input file supplied [-i filename]\n"

    try:
        if options.killfile is not None and os.path.exists(options.killfile):
            print "The purpose of the killfile option is to create that file"
            print "only when you want fitallb to stop"
            print "If the file already exists when you run fitallb it is"
            print "never going to get started"
            raise ValueError("Your killfile "+options.killfile+" already exists")
    except:
        pass

    #Read and check input
    far = check_input.parse_input(input_file=options.filename)  # Make instance of parse_input class
    far.read()              # read input file
    far.check()             # check validity of input
    if far.missing == True: # if problem exit
        logging.info('MISSING ITEMS')
        sys.exit()
    far.initialize()        # if ok initialize
    # mandatory farfield info 
    far.read_par(far.files['par_file']) # read detector.par file
    far.read_flt(far.files['flt_file']) # read peaks_t##.flt file
    far.read_res()                      # read paramters file to resume refinement                              NB! optional
    if far.files['res_file'] == None:
        far.read_log()                  # read grainspotter.log file
        far.read_rej(far.files['rej_file'])                      # read file containing rejected peaks to resume refinement   NB! optional
    far.set_start()                     # set values and errors for refinement start
    check_input.set_globals(far)
    
    # optional nearfield info
    if far.files['near_flt_file'] != None:
        check_input.interrupt(options.killfile)
        assert far.files['near_par_file'] != None, 'near_par_file parameter file for near-field detector is missing'
        near = deepcopy(far)
        # take special care of near-field keywords (eg copy near_dety_size to dety_size)
        for key in near.fit.keys():
            if key[0:5] == 'near_': 
                near.fit[key[5:len(key)]] = near.fit[key]
        near.fit['stem'] = far.fit['stem'] + '_near'
#        print 'far.files', far.files
#        print 'near.files', near.files
        near.read_par(near.files['near_par_file'])
        near.read_flt(near.files['near_flt_file'])
        keys = ['cell__a', 'cell__b', 'cell__c', 'cell_alpha', 'cell_beta', 'cell_gamma', 'cell_lattice_[P,A,B,C,I,F,R]', 'chi', 'omegasign', 'wavelength']
        for key in keys:
            assert near.param[key] == far.param[key], '%s is different for near- and far-field detectors' % key
        # in case of different wedge use farfield value and refine
        if near.param['wedge'] != far.param['wedge']:
            near.param['wedge'] = far.param['wedge']
        near.values = far.values
        near.errors = far.errors
        near.fit['skip'] = far.fit['skip']
        check_input.set_globals(near)
        # match peaks on nearfield detector and reject outliers
        from FitAllB import near_field
        check_input.interrupt(options.killfile)
        near_field.find_refl(near)
        check_input.interrupt(options.killfile)
        near_field.match(near)
        near.read_rej(near.files['near_rej_file'])                      # read file containing rejected peaks to resume refinement   NB! optional
        from FitAllB import error
        error.vars(near)
        from FitAllB import build_fcn
        build_fcn.FCN(near)
        check_input.interrupt(options.killfile)
        near.reject()
        near.write_rej()
        near.fit['reforder'] = ['start', 'rotpos', 'end']
        near.fit['goon'] = near.fit['reforder'][0]
        # nearfield refinement
        from FitAllB import fit
        fit.refine(near,options.killfile)
    
    
    #  Farfield outlier rejection
    #check_input.set_globals(far)
    if far.files['res_file'] != None:
        check_input.interrupt(options.killfile)
        from FitAllB import near_field
        near_field.find_refl(far)
        check_input.interrupt(options.killfile)
        near_field.match(far)
        far.read_rej(far.files['rej_file'])                      # read file containing rejected peaks to resume refinement   NB! optional
    from FitAllB import error
    error.vars(far)
    from FitAllB import build_fcn
    build_fcn.FCN(far)
    check_input.interrupt(options.killfile)
    far.reject()            
    far.write_rej()                 
    
    # farfield refinement
    far.fit['reforder'] = ['grain', 'final', 'end'] 
    far.fit['goon'] = far.fit['reforder'][0]
    from FitAllB import fit
    fit.refine(far,options.killfile)
    
        
    # program ends here after deleting fcn.py and fcn.pyc
    print '\nNormal termination of FitAllB'
    os.remove('%s/fcn.py' % far.fit['direc'])
    os.remove('%s/fcn.pyc' % far.fit['direc'])
