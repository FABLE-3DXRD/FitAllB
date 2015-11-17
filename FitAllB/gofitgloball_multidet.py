import sys
from FitAllB import check_input_multidet 
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

    return options

def run(options):
    if options.filename == None:
        #parser.print_help()
        print "\nNo input file supplied [-i filename]\n"
        sys.exit()
    
        
    #Read and check input
    far = check_input_multidet.parse_input(input_file=options.filename)  # Make instance of parse_input class
    far.read()              # read input file
    far.check()             # check validity of input
    if far.missing == True: # if problem exit
        logging.info('MISSING ITEMS')
        sys.exit()
    far.initialize()        # if ok initialize
    far.read_res()                      # read paramters file to resume refinement                            
    # mandatory farfield info 
    far.read_par(far.files['par_file']) # read detector.par file
    far.read_flt(far.files['flt_file']) # read peaks_t##.flt file
    keys = ['dety_size','detz_size']
    for key in keys:
        far.fit[key+'0'] = far.fit[key]
    keys = ['o11','o12','o21','o22']
    for key in keys:
        far.param[key+'0'] = far.param[key]
    # second detector
    if far.fit['no_det'] > 1:
        far1 = deepcopy(far)
        del(far1.labels)
        del(far1.tth)
        del(far1.eta)
        del(far1.nrefl)
        del(far1.id)
        del(far1.h)
        del(far1.k)
        del(far1.l)
        keys = ['dety_size','detz_size']
        for key in keys:
            far1.fit[key] = far.fit[key+'1']
        far1.read_par(far.files['par_file_1']) # read detector.par file
        far1.read_flt(far.files['flt_file_1']) # read peaks_t##.flt file
        #put info from far1 in far
        keys = ['cell__a','cell__b','cell__c','cell_alpha','cell_beta','cell_gamma','cell_lattice_[P,A,B,C,I,F,R]','chi','omegasign','wavelength','wedge']
        for key in keys:
            assert far1.param[key] == far.param[key], '%s is different for par_file and and par_file_1' %key
        keys = ['distance','o11','o12','o21','o22','tilt_x','tilt_y','tilt_z','y_center','z_center','y_size','z_size']
        for key in keys:
            far.param[key+'1'] = far1.param[key]
        print 'parameters of second detector read successfully'
        
    # third detector, NB more can be added in this way!!
    if far.fit['no_det'] > 2:
        far2 = deepcopy(far)
        del(far2.labels)
        del(far2.tth)
        del(far2.eta)
        del(far2.nrefl)
        del(far2.id)
        del(far2.h)
        del(far2.k)
        del(far2.l)
        keys = ['dety_size','detz_size']
        for key in keys:
            far2.fit[key] = far.fit[key+'2']
        far2.read_par(far.files['par_file_2']) # read detector.par file
        far2.read_flt(far.files['flt_file_2']) # read peaks_t##.flt file
        #put info from far1 in far
        keys = ['cell__a','cell__b','cell__c','cell_alpha','cell_beta','cell_gamma','cell_lattice_[P,A,B,C,I,F,R]','chi','omegasign','wavelength','wedge']
        for key in keys:
            assert far2.param[key] == far.param[key], '%s is different for par_file and and par_file_1' %key
        keys = ['distance','o11','o12','o21','o22','tilt_x','tilt_y','tilt_z','y_center','z_center','y_size','z_size']
        for key in keys:
            far.param[key+'2'] = far2.param[key]
        print 'parameters of third detector read successfully'
    
    # set values and errors for refinement start
    far.set_start()                     
    check_input_multidet.set_globals(far)
    if far.fit['no_det'] > 1:
        far1.set_start()   
        check_input_multidet.set_globals(far1)
    if far.fit['no_det'] > 2:
        far2.set_start()   
        check_input_multidet.set_globals(far2)

 
#    for key in far.param.keys():
#        print key,far.param[key]
    
    # forward projection to assign reflections
    if far.files['res_file'] != None and far.labels == None:
        from FitAllB import near_field
        near_field.find_refl(far)
        print '\nFirst detector'
        near_field.match(far)
    from FitAllB import error
    error.vars(far)
    from FitAllB import build_fcn
    build_fcn.FCN(far)
    far.reject()            
    far.nrefl = [far.nrefl]
    far.id = [far.id]
    far.h = [far.h]
    far.k = [far.k]
    far.l = [far.l]
    far.w = [far.w]
    far.dety = [far.dety]
    far.detz = [far.detz]
#    far.vars = [far.vars]
    far.param['total_refl'] = [far.param['total_refl']]

    if far.fit['no_det'] > 1:
        if far1.files['res_file'] != None and far1.labels == None:
            from FitAllB import near_field
            near_field.find_refl(far1)
            print '\nSecond detector'
            near_field.match(far1)
        error.vars(far1)
        from FitAllB import build_fcn
        build_fcn.FCN(far1)
        far1.reject()            
        far.nrefl.append(far1.nrefl)
        far.id.append(far1.id)
        far.h.append(far1.h)
        far.k.append(far1.k)
        far.l.append(far1.l)
        far.w.append(far1.w)
        far.dety.append(far1.dety)        
        far.detz.append(far1.detz) 
#        far.vars.append(far1.vars) 
        far.param['total_refl'].append(far1.param['total_refl'])
    if far.fit['no_det'] > 2:
        if far2.files['res_file'] != None and far2.labels == None:
            from FitAllB import near_field
            near_field.find_refl(far2)
            print '\nThird detector'
            near_field.match(far2)
        error.vars(far2)
        from FitAllB import build_fcn
        build_fcn.FCN(far2)
        far2.reject()            
        far.nrefl.append(far2.nrefl)
        far.id.append(far2.id)
        far.h.append(far2.h)
        far.k.append(far2.k)
        far.l.append(far2.l)
        far.w.append(far2.w)
        far.dety.append(far2.dety)        
        far.detz.append(far2.detz)        
#        far.vars.append(far2.vars) 
        far.param['total_refl'].append(far2.param['total_refl'])
         
#    print far.nrefl

    import error_multidet
    error_multidet.vars(far)
#    print far.vars
   
    
    # farfield refinement
    for k in range(far.fit['cycle']):
    #while len(far.fit['newreject_grain']) > 0:
        # refine grain paramters
        far.fit['reforder'] = ['start%s' %k,'rotpos%s' %k,'end'] 
 #       if k==0:
 #           far.fit['reforder'] = ['start%s' %k,'end'] 
        far.fit['goon'] = far.fit['reforder'][0]
        far.fit['newreject_grain'] = range(far.no_grains+1)
        # refine grains
        from FitAllB import fit_multidet
        far.residual = []
        far.volume = []
        far.mean_ia = []
        for i in range(far.no_grains):
            far.residual.append([])
            far.volume.append([])
            far.mean_ia.append([])
            for m in range(far.fit['no_det']):
                for j in range(far.nrefl[m][i]):
                    far.residual[i].append(1)
                    far.volume[i].append(1)
                    far.mean_ia[i].append(1)
        from FitAllB import build_fcn_multidet
        build_fcn_multidet.FCN(far)
        import fcn
        reload(fcn)
        fit_multidet.refine(far)
        # refine globals
        far.fit['reforder'] = ['globals%s' %k,'end'] 
        far.fit['goon'] = far.fit['reforder'][0]
        from FitAllB import fitga_multidet
        far.residual = []
        far.volume = []
        far.mean_ia = []
        for i in range(far.no_grains):
            far.residual.append([])
            far.volume.append([])
            far.mean_ia.append([])
            for m in range(far.fit['no_det']):
                for j in range(far.nrefl[m][i]):
                    far.residual[i].append(1)
                    far.volume[i].append(1)
                    far.mean_ia[i].append(1)
# ---- OK until here
        fitga_multidet.refine(far)
        from FitAllB import reject_multidet
        reject_multidet.residual(far,far.fit['rej_resmedian'])       
        reject_multidet.mean_ia(far,far.fit['rej_ia']*float(k+2)/float(k+1))       
#        reject.intensity(far)       
        from FitAllB import write_output_multidet
        write_output_multidet.write_rej(far,message='globals%i' %k)                 
        
        no_ref = []
        for i in range(far.no_grains):
            no_ref.append(0)
            for k in range(far.fit['no_det']):
                no_ref[i] = no_ref[i] + far.nrefl[k][i]

        for i in range(far.no_grains):
            if no_ref[i] < far.fit['min_refl'] and i+1 not in far.fit['skip']:
                far.fit['skip'].append(i+1)
        far.fit['skip'].sort()
     
    
# --- OK after this point    
    
        
    # program ends here after deleting fcn.py and fcn.pyc
    print '\nNormal termination of FitGlobAll_MultiDet'
    os.remove('%s/fcn.py' %far.fit['direc'])
    os.remove('%s/fcn.pyc' %far.fit['direc'])
