#!/usr/bin/env python

# Modules to import 
import sys
from xfab import tools
from FitAllB import check_input
import logging
from optparse import OptionParser

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
inp = check_input.parse_input(input_file=options.filename)  # Make instance of parse_input class
inp.read()              # read input file
inp.check()             # check validity of input
if inp.missing == True: # if problem exit
    logging.info('MISSING ITEMS')
    sys.exit()
inp.initialize()        # if ok initialize
inp.read_par()          # read detector.par file
inp.read_flt()          # read peaks_t##.flt file
inp.read_log()          # read grainspotter.log file
inp.read_res()          # read paramters file to resume refinement                              NB! optional
inp.read_rej()          # read file containing rejected peaks to resume refinement   NB! optional
inp.set_start()         # set values and errors for refinement start

if inp.fit['resume'] == None: # do outlier rejection
    inp.reject()            
else:                          # if refinement is resumed build residual and volume arrays 
    inp.residual = []
    inp.volume = []
    for i in range(inp.no_grains):
        inp.residual.append([])
        inp.volume.append([])
        for j in range(inp.nrefl[i]):
            inp.residual[i].append(1)
            inp.volume[i].append(1)
    from FitAllB import reject
    reject.intensity(inp)       # necessary to get correct volumes in output file, very few peaks actually rejected

inp.write_rej()         # write the files rejected during friedel and merge to rejection file


    
# refinement
while inp.fit['goon'] != 'end':
#    print '\nNumber of assigned reflections:'
#    print inp.nrefl
# calculate experimental errors using the present values 
    from FitAllB import error
    error.vars(inp)
# build functions to minimise
    from FitAllB import build_fcn
    build_fcn.FCN(inp)
# minuit fitting
    from FitAllB import fit
    lsqr = fit.fit_minuit(inp)
    lsqr.refine()


