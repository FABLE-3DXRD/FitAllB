#!/usr/bin/env python

# Modules to import 
from __future__ import print_function
from FitAllB import gofitglobalgrain

import logging
logging.basicConfig(level=logging.INFO,format='\n%(levelname)s: %(message)s')

if __name__=="__main__":
    options = None
    try:
        from optparse import OptionParser
        parser = OptionParser()
        options  = gofitglobalgrain.get_options(parser)
        print(options)
        gofitglobalgrain.run(options)
    except:
        if options != None:
            parser.print_help()
        raise 
