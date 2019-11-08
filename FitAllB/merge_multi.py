#!/usr/bin/env python



import os
from copy import deepcopy
from . import merge
import time



layer = list(range(28-1,-1,-1))
ranges = ['c','d','a','b']
thresholds = [200,500,1000,2000,5000,10000,20000]


output = deepcopy(ranges)

for l in layer:
        t1 = time.clock()
        for r in range(len(ranges)):
            stem = 'grwi_mapping_%0.3d_f%s'  %(l,ranges[r])
            output[r] = 'grwi_mapping_%0.3d_f%s.flt'  %(l,ranges[r])
            logfile = 'grwi_mapping_%0.3d_f%s.log'  %(l,ranges[r])
            print(output[r], thresholds)
            merge.merge(stem,thresholds,output[r],logfile)
	    if ranges[r] == 'c' or ranges[r] == 'd':
		    command = 'python omega_rotate.py %s %s' %(output[r],output[r])
		    os.system(command)
        merged = 'grwi_mapping_%0.3d_frelon.flt'  %(l)
        merge.two2one(output[0],output[1],merged)
        for r in range(2,len(ranges)):
            merge.two2one(merged,output[r],merged)
        merge.spot3d_id(merged)
        print('Created merged file', merged, 'spending ', int(time.clock()-t1), 's\n')
            
            
            
            
