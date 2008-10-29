import numpy as n
from xfab import tools
from polyxsim import reflections
from copy import deepcopy
import time
import minuit
import sys
import logging
logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')


def friedel(inp):
        """
        This function eliminates the peaks that have an intensity that differs significantly from that of the Friedel mate
    
        Jette Oddershede, July 24 2008
        """
    
        # remove peaks assigned more than once to the same grain
        double = 0
        for i in range(inp.no_grains):
            for j in range(inp.nrefl[i]-1,-1,-1): # loop backwards to make pop work
                if inp.id[i].count(inp.id[i][j]) > 1:
                    double = double + 1
                    reject(inp,i,j,'unique')
        print 'Removed', double, 'reflections which were assigned more than once to the same grain'
        
        # remove single Friedel mates with poor match
        delete = 0
        for i in range(inp.no_grains):
            if i+1 in inp.fit['skip']:
                pass
            else:   
                for peak1 in range(len(inp.id[i])-1,-1,-1):
                    h1 = inp.h[i][peak1]
                    k1 = inp.k[i][peak1]
                    l1 = inp.l[i][peak1]
                    ds1 = tools.sintl(inp.unit_cell,[h1,k1,l1])
                    iavg = 0
                    nn = 0
                    roi = []
                    for peak2 in range(len(inp.id[i])):
                        h2 = inp.h[i][peak2]
                        k2 = inp.k[i][peak2]
                        l2 = inp.l[i][peak2]
                        ds2 = tools.sintl(inp.unit_cell,[h2,k2,l2])
                        if ds2 < ds1:
                            pass
                        elif ds2 > ds1:
                            break
                        else:
                            if peak2 != peak1:
                                roi.append(peak2)
                            iavg = iavg + inp.F2vol[inp.id[i][peak2]]
                            nn = nn + 1
                    if nn !=0:
                        iavg = iavg/nn
                    for peak2 in roi:
                        h2 = inp.h[i][peak2]
                        k2 = inp.k[i][peak2]
                        l2 = inp.l[i][peak2]
                        ds2 = tools.sintl(inp.unit_cell,[h2,k2,l2])
                        if -h1 == h2 and -k1 == k2 and -l1 == l2:
                            i1 = inp.F2vol[inp.id[i][peak1]]
                            i2 = inp.F2vol[inp.id[i][peak2]]
                            if i2 > 0 and nn > 1:
                                if (i1/i2 < 0.5 or i1/i2 > 2):
                                    i_p1 = int((iavg*nn-i1)/(nn-1))
                                    # test if removing peak1 helps:
                                    if i2/i_p1 > 0.5 and i2/i_p1 < 2:
                                        delete = delete + 1
                                        reject(inp,i,j,'friedel')
                                        break

        print 'Removed', delete, 'reflections which did not match their Friedel mates in intensity'

        # remove Friedel pairs with poor match
        delete = 0
        for i in range(inp.no_grains):
            if i+1 in inp.fit['skip']:
                pass
            else:   
                for peak1 in range(len(inp.id[i])-1,-1,-1):
                    if peak1 < len(inp.id[i]):
                        h1 = inp.h[i][peak1]
                        k1 = inp.k[i][peak1]
                        l1 = inp.l[i][peak1]
                        ds1 = tools.sintl(inp.unit_cell,[h1,k1,l1])
                        iavg = 0
                        nn = 0
                        roi = []
                        for peak2 in range(len(inp.id[i])-1,-1,-1):
                            h2 = inp.h[i][peak2]
                            k2 = inp.k[i][peak2]
                            l2 = inp.l[i][peak2]
                            ds2 = tools.sintl(inp.unit_cell,[h2,k2,l2])
                            if ds2 > ds1:
                                pass
                            elif ds2 < ds1:
                                break
                            else:
                                if peak2 != peak1:
                                    roi.append(peak2)
                                iavg = iavg + inp.F2vol[inp.id[i][peak2]]
                                nn = nn + 1
                        if nn !=0:
                            iavg = iavg/nn
                        for peak2 in roi:
                            h2 = inp.h[i][peak2]
                            k2 = inp.k[i][peak2]
                            l2 = inp.l[i][peak2]
                            ds2 = tools.sintl(inp.unit_cell,[h2,k2,l2])
                            if -h1 == h2 and -k1 == k2 and -l1 == l2:
                                i1 = inp.F2vol[inp.id[i][peak1]]
                                i2 = inp.F2vol[inp.id[i][peak2]]
                                if (i1/i2 < 0.5 or i1/i2 > 2):
#                                    delete = delete + 1
#                                    reject(inp,i,peak1,'friedel')
#                                    reject(inp,i,peak2,'friedel')
                                    break

        print 'Removed', delete, 'entire Friedel pairs with non-matching intensities', delete*2, 'reflections'
        insignificant(inp)

    
def intensity(inp):
        """
        Reject peaks based on intensity
        
        Jette Oddershede, August 27 2008
        """
        
        if  inp.files['structure_file'] != None:
            inp.param['structure_phase_0'] = inp.files['structure_file']
            xtal_structure = reflections.open_structure(inp.param,0)
            hkl = reflections.gen_miller(inp.param,0)
            hkl = reflections.calc_intensity(hkl,xtal_structure)
#            print hkl
            
            for i in range(inp.no_grains):
                for j in range(inp.nrefl[i]):
                    h = inp.h[i][j]
                    k = inp.k[i][j]
                    l = inp.l[i][j]
                    value = False
                    for m in range(len(hkl)):
                        if hkl[m][0] == h and hkl[m][1] == k and hkl[m][2] == l:
                            inp.volume[i][j] = inp.F2vol[inp.id[i][j]]/hkl[m][3]
                            value = True
                            break
                        else:
                            pass
                    if value == False:
#                        print i+1,j+1,h,k,l,'should never go here, ask osho for help...'
                        inp.volume[i][j] = -1
                        
            data = deepcopy(inp.volume)
            minvol = []
            maxvol = []
            for i in range(inp.no_grains):
                if i+1 in inp.fit['skip']:
                    pass
                else:
                    rej = []
                    newreject = 1
                    while newreject > 0:
                        tmp = len(rej)
                        mad(data[i],rej,5)
                        newreject = len(rej) - tmp
                    avgdata = n.sum(data[i])/len(data[i])
                    sigdata = spread(data[i])
                    if len(rej) > 1:
                        avgrej = n.sum(rej)/len(rej)
                        sigrej = spread(rej)
                    elif len(rej) == 1:
                        avgrej = rej[0]
                        sigrej = 0                
                    else:
                        avgrej = 0
                        sigrej = 0
                if len(data[i]) > 0:
                    minvol.append(max(0,min(data[i])))
                    maxvol.append(max(data[i]))
                else: 
                    minvol.append(0)
                    maxvol.append(-1)
                    if i+1 not in inp.fit['skip']:
                        inp.fit['skip'].append(i+1)
                    
#                   print '\n',i, avgdata, sigdata, len(data[i]), minvol[i], maxvol[i],'\n   ',avgrej, sigrej, len(reject)#,'\n',data[i],'\n',reject

            delete = 0
            for i in range(inp.no_grains):
                if i+1 in inp.fit['skip']:
                    pass
                else:
                    for j in range(inp.nrefl[i]-1,-1,-1):
                        if inp.volume[i][j] < minvol[i] or inp.volume[i][j] > maxvol[i]:
                            reject(inp,i,j,'intensity')
                            delete = delete + 1
                    
            print 'Rejected', delete, 'peaks because of different intensity scales'
            insignificant(inp)
            
                    			
def residual(inp,limit,only=None):
        """
        Reject outliers peaks with a distance to the calculated peak position of
        more than limit times the mean distance for the given grain	
		
		Jette Oddershede, Risoe DTU, May 15 2008
        """
		
        # must update inp.vars because the order here is [i][j] in stead of [id[i][j]], the latter doesn't change when peaks are rejected, the former does.
        try:
            self.values = deepcopy(self.m.values)
        except:
            pass
        # calculate experimental errors using the present values 
        import error
        error.vars(inp)
        # build functions to minimise
        import build_fcn
        build_fcn.FCN(inp)
        #refinement update
        import fcn
        reload(fcn)
        for i in range(inp.no_grains):
            if i+1 in inp.fit['skip']:
                pass
            else:				
                for j in range(inp.nrefl[i]): 
                    inp.residual[i][j] = int(fcn.peak(inp.h[i][j],inp.k[i][j],inp.l[i][j],
                                              inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],inp.vars[i][j], 
                                              inp.values['wx'],inp.values['wy'],
                                              inp.values['tx'],inp.values['ty'],inp.values['tz'],
                                              inp.values['py'],inp.values['pz'],
                                              inp.values['cy'],inp.values['cz'],
                                              inp.values['L'],
                                              inp.values['x%s' %i],inp.values['y%s' %i],inp.values['z%s' %i], 
                                              inp.values['phia%s' %i],inp.values['PHI%s' %i],inp.values['phib%s' %i], 
                                              inp.values['epsaa%s' %i],inp.values['epsab%s' %i],inp.values['epsac%s' %i], 
                                              inp.values['epsbb%s' %i],inp.values['epsbc%s' %i],inp.values['epscc%s' %i])) 
                                                    
        data = deepcopy(inp.residual)
        maxres = [0]*inp.no_grains
        for i in range(inp.no_grains):
            data[i].sort()
            if i+1 in inp.fit['skip']:
                pass
            else:		
                mean = int(n.sum(data[i])/len(data[i]))
                medi = median(data[i])
#                print i, len(data[i]), medi, mean,'\n',data[i]
                while mean > limit*medi:
                    data[i].pop()
                    mean = int(n.sum(data[i])/inp.nrefl[i])
                    medi = median(data[i])
                maxres[i] = max(data[i])
#                print i, len(data[i]),medi,mean,'\n',data[i],'\n'
        
        delete = 0
        if only==None:
            only = range(1,1+inp.no_grains)        
        for i in range(inp.no_grains):
            if i+1 in inp.fit['skip'] or i+1 not in only:
                pass
            else:				
                for j in range(inp.nrefl[i]-1,-1,-1): # loop backwards to make pop work
                    if inp.residual[i][j] > maxres[i]:
                        delete = delete + 1
                        reject(inp,i,j,'residual')
        print 'Rejected', delete, 'reflection based on residuals'
        insignificant(inp)
                    

def multi(inp):
        """
        Peaks assigned to more than one grain are rejected for a grain if both the residual and the volume are off
        NB! Define off
        
        Jette Oddershede, 1 September 2008
        """

        # handling reflection assigned to more than one grain
        grain = []
        peak = []
        for k in range(inp.param['total_refl']):
            grain.append([])
            peak.append([])
            
        for i in range(inp.no_grains):
            if i+1 not in inp.fit['skip']:
                for j in range(inp.nrefl[i]):
                    grain[inp.id[i][j]].append(i)
                    peak[inp.id[i][j]].append(j)
                        
        multi = 0
        bad = []
        if  inp.files['structure_file'] != None:
            volavg = []
            for i in range(inp.no_grains):
                if len(inp.volume[i]) > 0:
                    volavg.append(sum(inp.volume[i])/len(inp.volume[i]))
                else:
                    volavg.append(0)
            
        for k in range(inp.param['total_refl']):
            if len(grain[k]) > 1:
                multi = multi + 1
                for m in range(len(grain[k])):
#                    print grain[k][m],peak[k][m],k, \
#                          inp.residual[grain[k][m]][peak[k][m]], sum(inp.residual[grain[k][m]])/len(inp.residual[grain[k][m]]),\
#                          inp.volume[grain[k][m]][peak[k][m]], sum(inp.volume[grain[k][m]])/len(inp.volume[grain[k][m]])
                    for o in range(len(grain[k])):
                        #reject if largest residual...
                        if inp.residual[grain[k][m]][peak[k][m]] > inp.residual[grain[k][o]][peak[k][o]]:
                            #... and largest distance from mean volume if structure file given
                            if  inp.files['structure_file'] != None:
                                if abs(inp.volume[grain[k][m]][peak[k][m]]-volavg[grain[k][m]]) > abs(inp.volume[grain[k][o]][peak[k][o]]-volavg[grain[k][o]]):
                                    bad.append([grain[k][m],peak[k][m]])
                            else:
                                bad.append([grain[k][m],peak[k][m]])
        unique_list(bad)                      
        #print bad
        print 'Number of reflections assigned to more than one grain',multi       

                                        

        # for peaks assigned to more than one grain remove all but the best assignment                
        for i in range(len(bad)-1,-1,-1):
            reject(inp,bad[i][0],bad[i][1],'multi')
        print 'Delete',len(bad), 'reflection because they are assigned to more than one grain' 
        insignificant(inp)
                                    
        
def merge(inp):
        """
        This function merges grain if the fraction of similar peaks exceeds the overlap parameter
    
        Jette Oddershede, August 20 2008
        """
        
#        volsig = []
#        for i in range(inp.no_grains):
#            if i+1 not in inp.fit['skip']:
#                volsig.append(spread(inp.volume[i]))
#            else:
#                volsig.append([])
#        resavg = []
#        for i in range(inp.no_grains):
#            if i+1 not in inp.fit['skip']:
#                resavg.append(sum(inp.residual[i])/len(inp.residual[i]))
#            else:
#                resavg.append([])
        # merge overlapping grains
        for gr1 in range(inp.no_grains-1):
            if gr1+1 not in inp.fit['skip']:
                for gr2 in range(gr1+1,inp.no_grains):
                    if gr2+1 not in inp.fit['skip']:
                        multilimit = inp.fit['overlap'] * min(inp.nrefl[gr1],inp.nrefl[gr2])
                        multi = 0
                        for peak1 in range(inp.nrefl[gr1]):
                            for peak2 in range(inp.nrefl[gr2]): 
                                if inp.id[gr1][peak1] == inp.id[gr2][peak2]:
                                    multi = multi + 1
                        if multi > multilimit:
                            print 'Equal grains:', gr1+1, 'and', gr2+1, 'with number of equal refl:', multi, '(total', inp.nrefl[gr1], 'and', inp.nrefl[gr2],')'
                            if inp.nrefl[gr1] < inp.nrefl[gr2]:# and resavg[gr1] > resavg[gr2] and volsig[gr1] > volsig[gr2]: #remove residual and volume criteria
                                inp.fit['skip'].append(gr1+1)
                                print 'Skip grain', gr1+1
                            elif inp.nrefl[gr1] > inp.nrefl[gr2]:# and resavg[gr1] < resavg[gr2] and volsig[gr1] < volsig[gr2]: #remove residual and volume criteria
                                inp.fit['skip'].append(gr2+1)
                                print 'Skip grain', gr2+1
                            else:                                
                                print 'Merge grains', gr1+1, 'and', gr2+1
                                inp.fit['skip'].append(gr2+1)
                                set1 = set(inp.id[gr1][:])
                                set2 = set(inp.id[gr2][:])
                                for peak in range(inp.nrefl[gr1]-1,-1,-1): # loop backwards to make pop work
                                    if inp.id[gr1][peak] not in inp.id[gr2]:
                                        reject(inp,gr1,peak,'merge')
        unique_list(inp.fit['skip'])
        insignificant(inp)

            
def insignificant(inp):
        """
        Remove grains with less than 12 peaks as being insignificant
        """
        for i in range(inp.no_grains):
            if inp.nrefl[i] < 12 and i+1 not in inp.fit['skip']:
                inp.fit['skip'].append(i+1)
        inp.fit['skip'].sort()
        
 
def poor(inp):
    """
    Remove poor grains
    """
    pass


 
# Helpful functions               
               
def reject(inp,i,j,message):
        """
        Reject peak j from grain i and print the message as rejectvalue
        """
        inp.fit['rejectgrain'].append(i)
        inp.fit['rejectid'].append(inp.id[i].pop(j))
        inp.fit['hh'].append(inp.h[i].pop(j))
        inp.fit['kk'].append(inp.k[i].pop(j))
        inp.fit['ll'].append(inp.l[i].pop(j))
        inp.fit['rejectvalue'].append(message)
        inp.nrefl[i] = inp.nrefl[i] - 1
        inp.fit['outliers'] = inp.fit['outliers'] + 1
        inp.newreject = inp.newreject + 1
        inp.fit['newreject_grain'].append(i+1)
        inp.residual[i].pop(j)
        inp.volume[i].pop(j)
               
               
def median(numbers):
   # Sort the list and take the middle element.
   nn = len(numbers)
   copy = numbers[:] # So that "numbers" keeps its original order
   copy.sort()
   if nn & 1:         # There is an odd number of elements
      return copy[nn // 2]
   else:
      return (copy[nn // 2 - 1] + copy[nn // 2]) / 2

def spread(data):
        data = n.array(data)
        vars = n.sum(data*data) - n.sum(data)**2/len(data)
        return n.sqrt(vars/(len(data)-1))
        
def unique_list(list):
        list.sort()
        for i in range(len(list)-1,0,-1):
            if list[i] == list[i-1]:
                list.pop(i)
        return list

def mad(data,reject,limit):
        """
        Perform outlier rejection based on median absolute deviation
        Move all data points more than 5 median absolute deviations to the reject list
        
        Jette Oddershede 28 August 2008
        """

        if len(data) > 1:
            data.sort()
            medi = median(data)
            maddata = []
            for j in range(len(data)):
                maddata.append(abs(data[j]-medi))
            maddata.sort()
            mad = limit*maddata[len(maddata)/2]
            for j in range(len(data)-1,-1,-1):
                if data[j] < medi-mad or data[j] > medi+mad:
                    reject.append(data.pop(j))
            reject.sort()
    
