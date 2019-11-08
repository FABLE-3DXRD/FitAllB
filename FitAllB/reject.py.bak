import numpy as n
from xfab import tools
from polyxsim import reflections
from copy import deepcopy
import time
try:
    from iminuit import Minuit
except ImportError:
    from minuit import Minuit
import sys
import logging
logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')


def edge(inp):

    delete = 0
    for i in range(inp.no_grains):
        for j in range(inp.nrefl[i]-1,-1,-1):
            if inp.Sww[inp.id[i][j]] < -50:
                reject(inp,i,j,'omega edge')
                delete = delete + 1
            elif inp.Syy[inp.id[i][j]] < -50:
                reject(inp,i,j,'y edge')
                delete = delete + 1
            elif inp.Szz[inp.id[i][j]] < -50:
                reject(inp,i,j,'z edge')
                delete = delete + 1
                    
    print 'Rejected', delete, 'peaks close to the edges of the detector or the omega ranges'
    insignificant(inp)
                
    
        
def intensity(inp):
        """
        Reject peaks based on intensity
        
        Jette Oddershede, August 27 2008
        Absorption correction added December 2011
        """
        
        if  inp.files['structure_file'] != None:
            inp.param['structure_phase_0'] = inp.files['structure_file']
            xtal_structure = reflections.open_structure(inp.param,0)
            hkl = reflections.gen_miller(inp.param,0)
            hkl = reflections.calc_intensity(hkl,xtal_structure)
#            print hkl
            
            if inp.fit['abs_mu'] > 0:
                xmin = min(inp.fit['abs_xlim'])*1e3
                xmax = max(inp.fit['abs_xlim'])*1e3
                ymin = min(inp.fit['abs_ylim'])*1e3
                ymax = max(inp.fit['abs_ylim'])*1e3
            
            for i in range(inp.no_grains):
                x = deepcopy(inp.values['x%s' %i])
                y = deepcopy(inp.values['y%s' %i])
                z = deepcopy(inp.values['z%s' %i])
                #Make sure that grain not outside of bounds for doing the correction
                if inp.fit['abs_mu'] > 0:
                    if x < xmin:
                        x = xmin
                    elif x > xmax:
                        x = xmax
                    if y < ymin:
                        y = ymin                
                    elif y > ymax:
                        y = ymax                
                for j in range(inp.nrefl[i]):
                    h = inp.h[i][j]
                    k = inp.k[i][j]
                    l = inp.l[i][j]
                    # Absorption correction start
                    if inp.fit['abs_mu'] > 0:
                        w = inp.w[inp.id[i][j]]
                        dety = inp.dety[inp.id[i][j]]
                        detz = inp.detz[inp.id[i][j]]
                        Omega = tools.form_omega_mat_general(w*n.pi/180.,inp.values['wx']*n.pi/180.,inp.values['wy']*n.pi/180.)
                        R = tools.detect_tilt(inp.values['tx'],inp.values['ty'],inp.values['tz'])
                        d_out = n.dot(R,n.array([[0],
                                                [(dety-inp.values['cy'])*inp.values['py']],
                                                [(detz-inp.values['cz'])*inp.values['pz']]])) 
                        d_out = d_out + n.array([[inp.values['L']],[0],[0]]) - n.dot(Omega,n.array([[x],[y],[z]]))
                        d_out = n.dot(n.transpose(Omega),d_out/n.sqrt(n.sum(d_out**2)))
                        d_in = n.dot(n.transpose(Omega),n.array([[-1],[0],[0]])) #from grain to source!
                        # distance to bounding planes assuming d_in and d_out unit vectors
                        xdist_in = 1e6
                        xdist_out = 1e6
                        ydist_in = 1e6
                        ydist_out = 1e6
                        if d_in[0] < 0:
                            xdist_in = (x-xmin)/abs(d_in[0,0])
                        elif d_in[0] > 0:
                            xdist_in = (xmax-x)/abs(d_in[0,0])
                        if d_in[1] < 0:
                            ydist_in = (y-ymin)/abs(d_in[1,0])
                        elif d_in[1] > 0:
                            ydist_in = (ymax-y)/abs(d_in[1,0])
                        dist_in = min(xdist_in,ydist_in)
                        if d_out[0] < 0:
                            xdist_out = (x-xmin)/abs(d_out[0,0])
                        elif d_out[0] > 0:
                            xdist_out = (xmax-x)/abs(d_out[0,0])
                        if d_out[1] < 0:
                            ydist_out = (y-ymin)/abs(d_out[1,0])
                        elif d_out[1] > 0:
                            ydist_out = (ymax-y)/abs(d_out[1,0])
                        dist_out = min(xdist_out,ydist_out)
                        dist = dist_in + dist_out
                        # abs_mu in mm-1 and dist in microns
                        exponential = n.exp(inp.fit['abs_mu']*dist*1e-3)
                        #print h,k,l,w,x,y,z,dist_in*1e-3,dist_out*1e-3,dist*1e-3,exponential
                    else:
                        exponential = 1.
                    
                    # Absorption correction end
                    value = False
                    for m in range(len(hkl)):
                        if hkl[m][0] == h and hkl[m][1] == k and hkl[m][2] == l:
                            inp.volume[i][j] = inp.F2vol[inp.id[i][j]]/hkl[m][3]*exponential
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
                        mad(data[i],rej,inp.fit['rej_vol'])
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
            
                                
def mean_ia(inp,limit,only=None):
        """
        Calculate the internal angle for each peak and store in inp.mena_ia[inp.no_grains][inp.nrefl[i]]
        Jette Oddershede Januar 2009
        """
        
        import build_fcn
        build_fcn.FCN(inp)
        import fcn
        reload(fcn)

        delete = 0
        for i in range(inp.no_grains):
            if i+1 in inp.fit['skip']:
                pass
            else:
                rod = n.array([inp.rod[i][0]+inp.values['rodx%s' %i],inp.rod[i][1]+inp.values['rody%s' %i],inp.rod[i][2]+inp.values['rodz%s' %i]])
                for j in range(inp.nrefl[i]-1,-1,-1):
                    Omega = tools.form_omega_mat_general(inp.w[inp.id[i][j]]*n.pi/180,inp.values['wx']*n.pi/180,inp.values['wy']*n.pi/180)
                    gexp = fcn.gexp(inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],
                                    inp.values['wx'],inp.values['wy'],inp.values['tx'],inp.values['ty'],inp.values['tz'],
                                    inp.values['py'],inp.values['pz'],inp.values['cy'],inp.values['cz'],inp.values['L'],
                                    inp.values['x%s' %i],inp.values['y%s' %i],inp.values['z%s' %i])
                    gcalc = fcn.gcalc(inp.values['a'],inp.values['b'],inp.values['c'],inp.values['alpha'],inp.values['beta'],inp.values['gamma'],
                                      inp.h[i][j],inp.k[i][j],inp.l[i][j],
                                      rod[0],rod[1],rod[2],
                                      inp.values['epsaa%s' %i],inp.values['epsab%s' %i],inp.values['epsac%s' %i],
                                      inp.values['epsbb%s' %i],inp.values['epsbc%s' %i],inp.values['epscc%s' %i])
#                    gexp = n.dot(Omega,gexp)
#                    gcalc = n.dot(Omega,gcalc)
                    inp.mean_ia[i][j] = IA(n.transpose(gexp)[0],n.transpose(gcalc)[0])
#                    print i+1,inp.mean_ia[i][j]
                    if inp.mean_ia[i][j] > limit:
                        delete = delete + 1
                        reject(inp,i,j,'ia')

        if only != []:
            print 'Rejected', delete, 'reflection based on internal angles'
        insignificant(inp)


def mean_ia_old(inp,limit,only=None):
        """
        Calculate the internal angle for each peak and store in inp.mena_ia[inp.no_grains][inp.nrefl[i]]
        Jette Oddershede Januar 2009
        """
        
        import build_fcn
        build_fcn.FCN(inp)
        import fcn
        reload(fcn)

        for i in range(inp.no_grains):
            if i+1 in inp.fit['skip']:
                pass
            else:
                rod = n.array([inp.rod[i][0]+inp.values['rodx%s' %i],inp.rod[i][1]+inp.values['rody%s' %i],inp.rod[i][2]+inp.values['rodz%s' %i]])
                for j in range(inp.nrefl[i]):
                    Omega = tools.form_omega_mat_general(inp.w[inp.id[i][j]]*n.pi/180,inp.values['wx']*n.pi/180,inp.values['wy']*n.pi/180)
                    gexp = fcn.gexp(inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],
                                    inp.values['wx'],inp.values['wy'],inp.values['tx'],inp.values['ty'],inp.values['tz'],
                                    inp.values['py'],inp.values['pz'],inp.values['cy'],inp.values['cz'],inp.values['L'],
                                    inp.values['x%s' %i],inp.values['y%s' %i],inp.values['z%s' %i])
                    gcalc = fcn.gcalc(inp.values['a'],inp.values['b'],inp.values['c'],inp.values['alpha'],inp.values['beta'],inp.values['gamma'],
                                      inp.h[i][j],inp.k[i][j],inp.l[i][j],
                                      rod[0],rod[1],rod[2],
                                      inp.values['epsaa%s' %i],inp.values['epsab%s' %i],inp.values['epsac%s' %i],
                                      inp.values['epsbb%s' %i],inp.values['epsbc%s' %i],inp.values['epscc%s' %i])
                    gexp = n.dot(Omega,gexp)
                    gcalc = n.dot(Omega,gcalc)
#                    print int(inp.h[i][j]), int(inp.k[i][j]), int(inp.l[i][j]),inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],'gexp', 2*n.pi*n.transpose(gexp)[0]/inp.param['wavelength']
#                    print int(inp.h[i][j]), int(inp.k[i][j]), int(inp.l[i][j]),inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],'gcalc', 2*n.pi*n.transpose(gcalc)[0]/inp.param['wavelength']
                    inp.mean_ia[i][j] = IA(n.transpose(gexp)[0],n.transpose(gcalc)[0])
#                    inp.mean_ia[i][j] = IAforrod(n.transpose(gexp)[0],n.transpose(gcalc)[0],rod)
#                    print inp.h[i][j], inp.k[i][j], inp.l[i][j], inp.id[i][j], inp.mean_ia[i][j]

        data = deepcopy(inp.mean_ia)
        maxia = [0]*inp.no_grains
        for i in range(inp.no_grains):
            data[i].sort()
            if i+1 in inp.fit['skip']:
                pass
            else:       
                mean = n.sum(data[i])/len(data[i])
                medi = median(data[i])
#                print i, len(data[i]), medi, mean,'\n',data[i]
                while mean > limit*medi:
                    data[i].pop()
                    mean = n.sum(data[i])/inp.nrefl[i]
                    medi = median(data[i])
                maxia[i] = max(data[i])
#                print i, len(data[i]),medi,mean,'\n',data[i],'\n'
        
        delete = 0
        if only==None:
            only = range(1,1+inp.no_grains)        
        for i in range(inp.no_grains):
            if i+1 in inp.fit['skip'] or i+1 not in only:
                pass
            else:               
                for j in range(inp.nrefl[i]-1,-1,-1): # loop backwards to make pop work
                    if inp.mean_ia[i][j] > maxia[i]:
                        delete = delete + 1
                        reject(inp,i,j,'ia')
        if only != []:
            print 'Rejected', delete, 'reflection based on internal angles'
        insignificant(inp)


                    
def merge(inp):
        """
        This function merges grain if the fraction of similar peaks exceeds the overlap parameter
    
        Jette Oddershede, August 20 2008
        """
        
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
        if  inp.files['structure_file'] != None and inp.fit['rej_multi'] > 0:
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
                    if inp.fit['rej_multi'] < 0: # reject all multiple assignments
                        bad.append([grain[k][m],peak[k][m]])
                    elif inp.fit['rej_multi'] > 0: # let grains fight for reflections
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
        if inp.fit['rej_multi'] != 0:
            for i in range(len(bad)-1,-1,-1):
                reject(inp,bad[i][0],bad[i][1],'multi')
            print 'Delete',len(bad), 'reflection because they are assigned to more than one grain' 
            insignificant(inp)
                                    
        
def overflow(inp):

    delete = 0
    for i in range(inp.no_grains):
        for j in range(inp.nrefl[i]-1,-1,-1):
            if inp.Sww[inp.id[i][j]] < 0 and inp.Sww[inp.id[i][j]] > -20:
                reject(inp,i,j,'overflow')
                delete = delete + 1
                    
    print 'Rejected', delete, 'peaks because of overflow'
    insignificant(inp)

    
def residual(inp,limit,only=None):
        """
        Reject outliers peaks based on residuals until mean<limit*median    
        
        Jette Oddershede, Risoe DTU, May 15 2008
        """
        
        # must update inp.vars because the order here is [i][j] in stead of [id[i][j]], the latter doesn't change when peaks are rejected, the former does.
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
                    inp.residual[i][j] = fcn.peak(inp.values['a'],inp.values['b'],inp.values['c'],inp.values['alpha'],inp.values['beta'],inp.values['gamma'],
                                              inp.h[i][j],inp.k[i][j],inp.l[i][j],
                                              inp.w[inp.id[i][j]],inp.dety[inp.id[i][j]],inp.detz[inp.id[i][j]],
                                              #n.array([inp.Syy[inp.id[i][j]],inp.Szz[inp.id[i][j]],inp.Sww[inp.id[i][j]]]),
                                              inp.vars[i][j], 
                                              inp.values['wx'],inp.values['wy'],
                                              inp.values['tx'],inp.values['ty'],inp.values['tz'],
                                              inp.values['py'],inp.values['pz'],
                                              inp.values['cy'],inp.values['cz'],
                                              inp.values['L'],
                                              inp.values['x%s' %i],inp.values['y%s' %i],inp.values['z%s' %i], 
                                              inp.rod[i][0]+inp.values['rodx%s' %i],
                                              inp.rod[i][1]+inp.values['rody%s' %i],
                                              inp.rod[i][2]+inp.values['rodz%s' %i],
                                              inp.values['epsaa%s' %i],inp.values['epsab%s' %i],inp.values['epsac%s' %i], 
                                              inp.values['epsbb%s' %i],inp.values['epsbc%s' %i],inp.values['epscc%s' %i]) 
                                                    
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
        if only != []:
            print 'Rejected', delete, 'reflection based on residuals'
        insignificant(inp)
                    

def peak_spread(inp):
        """
        Calculate peak spread per grain in eta and tth
        So far no peak rejection based on these values
        
        Jette Oddershede, May 18 2011
        """
        
        for i in range(inp.no_grains):
            for j in range(inp.nrefl[i]):
                inp.spr_eta[i][j] = inp.sig_eta[inp.id[i][j]]
                inp.spr_tth[i][j] = inp.sig_tth[inp.id[i][j]]
                        
          
                                

# Helpful functions               
               
def reject(inp,i,j,message):
        """
        Reject peak j from grain i and print the message as rejectvalue
        """
        inp.vars[i].pop(j)
        inp.fit['rejectgrain'].append(i)
        inp.fit['rejectdet'].append(0)
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
        inp.mean_ia[i].pop(j)
        inp.spr_eta[i].pop(j)
        inp.spr_tth[i].pop(j)
        
               
def insignificant(inp):
        """
        Remove grains with less than inp.fit['min_refl'] peaks as being insignificant
        """
        for i in range(inp.no_grains):
            if inp.nrefl[i] < inp.fit['min_refl'] and i+1 not in inp.fit['skip']:
                inp.fit['skip'].append(i+1)
        inp.fit['skip'].sort()
        
 
def unique_list(list):
        list.sort()
        for i in range(len(list)-1,0,-1):
            if list[i] == list[i-1]:
                list.pop(i)
        return list


def IAforrod(gv1,gv2,rod):
        """
        Calculates the internal angle ia between gvectors gv1 and gv2 relative to a 
        rotation axis given as a rodrigues vector rod.
        gv1,gv2,rod n.array(1x3)
        Returns ia in degrees
    
        Jette Oddershede, Jan 2009
        """
    
        # rotation axis must be projected onto the plane of gv1xgv2 and gv2+gv2
        gvcross = n.cross(gv1,gv2)
        gvcross = gvcross/n.linalg.norm(gvcross)
        gvadd = gv1+gv2
        gvadd = gvadd/n.linalg.norm(gvadd)
        # calculate normal of this plane
        gvnorm = n.cross(gvcross,gvadd)
        gvnorm = gvnorm/n.linalg.norm(gvnorm)
        
        # rotation vector in plane of gv1xgv2 and gv2+gv2 now given as
        rot = rod - n.dot(rod,gvnorm)*gvnorm
        rot = rot/n.linalg.norm(rot)
        
        #project gv1 and gv2 onto the plane perpendicular to rot
        gv1_proj = gv1 - n.dot(gv1,rot)*rot
        gv2_proj = gv2 - n.dot(gv2,rot)*rot
    
        # The desired angle is now the angle between gv1_proj and gv2_proj
        gv1_proj = gv1_proj/n.linalg.norm(gv1_proj)
        gv2_proj = gv2_proj/n.linalg.norm(gv2_proj)
        ia = n.arccos(n.dot(gv1_proj,gv2_proj))
        
        # ia2 angle between rod and rot around norm
        rod = rod/n.linalg.norm(rod)
        ia2 = n.arccos(n.dot(rod,rot))
        norm = n.cross(rod,rot)
        norm = norm/n.linalg.norm(norm)
    
#        return (ia*180./n.pi,rot,ia2*180./n.pi,norm)
        return ia*180./n.pi

    
def IA(gv1,gv2):
        """
        Calculates the internal angle ia between gvectors gv1 and gv2
        gv1,gv2 n.array(1x3)
        Returns ia in degrees
    
        Jette Oddershede, Jan 2009
        """
    
        gv1 = gv1/n.linalg.norm(gv1)
        gv2 = gv2/n.linalg.norm(gv2)

        return n.arccos(n.dot(gv1,gv2))*180./n.pi

        
def median(numbers):
   # Sort the list and take the middle element.
   nn = len(numbers)
   copy = numbers[:] # So that "numbers" keeps its original order
   copy.sort()
   if nn & 1:         # There is an odd number of elements
      return copy[nn // 2]
   else:
      return (copy[nn // 2 - 1] + copy[nn // 2]) / 2
      
      
def median_absolute_deviation(numbers):
    med = median(numbers)    
    ad = []
    for i in range(len(numbers)):
        ad.append(abs(numbers[i]-med))
    mad = median(ad)
    return mad
      
def spread(data):
        data = n.array(data)
        vars = n.sum(data*data) - n.sum(data)**2/len(data)
        return n.sqrt(vars/(len(data)-1))
      
      
def mad(data,reject,limit):
        """
        Perform outlier rejection based on median absolute deviation
        Move all data points more than limit median absolute deviations to the reject list
        
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
    
def volume_multi(inp):
        """
        Calculate the average and spread of the grain volumes
        disregarding multiple assigned reflections
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
                        
        multi = []
        for k in range(inp.param['total_refl']):
            if len(grain[k]) > 1:
                for m in range(len(grain[k])):
                    multi.append([grain[k][m],peak[k][m]])
                    
        unique_list(multi)
        volume = deepcopy(inp.volume)
        
        for i in range(len(volume)-1,-1,-1):
            for j in range(len(volume[i])-1,-1,-1):
                if [i,j] in multi:
                    volume(i).pop(j)
                    
        avg = []
        spr = []
        for i in range(len(volume)):
            avg.append(n.sum(volume[i]))
            spr.append(spread(volume[i]))
            
        return (avg,spr)

                                        

