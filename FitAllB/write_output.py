import numpy as n
from xfab import tools
from xfab import symmetry
from xfab import detector
import reject
#import fcn
import sys
import logging
import conversion
from copy import deepcopy
from string import split
logging.basicConfig(level=logging.DEBUG,format='%(levelname)s %(message)s')

		
def write_cov(lsqr,i):
    """
    Calculate and save the covariance matrix for grain i
    Jette Oddershede June 13th 2008
    """
    
    # clear cov file at first visit 
    filename = '%s/%s_%s_cov.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'],lsqr.inp.fit['goon'])
    if i == 0:
        f = open(filename,'w')
        f.close()
    # now open for appending
    f = open(filename,'a' )
    if i == 0:
        f.write('\n \t\n**********%s********** \n\t' %lsqr.inp.fit['goon'])
    lines = 0
    print 'grain', i+1
    for entry1 in lsqr.grains[i]:
        if lsqr.mg.fixed[entry1] == False:
            if lines == 0:
                string = '\n         '
                for entry2 in lsqr.grains[i]:
                    if lsqr.mg.fixed[entry2] == False:
                        string = string + '%13s  ' %entry2
                string = string + '\n'
                f.write(string)
                lines = 1					
            string = '%13s  ' %entry1
            for entry2 in lsqr.grains[i]:
                if lsqr.mg.fixed[entry2] == False:
                    string = string + '%8e  ' %lsqr.mg.covariance[('%s' %entry1, '%s' %entry2)] 
            string = string + '\n'
            f.write(string)
        
    f.close()		


def write_cor(lsqr,i):
    """
    Calculate and save the correlation matrix for grain i
    Jette Oddershede June 13th 2008
    """
    
    # clear cor file at first visit
    filename = '%s/%s_%s_cor.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'],lsqr.inp.fit['goon'])
    if i == 0:
        f = open(filename,'w')
        f.close()
    # now open for appending
    f = open(filename,'a' )
    if i == 0:
        f.write('\n \t\n**********%s********** \n\t' %lsqr.inp.fit['goon'])
    lines = 0
    for entry1 in lsqr.grains[i]:
        if lsqr.mg.fixed[entry1] == False:
            if lines == 0:
                string = '\n         '
                for entry2 in lsqr.grains[i]:
                    if lsqr.mg.fixed[entry2] == False:
                        string = string + '%8s  ' %entry2
                string = string + '\n'
                f.write(string)
                lines = 1					
            string = '%8s  ' %entry1
            for entry2 in lsqr.grains[i]:
                if lsqr.mg.fixed[entry2] == False:
                    string = string + '%8f  ' %(lsqr.mg.covariance[('%s' %entry1, '%s' %entry2)]/(n.sqrt(lsqr.mg.covariance[('%s' %entry1, '%s' %entry1)])*n.sqrt(lsqr.mg.covariance[('%s' %entry2, '%s' %entry2)]))) 
            string = string + '\n'
            f.write(string)
        
    f.close()		


def write_global(lsqr):
    """
    Calculate and save the correlation and covariance matrices for the global parameters
    Jette Oddershede August 20th 2008
    """
    
    # clear cor file at first visit
    filename = '%s/%s_global.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'])
    if lsqr.inp.fit['goon'] == 'start':
        f = open(filename,'w')
        f.write('NB! The beam center coordinates refer to the internal coordinate system \n')
        f.write('    unlike the values of cy and cz output in the log file \n')
        f.close()
    # now open for appending
    f = open(filename,'a' )
    f.write('\n\n\n**********%s********** \n' %lsqr.inp.fit['goon'])
    f.write('\n\n**********correlation matrix**********\n')
    lines = 0
    for entry1 in lsqr.globals:
        if lsqr.m.fixed[entry1] == False:
            if lines == 0:
                string = '\n         '
                for entry2 in lsqr.globals:
                    if lsqr.m.fixed[entry2] == False:
                        string = string + '%8s  ' %entry2
                string = string + '\n'
                f.write(string)
                lines = 1					
            string = '%8s  ' %entry1
            for entry2 in lsqr.globals:
                if lsqr.m.fixed[entry2] == False:
                    string = string + '%8f  ' %(lsqr.m.covariance[('%s' %entry1, '%s' %entry2)]/(n.sqrt(lsqr.m.covariance[('%s' %entry1, '%s' %entry1)])*n.sqrt(lsqr.m.covariance[('%s' %entry2, '%s' %entry2)]))) 
            string = string + '\n'
            f.write(string)
        
    f.write('\n\n**********covariance matrix********** \n')
    lines = 0
    for entry1 in lsqr.globals:
        if lsqr.m.fixed[entry1] == False:
            if lines == 0:
                string = '\n         '
                for entry2 in lsqr.globals:
                    if lsqr.m.fixed[entry2] == False:
                        string = string + '%8s  ' %entry2
                string = string + '\n'
                f.write(string)
                lines = 1					
            string = '%8s  ' %entry1
            for entry2 in lsqr.globals:
                if lsqr.m.fixed[entry2] == False:
                    string = string + '%8e  ' %(lsqr.m.covariance[('%s' %entry1, '%s' %entry2)]) 
            string = string + '\n'
            f.write(string)

    f.close()		


def write_values(lsqr):
    """
    Save the fitted grain parameters, pos, U and eps

    INPUT: The parameter set from the final fitting
    OUTPUT: grainno mean_IA grain_volume x y z rodx rody rodz 
            U11 U12 U13 U21 U22 U23 U31 U32 U33 
            eps11   eps22   eps33   eps23   eps13   eps12
            eps11_s eps22_s eps33_s eps23_s eps13_s eps12_s
            sig11   sig22   sig33   sig23   sig13   sig12
            sig11_s sig22_s sig33_s sig23_s sig13_s sig12_s
 
    Jette Oddershede, Risoe DTU, April 21 2008
    """


    filename = '%s/%s_%s.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'],lsqr.inp.fit['goon'])
    f = open(filename,'w')
    format = "%d "*1 + "%f "*8 + "%0.12f "*9 + "%e "*24 +"\n"
    out = "# grainno mean_IA grainvolume x y z rodx rody rodz U11 U12 U13 U21 U22 U23 U31 U32 U33 eps11 eps22 eps33 eps23 eps13 eps12 "
    out = out + "eps11_s eps22_s eps33_s eps23_s eps13_s eps12_s sig11 sig22 sig33 sig23 sig13 sig12 sig11_s sig22_s sig33_s sig23_s sig13_s sig12_s\n"
    f.write(out)
    for i in range(lsqr.inp.no_grains):
        if i+1 in lsqr.inp.fit['skip']:
            pass
        else:
            U = tools.rod_to_u([lsqr.inp.rod[i][0]+lsqr.m.values['rodx%s' %i],lsqr.inp.rod[i][1]+lsqr.m.values['rody%s' %i],lsqr.inp.rod[i][2]+lsqr.m.values['rodz%s' %i],])
            eps = n.array([[lsqr.m.values['epsaa%s' %i],lsqr.m.values['epsab%s' %i],lsqr.m.values['epsac%s' %i]],
                           [lsqr.m.values['epsab%s' %i],lsqr.m.values['epsbb%s' %i],lsqr.m.values['epsbc%s' %i]],
                           [lsqr.m.values['epsac%s' %i],lsqr.m.values['epsbc%s' %i],lsqr.m.values['epscc%s' %i]]])
            eps_s = conversion.grain2sample(eps,U)
            sig = conversion.strain2stress(eps,lsqr.inp.C)
            sig_s = conversion.grain2sample(sig,U)            
            out = format %(i+1,
    			           n.sum(lsqr.inp.mean_ia[i])/len(lsqr.inp.mean_ia[i]),#0,
	    				   sum(lsqr.inp.volume[i])/lsqr.inp.nrefl[i],
                           lsqr.m.values['x%s' %i]/1000,
                           lsqr.m.values['y%s' %i]/1000,
                           lsqr.m.values['z%s' %i]/1000,
                           lsqr.inp.rod[i][0]+lsqr.m.values['rodx%s' %i],
                           lsqr.inp.rod[i][1]+lsqr.m.values['rody%s' %i],
                           lsqr.inp.rod[i][2]+lsqr.m.values['rodz%s' %i],
                           U[0,0],
                           U[0,1],
                           U[0,2],
                           U[1,0],
                           U[1,1],
                           U[1,2],
                           U[2,0],
                           U[2,1],
                           U[2,2],
                           eps[0][0],
                           eps[1][1],
                           eps[2][2],
                           eps[1][2],
                           eps[0][2],
                           eps[0][1],
                           eps_s[0][0],
                           eps_s[1][1],
                           eps_s[2][2],
                           eps_s[1][2],
                           eps_s[0][2],
                           eps_s[0][1],
                           sig[0][0],
                           sig[1][1],
                           sig[2][2],
                           sig[1][2],
                           sig[0][2],
                           sig[0][1],
                           sig_s[0][0],
                           sig_s[1][1],
                           sig_s[2][2],
                           sig_s[1][2],
                           sig_s[0][2],
                           sig_s[0][1]
                          )
            f.write(out)
    f.close()   


def write_errors(lsqr,i):
    """
    Save the fitted grain error parameters, pos, U and eps

    INPUT: The error set from the final fitting
    OUTPUT: grainno mean_IA grainvolume x y z rodx rody rodz
            U11 U12 U13 U21 U22 U23 U31 U32 U33 
            eps11   eps22   eps33   eps23   eps13   eps12
            eps11_s eps22_s eps33_s eps23_s eps13_s eps12_s
            sig11   sig22   sig33   sig23   sig13   sig12
            sig11_s sig22_s sig33_s sig23_s sig13_s sig12_s
 
    Jette Oddershede, Risoe DTU, June 23 2008
    """
    

    # error transformation into other coordinate system and from strain to stress under construction
    U = tools.rod_to_u([lsqr.inp.rod[i][0]+lsqr.m.values['rodx%s' %i],lsqr.inp.rod[i][1]+lsqr.m.values['rody%s' %i],lsqr.inp.rod[i][2]+lsqr.m.values['rodz%s' %i],])
    eps_ref = False
    rodx_err = 0
    rody_err = 0
    rodz_err = 0
    for key in lsqr.mg.covariance.keys():
        if 'epsaa%s' %i in key:
            eps_ref = True
            break

    if eps_ref == False:
        cov_eps = n.zeros((6,6))
    else:
        free = lsqr.mg.fixed.values().count(False)
        assert free >= 6, 'wrong dimensions of covariance matrix'
        covariance = n.zeros((free,free))
        j1 = 0
        for entry1 in lsqr.grains[i]:
            if lsqr.mg.fixed[entry1] == False:
                j2 = 0
                for entry2 in lsqr.grains[i]:
                    if lsqr.mg.fixed[entry2] == False:
                        covariance[j1][j2] = lsqr.mg.covariance[('%s' %entry1, '%s' %entry2)]
                        j2 = j2 + 1
                j1 = j1 + 1
                
        if free == 6:
            cov_eps = covariance
        else:
            derivatives = n.linalg.inv(covariance)
            derivatives_rod = derivatives[len(derivatives)-9:len(derivatives)-6,len(derivatives)-9:len(derivatives)-6]
            cov_rod = n.linalg.inv(derivatives_rod)
            vars_rod = n.linalg.eig(cov_rod)[0]
            derivatives = derivatives[len(derivatives)-6:,len(derivatives)-6:]
            cov_eps = n.linalg.inv(derivatives)

    # old way of doing it including covariances with orientations in strain errors
    try:
        cov_eps = n.array([[lsqr.mg.covariance[('epsaa%s' %i, 'epsaa%s' %i)],lsqr.mg.covariance[('epsaa%s' %i, 'epsbb%s' %i)],lsqr.mg.covariance[('epsaa%s' %i, 'epscc%s' %i)],
                        lsqr.mg.covariance[('epsaa%s' %i, 'epsbc%s' %i)],lsqr.mg.covariance[('epsaa%s' %i, 'epsac%s' %i)],lsqr.mg.covariance[('epsaa%s' %i, 'epsab%s' %i)]],
                       [lsqr.mg.covariance[('epsbb%s' %i, 'epsaa%s' %i)],lsqr.mg.covariance[('epsbb%s' %i, 'epsbb%s' %i)],lsqr.mg.covariance[('epsbb%s' %i, 'epscc%s' %i)],
                        lsqr.mg.covariance[('epsbb%s' %i, 'epsbc%s' %i)],lsqr.mg.covariance[('epsbb%s' %i, 'epsac%s' %i)],lsqr.mg.covariance[('epsbb%s' %i, 'epsab%s' %i)]],
                       [lsqr.mg.covariance[('epscc%s' %i, 'epsaa%s' %i)],lsqr.mg.covariance[('epscc%s' %i, 'epsbb%s' %i)],lsqr.mg.covariance[('epscc%s' %i, 'epscc%s' %i)],
                        lsqr.mg.covariance[('epscc%s' %i, 'epsbc%s' %i)],lsqr.mg.covariance[('epscc%s' %i, 'epsac%s' %i)],lsqr.mg.covariance[('epscc%s' %i, 'epsab%s' %i)]],
                       [lsqr.mg.covariance[('epsbc%s' %i, 'epsaa%s' %i)],lsqr.mg.covariance[('epsbc%s' %i, 'epsbb%s' %i)],lsqr.mg.covariance[('epsbc%s' %i, 'epscc%s' %i)],
                        lsqr.mg.covariance[('epsbc%s' %i, 'epsbc%s' %i)],lsqr.mg.covariance[('epsbc%s' %i, 'epsac%s' %i)],lsqr.mg.covariance[('epsbc%s' %i, 'epsab%s' %i)]],
                       [lsqr.mg.covariance[('epsac%s' %i, 'epsaa%s' %i)],lsqr.mg.covariance[('epsac%s' %i, 'epsbb%s' %i)],lsqr.mg.covariance[('epsac%s' %i, 'epscc%s' %i)],
                        lsqr.mg.covariance[('epsac%s' %i, 'epsbc%s' %i)],lsqr.mg.covariance[('epsac%s' %i, 'epsac%s' %i)],lsqr.mg.covariance[('epsac%s' %i, 'epsab%s' %i)]],
                       [lsqr.mg.covariance[('epsab%s' %i, 'epsaa%s' %i)],lsqr.mg.covariance[('epsab%s' %i, 'epsbb%s' %i)],lsqr.mg.covariance[('epsab%s' %i, 'epscc%s' %i)],
                        lsqr.mg.covariance[('epsab%s' %i, 'epsbc%s' %i)],lsqr.mg.covariance[('epsab%s' %i, 'epsac%s' %i)],lsqr.mg.covariance[('epsab%s' %i, 'epsab%s' %i)]]])
    except:
        cov_eps = n.zeros((6,6))
            
            
            
            
            
    cov_eps_s = conversion.CovarianceRotation(cov_eps,U)
    cov_sig = conversion.CovarianceTransformation(cov_eps,lsqr.inp.C)
    cov_sig_s = conversion.CovarianceRotation(cov_sig,U)
    
    filename = '%s/%s_errors.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'])
    filename2 = '%s/%s_%s_errors.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'],lsqr.inp.fit['goon'])
    try:
        f = open(filename,'r')
        lines = f.readlines()
        f.close()
    except:
        f = open(filename,'w')
        out = "# grainno mean_IA grainvolume x y z rodx rody rodz U11 U12 U13 U21 U22 U23 U31 U32 U33 eps11 eps22 eps33 eps23 eps13 eps12 "
        out = out + "eps11_s eps22_s eps33_s eps23_s eps13_s eps12_s sig11 sig22 sig33 sig23 sig13 sig12 sig11_s sig22_s sig33_s sig23_s sig13_s sig12_s\n"
        f.write(out)
        f.close()   
        f = open(filename,'r')
        lines = f.readlines()
        f.close()

    index = 0
    for j in range(1,len(lines)):
        if i+1 == eval(split(lines[j])[0]):
            index = j
            break
    if index == 0:
        index = len(lines)
        lines.append('')
    
    format = "%d "*1 + "%f "*8 + "%0.1f "*9 + "%e "*24 +"\n"
    
    lines[index] = format %(i+1,
			       reject.spread(lsqr.inp.mean_ia[i]),#0
                   reject.spread(lsqr.inp.volume[i]),
                   lsqr.mg.errors['x%s' %i]/1000,
                   lsqr.mg.errors['y%s' %i]/1000,
                   lsqr.mg.errors['z%s' %i]/1000,
                   lsqr.mg.errors['rodx%s' %i],
                   lsqr.mg.errors['rody%s' %i],
                   lsqr.mg.errors['rodz%s' %i],
  			       0,
   			       0,
   			       0,
   			       0,
   			       0,
   			       0,
   			       0,
   			       0,
   			       0,
                   n.sqrt(cov_eps[0][0]),
                   n.sqrt(cov_eps[1][1]),
                   n.sqrt(cov_eps[2][2]),
                   n.sqrt(cov_eps[3][3]),
                   n.sqrt(cov_eps[4][4]),
                   n.sqrt(cov_eps[5][5]),
                   n.sqrt(cov_eps_s[0][0]),
                   n.sqrt(cov_eps_s[1][1]),
                   n.sqrt(cov_eps_s[2][2]),
                   n.sqrt(cov_eps_s[3][3]),
                   n.sqrt(cov_eps_s[4][4]),
                   n.sqrt(cov_eps_s[5][5]),
                   n.sqrt(cov_sig[0][0]),
                   n.sqrt(cov_sig[1][1]),
                   n.sqrt(cov_sig[2][2]),
                   n.sqrt(cov_sig[3][3]),
                   n.sqrt(cov_sig[4][4]),
                   n.sqrt(cov_sig[5][5]),
                   n.sqrt(cov_sig_s[0][0]),
                   n.sqrt(cov_sig_s[1][1]),
                   n.sqrt(cov_sig_s[2][2]),
                   n.sqrt(cov_sig_s[3][3]),
                   n.sqrt(cov_sig_s[4][4]),
                   n.sqrt(cov_sig_s[5][5])
                  )
    f = open(filename,'w')
    for j in range(len(lines)):
        f.write(lines[j])
    f.close()   
    f = open(filename2,'w')
    for j in range(len(lines)):
        f.write(lines[j])
    f.close()   

		
def write_log(lsqr):
    """
    Write fitting and rejection info in the log file
    """

    filename = '%s/%s_log.log' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'])
    # clear log file at first visit
    filename = '%s/%s_log.log' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'])
    if lsqr.inp.fit['goon'] == 'start':
        f = open(filename,'w')
        f.close()
    # now open for appending
    f = open(filename,'a' )
    # start writing output
    if lsqr.inp.fit['goon'] == 'start':
        f.write('%s\n' %lsqr.inp.fit['title'])
        f.write('Number of input grains %i \n' %lsqr.inp.no_grains)
        f.write('Number of fitted grains %i \n' %(lsqr.inp.no_grains-len(lsqr.inp.fit['skip'])))
        if len(lsqr.inp.fit['skip']) > 0:
            f.write('Skip the following grain numbers: %s\n' %lsqr.inp.fit['skip'])
    f.write('\n \t\n**********%s********** \n\t\n' %lsqr.inp.fit['goon'])
    if lsqr.ref == True:
        if 'start' in lsqr.inp.fit['goon']:
            f.write('Tolerance, %e \n' %lsqr.m.tol)
        else:
            f.write('Tolerance, %e \n' %lsqr.mg.tol)
    f.write('Time %i s \n' %lsqr.time)
    f.write('Value %e \n \t\n' %lsqr.fval)
    f.write('Grain data file: %s_%s.txt \n' %(lsqr.inp.fit['stem'],lsqr.inp.fit['goon']))
    # print values and errors of global parameters
    for entries in lsqr.globals:
        if 'c' not in entries: # skip centres, these must be converted to detector.par convention before output
            if lsqr.m.fixed[entries] == True or lsqr.m.fixed[entries] == True:
                f.write('%s %f\n' %(entries, lsqr.m.values[entries]))
            else:
                f.write('%s %f +- %f\n' %(entries, lsqr.m.values[entries], lsqr.m.errors[entries]))
                    
    # convert beam center to detector.par convention
    (z_center, y_center) = detector.detyz_to_xy([lsqr.m.values['cy'],lsqr.m.values['cz']],
                                                lsqr.inp.param['o11'],lsqr.inp.param['o12'],lsqr.inp.param['o21'],lsqr.inp.param['o22'],
                                                lsqr.inp.fit['dety_size'],lsqr.inp.fit['detz_size'])
    (z_error, y_error) = n.dot(n.array([[abs(lsqr.inp.param['o11']),abs(lsqr.inp.param['o12'])],[abs(lsqr.inp.param['o21']),abs(lsqr.inp.param['o22'])]]),
                               n.array([lsqr.m.errors['cy'],lsqr.m.errors['cz']]))
    if lsqr.m.fixed['cy'] == True:
        f.write('%s %f\n' %('cy', y_center))
        f.write('%s %f\n' %('cz', z_center))
    else:
        f.write('%s %f +- %f\n' %('cy', y_center, y_error))
        f.write('%s %f +- %f\n' %('cz', z_center, z_error))
       
	# print info on poor grains and rejected peaks	
    f.write('\nPoor grains: %s' %lsqr.inp.fit['poor'])
    f.write('\nPoor grains values: %s' %lsqr.poor_value)
    f.write('\nPoor grains nrefl: %s\n' %lsqr.poor_nrefl)
    f.write('\nSkip grains: %s\n' %lsqr.inp.fit['skip'])
    f.write('\nNumber of refined grains: %s\n' %(lsqr.inp.no_grains-len(lsqr.inp.fit['skip'])))
    f.write('Number of rejected outliers using limit %f (new, total): (%i , %i)' %(lsqr.inp.fit['limit'][1],lsqr.inp.newreject,lsqr.inp.fit['outliers']))
    observations = 0
    for i in range(lsqr.inp.no_grains):
        if i+1 in lsqr.inp.fit['skip']:
            pass
        else:
            observations = observations + lsqr.inp.nrefl[i]
    f.write('\nNumber of assigned peaks: %i\n' %observations)                
 
    f.close()

        
def write_rej(inp, message = None):
    """
    Write rejection info in the rej file
    """
          
    if inp.files['rej_file'] == None:
        filename = '%s/%s_rej.txt' %(inp.fit['direc'],inp.fit['stem'])
    else:
        filename = inp.files['rej_file']
    # open for appending and write message
    f = open(filename,'a')
    if message != None:
        f.write('\n\n'+message+'\n')
    f.write('Skip grains: %s' %inp.fit['skip'])
    # write rejected peaks
    for j in range(inp.fit['outliers']-inp.newreject,inp.fit['outliers']):
        f.write('\n%i Rejected peak id %i from grain %i (hkl: %i %i %i ): %s' 
                %(j+1,inp.fit['rejectid'][j],inp.fit['rejectgrain'][j]+1,inp.fit['hh'][j],inp.fit['kk'][j],inp.fit['ll'][j],inp.fit['rejectvalue'][j]))
 
    f.close()
    
def write_gvectors(lsqr):
    """
    Write all the model gvectors used in the fit (fcn.gcalc) in _gcalc.txt 
    and all the experimental gvectors (fcn.gexp) in _gexp.txt
    """
    
    file_gcalc = '%s/%s_gcalc.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'])
    file_gexp = '%s/%s_gexp.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'])
    
    fcalc = open(file_gcalc,'w')
    fexp = open(file_gexp,'w')
    
    fcalc.write('# peak omega dety detz gv1 gv2 gv3 h k l IA grainno\n')
    fexp.write('# peak omega dety detz gv1 gv2 gv3 h k l IA grainno\n')

    format = "%d "*1 + "%f "*3 + "%0.12f "*3 + "%d "*3 + "%0.12f " + "%d " + "\n"
    
    for i in range(lsqr.inp.no_grains):
        if i+1 in lsqr.inp.fit['skip']:
            pass
        else:
            for j in range(lsqr.inp.nrefl[i]):
                gcalc = fcn.gcalc(lsqr.inp.h[i][j],lsqr.inp.k[i][j],lsqr.inp.l[i][j],
                                  lsqr.inp.rod[i][0]+lsqr.m.values['rodx%s' %i],
                                  lsqr.inp.rod[i][1]+lsqr.m.values['rodx%s' %i],
                                  lsqr.inp.rod[i][2]+lsqr.m.values['rodx%s' %i],
                                  lsqr.m.values['epsaa%s' %i],lsqr.m.values['epsab%s' %i],lsqr.m.values['epsac%s' %i],
                                  lsqr.m.values['epsbb%s' %i],lsqr.m.values['epsbc%s' %i],lsqr.m.values['epscc%s' %i])
                gexp = fcn.gexp(lsqr.inp.w[lsqr.inp.id[i][j]],
                                lsqr.inp.dety[lsqr.inp.id[i][j]],
                                lsqr.inp.detz[lsqr.inp.id[i][j]],
                                lsqr.m.values['wx'],lsqr.m.values['wy'],
                                lsqr.m.values['tx'],lsqr.m.values['ty'],lsqr.m.values['tz'],
                                lsqr.m.values['py'],lsqr.m.values['pz'],
                                lsqr.m.values['cy'],lsqr.m.values['cz'],lsqr.m.values['L'],
                                lsqr.m.values['x%s' %i],lsqr.m.values['y%s' %i],lsqr.m.values['z%s' %i])
                outcalc = format  %(lsqr.inp.id[i][j],
                                    lsqr.inp.w[lsqr.inp.id[i][j]],
                                    lsqr.inp.dety[lsqr.inp.id[i][j]],
                                    lsqr.inp.detz[lsqr.inp.id[i][j]],
                                    gcalc[0][0],gcalc[1][0],gcalc[2][0],
                                    lsqr.inp.h[i][j],lsqr.inp.k[i][j],lsqr.inp.l[i][j],
                                    lsqr.inp.mean_ia[i][j],i+1)
                outexp = format   %(lsqr.inp.id[i][j],
                                    lsqr.inp.w[lsqr.inp.id[i][j]],
                                    lsqr.inp.dety[lsqr.inp.id[i][j]],
                                    lsqr.inp.detz[lsqr.inp.id[i][j]],
                                    gexp[0][0],gexp[1][0],gexp[2][0],
                                    lsqr.inp.h[i][j],lsqr.inp.k[i][j],lsqr.inp.l[i][j],
                                    lsqr.inp.mean_ia[i][j],i+1)
                fcalc.write(outcalc)       
                fexp.write(outexp)
                   
    fcalc.close()
    fexp.close()
                             
    pos_true = n.array([-421.679, -374.656, -4.720]) 
    U_true = n.array([[-0.990683685084, -0.136182028039, 0.000539766854],
                       [0.125499230349, -0.914494074861, -0.384643380583], 
                       [0.052875129230, -0.380992181394, 0.923065099777]]) 
    eps_true = n.array([-5.940498e-006, 7.406837e-004, 3.340944e-004, -6.852152e-004, 3.126998e-004, 2.355405e-004])
    B_true = tools.epsilon_to_b(eps_true,lsqr.inp.param['unit_cell'])
    UB_true = n.dot(U_true,B_true)

    pos = n.array([lsqr.m.values['x0'],lsqr.m.values['y0'],lsqr.m.values['z0']])
    U = tools.rod_to_u([lsqr.inp.rod[0][0]+lsqr.m.values['rodx0'],lsqr.inp.rod[0][1]+lsqr.m.values['rody0'],lsqr.inp.rod[0][2]+lsqr.m.values['rodz0']])
    eps = n.array([lsqr.m.values['epsaa0'],lsqr.m.values['epsab0'],lsqr.m.values['epsac0'],
                   lsqr.m.values['epsbb0'],lsqr.m.values['epsbc0'],lsqr.m.values['epscc0']])
    B = tools.epsilon_to_b(eps,lsqr.inp.param['unit_cell'])
    UB = n.dot(U,B)
    
    mislist = symmetry.Umis(U_true,U,7)
    mis = 1.
    for k in range(len(mislist)):
        if mislist[k][1] < mis:
            mis = mislist[k][1]

    file_ub = '%s/%s_ub.txt' %(lsqr.inp.fit['direc'],lsqr.inp.fit['stem'])
    fub = open(file_ub,'w')

    posformat = "%f "*3 + "\n"
    format = "%0.12f "*3 + "\n"
    fub.write('pos_true\n')
    fub.write(posformat %(pos_true[0],pos_true[1],pos_true[2]))
    fub.write('U_true\n')
    fub.write(format %(U_true[0][0],U_true[0][1],U_true[0][2]))
    fub.write(format %(U_true[1][0],U_true[1][1],U_true[1][2]))
    fub.write(format %(U_true[2][0],U_true[2][1],U_true[2][2]))    
    fub.write('B_true\n')
    fub.write(format %(B_true[0][0],B_true[0][1],B_true[0][2]))
    fub.write(format %(B_true[1][0],B_true[1][1],B_true[1][2]))
    fub.write(format %(B_true[2][0],B_true[2][1],B_true[2][2]))    
    fub.write('UB_true\n')
    fub.write(format %(UB_true[0][0],UB_true[0][1],UB_true[0][2]))
    fub.write(format %(UB_true[1][0],UB_true[1][1],UB_true[1][2]))
    fub.write(format %(UB_true[2][0],UB_true[2][1],UB_true[2][2]))    
    fub.write('pos\n')
    fub.write(posformat %(pos[0],pos[1],pos[2]))
    fub.write('U\n')
    fub.write(format %(U[0][0],U[0][1],U[0][2]))
    fub.write(format %(U[1][0],U[1][1],U[1][2]))
    fub.write(format %(U[2][0],U[2][1],U[2][2]))    
    fub.write('B\n')
    fub.write(format %(B[0][0],B[0][1],B[0][2]))
    fub.write(format %(B[1][0],B[1][1],B[1][2]))
    fub.write(format %(B[2][0],B[2][1],B[2][2]))    
    fub.write('UB\n')
    fub.write(format %(UB[0][0],UB[0][1],UB[0][2]))
    fub.write(format %(UB[1][0],UB[1][1],UB[1][2]))
    fub.write(format %(UB[2][0],UB[2][1],UB[2][2]))    
    fub.write('misorientation: %s\n' %mis)
    
    
    