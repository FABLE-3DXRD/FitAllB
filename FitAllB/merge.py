
import ImageD11.columnfile as ic
import numpy as n
from copy import deepcopy
import sys

def merge(stem,thresholds,output,logfile):
    """
    Example:    merge.merge('peaks',[200,500],'merged.flt')
                This merges files peaks_t200.flt and peaks_t500.flt 
                into merged.flt
    Algorithm:  For each peak at the lowest threshold fint out to 
                which peaks at higher thresholds this corresponds. 
                If splitting(=overlap) is encountered take each 
                contributing peak at the lowest threshold without 
                overlap. If intensity sum of contributions exceeds 
                parent both the parent and its possible contributors 
                are rejected. Finally the spot3d_id is made unique.
    
    Jette Oddershede, October 2008
    """

    flt = [0]*len(thresholds)       #flt[thresholds]
    data = [0]*len(thresholds)      #data[thresholds][peaks][titles]
    peaks = [0]*len(thresholds)     #data[thresholds]
    overflow = [0]*len(thresholds)     #data[thresholds]
    singles = [0]*len(thresholds)     #data[thresholds]
    
    log = open(logfile,'w')

    # read files and construct the data list sorted according to increasing omega
    for i in range(len(thresholds)):
        flt[i] = ic.columnfile('%s_t%s.flt' %(stem,thresholds[i])) 
        if i == 0:
            titles = flt[i].titles
        if i > 0:
            assert thresholds[i] > thresholds[i-1], 'wrong order of thresholds'
            flt[i].titles == titles, 'wrong order of titles in flt %s_t%i' %(stem,thresholds[i])
        temp = deepcopy(n.transpose(flt[i].bigarray))
        temp = deepcopy(temp[n.argsort(temp,0)[:,titles.index('IMax_int')],:])
        temp = temp[::-1] # reverse order so that strongest peaks come first
        data[i] = temp.tolist()
        # remove peaks that are overflown or peaks that are smaller than 4 pixels
        for j in range(len(data[i])-1,-1,-1):
                if data[i][j][titles.index('Min_o')] == 0 and data[i][j][titles.index('Min_s')] == 0 and data[i][j][titles.index('Min_f')] == 0:
                    data[i].pop(j)
                    overflow[i] = overflow[i] + 1
                elif data[i][j][titles.index('Number_of_pixels')] > 700:
                    data[i].pop(j)
                    overflow[i] = overflow[i] + 1
                elif data[i][j][titles.index('IMax_int')] > 65000:
                    data[i].pop(j)
                    overflow[i] = overflow[i] + 1
                elif data[i][j][titles.index('Number_of_pixels')] < 2:
                    data[i].pop(j)
                    singles[i] = singles[i] + 1
                    
        peaks[i] = len(data[i])
    
    
#    print 'Threshold   read   overflown     singles     peaks'
    log.write('\nThreshold   read   overflown     singles     peaks\n')
    for i in range(len(thresholds)):                  
#            print thresholds[i],'     ',peaks[i]+overflow[i]+singles[i],'     ',overflow[i],'   ',singles[i],'   ',peaks[i]
            log.write('%i     %i     %i     %i     %i\n' %(thresholds[i],peaks[i]+overflow[i]+singles[i],overflow[i],singles[i],len(data[i])))

            
            
    assign = []
    thresholds = thresholds + [2**16]
    jmax = len(data[0])-1
    for j in range(len(data[0])-1,-1,-1):
        assign.append([data[0][j]+[0]])
        for i in range(0,len(thresholds)-1):
            if i == 0: 
                if data[0][j][titles.index('IMax_int')] < thresholds[1]:
                    data[0].pop(j)   
                    break
            else:
                for jj in range(len(data[i])-1,-1,-1):
                    if data[i][jj][titles.index('IMax_int')] < data[0][j][titles.index('IMax_int')] - 1 :
                        pass
                    elif data[i][jj][titles.index('IMax_int')] > data[0][j][titles.index('IMax_int')] + 1:
                        break
                    else:
                        if abs(data[i][jj][titles.index('IMax_int')] - data[0][j][titles.index('IMax_int')]) <.0001 and\
                            data[i][jj][titles.index('IMax_o')] == data[0][j][titles.index('IMax_o')] and \
                            data[i][jj][titles.index('IMax_s')] == data[0][j][titles.index('IMax_s')] and \
                            data[i][jj][titles.index('IMax_f')] == data[0][j][titles.index('IMax_f')]:
                            assign[jmax-j].append(data[i].pop(jj)+[i])
                            break
                if len(assign[jmax-j])<i+1 and assign[jmax-j][-1][-1] != i:
                    log.write('%f %f %f %f %s %i\n' %(data[0][j][titles.index('sc')],data[0][j][titles.index('fc')],data[0][j][titles.index('omega')],data[0][j][titles.index('IMax_int')], 'no fit at t= ',thresholds[i]))
#                    print data[0][j][titles.index('sc')],data[0][j][titles.index('fc')],data[0][j][titles.index('omega')],data[0][j][titles.index('IMax_int')], 'no fit at t= ',thresholds[i]
                elif len(assign[jmax-j])<i+1 and assign[jmax-j][-1][-1] == i and assign[jmax-j][0][-1] == 0:
                    assign[jmax-j].pop(0)
                if data[0][j][titles.index('IMax_int')] < thresholds[i+1]:
                    data[0].pop(j)   
                    break
    
    
    
#    print 'Threshold   peaks   left     after first matching'
    log.write('\nThreshold   peaks   left     after first matching\n')
    for i in range(len(thresholds)-1):                  
        if i == 0:
#            print thresholds[i],'     ',peaks[i],'     ',len(assign)
            log.write('%i     %i     %i\n' %(thresholds[i],peaks[i],len(assign)))
        else:
#            print thresholds[i],'     ',peaks[i],'     ',len(data[i])
            log.write('%i     %i     %i\n' %(thresholds[i],peaks[i],len(data[i])))
    
    
    
    for i in range(1,len(thresholds)-1):
        for j in range(len(data[i])-1,-1,-1):
            temp = []
            if data[i][j][titles.index('IMax_int')] < thresholds[i+1]:
                temp.append(data[i].pop(j)+[i])
            else:
                temp = [data[i][j]+[i]]
                for ii in range(i+1,len(thresholds)-1):
                    for jj in range(len(data[ii])-1,-1,-1):
                        if data[ii][jj][titles.index('IMax_int')] < data[i][j][titles.index('IMax_int')] - 1 :
                            pass
                        elif data[ii][jj][titles.index('IMax_int')] > data[i][j][titles.index('IMax_int')] + 1:
                            break
                        else:
                            if abs(data[ii][jj][titles.index('IMax_int')] - data[i][j][titles.index('IMax_int')]) <.0001 and\
                                data[ii][jj][titles.index('IMax_o')] == data[i][j][titles.index('IMax_o')] and \
                                data[ii][jj][titles.index('IMax_s')] == data[i][j][titles.index('IMax_s')] and \
                                data[ii][jj][titles.index('IMax_f')] == data[i][j][titles.index('IMax_f')]:
                                temp.append(data[ii].pop(jj)+[ii])
                                break
                            
                    if data[i][j][titles.index('IMax_int')] < thresholds[ii+1]:
                        data[i].pop(j)
                        break

            for k in range(len(assign)):
                if temp[0][titles.index('f_raw')] >= assign[k][0][titles.index('Min_f')] and\
                   temp[0][titles.index('f_raw')] <= assign[k][0][titles.index('Max_f')] and\
                   temp[0][titles.index('s_raw')] >= assign[k][0][titles.index('Min_s')] and\
                   temp[0][titles.index('s_raw')] <= assign[k][0][titles.index('Max_s')] and\
                   temp[0][titles.index('omega')] >= assign[k][0][titles.index('Min_o')] and\
                   temp[0][titles.index('omega')] <= assign[k][0][titles.index('Max_o')] :
                    delete = 0
                    for m in range(len(assign[k])-1,0,-1):
                        if assign[k][m][-1] > temp[0][-1]:
                            pass
                        elif assign[k][m][-1] == temp[0][-1]:
                            intensity = temp[0][titles.index('sum_intensity')]+assign[k][m][titles.index('sum_intensity')]
                            if intensity <= assign[k][m-1][titles.index('sum_intensity')]:
                                delete = temp[0][-1]
                            elif assign[k][m][-1] < temp[0][-1] and delete != 0:
                                assign[k].pop(m)
                    if delete != 0:
                        if len(assign[k]) == 0:
                            assign.pop(k)
                        break
            assign.append(temp)

    final = []
    for j in range(len(assign)):
        final.append(assign[j][0])
    final = n.array(final)
    final = final[n.argsort(final,0)[:,titles.index('omega')],:]
    
    for j in range(len(final)):
        final[j][titles.index('spot3d_id')] = j
#        print data[0][j][titles.index('omega')]
    
    f = open(output,'w')
    out = '#'
    for i in range(len(titles)):
        if titles[i] in list(ic.FORMATS.keys()):      
            out = out + ' %s' %titles[i]
    for j in range(len(final)):
        out = out + '\n'
#        print '\r',j,len(data[0]),
        sys.stdout.flush()
        for i in range(len(titles)):
            if titles[i] in list(ic.FORMATS.keys()):      
                out = out + ' ' + ic.FORMATS[titles[i]] %(final[j][titles.index(titles[i])])
    f.write(out)        
    f.close()
    
    
    
    
#    print 'Threshold   peaks   left     final'
    log.write('\nThreshold   peaks   left     final\n')
    for i in range(len(thresholds)-1):                  
        if i == 0:
#            print thresholds[i],'     ',peaks[i],'     ',len(assign)
            log.write('%i     %i     %i\n' %(thresholds[i],peaks[i],len(assign)))
        else:
#            print thresholds[i],'     ',peaks[i],'     ',len(data[i])
            log.write('%i     %i     %i\n' %(thresholds[i],peaks[i],len(data[i])))
    log.close()
    
#    for i in range(len(data[len(thresholds)-2])):
#        print data[len(thresholds)-2][i]
    

def two2one(in1,in2,out):
    """
    Example:    merge.two2one('peaks_fa.flt','peaks_fb.flt','merged.flt')
                This merges files peaks_fa.flt and peaks_fb.flt 
                into merged.flt 
    
    Jette Oddershede, October 2008
    """
    
    f1 = ic.columnfile(in1) 
    f2 = ic.columnfile(in2) 
    print(in1,f1.nrows)
    print(in2,f2.nrows)
    assert f1.titles == f2.titles, 'Attempting to merge columnfiles with different column titles'

    f = open(in1,'r')
    flt1 = f.read()
    f.close()
    f = open(in2,'r')
    flt2lines = f.readlines()
    f.close()
    flt2 = ''
    for i in range(len(flt2lines)):
        if '#' not in flt2lines[i]:
            flt2 = flt2 + flt2lines[i]
    flt = open(out,'w')
    flt.write(flt1+'\n'+flt2)
    flt.close
    print(out,f1.nrows+f2.nrows)
    
    
def spot3d_id(filename):
    """
    Example:    merge.spot3d_id('peaks.flt')
                This makes unique spot3d_id in file peaks.flt
    
    Jette Oddershede, October 2008
    """
    flt = ic.columnfile(filename)
    flt.setcolumn(n.arange(flt.nrows),'spot3d_id')
    flt.writefile(filename)

