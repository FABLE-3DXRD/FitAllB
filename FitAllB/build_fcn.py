def FCN(inp):	
    """
    Function to build fcn.py, the module containing the functions to be 
    minimised by Minuit. This is necessary because Minuit must be called 
    with a function where all refinable parameters are given in the 
    function call, thus strings/arrays are not allowed. fcn.py contains 
    least squares expresions for refining the difference between 
    experimental and calculated g-vectors
	
    Helpful functions:
        gexp: 	experimental g-vector
        gcalc: 	calculated g-vector
        peak: 	single peak contribution to the least squares sum
                (gexp-gcalc)T.cov_inv.(gexp-gcalc) 
                where cov is diagonal and the elements are given by vars	
    Functions to be minimised:
        FCNpeak:	simple call of peak 
        FCNgrain:	peak summed over a grain
        FCN:		peak summed over entire sample
		
    These are all called with the following refinable parameters:
        globals: 		wx,wy,tx,ty,tz,py,pz,cy,cz,L
        for each grain: x,y,z,rodx,rody,rodz,
                        epsaa,epsab,epsac,epsbb,epsbc,epscc
		
    Additional needed parameters that are written in fcn.py
        wavelength:	wavelength (A)	
        unit_cell:	[a,b,c,alpha,beta,gamma] of unstrained cell (A and deg)
        O:			[[o11,o12],[o21,o22]] detector flip matrix
        no_grains:  number of grains
        nrefl[i]:	list of number of assigned reflections for each grain, 
                    i=range(no_refl)
        id[i][j]:	list of peak ids for each assigned reflection, 
                    i=range(no_refl) and j=range(nrefl[i])
        h[i][j]:    list of h,k,l for each assigned reflection, 
        k[i][j]:	i=range(no_refl) and j=range(nrefl[i])
        l[i][j]:
        w[k]:		omega for each measured reflection, 
                    id[i][j] used for identification, k=range(max(id[i][j])) 
        dety[k]:	dety for each measured reflection, 
                    id[i][j] used for identification, k=range(max(id[i][j])) 
        detz[k]:	detz for each measured reflection, 
                    id[i][j] used for identification, k=range(max(id[i][j])) 
        vars[k]:	error vector of gexp for each measuered reflection, 
                    id[i][j] used for identification, k=range(max(id[i][j])) 
		
    Jette Oddershede, July 2008
    """
	
	
    import numpy as n

    string = 'import numpy as n\n'
    string = string + 'from xfab import tools\n\n'

    # wavelength
    string = string + 'wavelength = %f \n' %inp.param['wavelength']
    # unit cell
    string = string + 'unit_cell = n.array([%f,%f,%f,%f,%f,%f]) \n' %(inp.param['cell__a'],inp.param['cell__b'],inp.param['cell__c'],inp.param['cell_alpha'],inp.param['cell_beta'],inp.param['cell_gamma'])
    # detector flip
    string = string + "O = n.array([[%i,%i],[%i,%i]])\n"     %(inp.param['o11'],inp.param['o12'],inp.param['o21'],inp.param['o22'])
    # no_grains
    string = string + "no_grains = %i \n" %inp.no_grains
    string = string + "skip = %s \n" %inp.fit['skip']
    # nrefl
    string = string + "nrefl = [" 
    for i in range(inp.no_grains):
        string = string + '%i,' %inp.nrefl[i]
    string = string + ']\n\n'
    # rodrigues vector from GrainSpotter since only deviations from this are refined
    string = string + "rod = [" 
    for i in range(inp.no_grains):
        string = string + '%s,' %inp.rod[i]
    string = string + ']\n\n'

    # id
    string = string + "id = [\n" 
    for i in range(inp.no_grains):
        string = string + '\t[' 
        for j in range(0,inp.nrefl[i]):
            string = string + '%i,' %inp.id[i][j]
        string = string + '],\n'
    string = string + '\t]\n'
	
    # h
    string = string + "h = [\n" 
    for i in range(inp.no_grains):
        string = string + '\t['
        for j in range(0,inp.nrefl[i]):
            string = string + '%i,' %inp.h[i][j]
        string = string + '],\n'
    string = string + '\t]\n'
	
    # k
    string = string + "k = [\n" 
    for i in range(inp.no_grains):
        string = string + '\t[' 
        for j in range(0,inp.nrefl[i]):
            string = string + '%i,' %inp.k[i][j]
        string = string + '],\n'
    string = string + '\t]\n'
	
    # l
    string = string + "l = [\n" 
    for i in range(inp.no_grains):
        string = string + '\t[' 
        for j in range(0,inp.nrefl[i]):
            string = string + '%i,' %inp.l[i][j]
        string = string + '],\n'
    string = string + '\t]\n'
	
    # w
    string = string + "w = [" 
    nn = 1
    for i in range(0,inp.param['total_refl']):
        nn = nn + 1
        string = string + '%f,' %inp.w[i]
        if nn == 20:
            nn = 0
            string = string + '\n\t'
    string = string + ']\n'
	
    # dety
    string = string + "dety = [" 
    nn = 1
    for i in range(0,inp.param['total_refl']):
        nn = nn + 1
        string = string + '%f,' %inp.dety[i]
        if nn == 20:
            nn = 0
            string = string + '\n\t'
    string = string + ']\n'
	
    # detz
    string = string + "detz = [" 
    nn = 1
    for i in range(0,inp.param['total_refl']):
        nn = nn + 1
        string = string + '%f,' %inp.detz[i]
        if nn == 20:
            nn = 0
            string = string + '\n\t'
    string = string + ']\n'

    # vars
    string = string + "vars = [\n" 
    for i in range(inp.no_grains):
        nn = 1
        string = string + '\t[' 
        for j in range(0,inp.nrefl[i]):
            nn = nn + 1
            string = string + '%s,' %inp.vars[i][j]
            if nn == 5:
                nn = 0
                string = string + '\n\t'
        string = string + '],\n'
    string = string + '\t]\n\n'
    
	
	
# helpful functions

    string = string + 'def gexp(w,dety,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z):\n' 
    string = string + '\t Omega = tools.quart2Omega(w,wx*n.pi/180,wy*n.pi/180)\n'
    string = string + '\t R = tools.detect_tilt(tx,ty,tz)\n'
    string = string + '\t d = n.dot(R,n.array([[0],[(dety-cy)*py],[(detz-cz)*pz]])) \n'
    string = string + '\t d = d + n.array([[L],[0],[0]]) - n.dot(Omega,n.array([[x],[y],[z]]))\n'
    string = string + '\t gexp = n.dot(n.linalg.inv(Omega),(d/n.sqrt(n.sum(d**2)) - n.array([[1],[0],[0]])))\n'
    string = string + '\t return gexp \n\n'

    string = string + 'def gcalc(h,k,l,w,dety,detz,wx,wy,rodx,rody,rodz,epsaa,epsab,epsac,epsbb,epsbc,epscc):\n' 
    string = string + '\t Omega = tools.quart2Omega(w,wx*n.pi/180,wy*n.pi/180)\n'
    string = string + "\t B = tools.epsilon2B(n.array([epsaa,epsab,epsac,epsbb,epsbc,epscc]),unit_cell)\n" 
    string = string + '\t U = tools.rod2U([rodx,rody,rodz])\n'
    string = string + '\t Bhkl = n.dot(B,n.array([[h],[k],[l]]))\n'
    string = string + '\t gcalc = (wavelength/(2*n.pi))*n.dot(U,Bhkl) \n'
    string = string + '\t return gcalc \n\n'
	
    string = string + 'def peak(h,k,l,w,dety,detz,vars,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z,rodx,rody,rodz,epsaa,epsab,epsac,epsbb,epsbc,epscc):\n' 
    string = string + '\t diff = gexp(w,dety,detz,wx,wy,tx,ty,tz,py,pz,cy,cz,L,x,y,z)-gcalc(h,k,l,w,dety,detz,wx,wy,rodx,rody,rodz,epsaa,epsab,epsac,epsbb,epsbc,epscc)\n'
    string = string + '\t result = n.sum(diff*diff/n.array([[vars[0]],[vars[1]],[vars[2]]]))\n'
    string = string + '\t return result \n\n'

	
# FCN function for all grains
	
    string = string + 'def FCN(wx,wy,tx,ty,tz,py,pz,cy,cz,L' 
    for i in range(inp.no_grains):
        string = string + ',\n \t x%s,y%s,z%s,rodx%s,rody%s,rodz%s,epsaa%s,epsab%s,epsac%s,epsbb%s,epsbc%s,epscc%s' \
	                          %(i,i,i,i,i,i,i,i,i,i,i,i)
    string = string + '):\n \n'

    # xyz 
    string = string + "\t x = ["
    for i in range(inp.no_grains):
        string = string + 'x%i,' %i
    string = string + ']\n'
    string = string + "\t y = ["
    for i in range(inp.no_grains):
        string = string + 'y%i,' %i
    string = string + ']\n'
    string = string + "\t z = ["
    for i in range(inp.no_grains):
        string = string + 'z%i,' %i
    string = string + ']\n'

    # Rodrigues vectors
    string = string + "\t rodx = ["
    for i in range(inp.no_grains):
        string = string + 'rodx%i,' %i
    string = string + ']\n'
    string = string + "\t rody = ["
    for i in range(inp.no_grains):
        string = string + 'rody%i,' %i
    string = string + ']\n'
    string = string + "\t rodz = ["
    for i in range(inp.no_grains):
        string = string + 'rodz%i,' %i
    string = string + ']\n'
	
    # strain tensor
    string = string + "\t epsaa = ["
    for i in range(inp.no_grains):
        string = string + 'epsaa%i,' %i
    string = string + ']\n'
    string = string + "\t epsab = ["
    for i in range(inp.no_grains):
        string = string + 'epsab%i,' %i
    string = string + ']\n'
    string = string + "\t epsac = ["
    for i in range(inp.no_grains):
        string = string + 'epsac%i,' %i
    string = string + ']\n'
    string = string + "\t epsbb = ["
    for i in range(inp.no_grains):
        string = string + 'epsbb%i,' %i
    string = string + ']\n'
    string = string + "\t epsbc = ["
    for i in range(inp.no_grains):
        string = string + 'epsbc%i,' %i
    string = string + ']\n'
    string = string + "\t epscc = ["
    for i in range(inp.no_grains):
        string = string + 'epscc%i,' %i
    string = string + ']\n\n'

    # initialise sum
    string = string + '\t sum = 0 \n \n'

    string = string + '\t for i in range(no_grains):\n'
    string = string + '\t\t if i+1 in skip:\n'
    string = string + '\t\t\t pass \n'
    string = string + '\t\t else:\n'
    string = string + '\t\t\t for j in range(nrefl[i]):\n'
    string = string + '\t\t\t\t sum = sum + peak(h[i][j],k[i][j],l[i][j],w[id[i][j]],dety[id[i][j]],detz[id[i][j]],vars[i][j], ' 
    string = string + 'wx,wy,tx,ty,tz,py,pz,cy,cz,L,x[i],y[i],z[i],rod[i][0]+rodx[i],rod[i][1]+rody[i],rod[i][2]+rodz[i],epsaa[i],epsab[i],epsac[i],epsbb[i],epsbc[i],epscc[i]) \n'
    string = string + '\n'
    string = string + '\t return sum \n\n\n'


	
# FCNgrain function for a single grain 
	
    string = string + 'def FCNgrain(i,wx,wy,tx,ty,tz,py,pz,cy,cz,L' 
    for i in range(inp.no_grains):
        string = string + ',\n \t x%s,y%s,z%s,rodx%s,rody%s,rodz%s,epsaa%s,epsab%s,epsac%s,epsbb%s,epsbc%s,epscc%s' \
	                          %(i,i,i,i,i,i,i,i,i,i,i,i)
    string = string + '):\n \n'
    string = string + '\t i=int(i)\n'	

    # xyz 
    string = string + "\t x = ["
    for i in range(inp.no_grains):
        string = string + 'x%i,' %i
    string = string + ']\n'
    string = string + "\t y = ["
    for i in range(inp.no_grains):
        string = string + 'y%i,' %i
    string = string + ']\n'
    string = string + "\t z = ["
    for i in range(inp.no_grains):
        string = string + 'z%i,' %i
    string = string + ']\n'

    # Rodrigues vectors
    string = string + "\t rodx = ["
    for i in range(inp.no_grains):
        string = string + 'rodx%i,' %i
    string = string + ']\n'
    string = string + "\t rody = ["
    for i in range(inp.no_grains):
        string = string + 'rody%i,' %i
    string = string + ']\n'
    string = string + "\t rodz = ["
    for i in range(inp.no_grains):
        string = string + 'rodz%i,' %i
    string = string + ']\n'
	
    # strain tensor
    string = string + "\t epsaa = ["
    for i in range(inp.no_grains):
        string = string + 'epsaa%i,' %i
    string = string + ']\n'
    string = string + "\t epsab = ["
    for i in range(inp.no_grains):
        string = string + 'epsab%i,' %i
    string = string + ']\n'
    string = string + "\t epsac = ["
    for i in range(inp.no_grains):
        string = string + 'epsac%i,' %i
    string = string + ']\n'
    string = string + "\t epsbb = ["
    for i in range(inp.no_grains):
        string = string + 'epsbb%i,' %i
    string = string + ']\n'
    string = string + "\t epsbc = ["
    for i in range(inp.no_grains):
        string = string + 'epsbc%i,' %i
    string = string + ']\n'
    string = string + "\t epscc = ["
    for i in range(inp.no_grains):
        string = string + 'epscc%i,' %i
    string = string + ']\n\n'

    # initialise sum
    string = string + '\t sum = 0 \n \n'

    string = string + '\t for j in range(nrefl[i]):\n'
    string = string + '\t\t sum = sum + peak(h[i][j],k[i][j],l[i][j],w[id[i][j]],dety[id[i][j]],detz[id[i][j]],vars[i][j], ' 
    string = string + 'wx,wy,tx,ty,tz,py,pz,cy,cz,L,x[i],y[i],z[i],rod[i][0]+rodx[i],rod[i][1]+rody[i],rod[i][2]+rodz[i],epsaa[i],epsab[i],epsac[i],epsbb[i],epsbc[i],epscc[i]) \n'
    string = string + '\n'
    string = string + '\t return sum \n\n\n'

	
	
# FCNpeak function calculating the contribution from a single peak including all variables	
	
    string = string + 'def FCNpeak(i,j,wx,wy,tx,ty,tz,py,pz,cy,cz,L' 
    for i in range(inp.no_grains):
        string = string + ',\n \t x%s,y%s,z%s,rodx%s,rody%s,rodz%s,epsaa%s,epsab%s,epsac%s,epsbb%s,epsbc%s,epscc%s' \
	                          %(i,i,i,i,i,i,i,i,i,i,i,i)
    string = string + '):\n \n'
    string = string + '\t i=int(i)\n'	
    string = string + '\t j=int(j)\n'	

    # xyz 
    string = string + "\t x = ["
    for i in range(inp.no_grains):
        string = string + 'x%i,' %i
    string = string + ']\n'
    string = string + "\t y = ["
    for i in range(inp.no_grains):
        string = string + 'y%i,' %i
    string = string + ']\n'
    string = string + "\t z = ["
    for i in range(inp.no_grains):
        string = string + 'z%i,' %i
    string = string + ']\n'

    # Rodrigues vectors
    string = string + "\t rodx = ["
    for i in range(inp.no_grains):
        string = string + 'rodx%i,' %i
    string = string + ']\n'
    string = string + "\t rody = ["
    for i in range(inp.no_grains):
        string = string + 'rody%i,' %i
    string = string + ']\n'
    string = string + "\t rodz = ["
    for i in range(inp.no_grains):
        string = string + 'rodz%i,' %i
    string = string + ']\n'
	
    # strain tensor
    string = string + "\t epsaa = ["
    for i in range(inp.no_grains):
        string = string + 'epsaa%i,' %i
    string = string + ']\n'
    string = string + "\t epsab = ["
    for i in range(inp.no_grains):
        string = string + 'epsab%i,' %i
    string = string + ']\n'
    string = string + "\t epsac = ["
    for i in range(inp.no_grains):
        string = string + 'epsac%i,' %i
    string = string + ']\n'
    string = string + "\t epsbb = ["
    for i in range(inp.no_grains):
        string = string + 'epsbb%i,' %i
    string = string + ']\n'
    string = string + "\t epsbc = ["
    for i in range(inp.no_grains):
        string = string + 'epsbc%i,' %i
    string = string + ']\n'
    string = string + "\t epscc = ["
    for i in range(inp.no_grains):
        string = string + 'epscc%i,' %i
    string = string + ']\n\n'

    string = string + '\t sum = peak(h[i][j],k[i][j],l[i][j],w[id[i][j]],dety[id[i][j]],detz[id[i][j]],vars[i][j], ' 
    string = string + 'wx,wy,tx,ty,tz,py,pz,cy,cz,L,x[i],y[i],z[i],rod[i][0]+rodx[i],rod[i][1]+rody[i],rod[i][2]+rodz[i],epsaa[i],epsab[i],epsac[i],epsbb[i],epsbc[i],epscc[i]) \n'
    string = string + '\t return sum \n\n\n'
	


    f = open('./%s/fcn.py' %inp.fit['direc'],'w')
    f.write(string)
    f.close()

