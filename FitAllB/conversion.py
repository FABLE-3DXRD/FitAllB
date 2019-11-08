#! /usr/bin/env python
# stuff for representations of tensor components: tens.py



from copy import deepcopy
from numpy import *


sqr6i = 1./sqrt(6.)
sqr3i = 1./sqrt(3.)
sqr2i = 1./sqrt(2.)
sqr2  = sqrt(2.)
sqr2b3 = sqrt(2./3.)


def symmToMVvec(A):
    """
    convert from symmetric matrix to Mandel-Voigt vector
    representation (JVB)
    """ 
    mvvec = zeros(6, dtype='float64')
    mvvec[0] = A[0,0]
    mvvec[1] = A[1,1]
    mvvec[2] = A[2,2]
    mvvec[3] = sqr2 * A[1,2]
    mvvec[4] = sqr2 * A[0,2]
    mvvec[5] = sqr2 * A[0,1] 
    return mvvec


def MVvecToSymm(A):
    """
    convert from Mandel-Voigt vector to symmetric matrix 
    representation (JVB) 
    """ 
    symm = zeros((3, 3), dtype='float64')
    symm[0, 0] = A[0]
    symm[1, 1] = A[1]
    symm[2, 2] = A[2]
    symm[1, 2] = A[3]*sqr2i
    symm[0, 2] = A[4]*sqr2i
    symm[0, 1] = A[5]*sqr2i 
    symm[2, 1] = A[3]*sqr2i
    symm[2, 0] = A[4]*sqr2i
    symm[1, 0] = A[5]*sqr2i 
    return symm


def MVCOBMatrix(R):
    """
    GenerateS array of 6 x 6 basis transformation matrices for the
    Mandel-Voigt tensor representation in 3-D given by: 
    
    [A] = [[A_11, A_12, A_13],
           [A_12, A_22, A_23],
           [A_13, A_23, A_33]]
        |
        |
        V
    {A} = [A_11, A_22, A_33, sqrt(2)*A_23, sqrt(2)*A_13, sqrt(2)*A_12]
              
    where the operation R * A *R.T (in tensor notation) is obtained by
    the matrix-vector product [T]*{A}.
    
    USAGE
    
        T = MVCOBMatrix(R)
    
    INPUTS
    
        1) R is (3, 3) an ndarray representing a change of basis matrix
    
    OUTPUTS
    
        1) T is (6, 6), an ndarray of transformation matrices as
           described above
    
    NOTES
    
        1) Compoments of symmetric 4th-rank tensors transform in a
           manner analogous to symmetric 2nd-rank tensors in full
           matrix notation. 


    SEE ALSO


    symmToMVvec, quatToMat
    """
    T = zeros((6, 6), dtype='float64')
    
    T[0, 0] = R[0, 0]**2
    T[0, 1] = R[0, 1]**2
    T[0, 2] = R[0, 2]**2
    T[0, 3] = sqr2 * R[0, 1] * R[0, 2]
    T[0, 4] = sqr2 * R[0, 0] * R[0, 2]
    T[0, 5] = sqr2 * R[0, 0] * R[0, 1]
    T[1, 0] = R[1, 0]**2
    T[1, 1] = R[1, 1]**2
    T[1, 2] = R[1, 2]**2
    T[1, 3] = sqr2 * R[1, 1] * R[1, 2]
    T[1, 4] = sqr2 * R[1, 0] * R[1, 2]
    T[1, 5] = sqr2 * R[1, 0] * R[1, 1]
    T[2, 0] = R[2, 0]**2
    T[2, 1] = R[2, 1]**2
    T[2, 2] = R[2, 2]**2
    T[2, 3] = sqr2 * R[2, 1] * R[2, 2]
    T[2, 4] = sqr2 * R[2, 0] * R[2, 2]
    T[2, 5] = sqr2 * R[2, 0] * R[2, 1]
    T[3, 0] = sqr2 * R[1, 0] * R[2, 0]
    T[3, 1] = sqr2 * R[1, 1] * R[2, 1]
    T[3, 2] = sqr2 * R[1, 2] * R[2, 2]
    T[3, 3] = R[1, 2] * R[2, 1] + R[1, 1] * R[2, 2]
    T[3, 4] = R[1, 2] * R[2, 0] + R[1, 0] * R[2, 2]
    T[3, 5] = R[1, 1] * R[2, 0] + R[1, 0] * R[2, 1]
    T[4, 0] = sqr2 * R[0, 0] * R[2, 0]
    T[4, 1] = sqr2 * R[0, 1] * R[2, 1]
    T[4, 2] = sqr2 * R[0, 2] * R[2, 2]
    T[4, 3] = R[0, 2] * R[2, 1] + R[0, 1] * R[2, 2]
    T[4, 4] = R[0, 2] * R[2, 0] + R[0, 0] * R[2, 2]
    T[4, 5] = R[0, 1] * R[2, 0] + R[0, 0] * R[2, 1]
    T[5, 0] = sqr2 * R[0, 0] * R[1, 0]
    T[5, 1] = sqr2 * R[0, 1] * R[1, 1]
    T[5, 2] = sqr2 * R[0, 2] * R[1, 2]
    T[5, 3] = R[0, 2] * R[1, 1] + R[0, 1] * R[1, 2]
    T[5, 4] = R[0, 0] * R[1, 2] + R[0, 2] * R[1, 0]
    T[5, 5] = R[0, 1] * R[1, 0] + R[0, 0] * R[1, 1]
    return T


def normalProjectionOfMV(vec):
    # 
    # To perform n' * A * n as [N]*{A}
    #


    # normalize in place... col vectors!
    v2 = vec**2
    n  = vec / sqrt(tile(v2.sum(0), (vec.shape[0], 1)))
    
    nmat = array([
        n[0, :]**2, 
        n[1, :]**2, 
        n[2, :]**2, 
        sqr2*n[1, :]*n[2, :], 
        sqr2*n[0, :]*n[2, :], 
        sqr2*n[0, :]*n[1, :]])
    
    nmat = nmat.T
    return nmat
    
    
def grain2sample(grain,U):
    """
    Conversion of symmetric tensor in cartesian grain system to sample system
    via the use of the grain orientation matrix U
    grain = 3x3 symmetric tensor in grain system
    U = 3x3 unitary orientation matrix
    sample = 3x3 symmetric tensor in sample system
    
    Jette Oddershede September 16 2008
    """
    
    grainMV = symmToMVvec(grain)
    UMV = MVCOBMatrix(U)
    sampleMV = dot(UMV,grainMV)
    sample = MVvecToSymm(sampleMV)
    return sample
    
    
def sample2grain(sample,U):
    """
    Conversion of symmetric tensor in cartesian sample system to grain system
    via the use of the grain orientation matrix U
    sample = 3x3 symmetric tensor in sample system
    U = 3x3 unitary orientation matrix
    grain = 3x3 symmetric tensor in grain system
    
    Jette Oddershede, September 16 2008
    """
    
    sampleMV = symmToMVvec(sample)
    UinvMV = MVCOBMatrix(transpose(U))
    grainMV = dot(UinvMV,sampleMV)
    grain = MVvecToSymm(grainMV)
    return grain
    
    
def strain2stress(epsilon,C):
    """
    Conversion from strain to stress tensor using the 6x6 stiffness tensor C
    which can be formed by formStiffnessMV
    epsilon = 3x3 symmetric strain tensor
    sigma = 3x3 symmetric stress tensor
    
    Jette Oddershede, September 16 2008
    """
    epsilonMV = symmToMVvec(epsilon)
    sigmaMV = dot(C,epsilonMV)
    sigma = MVvecToSymm(sigmaMV)
    return sigma
    
    
def stress2strain(sigma,S):
    """
    Conversion from stress to strain tensor using the 6x6 compliance tensor S
    which can be formed by formComplianceMV
    sigma = 3x3 symmetric stress tensor
    epsilon = 3x3 symmetric strain tensor
    
    Jette Oddershede, September 16 2008
    """
    sigmaMV = symmToMVvec(sigma)
    epsilonMV = dot(S,sigmaMV)
    epsilon = MVvecToSymm(epsilonMV)
    return epsilon
    
    
def formStiffnessMV(crystal_system,c11=None,c12=None,c13=None,c14=None,c15=None,c16=None,
                                   c22=None,c23=None,c24=None,c25=None,c26=None,
                                   c33=None,c34=None,c35=None,c36=None,
                                   c44=None,c45=None,c46=None,
                                   c55=None,c56=None,
                                   c66=None):
    """
    Form the stiffness matrix to convert from strain to stress in the Mandel-Voigt notation
    for a given crystal system using the unique input stiffness constants (Voigt convention)
    
    Jette Oddershede, September 16 2008 after MATLAB routine by Joel Bernier
    """
    
    if crystal_system == 'isotropic':
        unique_list = 'c11,c12'
        unique = [c11,c12]
        full = [c11,c12,c12,0,0,0,c11,c12,0,0,0,c11,0,0,0,(c11-c12)/2,0,0,(c11-c12)/2,0,(c11-c12)/2]
    elif crystal_system == 'cubic':
        unique_list = 'c11,c12,c44'
        unique = [c11,c12,c44]
        full = [c11,c12,c12,0,0,0,c11,c12,0,0,0,c11,0,0,0,c44,0,0,c44,0,c44]
    elif crystal_system == 'hexagonal':
        unique_list = 'c11,c12,c13,c33,c44'
        unique = [c11,c12,c13,c33,c44]
        full = [c11,c12,c13,0,0,0,c11,c13,0,0,0,c33,0,0,0,c44,0,0,c44,0,(c11-c12)/2]
    elif crystal_system == 'trigonal_high':
        unique_list = 'c11,c12,c13,c14,c33,c44'
        unique = [c11,c12,c13,c14,c33,c44]
        full = [c11,c12,c13,c14,0,0,c11,c13,-c14,0,0,c33,0,0,0,c44,0,0,c44,c14/2,(c11-c12)/2]
    elif crystal_system == 'trigonal_low':
        unique_list = 'c11,c12,c13,c14,c25,c33,c44'
        unique = [c11,c12,c13,c14,c25,c33,c44]
        full = [c11,c12,c13,c14,-c25,0,c11,c13,-c14,c25,0,c33,0,0,0,c44,0,c25/2,c44,c14/2,(c11-c12)/2]
    elif crystal_system == 'tetragonal_high':
        unique_list = 'c11,c12,c13,c33,c44,c66'
        unique = [c11,c12,c13,c33,c44,c66]
        full = [c11,c12,c13,0,0,0,c11,c13,0,0,0,c33,0,0,0,c44,0,0,c44,0,c66]
    elif crystal_system == 'tetragonal_low':
        unique_list = 'c11,c12,c13,c16c33,c44,c66'
        unique = [c11,c12,c13,c16c33,c44,c66]
        full = [c11,c12,c13,0,0,c16,c11,c13,0,0,-c16,c33,0,0,0,c44,0,0,c44,0,c66]
    elif crystal_system == 'orthorhombic':
        unique_list = 'c11,c12,c13,c22,c23,c33,c44,c55,c66'
        unique = [c11,c12,c13,c22,c23,c33,c44,c55,c66]
        full = [c11,c12,c13,0,0,0,c22,c23,0,0,0,c33,0,0,0,c44,0,0,c55,0,c66]
    elif crystal_system == 'monoclinic':
        unique_list = 'c11,c12,c13,c15,c22,c23,c25,c33,c35,c44,c46,c55,c66'
        unique = [c11,c12,c13,c15,c22,c23,c25,c33,c35,c44,c46,c55,c66]
        full = [c11,c12,c13,0,c15,0,c22,c23,0,c25,0,c33,0,c35,0,c44,0,c46,c55,0,c66]
    elif crystal_system == 'triclinic':
        unique_list = 'c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66'
        unique = [c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66]
        full = unique
    else:
        print('crystal system', crystal_system, 'not supported')
        return
 
    assert None not in unique, 'For crytal_system %s, the following must be given:\n %s' %(crystal_system,unique_list)
    
    full = array(full)
    stiffness = zeros((6,6))
    stiffness[0][0:3] = full[0:3]
    stiffness[0][3:6] = full[3:6]*sqr2
    stiffness[1][1:3] = full[6:8]
    stiffness[1][3:6] = full[8:11]*sqr2
    stiffness[2][2] = full[11]
    stiffness[2][3:6] = full[12:15]*sqr2
    stiffness[3][3:6] = full[15:18]*2.
    stiffness[4][4:6] = full[18:20]*2.
    stiffness[5][5] = full[20]*2.
    for i in range(1,6):
        for j in range(0,i):
            stiffness[i][j] = stiffness[j][i]
            
    return stiffness
        
        
def formComplianceMV(crystal_system,s11=None,s12=None,s13=None,s14=None,s15=None,s16=None,
                                    s22=None,s23=None,s24=None,s25=None,s26=None,
                                    s33=None,s34=None,s35=None,s36=None,
                                    s44=None,s45=None,s46=None,
                                    s55=None,s56=None,
                                    s66=None):
    """
    Form the compliance matrix to convert from stress to strain in the Mandel-Voigt notation
    for a given crystal system using the unique input compliance constants (Voigt convention)
    
    Jette Oddershede, September 16 2008 after MATLAB routine by Joel Bernier
    """
    
    if crystal_system == 'isotropic':
        unique = [s11,s12]
        full = [s11,s12,s12,0,0,0,s11,s12,0,0,0,s11,0,0,0,(s11-s12)*2,0,0,(s11-s12)*2,0,(s11-s12)*2]
    elif crystal_system == 'cubic':
        unique = [s11,s12,s44]
        full = [s11,s12,s12,0,0,0,s11,s12,0,0,0,s11,0,0,0,s44,0,0,s44,0,s44]
    elif crystal_system == 'hexagonal':
        unique = [s11,s12,s13,s33,s44]
        full = [s11,s12,s13,0,0,0,s11,s13,0,0,0,s33,0,0,0,s44,0,0,s44,0,(s11-s12)*2]
    elif crystal_system == 'orthorhombic':
        unique = [s11,s12,s13,s22,s23,s33,s44,s55,s66]
        full = [s11,s12,s13,0,0,0,s22,s23,0,0,0,s33,0,0,0,s44,0,0,s55,0,s66]
    elif crystal_system == 'monoclinic':
        unique = [s11,s12,s13,s15,s22,s23,s25,s33,s35,s44,s46,s55,s66]
        full = [s11,s12,s13,0,s15,0,s22,s23,0,s25,0,s33,0,s35,0,s44,0,s46,s55,0,s66]
    elif crystal_system == 'triclinic':
        unique = [s11,s12,s13,s14,s15,s16,s22,s23,s24,s25,s26,s33,s34,s35,s36,s44,s45,s46,s55,s56,s66]
        full = unique
    else:
        print('crystal system', crystal_system, 'not supported')
        return
 
    assert None not in unique, 'Missing constant for %s' %crystal_system
    
    full = array(full)
    compliance = zeros((6,6))
    compliance[0][0:3] = full[0:3]
    compliance[0][3:6] = full[3:6]*sqr2i
    compliance[1][1:3] = full[6:8]
    compliance[1][3:6] = full[8:11]*sqr2i
    compliance[2][2] = full[11]
    compliance[2][3:6] = full[12:15]*sqr2i
    compliance[3][3:6] = full[15:18]*.5
    compliance[4][4:6] = full[18:20]*.5
    compliance[5][5] = full[20]*.5
    for i in range(1,6):
        for j in range(0,i):
            compliance[i][j] = compliance[j][i]
            
    return compliance
        
 
def covariance2MV(cov):
    """
    Input:  6x6 covariance matrix cov(a) of the vector a = (e11,e22,e33,e23,e13,e12) 
            whose components belong to a symmetric 2nd rank tensor
    Output: 6x6 covariance matrix of the vector aMV = (e11,e22,e33,sqrt(2)*e23,sqrt(2)*e13,sqrt(2)*e12),
            the Mandel-Voigt analogue to the input vector
    Code:   if i in {1,2,3} and j in {1,2,3} cov(aMV)ij = cov(a) 
            if i in {1,2,3} and j in {4,5,6} cov(aMV)ij = cov(a)*sqrt(2) 
            if i in {4,5,6} and j in {4,5,6} cov(aMV)ij = cov(a)*2 
            NB! covariance matrices are symmetric
            
    Jette Oddershede September 18 2008
    """
     
    covMV = deepcopy(cov)
    for i in range(6):
        for j in range(6):
            if i<3 and j<3:
                pass
            elif i>2 and j>2:
                covMV[i][j] = cov[i][j]*2.
            else:
                covMV[i][j] = cov[i][j]*sqr2
    
    return covMV
     

def MV2covariance(covMV):
    """
    Input:  6x6 covariance matrix of the vector aMV = (e11,e22,e33,sqrt(2)*e23,sqrt(2)*e13,sqrt(2)*e12),
            the Mandel-Voigt analogue to the input vector
    Output: 6x6 covariance matrix cov(a) of the vector a = (e11,e22,e33,e23,e13,e12) 
            whose components belong to a symmetric 2nd rank tensor
    Code:   if i in {1,2,3} and j in {1,2,3} cov(a)ij = cov(aMV) 
            if i in {1,2,3} and j in {4,5,6} cov(a)ij = cov(aMV)/sqrt(2) 
            if i in {4,5,6} and j in {4,5,6} cov(a)ij = cov(aMV)/2 
            NB! covariance matrices are symmetric
            
    Jette Oddershede September 18 2008
    """
     
    cov = deepcopy(covMV)
    for i in range(6):
        for j in range(6):
            if i<3 and j<3:
                pass
            elif i>2 and j>2:
                cov[i][j] = covMV[i][j]*.5
            else:
                cov[i][j] = covMV[i][j]*sqr2i
    
    return cov
     


     
def CovarianceTransformation(cov,T):
    """
    Input:  6x6 covariance matrix cov(a) of the vector a = (e11,e22,e33,e23,e13,e12) 
            whose components belong to a symmetric 2nd rank tensor
            and
            SMV Mandel-Voigt Compliance or CMV Mandel-Voigt Stiffness
    Jette Oddershede September 22 2008
    """
    
    covMV = covariance2MV(cov)
    covMV_trans = dot(T,dot(covMV,transpose(T)))
    cov_trans = MV2covariance(covMV_trans)
    return cov_trans
    

def CovarianceRotation(cov,U):
    """
    Input:  6x6 covariance matrix cov(a) of the vector a = (e11,e22,e33,e23,e13,e12) 
            whose components belong to a symmetric 2nd rank tensor
            and
            3x3 rotation matrix
    Jette Oddershede September 22 2008
    """
    
    covMV = covariance2MV(cov)
    UMV = MVCOBMatrix(U)
    covMV_trans = dot(UMV,dot(covMV,transpose(UMV)))
    cov_trans = MV2covariance(covMV_trans)
    return cov_trans
