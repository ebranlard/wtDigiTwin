import numpy as np
'''
This code generates a complete mass matrix using turbine bending and mass data

Reference:
     [1]: Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex
'''

''' 
Bending mode function which requires inputs of array x, coefficients of bending, and exponents associated with each coeff
'''

def polymode_dylan(x,coeff,exp):
    # --- Alternative way of doing it
    mode   = np.zeros(x.shape)
    ddmode = np.zeros(x.shape)
    for i in range(0,len(coeff)):
        mode += coeff[i]*x**exp[i]
        ddmode += coeff[i]*x**exp[i-2]*[i-1]*i
    return mode/mode[-1], ddmode/mode[-1]


def polymode(x,coeff,exp):
    """ 
    Computes a shape function described as a polynomial expression y = a_i x^e_i
        where the a_i are given by `coeff`
              the e_i are given by `exp`
    The shape function is normalized such as to have a unitary tip deflection

    INPUTS:
        x : spanwise dimension, from 0 to L, not dimensionless!

    Returns:
        U, dU, ddU the shape, slope and curvature
    """
    mode   = np.zeros(x.shape)
    dmode  = np.zeros(x.shape)
    ddmode = np.zeros(x.shape)
    # Polynomials assume x to be dimensionless
    x_max= x[-1] 
    x_bar=x/x[-1] 
    for i in range(0,len(coeff)):
        mode += coeff[i]*x_bar**exp[i]
        if exp[i]-1>=0:
            dmode += coeff[i]*exp[i]* x_bar**(exp[i]-1)
        if exp[i]-2>=0:
            ddmode += coeff[i]*exp[i]*(exp[i]-1) * x_bar**(exp[i]-2)
    # Scaling by the tip deflection, and include x_max for derivatives since derivates were computed w.r.t. x_bar not x
    scale= mode[-1]
    return mode/scale, dmode/(scale*x_max), ddmode/(scale*x_max*x_max)

def GKBeam(s_span, EI, ddU, bOrth=False):    #Compute Kgg
    """ 
       Computes generalized stiffness matrix for a beam
       Eq.(20) from [1]
       TODO torsion

    OPTIONAL INPUTS:
     - bOrth : if true, enforce orthogonality of modes
    """
    nU = len(ddU)
    KK0 = np.zeros((6+nU,6+nU))
    Kgg = np.zeros((nU,nU))
    Kgg[:,:] = np.nan
    for i in range(0,nU):
        for j in range(0,nU):
            Kgg[i,j] = np.trapz(EI[0,:]*ddU[i][0,:]*ddU[j][0,:] + EI[1,:]*ddU[i][1,:]*ddU[j][1,:] + EI[2,:]*ddU[i][2,:]*ddU[j][2,:],s_span)
    if bOrth:
        Kgg=Kgg*np.eye(nU)
    #print('Kgg\n',Kgg)
    KK0[6:,6:] = Kgg
    return KK0
    
def GMBeam(s_G, s_span, m, U=None, bOrth=False):
    """
    Computes generalized mass matrix for a beam.
    Eq.(2) from [1]

    Performing full integration of mass matrix without shape integral functions
    NOTE: Beam assumed to be along x for now (only because of Jxx)
    
    INPUTS
     - s_G    : [m] 3 x nSpan , location of cross sections COG
     - s_span : [m] span integration variable (e.g. s_G(1,:))
     - m      : [kg/m] cross section mass along the beam
     - jxxG   : [kg.m] second moment of inertia of cross section # TODO


    OPTIONAL INPUTS:
     - bOrth : if true, enforce orthogonality of modes
     - JxxG, if omitted, assumed to be 0 # TODO
     - U , if omitted, then rigid body (6x6) mass matrix is returned
    
    """

    if U is not None:
        nU = len(U)
    else:
        nU=0
    # --- Mxx
    M = np.trapz(m,s_span)
    Mxx = np.identity(3)*M
    #print('Mxx\n',Mxx)

    # --- Mxt
    C_x = np.trapz(s_G[0,:]*m,s_span)
    C_y = np.trapz(s_G[1,:]*m,s_span)
    C_z = np.trapz(s_G[2,:]*m,s_span)
    Mxt = np.array([[0, C_z, -C_y],[-C_z, 0, C_x],[C_y, -C_x, 0]])
    #print('Mxt\n',Mxt)

    # --- Mxg
    Mxg      = np.zeros((3,nU))
    Mxg[:,:] = np.nan
    for j in range(nU):
        Mxg[0,j] = np.trapz(U[j][0,:]*m,s_span)
        Mxg[1,j] = np.trapz(U[j][1,:]*m,s_span)
        Mxg[2,j] = np.trapz(U[j][2,:]*m,s_span)
    #print('Mxg\n',Mxg)
        
    # --- Mtt
    s00 = np.trapz(s_G[0,:]*s_G[0,:]*m,s_span)
    s01 = np.trapz(s_G[0,:]*s_G[1,:]*m,s_span)
    s02 = np.trapz(s_G[0,:]*s_G[2,:]*m,s_span)
    s11 = np.trapz(s_G[1,:]*s_G[1,:]*m,s_span)
    s12 = np.trapz(s_G[1,:]*s_G[2,:]*m,s_span)
    s22 = np.trapz(s_G[2,:]*s_G[2,:]*m,s_span)

    Mtt = np.zeros((3,3))
    Mtt[0,0] = s11 + s22;   Mtt[0,1] = -s01;       Mtt[0,2] = -s02
    Mtt[1,0] = -s01;        Mtt[1,1] = s00 + s22;  Mtt[1,2] = -s12
    Mtt[2,0] = -s02;        Mtt[2,1] = -s12;       Mtt[2,2] = s00+s11
    #print('Mtt\n',Mtt)

    # --- Mtg
    Mtg      = np.zeros((3,nU))
    Mtg[:,:] = np.nan
    for j in range(nU):
        Mtg[0,j] = np.trapz((-s_G[2,:]*U[j][1,:] + s_G[1,:]*U[j][2,:])*m,s_span)
        Mtg[1,j] = np.trapz(( s_G[2,:]*U[j][0,:] - s_G[0,:]*U[j][2,:])*m,s_span)
        Mtg[2,j] = np.trapz((-s_G[1,:]*U[j][0,:] + s_G[0,:]*U[j][1,:])*m,s_span)
    #print('Mtg\n',Mtg)
        
    # --- Mgg
    Mgg = np.zeros((nU,nU))
    for i in range(nU):
        for j in range(nU):
            Mgg[i,j] = np.trapz((U[i][0,:]*U[j][0,:] + U[i][1,:]*U[j][1,:] + U[i][2,:]*U[j][2,:])*m,s_span)
    if bOrth:
        Mgg=Mgg*np.eye(nU)
    #print('Mgg\n',Mgg)

    # --- Build complete mass matrix
    MM = np.zeros((6+nU,6+nU))
    MM[:3,:3]   = Mxx; MM[:3,3:6] = Mxt; MM[:3,6:] = Mxg
    MM[3:6,3:6] = Mtt; MM[3:6,6:] = Mtg
    MM[6:,6:]   = Mgg

    i_lower     = np.tril_indices(len(MM), -1)
    MM[i_lower] = MM.T[i_lower]
    return MM
