"""
Reference:
     [1]: Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex
"""

import numpy as np
import unittest


# --------------------------------------------------------------------------------}
# --- Bodies 
# --------------------------------------------------------------------------------{
class body():
    def __repr__(self):
        pass
    pass

def fCreateBodyRigid(Name,Mass,J_G,rho_G):
    """
    Creates a rigid body 
    For now, based on a fake class
    """
    B = body()
    B.s_G_inB = rho_G
    B.J_G_inB = J_G
    B.Mass    = Mass
    B.J_O_inB = fTranslateInertiaMatrixFromCOG(B.J_G_inB, B.Mass, B.s_G_inB)
    B.KK=np.zeros((6,6))
    B.MM = fGMRigidBody(B.Mass,B.J_O_inB,B.s_G_inB)
    B.nf=0
    return B

def fGMRigidBody(Mass,J,rho):
    """ Generalized mass matrix for a rigid body (i.e. mass matrix) Eq.(15) of [1] """
    S=Mass*fSkew(rho)
    MM=np.zeros((6,6))
    MM[0:3,0:3] = Mass*np.eye(3);
    MM[0:3,3:6] = -S;
    MM[3:6,0:3] = S ; # transpose(S)=-S;
    MM[3:6,3:6] = J ;
    return MM




# --------------------------------------------------------------------------------}
# --- B Matrices 
# --------------------------------------------------------------------------------{
def fB_inB(R_EI, B_I):
    """ Transfer a global B_I matrix (body I at point I) into a matrix in it's own coordinate.
    Simply multiply the top part and bottom part of the B matrix by the 3x3 rotation matrix R_EI
    e.g.
         B_N_inN = [R_EN' * B_N(1:3,:);  R_EN' * B_N(4:6,:)];
    """ 
    if len(B_I)==0:
        B_I_inI = np.array([])
    else:
        B_I_inI = np.vstack(( np.dot(R_EI.T, B_I[:3,:]),  np.dot(R_EI.T, B_I[3:,:]) ))
    return B_I_inI

def fB_aug(B_I_inI, nf_I, nf_Curr=None, nf_Prev=None):
    """
    Augments the B_I_inI matrix, to include nf_I flexible degrees of freedom.
    This returns the full B matrix on the left side of Eq.(11) from [1], 
    based on the Bx and Bt matrices on the right side of this equation
    """
    if len(B_I_inI)==0:
        if nf_I>0:
            BB_I_inI = np.vstack( (np.zeros((6,nf_I)), np.eye(nf_I)))
        else:
            BB_I_inI= np.zeros((6,0))
    else:
        if nf_Curr is not None:
            # Case of several flexible bodies connected to one point (i.e. blades)
#             if nf_Curr==0:
#                 raise NotImplementedError()
            nf_After=nf_I-nf_Prev-nf_Curr
            I = np.block( [np.zeros((nf_Curr,nf_Prev)), np.eye(nf_Curr), np.zeros((nf_Curr,nf_After))] )
        else:
            nf_Curr=nf_I
            I=np.eye(nf_I)

        BB_I_inI = np.block([ [B_I_inI, np.zeros((6,nf_I))], [np.zeros((nf_Curr,B_I_inI.shape[1])), I]]);

    return BB_I_inI



def fBMatRecursion(Bp,Bhat,R0p,r_pi):
    """ Recursive formulae for B' and Bhat 
    See discussion after Eq.(12) and (15) from [1]
    """
    # --- Safety checks
    if len(Bp)==0:
        n_p = 0
    elif Bp.ndim==2:
        n_p = Bp.shape[1]
    else:
        raise Exception('Bp needs to be empty or a 2d array')
    if len(Bhat)==0:
        ni = 0
    elif Bhat.ndim==2:
        ni = Bhat.shape[1]
    else:
        raise Exception('Bi needs to be empty or a 2d array')

    r_pi=r_pi.ravel().reshape((3,1))

    # TODO use Translate here
    Bi = np.zeros((6,ni+n_p))
    for j in range(n_p):
        Bi[:3,j] = Bp[:3,j]+np.cross(Bp[3:,j],r_pi.ravel()) # Recursive formula for Bt mentioned after Eq.(15)
        Bi[3:,j] = Bp[3:,j] # Recursive formula for Bx mentioned after Eq.(12)
    if ni>0:
        Bi[:3,n_p:] = np.dot(R0p,Bhat[:3,:]) # Recursive formula for Bt mentioned after Eq.(15)
        Bi[3:,n_p:] = np.dot(R0p,Bhat[3:,:]) # Recursive formula for Bx mentioned after Eq.(12)
    return Bi

def fBMatTranslate(Bp,r_pi):
    """
    Rigid translation of a B matrix to another point, i.e. transfer the velocities from a point to another: 
      - translational velocity:  v@J = v@I + om@I x r@IJ
      - rotational velocity   : om@J = om@I
    """
    Bi=np.zeros(Bp.shape)
    if Bp.ndim==1:
        raise NotImplementedError

    for j in range(Bp.shape[1]):
        Bi[0:3,j] = Bp[0:3,j]+np.cross(Bp[3:6,j],r_pi.ravel());
        Bi[3:6,j] = Bp[3:6,j]
    return Bi


def fBMB(BB_I_inI,MM):
    """ Computes the body generalized matrix: B'^t M' B 
    See Eq.(8) of [1] 
    """
    MM_I = np.dot(np.transpose(BB_I_inI), MM).dot(BB_I_inI)
    return MM_I

# --------------------------------------------------------------------------------}
# --- Rotation 
# --------------------------------------------------------------------------------{
def fRotx(t):
    Ax = np.array([[1,0,0],[0,np.cos(t),-np.sin(t)],[0,np.sin(t),np.cos(t)]])
    return Ax

def fRoty(t):
    Ay = np.array([[np.cos(t),0,np.sin(t)],[0,1,0],[-np.sin(t),0,np.cos(t)]])
    return Ay

def fRotz(t):
    Az = np.array([[np.cos(t),-np.sin(t),0],[np.sin(t),np.cos(t),0],[0,0,1]])
    return Az

# --------------------------------------------------------------------------------}
# --- Inertia functions 
# --------------------------------------------------------------------------------{
def fTranslateInertiaMatrix(I_A, Mass, r_BG, r_AG = np.array([0,0,0])):
    """
    Transform inertia matrix with respect to point A to the inertia matrix with respect to point B
    NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system. 
    NOTE: one of the vector r_BG or r_AG may be empty or 0 instead of [0,0,0];
    NOTE: if r_AG is not provided it is assumed to be 0, i.e. A=G
    To avoid this confusion you can use fTranslateInertiaMatrixFromCOG  and fTranslateInertiaMatrixToCOG
    
    INPUTS:
       I_A  : Inertia matrix 3x3 in the coordinate system A
       Mass : Mass of the body
       r_BG: vector from point B to COG of the body
    
    OPTIONAL INPUTS:
       r_AG: vector from point A to point G
    """
    if len(r_AG) < 3:
        r_AG = np.array([0,0,0])
    if len(r_BG) < 3:
        r_BG = np.array([0,0,0])   
    return I_A - Mass*(np.dot(fSkew(r_BG), fSkew(r_BG))-np.dot(fSkew(r_AG),fSkew(r_AG)))

def fTranslateInertiaMatrixToCOG(I_B = None,Mass = None,r_GB = None): 
    """ Transform inertia matrix with respect to point B to the inertia matrix with respect to the COG
    NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system.
    
    INPUTS:
      I_G  : Inertia matrix 3x3 with respect to COG
      Mass : Mass of the body
      r_GB: vector from COG to point B
    """
    I_G = I_B + Mass * np.dot(fSkew(r_GB), fSkew(r_GB))
    return I_G

def fTranslateInertiaMatrixFromCOG(I_G = None,Mass = None,r_BG = None): 
    """
    Transform inertia matrix with respect to COG to the inertia matrix with respect to point B
    NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system.
    INPUTS:
       I_G  : Inertia matrix 3x3 with respect to COG
       Mass : Mass of the body
       r_BG: vector from point B to COG of the body
    """
    I_B = I_G - Mass * np.dot(fSkew(r_BG),fSkew(r_BG))
    return I_B
    
def fSkew(x):
    """ Returns the skew symmetric matrix M, such that: cross(x,v) = M v """
    return np.array([[0, -x[2], x[1]],[x[2],0,-x[0]],[-x[1],x[0],0]])


# --------------------------------------------------------------------------------}
# --- TEST  
# --------------------------------------------------------------------------------{
class TesT(unittest.TestCase):
    def test_rot(self):
        # --- Identity matrix for 0 rotation
        np.testing.assert_almost_equal(fRotx(0),np.eye(3))

    def test_skew(self):
        # --- Testing  \tilde{u} . v  ==  u x v
        u=np.array([1,2,3])
        v=np.array([-1,0,-3])
        np.testing.assert_equal(np.cross(u,v), np.dot(fSkew(u),v))

    def test_inertia(self):
        # --- Transferring inertia at G to point A and B and then to each other
        I_G=np.diag([1,2,3]); 
        M=2;
        r_OG = np.array([ 1, 2, 10 ])
        r_OA = r_OG + np.array([5, 8 , 2] )
        r_OB = r_OG + np.array([4, -6, -3])
        r_AG = r_OG-r_OA
        r_BG = r_OG-r_OB
        I_A  = fTranslateInertiaMatrix(I_G,M,r_AG)        # I@ A
        I_B  = fTranslateInertiaMatrix(I_G,M,r_BG)        # I@ B
        I_B2 = fTranslateInertiaMatrix(I_A,M,r_BG,r_AG   )# I@B from A
        I_A2 = fTranslateInertiaMatrix(I_B,M,r_AG,r_BG   )# I@A from B
        np.testing.assert_equal(I_A,I_A2)
        np.testing.assert_equal(I_B,I_B2)

        # --- Transfer of inertia at A to COG then back at A
        I_A = np.eye(3)
        M = 12
        r_GA = np.array([3,4,5])
        I_G  = fTranslateInertiaMatrixToCOG  (I_A,M,r_GA)  # I@G from A
        I_A2 = fTranslateInertiaMatrixFromCOG(I_G,M,-r_GA) # I@A from G
        np.testing.assert_equal(I_A,I_A2)


    def test_BMat(self):
        # --- Example taken from Main_TNSB from YAMS
        alpha_y  = -1.466329361065083e-02
        B_N_ref = np.array([ 0.0000E+00, 0.0000E+00, 1.0000E+00, 0.0000E+00, alpha_y, 0.0000E+00]).reshape((6,1))
        B_S_ref = np.array([
          [ 1.466171724553691e-01 ,  0.000000000000000e+00],
          [ 0.000000000000000e+00 ,  0.000000000000000e+00],
          [ 1.002150044745554e+00 ,  0.000000000000000e+00],
          [ 0.000000000000000e+00 , -1.466276815184685e-02],
          [-1.466329361065083e-02 ,  0.000000000000000e+00],
          [ 0.000000000000000e+00 ,  9.998924958364900e-01 ]])
        r_NS_ref=np.array([1.466276815184686e-01, 0.000000000000000e+00, -9.998924958364899e+00]).reshape((3,1))

        R_TN     = fRoty(alpha_y)        ;
        q_psi    = 1
        z_NS     = - 10
        r_NS_inN = np.array([0, 0, z_NS]).reshape((3,1))
        # link E-T
        R_ET     = np.eye(3)
        # ---------------------------------------------
        # Link T-N
        r_TN     = np.zeros((3,1))
        r_TN[0]  = 1.0000E+02
        Bx_TN    = np.zeros((3,1))
        Bt_TN    = np.zeros((3,1))
        Bx_TN[2] = 1
        Bt_TN[1] = alpha_y
        B_TN     = np.vstack((Bx_TN,Bt_TN))
        B_T      = np.array([])
        B_N      = fBMatRecursion(B_T,B_TN,R_ET,r_TN)
        np.testing.assert_equal(B_N,B_N_ref)
        R_TN=fRoty(alpha_y);
        R_EN=np.dot(R_ET,R_TN)
        # ---------------------------------------------
        # Link N-S
        R_NS = fRotz(q_psi+np.pi) # Adding pi here , blade down
        R_ES = np.dot(R_EN, R_NS   )
        r_NS = np.dot(R_EN, r_NS_inN )
        np.testing.assert_almost_equal(r_NS,r_NS_ref)

        Bx_NS=np.array([0,0,0]).reshape((3,1))
        Bt_NS=np.array([0,0,1]).reshape((3,1))
        B_NS =np.vstack((Bx_NS,Bt_NS))

        B_S = fBMatRecursion(B_N,B_NS,R_EN,r_NS);
        np.testing.assert_almost_equal(B_S,B_S_ref)




if __name__=='__main__':
    unittest.main()
