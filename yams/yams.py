"""
Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

import numpy as np
import unittest
from .flexibility import GMBeam, GKBeam, GKBeamStiffnening, polymode

# --------------------------------------------------------------------------------}
# --- Bodies 
# --------------------------------------------------------------------------------{
class Body(object):
    def __init__(B):
        B.MM     = None
        B.Name   = ''

    def updateKinematics(o,x_0,R_0b,gz,v_0,a_v_0):
        # Updating position of body origin in global coordinates
        o.r_O = x_0[0:3]
        o.gzf = gz
        # Updating Transformation matrix
        o.R_0b=R_0b
        # Updating rigid body velocity and acceleration
        o.v_O_inB     = np.dot(R_0b, v_0[0:3])
        o.om_O_inB    = np.dot(R_0b, v_0[3:6])
        o.a_O_v_inB   = np.dot(R_0b, a_v_0[0:3])
        o.omp_O_v_inB = np.dot(R_0b, a_v_0[3:6])
        
    def __repr__(B):
        pass
    
    @property
    def Mass(B):
        if B.MM is None:
            return 0
        return B.MM[0,0]

    @property
    def nf(B):
        if hasattr(B,'PhiU'):
            return len(B.PhiU)
        else:
            return 0

# --------------------------------------------------------------------------------}
# --- Rigid Body 
# --------------------------------------------------------------------------------{
class RigidBody(Body):
    def __init__(B,Name,Mass,J_G,rho_G):
        """
        Creates a rigid body 
        For now, based on a fake class
        """
        super(Body,B).__init__()
        B.s_G_inB = rho_G
        B.J_G_inB = J_G
        B.J_O_inB = fTranslateInertiaMatrixFromCOG(B.J_G_inB, Mass, B.s_G_inB)
        B.KK=np.zeros((6,6))
        B.MM = fGMRigidBody(Mass,B.J_O_inB,B.s_G_inB)


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
# --- Beam Body 
# --------------------------------------------------------------------------------{
class BeamBody(Body):
    def __init__(B, s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=None, s_G0=None, bAxialCorr=False, bOrth=False, Mtop=0, bStiffening=True, gravity=None):
        """ 
          Points P0 - Undeformed mean line of the body
        """
        super(BeamBody,B).__init__()
        B.s_span = s_span
        B.m      = m
        B.s_G0   = s_G0
        B.PhiU   = PhiU
        B.PhiV   = PhiV
        B.PhiK   = PhiK
        B.jxxG   = jxxG
        B.s_P0   = s_P0
        B.EI     = EI
        if jxxG is None:
            B.jxxG   = 0*m
        if B.s_G0 is None:
            B.s_G0=B.s_P0
    
        B.s_G    = B.s_G0
        B.bAxialCorr = bAxialCorr
        B.bOrth      = bOrth
        B.Mtop       = Mtop

        B.computeMassMatrix()
        B.KK = GKBeam(B.s_span, B.EI, B.PhiK, bOrth=B.bOrth)
        if bStiffening:
            pass
            #print('>>>>>>>>>>>>>> TODO TODO TODO Geometrical stiffnening')
            #KKg = GKBeamStiffnening(B.s_span, B.PhiV, gravity, m, Mtop)
            #print(KKg)
        B.DD = np.zeros((6+B.nf,6+B.nf))

        # TODO
        B.V0         = np.zeros((3,B.nSpan))
        B.K0         = np.zeros((3,B.nSpan))
        B.rho_G0_inS = np.zeros((3,B.nSpan)) # location of COG in each cross section
        #[o.PhiV,o.PhiK] = fBeamSlopeCurvature(o.s_span,o.PhiU,o.PhiV,o.PhiK,1e-2);
        #[o.V0,o.K0]     = fBeamSlopeCurvature(o.s_span,o.s_P0,o.V0,o.K0,1e-2)    ;
        #if isempty(o.s_G0); o.s_G0=o.s_P0; end;
        #if isempty(o.rho_G0_inS); o.rho_G0_inS=np.zeros(3,o.nSpan); end;
        #if isempty(o.rho_G0    ); 
        #    o.rho_G0 =np.zeros(3,o.nSpan);
        #    for i=1:o.nSpan
        #        o.rho_G0(1:3,i) =fRotx(o.V0(1,i))*o.rho_G0_inS(:,i);

    def computeMassMatrix(B):
        B.MM = GMBeam(B.s_G, B.s_span, B.m, B.PhiU, jxxG=B.jxxG, bUseIW=True, main_axis='x', bAxialCorr=B.bAxialCorr, bOrth=B.bOrth)

    def updateKinematics(o,x_0,R_0b,gz,v_0,a_v_0):
        super(BeamBody,o).updateKinematics(x_0,R_0b,gz,v_0,a_v_0)
        # --- Calculation of deformations wrt straight beam axis, curvature (K) and velocities (UP)
        if o.nf>0:
            o.gzpf  = v_0[6:]
            o.gzppf = a_v_0[6:]
            # Deflections shape
            o.U  = np.zeros((3,o.nSpan));
            o.V  = np.zeros((3,o.nSpan));
            o.K  = np.zeros((3,o.nSpan));
            #o.U(1,:) = o.s_span; 
            o.UP = np.zeros((3,o.nSpan));
            for j in range(o.nf):
                o.U [0:3,:] = o.U [0:3,:] + o.gzf[j]  * o.PhiU[j][0:3,:]
                o.UP[0:3,:] = o.UP[0:3,:] + o.gzpf[j] * o.PhiU[j][0:3,:]
                o.V [0:3,:] = o.V [0:3,:] + o.gzf[j]  * o.PhiV[j][0:3,:]
                o.K [0:3,:] = o.K [0:3,:] + o.gzf[j]  * o.PhiK[j][0:3,:]
            o.V_tot=o.V+o.V0;
            o.K_tot=o.K+o.K0;

            # Position of mean line
            o.s_P=o.s_P0+o.U;

            # Position of deflected COG
            # TODO TODO TODO mean_axis not x
            o.rho_G      = np.zeros((3,o.nSpan))
            o.rho_G[1,:] = o.rho_G0_inS[1,:]*np.cos(o.V_tot[0,:])-o.rho_G0_inS[2,:]*np.sin(o.V_tot[0,:]);
            o.rho_G[2,:] = o.rho_G0_inS[1,:]*np.sin(o.V_tot[0,:])+o.rho_G0_inS[2,:]*np.cos(o.V_tot[0,:]);
            o.s_G = o.s_P+o.rho_G;
            # Alternative:
            #rho_G2     = zeros(3,o.nSpan);
            #rho_G2(2,:) = o.rho_G0(2,:).*cos(o.V(1,:))-o.rho_G0(3,:).*sin(o.V(1,:));
            #rho_G2(3,:) = o.rho_G0(2,:).*sin(o.V(1,:))+o.rho_G0(3,:).*cos(o.V(1,:));
            #compare(o.rho_G,rho_G2,'rho_G');
            # Position of connection point
            print('TODO connection points')
            #for ic=1:length(o.Connections)
            #    iNode=o.Connections{ic}.ParentNode;
            #    %o.Connections{ic}.s_C_inB = o.U(1:3,iNode);
            #    o.Connections{ic}.s_C_inB = o.s_P(1:3,iNode);

    @property
    def nSpan(B):
        return len(B.s_span)


# --------------------------------------------------------------------------------}
# --- Uniform Beam Body 
# --------------------------------------------------------------------------------{
class UniformBeamBody(BeamBody):
    def __init__(B, Name, nShapes, nSpan, L, EI0, m, Mtop=0, jxxG=None, GKt=None, bAxialCorr=True, bCompatibility=False, bStiffnessFromGM=False):

        import beams.theory as bt
        if jxxG is None:
            jxxG=0
        if GKt is None:
            GKt=0

        A=1; rho=A*m;
        x=np.linspace(0,L,nSpan);
        # Mode shapes
        freq,s_span,U,V,K = bt.UniformBeamBendingModes('unloaded-topmass-clamped-free',EI0,rho,A,L,x=x,Mtop=Mtop)
        PhiU = np.zeros((nShapes,3,nSpan)) # Shape
        PhiV = np.zeros((nShapes,3,nSpan)) # Slope
        PhiK = np.zeros((nShapes,3,nSpan)) # Curvature
        for j in np.arange(nShapes):  
            PhiU[j][2,:] = U[j,:] # Setting modes along z
            PhiV[j][2,:] = V[j,:]
            PhiK[j][2,:] = K[j,:]
        m       = m    * np.ones(nSpan)
        jxxG    = jxxG * np.ones(nSpan)
        EI      = np.zeros((3,nSpan))
        EI[1,:] = EI0
        EI[2,:] = EI0
        GKt     = GKt  * np.ones(nSpan)
        
        # --- Straight undeflected shape (and COG)
        s_P0      = np.zeros((3,nSpan))
        s_P0[0,:] = x

	# Create a beam body
        super(UniformBeamBody,B).__init__(s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=jxxG, bAxialCorr=bAxialCorr, Mtop=Mtop)


# --------------------------------------------------------------------------------}
# --- FAST Beam body 
# --------------------------------------------------------------------------------{
class FASTBeamBody(BeamBody):
    def __init__(B,body_type,ED,inp,Mtop,nShapes=2,main_axis='x',nSpan=40,bAxialCorr=False):
        """ 
        INPUTS:
           nSpan: number of spanwise station used (interpolated from input)
                  Use -1 or None to use number of stations from input file
        """
        # --- Reading coefficients
        exp = np.arange(2,7)
        if body_type.lower()=='blade':
            coeff = np.array([[ inp['BldFl1Sh(2)'], inp['BldFl2Sh(2)'], inp['BldEdgSh(2)']],
                              [ inp['BldFl1Sh(3)'], inp['BldFl2Sh(3)'], inp['BldEdgSh(3)']],
                              [ inp['BldFl1Sh(4)'], inp['BldFl2Sh(4)'], inp['BldEdgSh(4)']],
                              [ inp['BldFl1Sh(5)'], inp['BldFl2Sh(5)'], inp['BldEdgSh(5)']],
                              [ inp['BldFl1Sh(6)'], inp['BldFl2Sh(6)'], inp['BldEdgSh(6)']]])

            damp_zeta  = np.array([ inp['BldFlDmp(1)'], inp['BldFlDmp(2)'], inp['BldEdDmp(1)']])/100
            mass_fact = inp['AdjBlMs']                                              # Factor to adjust blade mass density (-)
          
            prop     = inp['BldProp']
            span_max = ED['TipRad']   # TODO TODO do somthing about hub rad
            s_bar, m, EIFlp, EIEdg  =prop[:,0], prop[:,3], prop[:,4], prop[:,5]

        elif body_type.lower()=='tower':
            coeff = np.array([[ inp['TwFAM1Sh(2)'], inp['TwFAM2Sh(2)'], inp['TwSSM1Sh(2)'], inp['TwSSM2Sh(2)']],
                              [ inp['TwFAM1Sh(3)'], inp['TwFAM2Sh(3)'], inp['TwSSM1Sh(3)'], inp['TwSSM2Sh(3)']],
                              [ inp['TwFAM1Sh(4)'], inp['TwFAM2Sh(4)'], inp['TwSSM1Sh(4)'], inp['TwSSM2Sh(4)']],
                              [ inp['TwFAM1Sh(5)'], inp['TwFAM2Sh(5)'], inp['TwSSM1Sh(5)'], inp['TwSSM2Sh(5)']],
                              [ inp['TwFAM1Sh(6)'], inp['TwFAM2Sh(6)'], inp['TwSSM1Sh(6)'], inp['TwSSM2Sh(6)']]])
            damp_zeta = np.array([inp['TwrFADmp(1)'], inp['TwrFADmp(2)'], inp['TwrSSDmp(1)'], inp['TwrSSDmp(2)']])/100 # structural damping ratio 
            mass_fact = inp['AdjTwMa']                                              # Factor to adjust tower mass density (-)

            prop     = inp['TowProp']
            span_max = ED['TowerHt']-ED['TowerBsHt']
            s_bar, m, EIFlp, EIEdg  = prop[:,0], prop[:,1], prop[:,2], prop[:,3]


        else:
            raise Exception('Body type not supported {}'.format(body_type))
        nShpMax=coeff.shape[1]
        if nShapes>nShpMax:
            raise Exception('A maximum of {} shapes function possible with FAST {} body'.format(nShpMax,body_type))

        gravity=ED['Gravity']

        # --- Interpolating structural properties
        m *= mass_fact
        if nSpan is None or nSpan<0:
            nSpan  = len(prop[:,0])
        s_span     = np.linspace(0,span_max,nSpan)
        m          = np.interp(s_span/span_max,s_bar,m)
        EIFlp      = np.interp(s_span/span_max,s_bar,EIFlp)
        EIEdg      = np.interp(s_span/span_max,s_bar,EIEdg)

        # --- Definition of main directions
        ShapeDir=np.zeros(nShpMax).astype(int)
        DirNames=['x','y','z']
        EI =np.zeros((3,nSpan))
        if main_axis=='x':
            iMain = 0  # longitudinal axis along x
            ShapeDir[0:2] = 2 # First two shapes are along z (flapwise/Fore-Aft)
            ShapeDir[2:]  = 1 # Third shape along along y (edgewise/Side-Side) Sign...
            EI[2,:] = EIFlp
            EI[1,:] = EIEdg
        elif main_axis=='z':
            iMain = 2  # longitudinal axis along z
            ShapeDir[0:2] = 0 # First two shapes are along x (flapwise/Fore-Aft)
            ShapeDir[2:]  = 1 # Third shape along along y (edgewise/Side-Side) Sign...
            EI[0,:] = EIFlp
            EI[1,:] = EIEdg
        else:
            raise NotImplementedError()

        # --- Undeflected shape
        s_P0          = np.zeros((3,nSpan))
        s_P0[iMain,:] = s_span 
        # TODO blade COG
        jxxG=m*0

        # --- Shape functions
        PhiU = np.zeros((nShapes,3,nSpan)) # Shape
        PhiV = np.zeros((nShapes,3,nSpan)) # Slope
        PhiK = np.zeros((nShapes,3,nSpan)) # Curvature
        for j in np.arange(nShapes):
            iAxis = ShapeDir[j]
            PhiU[j][iAxis,:], PhiV[j][iAxis,:], PhiK[j][iAxis,:] = polymode(s_span,coeff[:,j],exp)


        super(FASTBeamBody,B).__init__(s_span, s_P0, m, PhiU, PhiV, PhiK, EI, jxxG=jxxG, bAxialCorr=bAxialCorr, bOrth=body_type=='blade', gravity=gravity,Mtop=Mtop)

        # Damping
        B.DD=np.zeros((6+nShapes,6+nShapes))

        # Using diagonal damping
        for j in range(nShapes):
            m             = B.MM[6+j,6+j]
            k             = B.KK[6+j,6+j]
            om            = np.sqrt(k/m)
            xi            = damp_zeta[j]*2*np.pi
            c             = xi * m * om / np.pi
            B.DD[6+j,6+j] = c


        # --- Storing data in object
        B.main_axis = main_axis # TODO







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
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):
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
