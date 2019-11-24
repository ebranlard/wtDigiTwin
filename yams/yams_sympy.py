"""
Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""
import numpy as np
import sympy
from sympy import Symbol
from sympy import Matrix
from sympy.printing import lambdarepr
from sympy import init_printing
from sympy import lambdify
from sympy.abc import *
from sympy import trigsimp
from sympy import cos,sin
display=lambda x: sympy.pprint(x, use_unicode=False,wrap_line=False)

def colvec(v): 
    return Matrix([[v[0]],[v[1]],[v[2]]])

def cross(V1,V2):
    return [V1[1]*V2[2]-V1[2]*V2[1], V1[2]*V2[0]-V1[0]*V2[2], (V1[0]*V2[1]-V1[1]*V2[0]) ]
def eye(n): 
    return Matrix( np.eye(n).astype(int) )



# --------------------------------------------------------------------------------}
# --- Connections 
# --------------------------------------------------------------------------------{
class Connection():
    def __init__(self,Type,RelPoint=None,RelOrientation=None,JointRotations=None):
        if RelOrientation is None:
            RelOrientation=eye(3)
        if RelPoint is None:
            RelPoint=colvec([0,0,0])

        self.Type=Type
        
        self.s_C_0_inB = RelPoint
        self.s_C_inB   = self.s_C_0_inB
        self.R_ci_0    = RelOrientation
        self.R_ci      = self.R_ci_0     

        if self.Type=='Rigid':
            self.nj=0
        elif self.Type=='SphericalJoint':
            self.JointRotations=JointRotations;
            self.nj=len(self.JointRotations);
        else:
            raise NotImplementedError()

    def updateKinematics(j,q):
        j.B_ci=Matrix(np.zeros((6,j.nj)))
        if j.Type=='Rigid':
            j.R_ci=j.R_ci_0
        elif j.Type=='SphericalJoint':
            R=eye(3)
            myq    = q   [j.I_DOF,0];
            #myqdot = qdot[j.I_DOF];

            for ir,rot in enumerate(j.JointRotations):
                if rot=='x':
                    I=np.array([1,0,0])
                    Rj=R_x( myq[ir] )
                elif rot=='y':
                    I=np.array([0,1,0])
                    Rj=R_y( myq[ir] )
                elif rot=='z':
                    I=np.array([0,0,1])
                    Rj=R_z( myq[ir] )
                else:
                    raise Exception()
                # Setting Bhat column by column
                j.B_ci[3:,ir] = np.dot(R,I) # NOTE: needs to be done before R updates
                # Updating rotation matrix
                R      = np.dot(R , Rj )
                j.R_ci = Matrix(np.dot(R, j.R_ci_0 ))


# --------------------------------------------------------------------------------}
# --- Bodies 
# --------------------------------------------------------------------------------{
class Body(object):
    def __init__(B,Name=''):
        B.Children    = []
        B.Connections = []
        B.Name        = Name
        B.MM     = None
        B.B           = [] # Velocity transformation matrix
        B.updatePosOrientation(colvec([0,0,0]), eye(3))

    def updatePosOrientation(o,x_0,R_0b):
        o.r_O = x_0      # position of body origin in global coordinates
        o.R_0b=R_0b      # transformation matrix from body to global

    def connectTo(self, Child, Point=None, Type=None, RelOrientation=None, JointRotations=None):
        if Type =='Rigid':
            c=Connection(Type, RelPoint=Point, RelOrientation = RelOrientation)
        else: # TODO first node, last node
            c=Connection(Type, RelPoint=Point, RelOrientation=RelOrientation, JointRotations=JointRotations)
        self.Children.append(Child)
        self.Connections.append(c)

    def setupDOFIndex(o,n):
        nForMe=o.nf
        # Setting my dof index
        o.I_DOF=n+ np.arange(nForMe) 
        # Update
        n=n+nForMe
        for child,conn in zip(o.Children,o.Connections):
            # Connection first
            nForConn=conn.nj;
            conn.I_DOF=n+np.arange(nForConn)
            # Update
            n=n+nForConn;
            # Then Children
            n=child.setupDOFIndex(n)
        return n

    def __repr__(B):
        pass

    @property
    def R_bc(self):
        return eye(3);
    @property
    def Bhat_x_bc(self):
        return Matrix(np.zeros((3,0)))
    @property
    def Bhat_t_bc(self):
        return Matrix(np.zeros((3,0)))

    def updateChildrenKinematicsNonRecursive(p,q):
        # At this stage all the kinematics of the body p are known
        # Useful variables
        R_0p =  p.R_0b
        B_p  =  p.B
        r_0p  = p.r_O

        nf_all_children=sum([child.nf for child in p.Children])

        for ic,(body_i,conn_pi) in enumerate(zip(p.Children,p.Connections)):
            # Flexible influence to connection point
            R_pc  = p.R_bc
            Bx_pc = p.Bhat_x_bc
            Bt_pc = p.Bhat_t_bc
            # Joint influence to next body (R_ci, B_ci)
            conn_pi.updateKinematics(q) # TODO

            # Full connection p and j
            R_pi   = R_pc*conn_pi.R_ci  
            if conn_pi.B_ci.shape[1]>0:
                Bx_pi  = Matrix(np.column_stack((Bx_pc, np.dot(R_pc,conn_pi.B_ci[:3,:]))))
                Bt_pi  = Matrix(np.column_stack((Bt_pc, np.dot(R_pc,conn_pi.B_ci[3:,:]))))
            else:
                Bx_pi  = Bx_pc
                Bt_pi  = Bt_pc
              
            # Rotation of body i is rotation due to p and j
            R_0i = R_0p * R_pi

            # Position of connection point in P and 0 system
            r_pi_inP= conn_pi.s_C_inB
            r_pi    = R_0p * r_pi_inP
            B_i      = fBMatRecursion(B_p, Bx_pi, Bt_pi, R_0p, r_pi)
            B_i_inI  = fB_inB(R_0i, B_i)
            BB_i_inI = fB_aug(B_i_inI, body_i.nf)

            body_i.B      = B_i    
            body_i.B_inB  = B_i_inI
            body_i.BB_inB = BB_i_inI

            # --- Updating Position and orientation of child body 
            r_0i = r_0p + r_pi  # % in 0 system
            body_i.R_pb = R_pi 
            body_i.updatePosOrientation(r_0i,R_0i)

            # TODO flexible dofs and velocities/acceleration
            body_i.gzf  = q[body_i.I_DOF,0] # TODO use updateKinematics

    def getFullM(o,M):
        if not isinstance(o,GroundBody):
            MqB      = fBMB(o.BB_inB,o.MM)
            n        = MqB.shape[0]
            M[:n,:n] = M[:n,:n]+MqB     
        for c in o.Children:
            M=c.getFullM(M)
        return M
        
    def getFullK(o,K):
        if not isinstance(o,GroundBody):
            KqB      = fBMB(o.BB_inB,o.KK)
            n        = KqB.shape[0]
            K[:n,:n] = K[:n,:n]+KqB     
        for c in o.Children:
            K=c.getFullK(K)
        return K
        
    def getFullD(o,D):
        if not isinstance(o,GroundBody):
            DqB      = fBMB(o.BB_inB,o.DD)
            n        = DqB.shape[0]
            D[:n,:n] = D[:n,:n]+DqB     
        for c in o.Children:
            D=c.getFullD(D)
        return D


# --------------------------------------------------------------------------------}
# --- Ground Body 
# --------------------------------------------------------------------------------{
class GroundBody(Body):
    def __init__(B):
        super(GroundBody,B).__init__('Grd')
        B.nf   = 0

# --------------------------------------------------------------------------------}
# --- Rigid Body 
# --------------------------------------------------------------------------------{
class RigidBody(Body):
    def __init__(B, Name, Mass, J_G, rho_G):
        """
        Creates a rigid body 
        """
        super(RigidBody,B).__init__(Name)
        B.nf  = 0
        B.s_G_inB = rho_G
        B.J_G_inB = J_G  
        B.Mass    = Mass 

# --------------------------------------------------------------------------------}
# --- Beam Body 
# --------------------------------------------------------------------------------{
class BeamBody(Body):
    def __init__(B,Name,nf,main_axis='z',nD=2):
        super(BeamBody,B).__init__(Name)
        B.nf  = nf
        B.nD  = nD
        B.main_axis = main_axis
    @property
    def alpha_couplings(self):
        return  np.dot(self.Bhat_t_bc , self.gzf)

    @property
    def R_bc(self):
        if self.main_axis=='x':
            alpha_y= symbols('alpha_y') #-p.V(3,iNode);
            alpha_z= symbols('alpha_z') # p.V(2,iNode);
            return R_y(alpha_y)*R_z(alpha_z)

        elif self.main_axis=='z':
            alpha_x= symbols('alpha_x') #-p.V(2,iNode);
            alpha_y= symbols('alpha_y') # p.V(1,iNode);
            return R_x(alpha_x)*R_y(alpha_y)
        else:
            raise NotImplementedError()

    @property
    def Bhat_x_bc(self):
        #      Bx_pc(:,j)=p.PhiU{j}(:,iNode);
        Bhat_x_bc = Matrix(np.zeros((3,self.nf)).astype(int))
        if self.main_axis=='z':
            for j in np.arange(self.nf):
                if j<self.nf/2 or self.nD==1:
                    Bhat_x_bc[0,j]=symbols('ux{:d}c'.format(j+1)) # p.PhiU{j}(:,iNode);  along x
                else:
                    Bhat_x_bc[1,j]=symbols('uy{:d}c'.format(j+1)) # p.PhiU{j}(:,iNode);  along y
        elif self.main_axis=='x':
            for j in np.arange(self.nf):
                if j<self.nf/2 or self.nD==1:
                    Bhat_x_bc[2,j]=symbols('uz{:d}c'.format(j+1)) # p.PhiU{j}(:,iNode);  along z
                else:
                    Bhat_x_bc[1,j]=symbols('uy{:d}c'.format(j+1)) # p.PhiU{j}(:,iNode);  along y
        return Bhat_x_bc

    @property
    def Bhat_t_bc(self):
        #      Bt_pc(:,j)=[0; -p.PhiV{j}(3,iNode); p.PhiV{j}(2,iNode)];
        Bhat_t_bc = Matrix(np.zeros((3,self.nf)).astype(int))
        if self.main_axis=='z':
            for j in np.arange(self.nf):
                if j<self.nf/2  or self.nD==1:
                    Bhat_t_bc[1,j]=symbols('vy{:d}c'.format(j+1))
                else:
                    Bhat_t_bc[0,j]=-symbols('vx{:d}c'.format(j+1))
        elif self.main_axis=='x':
            for j in np.arange(self.nf):
                if j<self.nf/2 or self.nD==1:
                    Bhat_t_bc[1,j]=-symbols('vz{:d}c'.format(j+1))
                else:
                    Bhat_t_bc[2,j]=symbols('vy{:d}c'.format(j+1))
        return Bhat_t_bc




# --------------------------------------------------------------------------------}
# --- Rotation 
# --------------------------------------------------------------------------------{
def R_x(t):
    return Matrix( [[1,0,0], [0,cos(t),-sin(t)], [0,sin(t),cos(t)]])
def R_y(t):
    return Matrix( [[cos(t),0,sin(t)], [0,1,0], [-sin(t),0,cos(t)] ])
def R_z(t): 
    return Matrix( [[cos(t),-sin(t),0], [sin(t),cos(t),0], [0,0,1]])
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
        B_I_inI = Matrix(np.array([]))
    else:
        B_I_inI = Matrix(np.vstack(( R_EI.T* B_I[:3,:],  R_EI.T * B_I[3:,:])))
    return B_I_inI

def fB_aug(B_I_inI, nf_I, nf_Curr=None, nf_Prev=None):
    """
    Augments the B_I_inI matrix, to include nf_I flexible degrees of freedom.
    This returns the full B matrix on the left side of Eq.(11) from [1], 
    based on the Bx and Bt matrices on the right side of this equation
    """
    if len(B_I_inI)==0:
        if nf_I>0:
            BB_I_inI = Matrix(np.vstack( (np.zeros((6,nf_I)).astype(int), np.eye(nf_I).astype(int))) )
        else:
            BB_I_inI= Matrix(np.zeros((6,0)).astype(int))
    else:
        if nf_Curr is not None:
            # Case of several flexible bodies connected to one point (i.e. blades)
            nf_After=nf_I-nf_Prev-nf_Curr
            I = np.block( [np.zeros((nf_Curr,nf_Prev)), np.eye(nf_Curr), np.zeros((nf_Curr,nf_After))] )
        else:
            nf_Curr=nf_I
            I=np.eye(nf_I)

        BB_I_inI = np.block([ [B_I_inI, np.zeros((6,nf_I))], [np.zeros((nf_Curr,B_I_inI.shape[1])), I]]);

    return Matrix(BB_I_inI)


def fBMatRecursion(Bp, Bhat_x, Bhat_t, R0p, r_pi):
    """ Recursive formulae for B' and Bhat 
    See discussion after Eq.(12) and (15) from [1]
    """
    # --- Safety checks
    if len(Bp)==0:
        n_p = 0
    elif len(Bp.shape)==2:
        n_p = Bp.shape[1]
    else:
        raise Exception('Bp needs to be empty or a 2d array')
    if len(Bhat_x)==0:
        ni = 0
    elif len(Bhat_x.shape)==2:
        ni = Bhat_x.shape[1]
    else:
        raise Exception('Bi needs to be empty or a 2d array')

    r_pi=colvec(r_pi)

    # TODO use Translate here
    Bi = Matrix(np.zeros((6,ni+n_p)))
    for j in range(n_p):
        Bi[:3,j] = Bp[:3,j]+cross(Bp[3:,j],r_pi) # Recursive formula for Bt mentioned after Eq.(15)
        Bi[3:,j] = Bp[3:,j] # Recursive formula for Bx mentioned after Eq.(12)
    if ni>0:
        Bi[:3,n_p:] = R0p*Bhat_x[:,:] # Recursive formula for Bx mentioned after Eq.(15)
        Bi[3:,n_p:] = R0p*Bhat_t[:,:] # Recursive formula for Bt mentioned after Eq.(12)
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

