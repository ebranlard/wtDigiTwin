import numpy as np
import unittest
from scipy.linalg import expm

def fEstimateKFTimeStep(u1,y1,z0,A,B,G,J,P0,Q,R): 
    """ Performs one time step of Kalman filter estimation 

    INPUTS:
     u1: inputs at time n
     y1: measurements at time n
     z0: Kalman state estimate at time n-1

     Equations number are compared to the following reference:
     [1] Lourens"""
        
    # estimate next step
    z1m   = A.dot(z0)  + B.dot(u1)
    y1hat = G.dot(z1m) + J.dot(u1)
    P1m   = (A.dot(P0)).dot(A.T) + Q
    
    # Calculate Kalman gain
    # same as Lk from [1] - And their Rtilde_k is G*P1m*G'+R
    Kk =  np.dot(P1m,G.T).dot( np.linalg.inv(((G.dot(P1m)).dot(G.T) + R)))   # TODO check me
    # update estimate with measurement
    z1 = z1m + Kk.dot(y1 - y1hat)
    
    P1 = (np.eye(A.shape[0]) - Kk.dot(G) ).dot(P1m)
    return z1,P1,Kk


def fKFDiscretize(Ac,Bc,dt,method='exponential'):
    """ Discretize the continuous states matrices Ac, Bc
    
    "Real" system:
        zdot = Ac.z + Bc.u + wd
          y  = Gc.z + Jc.u + wn
    "Discretized" system:
        z_{k+1} = Ad.z_k + Bd.u_k + wd
          y_k   = Gd.z_k + Jd.u_k + wn 

    INPUTS:
        methods: 'exponential', 'eigenvalues', 'forward_euler'

    OUTPUTS:
       Ad,Bd
    """
    # --- A matrix
    if method=='exponential':
        # Using matrix exponential directly
        Ad = expm(Ac * dt)
        if np.linalg.det(Ac) == 0:
            # TODO TODO missing square of dt for A?
            print('[WARN] Matrix A is singular, using forward euler to discretize B matrix\n' % ())
            # Forward euler
            Bd = dt * Bc
        else:
            mA_B = np.linalg.solve(Ac,Bc)
            Bd = np.dot( (Ac - np.eye(Ac.shape[0])) ,  mA_B)
    elif method=='eigenvalues':
        raise NotImplementedError()
        # Using eigenvalues
        #Q,Lambda = eig(Ac)
        #Ad = real(Q * expm(Lambda * dt) / Q)
        #Bd = Bc * dt
    elif method=='forward_euler':
        # Using forward Euler
        Ad = np.eye(Ac.shape[0]) + Ac * dt # TODO TODO missing square of dt?
        Bd = Bc * dt
    else:
        raise Exception('Unknown discretization method: %s',method)
    
    # --- G,J matrices
    #Gd = Gc
    #Jd = Jc
    return Ad,Bd

def fBuildSystem_Linear(M,C,K,Sa,Sv,Sd,Sp=None,Rp=None,Yp=None,Method='default'):
    """ Takes system matrices of a mechanical system, returns a state matrix.
    The state matrix may be an "augmented matrix", in which case Sp, Rp, should be provided

    - Mechanical equation:
       M qddot + Cqdot + Kq = f
                                      (f = Sp.p, for augmented)
    - Output equation:
       y = Sa.qddot + Sv.qdot + Sd.q 
                                      (+ Yp.p, for some augmented system)
    - (Augmented load evolution:
         pdot = Rp.p , only for augmented system)

    State Equation
        zdot = Ac.z + Bc.u + wd
    Measurement Equation
          y  = Gc.z + Jc.u + wn
    """
    nDOF = M.shape[0]
    nY   = Sd.shape[0]
    if 'default' == Method:
        Z=np.zeros((nDOF,nDOF))
        I=np.eye(nDOF)
        Ac = np.block( [ [Z , I ], [ mM_K, mM_C] ])
        Bc = 0
        Gc = np.block( [ Sd + np.dot(Sa,mM_K),  Sv + np.dot(Sa, mM_C) ] )
        Jc = 0
    else:
        if 'augmented_first_order' == Method:
            # Needs Sp and Rp to be defined!
            if Sp is None or Rp is None:
                raise Exception('Both Sp and Rp needs to be set with augmented first order method')
            nP = Sp.shape[1]
            if Yp is None:
                Yp=np.zeros((nY,nP))

            Z    = np.zeros((nDOF,nDOF))
            Znnp = np.zeros((nDOF,nP  ))
            Znpn = np.zeros((nP  ,nDOF))
            I    = np.eye(nDOF)
            mM_K = np.linalg.solve(-M,K)
            mM_C = np.linalg.solve(-M,C)
            M_Sp  = np.linalg.solve(M,Sp)
            Ac = np.block( [ [Z, I ,Znnp] , [mM_K, mM_C, M_Sp], [Znpn, Znpn, Rp] ])
            Bc = 0
            Gc = np.block( [Sd + np.dot(Sa,mM_K), Sv + np.dot(Sa,mM_C), Yp+np.dot(Sa,M_Sp) ])
#             print('Sd..:\n', Sd + np.dot(Sa,mM_K))
#             print('Sv..:\n', Sv + np.dot(Sa,mM_C))
#             print('Sp..:\n', Yp+np.dot(Sa,M_Sp) )
            Jc = 0
        else:
            raise Exception('Method %s not implemented')
    
    return Ac,Bc,Gc,Jc



def EmptyStateMat(nX,nU,nY):
    """ Returns state matrices with proper dimensions, filled with 0 """
    Ac = np.zeros((nX,nX))
    Gc = np.zeros((nY,nX))
    Bc = np.zeros((nX,nU))
    Jc = np.zeros((nY,nU))
    return Ac,Bc,Gc,Jc

def EmptySystemMat(nDOF_2nd, nY, nP=None):
    """ Returns matrices with proper dimensions, filled with 0
    INPUTS:
       - nDOF_2nd: number of "mechanical" degrees of freedoms, when the equations are written 
                   as a secnod order system, e.g. M qddot = F, then nDOF_2nd is the size of qddot
       - nY: Number of outputs
       - nP: Number of extended loads, if states are to be augmented
    NOTE:
        full augmented state vector has size 2*nDOF_2nd + nP
    """
    M=np.zeros((nDOF_2nd,nDOF_2nd))
    C=np.zeros((nDOF_2nd,nDOF_2nd))
    K=np.zeros((nDOF_2nd,nDOF_2nd))
    Sa = np.zeros((nY,nDOF_2nd))
    Sv = np.zeros((nY,nDOF_2nd))
    Sd = np.zeros((nY,nDOF_2nd))
    if nP is not None:
        Yp = np.zeros((nY,nP))
        Sp = np.zeros((nDOF_2nd,nP))
        Rp = np.zeros((nP,nP))
        return M,C,K,Sa,Sv,Sd,Sp,Rp,Yp
    else:
        return M,C,K,Sa,Sv,Sd




# --------------------------------------------------------------------------------}
# --- TESTS
# --------------------------------------------------------------------------------{
class Test(unittest.TestCase):

    def test_discretize_exp(self):
        # Discretize a continuous system using the exponential matrix method
        nX=3
        nU=1
        nY=2
        Ac, Bc, Gc, Jc = EmptyStateMat(nX,nU,nY)
        Ac[0,0]=1
        Ac[0,1]=4
        Ac[0,2]=4
        Ac[1,1]=2
        Ac[1,2]=5
        Ac[2,2]=3
        Bc[0,0]=1/2
        Bc[2,0]=1
        dt=0.5
        Ad,Bd = fKFDiscretize(Ac,Bc,dt,method='exponential')
        #print('Ac\n',Ad)
        #print('Bc\n',Bd)
        Ad_ref=np.array([
               [1.64872, 4.27824,12.60440],
               [0.00000, 2.71828, 8.81704],
               [0.00000, 0.00000, 4.48169]])
        Bd_ref = np.array([ [-2.00000], [ 0.83333], [ 0.66667]])
        np.testing.assert_almost_equal(Ad,Ad_ref,decimal=5)
        np.testing.assert_almost_equal(Bd,Bd_ref,decimal=5)

    def test_discretize_forward(self):
        # Discretize a continuous system using the forward euler method
        nX=3
        nU=1
        nY=2
        Ac, Bc, Gc, Jc = EmptyStateMat(nX,nU,nY)
        Ac[0,0]=0
        Ac[0,1]=1
        Ac[0,2]=0
        Ac[1,1]=0
        Ac[1,2]=1
        Ac[2,2]=0
        Bc[0,0]=1/2
        Bc[2,0]=1
        dt=0.5
        Ad_xp,Bd_xp = fKFDiscretize(Ac,Bc,dt,method = 'exponential'  )
        Ad_fe,Bd_fe = fKFDiscretize(Ac,Bc,dt,method = 'forward_euler')
        np.testing.assert_almost_equal(Bd_fe,Bd_xp,decimal=5)

        # TODO TODO, what to do of A
        #np.testing.assert_almost_equal(Ad_fe,Ad_xp,decimal=5)



    def test_build_shaftonly(self):
        # Build continuous matrices for the "Shaft only" case
        np.set_printoptions(linewidth=500)
        nDOF_2nd = 1           # Mech DOFs     :  q  = [psi]
        nDOF_1st = 2*nDOF_2nd  # Mech state    :  x  = [psi,psi_dot]
        nP       = 2           # Extended loads:  p  = [Ma ,Mg]
        nDOF_ext = nDOF_1st+nP # Ext. state    : x~  = [psi,psi_dot ,Ma,Mg]
        nY       = 2           # Outputs       :  y  = [psi_dot,Mg]
        nU       = 0           # Inputs
        J_LSS = 2.0
        # General init
        M,C,K,Sa,Sv,Sd,Sp,Rp,Yp = EmptySystemMat(nDOF_2nd, nY, nP)
        # Setting values
        M[0,0]  = J_LSS
        Sv[0,0] = 1
        Sp[0,0] = 1
        Sp[0,1] = -1
        Yp[1,1] = 1 # Direct feed through of force

        Ac,Bc,Gc,Jc = fBuildSystem_Linear(M,C,K,Sa,Sv,Sd,Sp=Sp,Rp=Rp,Yp=Yp,Method='augmented_first_order')
        # Reference values for test
        Ac_ref, Bc_ref, Gc_ref, Jc_ref = EmptyStateMat(nDOF_ext,nU,nY)
        Ac_ref[0,1] = 1
        Ac_ref[1,2] = 1/J_LSS
        Ac_ref[1,3] = -1/J_LSS
        Ac_ref[1,3] = -1/J_LSS
        Gc_ref[0,1] = 1
        Gc_ref[1,3] = 1

        np.testing.assert_equal(Ac,Ac_ref)
        np.testing.assert_equal(Bc,Bc_ref)
        np.testing.assert_equal(Gc,Gc_ref)
        np.testing.assert_equal(Jc,Jc_ref)

    def test_build_tower1shaft(self):
        # Build continuous matrices for the "Tower (1 mode) + Shaft " case
        nDOF_2nd = 2           # Mech DOFs     :  q  = [u, psi]
        nDOF_1st = 2*nDOF_2nd  # Mech state    :  x  = [u, psi, u_dot, psi_dot]
        nP       = 3           # Extended loads:  p  = [T, Ma ,Mg]
        nDOF_ext = nDOF_1st+nP # Ext. state    : x~  = [u, psi, udot, psi_dot , T, Qa,Qg]
        nY       = 3           # Outputs       :  y  = [u_ddot,psi_dot,Mg]
        nU       = 0           # Inputs
        M_twr = 3.0
        K_twr = 5.0
        C_twr = 10.0
        J_LSS = 2.0
        # General init of zero matrices
        M,C,K,Sa,Sv,Sd,Sp,Rp,Yp = EmptySystemMat(nDOF_2nd, nY, nP)
        # Setting values
        M[0,0]  = M_twr
        M[1,1]  = J_LSS
        K[0,0]  = K_twr
        C[0,0]  = C_twr
        Sa[0,0] = 1  # uddot = qddot[0]
        Sv[1,1] = 1  # psidot = qdot[1]
        Yp[2,2] = 1  # Direct feed-through of Mg
        Sp[0,0] = 1  # T = p[0]
        Sp[1,1] = 1  # dQ = p[1] -p[2]
        Sp[1,2] = -1 # dQ = p[1] -p[2]

        Ac,Bc,Gc,Jc = fBuildSystem_Linear(M,C,K,Sa,Sv,Sd,Sp=Sp,Rp=Rp,Yp=Yp,Method='augmented_first_order')

        # Reference values for test, NOTE: only true because no couplings assumed here
        Ac_ref, Bc_ref, Gc_ref, Jc_ref = EmptyStateMat(nDOF_ext,nU,nY)
        Gc_ref[0,0] = -K_twr/M_twr
        Gc_ref[0,2] = -C_twr/M_twr
        Gc_ref[0,4] =  1/M_twr
        Gc_ref[1,3] =  1           # psi_dot = x~[3]
        Gc_ref[2,6] =  1           # Mq      = x~[6]

        np.testing.assert_almost_equal(Gc,Gc_ref)

#         print('Ac\n',Ac)
#         print('Bc\n',Bc)
#         print('Gc\n',Gc)
#         print('Jc\n',Jc)

if __name__=='__main__':
    unittest.main()
