
##
import numpy as np
import copy
import unittest
try:
    from .yams import *
except:
    from yams import *

def main():

    # Main Parameters
    nSpan_twr   = 101
    nSpan_bld   = 61
    nShapes_twr = 1 # 0,1,2
    nShapes_bld = 1 # 0,1,2

    bCompat =False

    bSft     = 1
    bHubMass = 1
    bNacMass = 1
    bBldMass = 1
    bInit    = 1
    nB       = 2 # 2 or 3
    main_axis ='x'

    nDOF = 1 + nShapes_twr + nShapes_bld * nB

    q = np.zeros((nDOF,1))
    if bInit:
        q[:] = 1 # Define some kind of initial conditions

    ## --- Strucural and geometrical Inputs
    L_twr   = 100
    EI_twr  = 2*10**12
    m_twr   = 9*10**3
    L_bld   = 60
    EI_bld  = 2*10**10
    m_bld   = 5*10**2
    GKt_bld = 7*10**11
    jxx_bld = 10*5

    r_ET        = np.array([[0]    ,[0],[0]]  )
    r_TN        = np.array([[L_twr],[0],[0]]  )
    r_NGnac_inN = np.array([[0]    ,[0],[2.0]])
    r_NS_inN    = np.array([[0]    ,[0],[-10]])
    r_SGhub_inS = np.array([[0]    ,[0],[0]]  )
    r_SR_inS    = np.array([[0]    ,[0],[0]]  )
    r_RGhub_inS = np.array([[0]    ,[0],[0]]  )

    M_hub=10**5*bHubMass
    IR_hub = np.zeros((3,3))
    IR_hub[0,0] = 2*10**5
    IR_hub[1,1] = 2*10**5
    IR_hub[2,2] = 3*10**5 
    IR_hub = IR_hub * bHubMass

    M_nac   = 4*10**5 
    I0_nac=np.zeros((3,3)) 
    I0_nac[0,0]=7*10**6
    I0_nac[1,1]=3*10**6
    I0_nac[2,2]=1*10**6 
    # Inertias not at COG...
    IG_hub = fTranslateInertiaMatrix(IR_hub, M_hub, np.array([0,0,0]), r_RGhub_inS)
    IG_nac = fTranslateInertiaMatrixToCOG(I0_nac,M_nac, -r_NGnac_inN)
    Psi=nShapes_twr+1

    ## Derived parameters
    iPsi = nShapes_twr # Index of DOF corresponding to azimuth
    # --------------------------------------------------------------------------------}
    ## --- Creating bodies
    # --------------------------------------------------------------------------------{

    # Bld
    # TODO
    # TODO - THIS HAS SOME INITIAL CONDITION IN IT
    #Twr=UniformBeamBody('Tower', nShapes_twr, nSpan_twr, L_twr, EI_twr , m_twr, Mtop=0, jxxG=None, GKt=None, bCompatibility=bCompat)
    Bld1=Body()
    Bld1.MM = np.array([
     [  3.0000E+04,   0.0000E+00,   0.0000E+00,   0.0000E+00,   5.2444E+03,   0.0000E+00,  -2.4905E+02,  -1.1333E+03],
     [  0.0000E+00,   3.0000E+04,   0.0000E+00,  -5.2401E+03,   0.0000E+00,   9.0000E+05,   0.0000E+00,   0.0000E+00],
     [  0.0000E+00,   0.0000E+00,   3.0000E+04,   0.0000E+00,  -9.0000E+05,   0.0000E+00,   1.1746E+04,  -6.5057E+03],
     [  0.0000E+00,  -5.2401E+03,   0.0000E+00,   6.0150E+06,   0.0000E+00,  -4.3043E+05,   0.0000E+00,   0.0000E+00],
     [  5.2444E+03,   0.0000E+00,  -9.0000E+05,   0.0000E+00,   3.6015E+07,   0.0000E+00,  -5.1196E+05,   8.1533E+04],
     [  0.0000E+00,   9.0000E+05,   0.0000E+00,  -4.3043E+05,   0.0000E+00,   3.6000E+07,   0.0000E+00,   0.0000E+00],
     [ -2.4905E+02,   0.0000E+00,   1.1746E+04,   0.0000E+00,  -5.1196E+05,   0.0000E+00,   7.5019E+03,   4.2759E+00],
     [ -1.1333E+03,   0.0000E+00,  -6.5057E+03,   0.0000E+00,   8.1533E+04,   0.0000E+00,   4.2759E+00,   7.5066E+03]])
    Bld1.PhiU =np.zeros(nShapes_bld) # HACK since cannot access nf
    Bld1.KK = np.zeros((8, 8))
    Bld1.KK[6,6:]= np.array([ 2.8624E+05, -1.0224E+03])
    Bld1.KK[7,6:]= np.array([-1.0224E+03,  1.1249E+07])
    Bld1.MM = Bld1.MM[:6+nShapes_bld,:6+nShapes_bld]
    Bld1.KK = Bld1.KK[:6+nShapes_bld,:6+nShapes_bld]
    Bld2= copy.deepcopy(Bld1)
    Bld3= copy.deepcopy(Bld1)
    #print('Bld1.KK')
    #print(Bld1.KK)
    #print('Bld1.MM')
    #print(Bld1.MM)

    # ShaftHub Body 
    Sft=RigidBody('ShaftHub',M_hub,IG_hub,r_SGhub_inS);
    #print('Sft.MM')
    #print(Sft.MM)

    # Nacelle Body
    Nac=RigidBody('Nacelle',M_nac,IG_nac,r_NGnac_inN);
    #print('Nac.MM')
    #print(Nac.MM)

    # Tower Body
    # TODO
    if nB==2:
        Mtop=(Bld1.Mass+Bld2.Mass)*bBldMass + Sft.Mass + Nac.Mass;
    else:
        Mtop=(Bld1.Mass+Bld2.Mass+Bld3.Mass)*bBldMass + Sft.Mass + Nac.Mass;
    # TODO - THIS HAS SOME INITIAL CONDITION IN IT

    Twr=UniformBeamBody('Tower', nShapes_twr, nSpan_twr, L_twr, EI_twr , m_twr, Mtop=Mtop, bAxialCorr=False)
    #  Temporary
    x_0=np.array([[0],[0],[0]])
    R_0b=np.eye(3)
    gz=q[0:nShapes_twr,0]
    v_0   = np.zeros(6+nShapes_twr)
    a_v_0 = np.zeros(6+nShapes_twr)
    # TODO: to fully match matlab code, need "UseShapeIntegral" implemented
    Twr.updateKinematics(x_0,R_0b,gz,v_0,a_v_0)
    Twr.computeMassMatrix()
    #print('Twr.KK')
    #print(Twr.KK)
    #print('Twr.MM')
    #print(Twr.MM)

#     CyT=- np.array([ Twr.PhiV[0][2,-1],  1.5065E-01, 0, 0]) # End value of shapes functions in y direction
    CxT=  np.zeros(Twr.nf)
    CyT=  np.zeros(Twr.nf)
    CzT=  np.zeros(Twr.nf)
    for j,v in enumerate(Twr.PhiV):
        if main_axis=='x':
            iMainAxis=2 # TODO
            CyT[j]=-v[2,-1] # A deflection along z gives a negative angle around y
            CzT[j]= v[0,-1] # A deflection along y gives a positive angle around z
            #print('Alpha y - mode {}:'.format(j+1),CyT[j])
        elif main_axis=='z':
            CyT[j]= v[0,-1] # A deflection along x gives a positive angle around y
            CxT[j]=-v[0,-1] # A deflection along y gives a positive angle around z

    # --------------------------------------------------------------------------------}
    ## --- "Manual connection"
    # --------------------------------------------------------------------------------{
    # link E-T
    nf_T=Twr.nf;
    R_ET     = np.identity(3)
    B_T      = np.array([])
    # B_T      = fBMatRecursion(,np.vstack((Bx_TN,Bt_TN)),R_ET,r_ET)
    B_T_inT  = fB_inB(R_ET, B_T)
    BB_T_inT = fB_aug(B_T_inT, nf_T)
    MM_T     = fBMB(BB_T_inT,Twr.MM)
    KK_T     = fBMB(BB_T_inT,Twr.KK)
    #print('BB_T_inT')
    #print(BB_T_inT)
    #print('MM_T',MM_T)
    #print('KK_T',KK_T)

    CyT=CyT[:Twr.nf]

    # ---------------------------------------------
    # Link T-N
    nf_N = 0 
    if nShapes_twr == 0:
        Bx_TN = np.array([])
        Bt_TN = np.array([])
        alpha_y=0
    elif nShapes_twr == 1:
        Bx_TN = np.array([[0],[0],[1]])
        Bt_TN = np.array([[0],[CyT[0]],[0]])
        alpha_y = np.dot(CyT.ravel(), q[0,0].ravel())
    elif nShapes_twr == 2:
        Bx_TN = np.array([[0,0],[0,0],[1,1]])
        Bt_TN = np.array([[0,0],[CyT[0],CyT[1]],[0,0]])
        alpha_y = np.dot(CyT.ravel() , q[:2,0].ravel())
    else:
        # TODO use CzT
        raise NotImplementedError()
    #print('alpha_y',alpha_y)
    R_TN     = fRoty(alpha_y)
    R_EN     = np.dot(R_ET, R_TN)
    B_N      = fBMatRecursion(B_T,np.vstack((Bx_TN,Bt_TN)),R_ET,r_TN)
    B_N_inN  = fB_inB(R_EN, B_N)
    BB_N_inN = fB_aug(B_N_inN, nf_N)
    MM_N     = fBMB(BB_N_inN,Nac.MM)
    KK_N     = fBMB(BB_N_inN,Nac.KK)
    #print('BB_N_inN')
    #print(BB_N_inN)

    # ---------------------------------------------
    # Link N-S
    nf_S = Sft.nf
    q_psi = q[iPsi,0]

    R_NS     = fRotz(q_psi + np.pi)
    R_ES     = np.dot(R_EN, R_NS)
    r_NS     = np.dot(R_EN, r_NS_inN)
    Bx_NS    = np.array([[0],[0],[0]])
    Bt_NS    = np.array([[0],[0],[1]])
    B_S      = fBMatRecursion(B_N,np.vstack((Bx_NS,Bt_NS)),R_EN,r_NS)
    B_S_inS  = fB_inB(R_ES, B_S)
    BB_S_inS = fB_aug(B_S_inS, nf_S)
    MM_S     = fBMB(BB_S_inS,Sft.MM)
    KK_S     = fBMB(BB_S_inS,Sft.KK)
    #print('BB_S_inS')
    #print(BB_S_inS)
    #print('MM_S')
    #print(MM_S)

    # ---------------------------------------------
    # Link S-B1
    nf_B1 = Bld1.nf
    nf_B2 = Bld2.nf
    if nB == 2:
        R_SB1 = fRotz(0 * np.pi + 0)
        R_SB2 = fRotz(0 * np.pi + np.pi)
        R_SB3 = np.zeros((3,3))
        nf_B3 = 0
    elif nB == 3:
        R_SB1 = fRotz(0 * np.pi + 0)
        R_SB2 = fRotz(0 * np.pi - 2 * np.pi / 3)
        R_SB3 = fRotz(0 * np.pi + 2 * np.pi / 3)
        nf_B3 = Bld3.nf
    else:
        raise NotImplementedError()
    nf_B = nf_B1 + nf_B2 + nf_B3

    R_EB1 = np.dot(R_ES, R_SB1)
    R_EB2 = np.dot(R_ES, R_SB2)
    R_EB3 = np.dot(R_ES, R_SB3)
    r_SR  = np.dot(R_ES, r_SR_inS)
    B_R = fBMatRecursion(B_S,[],R_ES,r_SR)
    B_R_bis = fBMatTranslate(B_S,r_SR)
    B_B1_inB1 = fB_inB(R_EB1, B_R)
    B_B2_inB2 = fB_inB(R_EB2, B_R)
    B_B3_inB3 = fB_inB(R_EB3, B_R)
    BB_B1_inB1 = fB_aug(B_B1_inB1, nf_B, nf_B1, 0          )
    BB_B2_inB2 = fB_aug(B_B2_inB2, nf_B, nf_B2, nf_B1      )
    BB_B3_inB3 = fB_aug(B_B3_inB3, nf_B, nf_B3, nf_B1+nf_B2)

    MM_B1 = fBMB(BB_B1_inB1,Bld1.MM)
    KK_B1 = fBMB(BB_B1_inB1,Bld1.KK)
    MM_B2 = fBMB(BB_B2_inB2,Bld2.MM)
    KK_B2 = fBMB(BB_B2_inB2,Bld2.KK)
    if nB == 3:
        MM_B3 = fBMB(BB_B3_inB3,Bld3.MM)
        KK_B3 = fBMB(BB_B3_inB3,Bld3.KK)
    else:
        MM_B3 = MM_B2 * 0
        KK_B3 = KK_B2 * 0
     
    #print('BB_B1_inB1')
    #print(BB_B1_inB1)
    # print('MM_B2')
    # print(MM_B2)

    # --- Final assembly
    MM = (MM_B1 + MM_B2 + MM_B3) * bBldMass
    MM[:iPsi+1,:iPsi+1] += MM_S
    MM[:nShapes_twr,:nShapes_twr] += MM_T + MM_N

    KK = (KK_B1 + KK_B2 + KK_B3)
    KK[:iPsi+1,:iPsi+1] += KK_S
    KK[:nShapes_twr,:nShapes_twr] += KK_T + KK_N

    ## Display to screen
    MM[np.abs(MM)< 1e-09] = 0
    #print('M ("manually" built)')
    #print(MM)
    #print('K ("manually" build)')
    #print(KK)
    return MM,KK


class Test(unittest.TestCase):
    def test_TNSB_article(self):
        MM,KK=main()
        np.testing.assert_allclose(MM[0,0],7.86e5 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[1,1],7.23e7 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[2,2],7.50e3 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[3,3],7.50e3 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[0,2],7.71e3 ,rtol  = 1e-3)
        np.testing.assert_allclose(MM[0,3],1.578e4,rtol = 1e-3 )
        np.testing.assert_allclose(KK[0,0],6.01e6 ,rtol  = 1e-3)
        np.testing.assert_allclose(KK[1,1],0.00e0 ,rtol  = 1e-3)
        np.testing.assert_allclose(KK[2,2],2.86e5 ,rtol  = 1e-3)
        np.testing.assert_allclose(KK[2,2],2.86e5 ,rtol  = 1e-3)


if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    unittest.main()
