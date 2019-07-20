##
import numpy as np
import copy
import matplotlib.pyplot as plt

import weio
from yams import *
from flexibility import *

np.set_printoptions(linewidth=500)

# --- Input data from ED file
ED  = weio.read('../_DLC120_ws13_ye000_s1_r1_ED.dat')
twr = weio.read('../_NREL5MW_ED_Tower_Offshore.dat')
bld = weio.read('../_NREL5MW_ED_Blade.dat')


# Main Parameters
nSpan_twr   = 101
nSpan_bld   = 61
nShapes_twr = 2 # 0,1,2
nShapes_bld = 3 # 0,1,2

bHubMass = 1
bNacMass = 1
bBldMass = 1
bInit    = 0
nB       = 3 # 2 or 3
main_axis ='x'

nDOF = 1 + nShapes_twr + nShapes_bld * nB

q = np.zeros((nDOF,1))
if bInit:
    q[:] = 1 # Define some kind of initial conditions

tilt_up=-ED['ShftTilt']*np.pi/180*0
## --- Strucural and geometrical Inputs
if main_axis=='x':
    r_ET        = np.array([[ED['TowerBsHt']]              ,[0],[0]]) # NOTE: could be used to get hub height
    r_TN        = np.array([[ED['TowerHt']-ED['TowerBsHt']],[0],[0]])
    r_NGnac_inN = np.array([[ED['NacCMzn']]                ,[0],[ED['NacCMxn']]] )
    r_NS_inN    = np.array([[ED['Twr2Shft']]               ,[0],[0]]) # S on tower axis
    r_SR_inS    = np.array([[0]                            ,[0],[ED['OverHang']]] ) # S and R 
    r_SGhub_inS = np.array([[0]                            ,[0],[ED['OverHang']+ED['HubCM']]]   ) # 
else:
    raise NotImplementedError()
r_RGhub_inS = - r_SR_inS + r_SGhub_inS

print('r_ET       ',r_ET       .T)
print('r_TN       ',r_TN       .T)
print('r_NGnac_inN',r_NGnac_inN.T)
print('r_NS_inN   ',r_NS_inN   .T)
print('r_SR_inS   ',r_SR_inS   .T)
print('r_SGhub_inS',r_SGhub_inS.T)
print('r_RGhub_inS',r_RGhub_inS.T)


if main_axis=='x':
    M_hub=ED['HubMass']*bHubMass
    IR_hub = np.zeros((3,3))
    IR_hub[0,0] = 0
    IR_hub[1,1] = 0
    IR_hub[2,2] = ED['HubIner'] + ED['GenIner']*ED['GBRatio']**2
    IR_hub = IR_hub * bHubMass

    M_nac   = ED['NacMass'] *bNacMass
    I0_nac=np.zeros((3,3)) 
    I0_nac[0,0]= ED['NacYIner']
    I0_nac[1,1]=0
    I0_nac[2,2]=0
    I0_nac = I0_nac*bNacMass
else:
    raise NotImplementedError()

# Inertias not at COG...
IG_hub = fTranslateInertiaMatrix(IR_hub, M_hub, np.array([0,0,0]), r_RGhub_inS)
IG_nac = fTranslateInertiaMatrixToCOG(I0_nac,M_nac, -r_NGnac_inN)
print('IG_hub')
print(IG_hub)
print('IG_nac')
print(IG_nac)

## Derived parameters
iPsi = nShapes_twr # Index of DOF corresponding to azimuth
# --------------------------------------------------------------------------------}
## --- Creating bodies
# --------------------------------------------------------------------------------{
def fast_flexible_body(body_type,ED,inp,nShapes=2,main_axis='x',nSpan=40):
    """ 
    INPUTS:
       nSpan: number of spanwise station used (interpolated from input)
              Use -1 or None to use number of stations from input file
    """
    B=body()
    B.nf = nShapes

    # --- Reading coefficients
    exp = np.arange(2,7)
    if body_type.lower()=='blade':
        coeff = np.array([[ inp['BldFl1Sh(2)'], inp['BldFl2Sh(2)'], inp['BldEdgSh(2)']],
                          [ inp['BldFl1Sh(3)'], inp['BldFl2Sh(3)'], inp['BldEdgSh(3)']],
                          [ inp['BldFl1Sh(4)'], inp['BldFl2Sh(4)'], inp['BldEdgSh(4)']],
                          [ inp['BldFl1Sh(5)'], inp['BldFl2Sh(5)'], inp['BldEdgSh(5)']],
                          [ inp['BldFl1Sh(6)'], inp['BldFl2Sh(6)'], inp['BldEdgSh(6)']]])

    elif body_type.lower()=='tower':
        coeff = np.array([[ inp['TwFAM1Sh(2)'], inp['TwFAM2Sh(2)'], inp['TwSSM1Sh(2)'], inp['TwSSM2Sh(2)']],
                          [ inp['TwFAM1Sh(3)'], inp['TwFAM2Sh(3)'], inp['TwSSM1Sh(3)'], inp['TwSSM2Sh(3)']],
                          [ inp['TwFAM1Sh(4)'], inp['TwFAM2Sh(4)'], inp['TwSSM1Sh(4)'], inp['TwSSM2Sh(4)']],
                          [ inp['TwFAM1Sh(5)'], inp['TwFAM2Sh(5)'], inp['TwSSM1Sh(5)'], inp['TwSSM2Sh(5)']],
                          [ inp['TwFAM1Sh(6)'], inp['TwFAM2Sh(6)'], inp['TwSSM1Sh(6)'], inp['TwSSM2Sh(6)']]])
    else:
        raise Exception('Body type not supported {}'.format(body_type))
    nShpMax=coeff.shape[1]
    if nShapes>nShpMax:
        raise Exception('A maximum of {} shapes function possible with FAST {} body'.format(nShpMax,body_type))

    # --- Structural properties
    if body_type.lower()=='blade':
        prop     = bld['BldProp']
        span_max = ED['TipRad']   # TODO TODO do somthing about hub rad
        s_bar, m, EIFlp, EIEdg  =prop[:,0], prop[:,3], prop[:,4], prop[:,5]

    elif body_type.lower()=='tower':
        prop     = twr['TowProp']
        span_max = ED['TowerHt']-ED['TowerBsHt']
        s_bar, m, EIFlp, EIEdg  = prop[:,0], prop[:,1], prop[:,2], prop[:,3]

    # --- Interpolating structural properties
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

    # --- Center of mass location
    cG      = np.zeros((3,nSpan))
    cG[iMain,:] = s_span # TODO blade COG

    # --- Shape functions
    U   = np.zeros((nShapes,3,nSpan)) # Shape
    dU  = np.zeros((nShapes,3,nSpan)) # Slope
    ddU = np.zeros((nShapes,3,nSpan)) # Curvature
    for j in np.arange(nShapes):
        iAxis = ShapeDir[j]
        U[j][iAxis,:], dU[j][iAxis,:], ddU[j][iAxis,:] = polymode(s_span,coeff[:,j],exp)

    B.MM = GMBeam(cG, s_span, m, U, bOrth=body_type=='blade')
    B.KK = GKBeam(s_span, EI, ddU , bOrth=body_type=='blade')

    # --- Storing data in object
    B.U      = U
    B.dU     = dU
    B.ddU    = dU
    B.s_span = s_span
    B.m = m
    B.EI = EI
    B.cG = cG
    B.main_axis = main_axis

    return B

# Bld
Bld1 = fast_flexible_body('blade',ED,bld,nShapes=nShapes_bld, nSpan=nSpan_bld, main_axis=main_axis)
Bld2 = copy.deepcopy(Bld1)
Bld3 = copy.deepcopy(Bld1)
print('Bld1.KK')
print(Bld1.KK)
print('Bld1.MM')
print(Bld1.MM)

print('I_gen_LSS', ED['GenIner']*ED['GBRatio']**2)
print('I_hub_LSS', ED['hubIner'])
print('I_rot_LSS', 3*Bld1.MM[5,5])
print('I_tot_LSS', 3*Bld1.MM[5,5]+ED['hubIner']+ED['GenIner']*ED['GBRatio']**2) 


# ShaftHub Body 
Sft=fCreateBodyRigid('ShaftHub',M_hub,IG_hub,r_SGhub_inS);
print('Sft.MM')
print(Sft.MM)

# Nacelle Body
Nac=fCreateBodyRigid('Nacelle',M_nac,IG_nac,r_NGnac_inN);
print('Nac.MM')
print(Nac.MM)

# Tower Body
Twr = fast_flexible_body('tower',ED,twr,nShapes=nShapes_twr, nSpan=nSpan_twr, main_axis=main_axis)
print('Twr.KK')
print(Twr.KK)
print('Twr.MM')
print(Twr.MM)
#     CyT=   np.array([ 0.02180798901   ,  0.41602854831]);
CxT=  np.zeros(Twr.nf)
CzT=  np.zeros(Twr.nf)
CzT=  np.zeros(Twr.nf)
for j,v in enumerate(Twr.dU):
    if main_axis=='x':
        iMainAxis=2 # TODO
        print('Alpha y - mode {}:'.format(j+1),v[iMainAxis,-1])
        CyT[j]=-v[2,-1] # A deflection along z gives a negative angle around y
        CzT[j]= v[0,-1] # A deflection along y gives a positive angle around z
    elif main_axis=='z':
        CyT[j]= v[0,-1] # A deflection along x gives a positive angle around y
        CxT[j]=-v[0,-1] # A deflection along y gives a positive angle around z

#                     Bt_pc=zeros(3,p.nf);
#                     for j=1:p.nf
#                         Bx_pc(:,j)=p.PhiU{j}(:,iNode);
#                         Bt_pc(:,j)=[0; -p.PhiV{j}(3,iNode); p.PhiV{j}(2,iNode)];
#                     end

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
print('BB_T_inT')
print(BB_T_inT)
print('MM_T',MM_T)
print('KK_T',KK_T)

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
print('alpha_y',alpha_y)
R_TN     = fRoty(alpha_y)
R_EN     = np.dot(R_ET, R_TN)
B_N      = fBMatRecursion(B_T,np.vstack((Bx_TN,Bt_TN)),R_ET,r_TN)
B_N_inN  = fB_inB(R_EN, B_N)
BB_N_inN = fB_aug(B_N_inN, nf_N)
MM_N     = fBMB(BB_N_inN,Nac.MM)
KK_N     = fBMB(BB_N_inN,Nac.KK)
# print('BB_N_inN')
# print(BB_N_inN)

# ---------------------------------------------
# Link N-S
nf_S = Sft.nf
q_psi = q[iPsi,0]

R_NS     = fRoty(-tilt_up)*fRotz(q_psi + np.pi) # << tilt 
R_ES     = np.dot(R_EN, R_NS)
r_NS     = np.dot(R_EN, r_NS_inN)
Bx_NS    = np.array([[0],[0],[0]])
Bt_NS    = np.array([[0],[0],[1]])
B_S      = fBMatRecursion(B_N,np.vstack((Bx_NS,Bt_NS)),R_EN,r_NS)
B_S_inS  = fB_inB(R_ES, B_S)
BB_S_inS = fB_aug(B_S_inS, nf_S)
MM_S     = fBMB(BB_S_inS,Sft.MM)
KK_S     = fBMB(BB_S_inS,Sft.KK)
print('BB_S_inS')
print(BB_S_inS)
print('MM_S')
print(MM_S)

# ---------------------------------------------
# Link S-B1
nf_B1 = Bld1.nf
nf_B2 = Bld2.nf
if nB == 2:
    R_SB1 = fRotz(1*np.pi + 0)
    R_SB2 = fRotz(1*np.pi + np.pi)
    R_SB3 = np.zeros((3,3))
    nf_B3 = 0
elif nB == 3:
    R_SB1 = fRotz(1*np.pi + 0)
    R_SB2 = fRotz(1*np.pi - 2*np.pi/3)
    R_SB3 = fRotz(1*np.pi + 2*np.pi/3)
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
 
# print('MM_B2')
# print(MM_B2)

# --- Final assembly
MM = (MM_B1 + MM_B2 + MM_B3)* bBldMass
MM[:iPsi+1,:iPsi+1] += MM_S 
MM[:nShapes_twr,:nShapes_twr] += MM_T + MM_N

KK = (KK_B1 + KK_B2 + KK_B3)
KK[:iPsi+1,:iPsi+1] += KK_S
KK[:nShapes_twr,:nShapes_twr] += KK_T + KK_N

## Display to screen
MM[np.abs(MM)< 1e-09] = 0
print('M ("manually" built)')
print(MM)
print('K ("manually" build)')
print(KK)
