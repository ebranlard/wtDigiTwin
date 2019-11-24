"""
TNSB refer to : Tower Nacelle, Shaft, Blades

This scripts provides some helper functions to simulates such an assembly of bodies using the Rayleigh-Rizt approximation and joint coordinates. 

The theory is provided in the reference below. The article also contains an example in the its section, which is reproduced in the file test_TNSB.py.


Reference:
     [1]: Branlard, Flexible multibody dynamics using joint coordinates and the Rayleigh-Ritz approximation: the general framework behind and beyond Flex, Wind Energy, 2019
"""

##
import numpy as np
import copy
import os
import unittest

try:
    from .yams import *
except:
    from yams import *

# --------------------------------------------------------------------------------}
# --- Manual assembly of a TNSB model 
# --------------------------------------------------------------------------------{
def manual_assembly(Twr,Nac,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis='x',tilt_up=0,DEBUG=False):

    # Main Parameters
    nDOF = len(q)
    iPsi = Twr.nf # Index of DOF corresponding to azimuth

    if main_axis=='z':
        raise NotImplementedError()

#     CyT=- np.array([ Twr.PhiV[0][2,-1],  1.5065E-01, 0, 0]) # End value of shapes functions in y direction
    CxT=  np.zeros(Twr.nf)
    CyT=  np.zeros(Twr.nf)
    CzT=  np.zeros(Twr.nf)
    for j,v in enumerate(Twr.PhiV):
        if main_axis=='x':
            iMainAxis=2 # TODO
            CyT[j]=-v[2,-1] # A deflection along z gives a negative angle around y
            CzT[j]= v[1,-1] # A deflection along y gives a positive angle around z # TODO TODO CHECK ME
            #print('Alpha y - mode {}:'.format(j+1),CyT[j])
        elif main_axis=='z':
            CxT[j]=-v[1,-1] # A deflection along y gives a negative angle around x # TODO TODO CHECK ME
            CyT[j]= v[0,-1] # A deflection along x gives a positive angle around y
    CyT=CyT[:Twr.nf]
    #                     Bt_pc=zeros(3,p.nf);
    #                     for j=1:p.nf
    #                         Bx_pc(:,j)=p.PhiU{j}(:,iNode);
    #                         Bt_pc(:,j)=[0; -p.PhiV{j}(3,iNode); p.PhiV{j}(2,iNode)];
    #                     end

    # --------------------------------------------------------------------------------}
    ## --- "Manual connection"
    # --------------------------------------------------------------------------------{
    # link E-T
    R_ET     = np.identity(3)
    B_T      = np.array([])
    # B_T      = fBMatRecursion(,np.vstack((Bx_ET,Bt_ET)),R_ET,r_ET)
    B_T_inT  = fB_inB(R_ET, B_T)
    BB_T_inT = fB_aug(B_T_inT, Twr.nf)
    MM_T     = fBMB(BB_T_inT,Twr.MM)
    KK_T     = fBMB(BB_T_inT,Twr.KK)
    DD_T     = fBMB(BB_T_inT,Twr.DD)


    # ---------------------------------------------
    # Link T-N
    # TODO
    if Twr.nf == 0:
        Bx_TN = np.array([])
        Bt_TN = np.array([])
        alpha_y=0
    elif Twr.nf == 1:
        if main_axis=='x':
            Bx_TN = np.array([[0],[0],[1]])
            Bt_TN = np.array([[0],[CyT[0]],[0]])
            alpha_y = np.dot(CyT.ravel(), q[0,0].ravel())
    elif Twr.nf == 2:
        if main_axis=='x':
            Bx_TN = np.array([[0,0],[0,0],[1,1]])
            Bt_TN = np.array([[0,0],[CyT[0],CyT[1]],[0,0]])
            alpha_y = np.dot(CyT.ravel() , q[:2,0].ravel())
    else:
        # TODO use CzT
        raise NotImplementedError()
    #print('alpha_y',alpha_y)
    R_TN     = R_y(alpha_y)
    R_EN     = np.dot(R_ET, R_TN)
    B_N      = fBMatRecursion(B_T,Bx_TN,Bt_TN,R_ET,r_TN_inT)
    B_N_inN  = fB_inB(R_EN, B_N)
    BB_N_inN = fB_aug(B_N_inN, Nac.nf)
    MM_N     = fBMB(BB_N_inN,Nac.MM)
    KK_N     = fBMB(BB_N_inN,Nac.KK)

    # ---------------------------------------------
    # Link N-S
    q_psi = q[iPsi,0]
    R_NS     = np.dot(R_y(-tilt_up),R_z(q_psi + np.pi)) # << tilt 
    R_ES     = np.dot(R_EN, R_NS)
    r_NS     = np.dot(R_EN, r_NS_inN)
    Bx_NS    = np.array([[0],[0],[0]])
    Bt_NS    = np.array([[0],[0],[1]])
    B_S      = fBMatRecursion(B_N,Bx_NS,Bt_NS,R_EN,r_NS)
    B_S_inS  = fB_inB(R_ES, B_S)
    BB_S_inS = fB_aug(B_S_inS, Sft.nf)
    MM_S     = fBMB(BB_S_inS,Sft.MM)
    KK_S     = fBMB(BB_S_inS,Sft.KK)

    # ---------------------------------------------
    # Link S-B1
    nB   = len(Blds)
    # Point R
    r_SR  = np.dot(R_ES, r_SR_inS)
    B_R = fBMatRecursion(B_S,[],[],R_ES,r_SR)
    B_R_bis = fBMatTranslate(B_S,r_SR)
    # Points B1, B2, B3
    MM_B      = np.zeros((nDOF,nDOF))
    KK_B      = np.zeros((nDOF,nDOF))
    DD_B      = np.zeros((nDOF,nDOF))
    nf_done=0
    nf_tot = sum([B.nf for B in Blds])
    for i,B in enumerate(Blds):
        psi_B= -i*2*np.pi/nB # 0 -2pi/2 2pi/3  or 0 pi
        if main_axis=='x':
            R_SB = R_z(0*np.pi + psi_B)
        elif main_axis=='z':
            R_SB = R_x(0*np.pi + psi_B)
        R_EB       = np.dot(R_ES, R_SB)
        B_B_inB    = fB_inB(R_EB, B_R)
        BB_B_inB   = fB_aug(B_B_inB, nf_tot, B.nf, nf_done)

        nf_done   += B.nf
        # Full matrices 
        MM_B +=     fBMB(BB_B_inB,B.MM)
        KK_B +=     fBMB(BB_B_inB,B.KK)
        DD_B +=     fBMB(BB_B_inB,B.DD)

     

    # --- Final assembly
    MM = MM_B 
    MM[:iPsi+1,:iPsi+1] += MM_S
    MM[:Twr.nf,:Twr.nf] += MM_T + MM_N

    KK = KK_B
    KK[:iPsi+1,:iPsi+1] += KK_S
    KK[:Twr.nf,:Twr.nf] += KK_T + KK_N

    DD = DD_B 
    DD[:Twr.nf,:Twr.nf] += DD_T
    ## Display to screen
    MM[np.abs(MM)< 1e-09] = 0
    if DEBUG:
        print('--------------------- Geom ---------------------')
        print('r_ET_inE   ',r_ET_inE   .T)
        print('r_TN_inT   ',r_TN_inT   .T)
        print('r_NS_inN   ',r_NS_inN   .T)
        print('r_SR_inS   ',r_SR_inS   .T)
        print('-------------------- Tower ---------------------')
        print('CyT\n',CyT)
        print('alpha_y',alpha_y)
        print('B_T\n',B_T)
        print('B_T_inT\n',B_T_inT)
        print('BB_T_inT\n',BB_T_inT)
#         print('MM_T\n',MM_T)
#         print('KK_T\n',KK_T)
#         print('DD_T\n',DD_T)
        print('------------------- Nacelle --------------------')
        print('B_N\n',B_N)
        print('B_N_inN\n',B_N_inN)
        print('BB_N_inN\n',BB_N_inN)
        print('MM_N\n',MM_N)
#         print('-------------------- Shaft ---------------------')
#         print('R_NS\n',R_NS)
#         print('BB_S_inS\n',BB_S_inS)
#         print('MM_S\n',MM_S)
#         print('------------------- Blades ---------------------')
#         #print('BB_B1_inB1')
#         #print(BB_B1_inB1)
#         print('MM_B\n',MM_B)
#         print('KK_B\n',KK_B)
#         print('DD_B\n',DD_B)
#         print('-------------------- Full ----------------------')
#         print('M ("manually" built)')
#         print(MM)
#         print('K ("manually" build)')
#         print(KK)

    ## Eigenvalue analysis
    #[Q,Lambda]=eig(K,M);
    #Omega2=diag(Lambda);
    #[Omega2,Isort]=sort(Omega2);
    #Q=Q(:,Isort);
    #f_eva= sqrt(Omega2)/(2*pi);
    #for i=1:length(f_eva);
    #    fprintf('f%d = %.3f \n',i,f_eva(i))


    # --- returning everthin in a structure class
    class Structure():
        pass
    Struct      = Structure()
    Struct.Twr  = Twr
    Struct.Nac  = Nac
    Struct.Sft  = Sft
    Struct.Blds = Blds
    Struct.MM   = MM
    Struct.KK   = KK
    Struct.DD   = DD

    Struct.iPsi = iPsi # Index 

    return Struct



if __name__=='__main__':
    np.set_printoptions(linewidth=500)
