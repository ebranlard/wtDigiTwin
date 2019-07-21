##
import numpy as np
import copy
import matplotlib.pyplot as plt
import os

from yams import *
from TNSB import *
import weio

# --------------------------------------------------------------------------------}
# --- Creating a TNSB model from a FAST model
# --------------------------------------------------------------------------------{
def FASTmodel2TNSB(EDfile,nB=3,nShapes_twr=2, nShapes_bld=0,nSpan_twr=101,nSpan_bld=61,bHubMass=1,bNacMass=1,bBldMass=1,DEBUG=False,main_axis ='x'):
    
    nDOF = 1 + nShapes_twr + nShapes_bld * nB # +1 for Shaft
    q = np.zeros((nDOF,1)) # TODO, full account of q not done

    # --- Input data from ED file
    ED      = weio.read(EDfile)
    rootdir = os.path.dirname(EDfile)
    bldfile = os.path.join(rootdir,ED['BldFile(1)'].strip('"')).replace('\\','/')
    twrfile = os.path.join(rootdir,ED['TwrFile'].strip('"')).replace('\\','/')
    twr     = weio.read(twrfile)
    bld     = weio.read(bldfile)

    tilt_up=-ED['ShftTilt']*np.pi/180*0
    ## --- Strucural and geometrical Inputs
    if main_axis=='x':
        r_ET_inE    = np.array([[ED['TowerBsHt']]              ,[0],[0]]) # NOTE: could be used to get hub height
        r_TN_inT    = np.array([[ED['TowerHt']-ED['TowerBsHt']],[0],[0]])
        r_NGnac_inN = np.array([[ED['NacCMzn']]                ,[0],[ED['NacCMxn']]] )
        r_NS_inN    = np.array([[ED['Twr2Shft']]               ,[0],[0]]) # S on tower axis
        r_SR_inS    = np.array([[0]                            ,[0],[ED['OverHang']]] ) # S and R 
        r_SGhub_inS = np.array([[0]                            ,[0],[ED['OverHang']+ED['HubCM']]]   ) # 
    else:
        raise NotImplementedError()
    r_RGhub_inS = - r_SR_inS + r_SGhub_inS

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

    # --------------------------------------------------------------------------------}
    ## --- Creating bodies
    # --------------------------------------------------------------------------------{
    # Bld
    Blds=[]
    Blds.append(FASTBeamBody('blade',ED,bld,nShapes=nShapes_bld, nSpan=nSpan_bld, main_axis=main_axis))
    Blds[0].MM *=bBldMass
    for iB in range(nB-1):
        Blds.append(copy.deepcopy(Blds[0]))
    # ShaftHub Body 
    Sft=RigidBody('ShaftHub',M_hub,IG_hub,r_SGhub_inS);
    # Nacelle Body
    Nac=RigidBody('Nacelle',M_nac,IG_nac,r_NGnac_inN);
    # Tower Body
    Twr = FASTBeamBody('tower',ED,twr,nShapes=nShapes_twr, nSpan=nSpan_twr, main_axis=main_axis)
    if DEBUG:
        print('IG_hub')
        print(IG_hub)
        print('IG_nac')
        print(IG_nac)
        print('I_gen_LSS', ED['GenIner']*ED['GBRatio']**2)
        print('I_hub_LSS', ED['hubIner'])
        print('I_rot_LSS', nB*Blds[0].MM[5,5])
        print('I_tot_LSS', nB*Blds[0].MM[5,5]+ED['hubIner']+ED['GenIner']*ED['GBRatio']**2) 
        print('r_NGnac_inN',r_NGnac_inN.T)
        print('r_SGhub_inS',r_SGhub_inS.T)
    # --------------------------------------------------------------------------------}
    # --- Manual assembly 
    # --------------------------------------------------------------------------------{
    Struct = manual_assembly(Twr,Nac,Sft,Blds,q,r_ET_inE,r_TN_inT,r_NS_inN,r_SR_inS,main_axis='x',tilt_up=tilt_up,DEBUG=DEBUG)
    return Struct


if __name__=='__main__':
    np.set_printoptions(linewidth=500)
    Struct= FASTmodel2TNSB(EDfile='../_data/NREL5MW_ED.dat', DEBUG=False)
