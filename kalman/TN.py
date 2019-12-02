import numpy as np
from .kalman import *
from .kalmanfilter import KalmanFilter
from .filters import moving_average
from ws_estimator.tabulated import TabulatedWSEstimator
import yams
from yams.TNSB_FAST import FASTmodel2TNSB

# --- External dependencies!
import welib.fastlib as fastlib
import weio


class KalmanFilterTN(KalmanFilter):
    def __init__(KF, FstFile, base,  bThrustInStates=True, nShapes_twr=1):
        """

        """

        nShapes_bld   = 0 # Hard coded for TN
        nDOF_2nd      = nShapes_twr+1 # Mech DOFs     :  q = [u, psi]
        
        if nShapes_twr>1:
            raise NotImplementedError()
        if bThrustInStates:
            sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega'] )
            sAug        = np.array(['Thrust' ,'Qaero'  ,'Qgen','WS'] )
            sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            sInp        = np.array(['pitch'])
        else:
            sStates     = np.array(['ut1'  ,'psi'  ,'ut1dot','omega','Qaero','Qgen'] )
            sAug        = np.array(['Qaero','Qgen','WS'] )
            sMeas       = np.array(['TTacc','omega','Qgen','pitch'])
            sInp        = np.array(['Thrust','pitch'])

        super(KalmanFilterTN, KF).__init__(sX0=sStates,sXa=sAug,sU=sInp,sY=sMeas)

        # --- Building state/outputs connection matrices
        M,C,K,Ya,Yv,Yq,Yp,Yu,Fp,Fu,Pp,Pq,Pv = EmptySystemMat (int(KF.nX0/2), KF.nY, KF.nP, KF.nU)

        # This below is problem specific
        if nShapes_twr==1 and bThrustInStates:
            Ya[0,0] = 1    # uddot                     = qddot[0]
            Yv[1,1] = 1    # psidot                    = qdot[1]
            Yp[2,2] = 1    # Direct feed-through of Mg
            Fp[0,0] = 1    # T                         = p[0]
            Fp[1,1] = 1    # dQ                        = p[1] -p[2]
            Fp[1,2] = -1   # dQ                        = p[1] -p[2]
            Yu[3,0] = 1    # pitch direct feedthrough
        else:
            raise NotImplementedError()


        # --- Mechanical system and turbine data
        WT = FASTmodel2TNSB(FstFile , nShapes_twr=nShapes_twr,nShapes_bld=nShapes_bld, DEBUG=False, bStiffening=True, main_axis='z')
#         nGear  = WT.ED['GBRatio']
        if nShapes_twr==1:
            # TODO aerodamping
            WT.DD      = WT.DD*3.5 # increased damping to account for aero damping
        print(WT)
        KF.WT=WT

        # --- Creating a wind speed estimator (reads tabulated aerodynamic data)
        KF.wse = TabulatedWSEstimator(fst_file=FstFile)
        KF.wse.load_files(base=base,suffix='')
        #print(wse)
        # --- Building continuous and discrete state matrices
        M,C,K = WT.MM, WT.DD, WT.KK
        Xx,Xu,Yx,Yu = BuildSystem_Linear(M,C,K,Ya,Yv,Yq,Fp=Fp,Pp=Pp,Yp=Yp,Yu=Yu,Method='augmented_first_order')
        KF.setMat(Xx,Xu,Yx,Yu)



    def loadMeasurements(KF, MeasFile, nUnderSamp=1, tRange=None):
        # --- Loading "Measurements"
        ColMap = {'ut1':'TTDspFA', 'psi':'Azimuth','ut1dot':'NcIMUTVxs','omega':'RotSpeed',
                 'Thrust':'RtAeroFxh','Qaero':'RtAeroMxh','Qgen':'GenTq',
                 'WS':'RtVAvgxh', 'pitch':'BldPitch1','TTacc':'NcIMUTAxs'}
        #          'WS':'Wind1VelX', 'pitch':'BldPitch1','TTacc':'NcIMUTAxs'}
        #          'Thrust':'RotThrust','Qaero':'RtAeroMxh','Qgen':'GenTq',
        # NOTE: RotThrust contain gravity and inertia

        # TODO 
        nGear  = KF.WT.ED['GBRatio']
        df=weio.read(MeasFile).toDataFrame()
        df.columns = [  v.split('_[')[0] for v in df.columns.values] 
        if tRange is not None:
            df=df[(df['Time']>= tRange[0]) & (df['Time']<= tRange[1])] # reducing time range
        df=df.iloc[::nUnderSamp,:]                      # reducing sampling
        time = df['Time'].values
        dt   = (time[-1] - time[0])/(len(time)-1)
        df['GenSpeed'] *= 2*np.pi/60 # converted to rad/s
        df['RotSpeed'] *= 2*np.pi/60 # converted to rad/s
        df['GenTq']    *= 1000*nGear # Convert to rot torque
        df['Azimuth']  *= np.pi/180  # rad
        df['RotTorq']  *= 1000 # [kNm]->[Nm]
        df['RotThrust']*= 1000 # [kN]->[N]
        KF.df=df

        # --- 
        KF.discretize(dt, method='exponential')
        KF.setTimeVec(time)
        KF.setCleanValues(df,ColMap)

        # --- Estimate sigmas from measurements
        sigX_c,sigY_c = KF.sigmasFromClean(factor=1)

    def prepareTimeStepping(KF):
        # --- Process and measurement covariances
        KF.P, KF.Q, KF.R = KF.covariancesFromSig()
        # --- Storage for plot
        KF.initTimeStorage()

    def prepareMeasurements(KF, NoiseRFactor=0, bFilterAcc=False, nFilt=15):
        # --- Creating noise measuremnts
        KF.setYFromClean(R=KF.R, NoiseRFactor=NoiseRFactor)
        if bFilterAcc:
            KF.set_vY('TTacc',  moving_average(KF.get_vY('TTacc'),n=nFilt) )



    def timeLoop(KF):
        # --- Initial conditions
        x = KF.initFromClean()
        P = KF.P

        iY={lab: i   for i,lab in enumerate(KF.sY)}
        for it in range(0,KF.nt-1):    
            if np.mod(it,500) == 0:
                print('Time step %8.0f t=%10.3f  WS=%4.1f Thrust=%.1f' % (it,KF.time[it],x[7],x[4]))
            # --- "Measurements"
            y  = KF.Y[:,it]

            # --- KF predictions
            u=KF.U_clean[:,it]
            x,P,_ = KF.estimateTimeStep(u,y,x,P,KF.Q,KF.R)

            # --- Estimate thrust and WS - Non generic code
            WS0=x[KF.iX['WS']]
            pitch     = y[KF.iY['pitch']]
            Qaero_hat = x[KF.iX['Qaero']]
            omega     = x[KF.iX['omega']]
            WS_hat = KF.wse.estimate(Qaero_hat, pitch, omega, WS0, relaxation = 0)
            Thrust = KF.wse.Thrust(WS_hat, pitch, omega)
            x[KF.iX['Thrust']] = Thrust
            x[KF.iX['WS']]     = WS_hat
            x[KF.iX['psi']]    = np.mod(x[KF.iX['psi']], 2*np.pi)

            # --- Store
            KF.X_hat[:,it+1] = x
            KF.Y_hat[:,it+1] = np.dot(KF.Yx,x) + np.dot(KF.Yu,u)

        KF.P = P


    def moments(KF):
        WT=KF.WT
        z_test = fastlib.ED_TwrGag(WT.ED) - WT.ED['TowerBsHt']
        EI     = np.interp(z_test, WT.Twr.s_span, WT.Twr.EI[0,:])
        kappa  = np.interp(z_test, WT.Twr.s_span, WT.Twr.PhiK[0][0,:])
        qx    = KF.X_hat[KF.iX['ut1']]
        KF.M_sim = [qx*EI[i]*kappa[i]/1000 for i in range(len(z_test))]                 # in [kNm]
        KF.M_ref = [KF.df['TwHt{:d}MLyt'.format(i+1)].values for i in range(len(z_test)) ] # in [kNm]
        return KF.M_sim, KF.M_ref

    def export(KF,OutputFile):
        M=np.column_stack([KF.time]+[KF.X_clean[j,:] for j,_ in enumerate(KF.sX)])
        M=np.column_stack([M]+[KF.X_hat  [j,:] for j,_ in enumerate(KF.sX)])
        M=np.column_stack([M]+[KF.Y      [j,:] for j,_ in enumerate(KF.sY)])
        M=np.column_stack([M]+[KF.Y_hat  [j,:] for j,_ in enumerate(KF.sY)])
        M=np.column_stack([M]+KF.M_ref)
        M=np.column_stack([M]+KF.M_sim)
        header='time'+','
        header+=','.join([s+'_ref' for s in KF.sX])+','
        header+=','.join([s+'_est' for s in KF.sX])+','
        header+=','.join([s+'_ref' for s in KF.sY])+','
        header+=','.join([s+'_est' for s in KF.sY])+','
        header+=','.join(['My_ref{:d}'.format(j) for j,_ in enumerate(KF.M_ref)])+','
        header+=','.join(['My_est{:d}'.format(j) for j,_ in enumerate(KF.M_sim)])
        np.savetxt(OutputFile,M,delimiter=',',header=header)



