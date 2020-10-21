import numpy as np
import pickle
import glob
import os
import weio



class FASTPeriodicOP(object):
    """ Class for a set of *.lin files, all assumed to be for the same periodic operating point"""
    def __init__(self,prefix,nLin=None):
        if nLin is None:
            linfiles= glob.glob(prefix + '*.*.lin')
            self.nLinTimes = len(linfiles)
        else:
            self.nLinTimes = nLin

        print(prefix, self.nLinTimes)

        self.prefix   = prefix
        self.Data     = []
        self.vAzim    = []
        self.vWS       = []
        self.vPitch    = []
        self.vRotSpeed = []
        self.vBu = []
        for i in np.arange(self.nLinTimes):
            linfilename= prefix+'.'+str(i+1)+'.lin'
            print(linfilename)
            if not os.path.exists(linfilename):
                print('Linearization file missing: ',linfilename)
            linfile=weio.read(linfilename)
            df=linfile.toDataFrame()
            self.Data.append(linfile)
            #self.A=lin['A']
            #B=linfile['B']
            #u=linfile['u']
            #self.C=lin['C']
            #self.D=lin['D']
            try:
                self.vWS.append(df['u']['WS_[m/s]'][0])
            except:
                print('Wind speed not found in input, assuming 0m/s')
                self.vWS.append(0)
            self.vRotSpeed.append(linfile.RotSpeed)
            self.vAzim.append(linfile.Azimuth)
            self.vPitch.append(df['u']['B1pitch_[rad]'][0]*180/np.pi)

        self.WS       = np.mean(self.vWS)
        self.Pitch    = np.mean(self.vPitch)
        self.RotSpeed = np.mean(self.vRotSpeed)

        self.x = df['x']
        self.y = df['y']
        self.u = df['u']
        try:
            self.EDdescr = linfile['EDDOF']
        except:
            print('EDDOF not available. A special version of OpenFAST is required.')


class FASTLin(object):
    """ Class for linearization data for different operating points (typically Campbell) """
    def __init__(self,folder='./', prefix='',nLin=None):

        fstfiles= glob.glob(folder + prefix + '*.*.lin')
        Sim_Prefix=np.unique(['.'.join(f.split('.')[:-2]) for f in fstfiles])
        nSim      = len(Sim_Prefix)
        # --- Read period operating points
        print('Reading linearizations for {} operating points'.format(nSim))
        self.OP_Data=[FASTPeriodicOP(pref,nLin=nLin) for pref in Sim_Prefix]
        # --- Sort by wind speed
        Isort = np.argsort(self.WS)
        self.OP_Data  = [self.OP_Data[i] for i in Isort]

        if self.MaxNLinTimes>1:
            IBad = [i for i in np.arange(nSim) if self.nLinTimes[i]<self.MaxNLinTimes and self.OP_Data[i].WS>0]
            if len(IBad)>0: 
                print('>>> The following simulations have insufficient number of data points:')
                for i in IBad:
                    print(self.OP_Data[i].prefix, self.OP_Data[i].nLinTimes)
            self.OP_Data = [self.OP_Data[i] for i in np.arange(nSim) if i not in IBad]

    @property
    def WS(self):
        return np.array([sim.WS for sim in self.OP_Data])

    @property
    def nLinTimes(self):
        return np.array([sim.nLinTimes for sim in self.OP_Data])

    @property
    def MaxNLinTimes(self):
        return np.max(self.nLinTimes)

    @property
    def nOP(self):
        return len(self.OP_Data)

    @property
    def xdescr(self):
        return self.OP_Data[0].x.columns.values
    @property
    def ydescr(self):
        return self.OP_Data[0].y.columns.values
    @property
    def EDdescr(self):
        return self.OP_Data[0].EDdescr
    @property
    def udescr(self):
        return self.OP_Data[0].u.columns.values
    @property
    def xop_mean(self):
        return np.mean(np.abs(np.array([op.x.values for op in self.OP_Data])),axis=0)
    @property
    def uop_mean(self):
        return np.mean(np.abs(np.array([op.u.values for op in self.OP_Data])),axis=0)
    @property
    def uop_mean(self):
        return np.mean(np.abs(np.array([op.u.values for op in self.OP_Data])),axis=0)

    @property
    def yop_mean(self):
        return np.mean(np.abs(np.array([op.y.values for op in self.OP_Data])),axis=0)

    def stats(self,matName,WS=None):
        if WS is None:
            WS = self.WS
            nOP=self.nOP
        else:
            nOP=len(WS)
        print('Returning stats for WS:',WS)
        M_mean=[]

        shape= self.OP_Data[0].Data[0][matName].shape

        M_all       = np.zeros( (nOP, self.MaxNLinTimes, shape[0],shape[1]))
        M_mean_perWS= np.zeros( (nOP, shape[0],shape[1]))
        M_std_perWS = np.zeros( (nOP, shape[0],shape[1]))

        # loop on operating points (e.g. WS)
        ii=0
        for iop, op in enumerate(self.OP_Data):
            if op.WS in WS:
                # Loop on linearization times (e.g. Azimuth)
                for iTimes in np.arange(self.MaxNLinTimes):
                    if op.nLinTimes==1:
                        M_all[ii,iTimes,:,:]=op.Data[0][matName]
                    else:
                        M_all[ii,iTimes,:,:]=op.Data[iTimes][matName]

                M_mean_perWS[ii,:,:] = np.mean(M_all[ii,:,:,:],axis=0)
                M_std_perWS [ii,:,:]  = np.std(M_all[ii,:,:,:],axis=0)
                ii+=1

        M_mean    = np.mean( M_mean_perWS, axis=0 )
        M_stdWS   = np.std ( M_mean_perWS, axis=0 ) # How much elements vary with wind speed
        M_stdAzim = np.mean( M_std_perWS , axis=0)  # How much elements vary due to azimuth

        return M_mean, M_mean_perWS, M_stdAzim, M_stdWS, M_all


    def save(self,filename):
        with open(filename,'wb') as f:
            pickle.dump(self,f)


#     def full_linear_model

