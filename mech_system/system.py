import numpy as np
from numpy.linalg import inv
from scipy.integrate import  solve_ivp #odeint


# --------------------------------------------------------------------------------}
# --- Functions independent of class 
# --------------------------------------------------------------------------------{
def StateMatrix(M,C,K):
    A = np.zeros((2*len(M),2*len(M)))
    A[:len(M),len(M):]  =  np.identity(len(M))
    A[len(M):,:len(M)]  = -np.dot(inv(M),K)
    A[len(M):,len(M):]  = -np.dot(inv(M),C)
    return A


def vec_interp(t,vTime,vF):
    """ Interpolate vector known at discrete values (vTime, vF) to a given time `t` """
    F    = np.zeros(vF.shape[0])
    for iDOF,F_DOF in enumerate(vF):
        F[iDOF] = np.interp(t,vTime,F_DOF)
    return F

def B_interp(t,M,vTime,vF):
    """ Interpolate B-vector from loads known at discrete values (vTime, vF) at a given time `t` """
    nDOF=len(vF)
    B = np.zeros((2*nDOF,1))
    F = vec_interp(t,vTime,vF)
    B[nDOF:,0] = np.dot(inv(M),F)
    return B

def B_reg(t,M,F):
    """ Return B vector from loads at time t and mass matrix """
    nDOF=len(F)
    B = np.zeros((2*nDOF,1))
    B[nDOF:,0] = np.dot(inv(M),F)
    return B

def dxdt(q, t, A, M, vTime, vF): 
    B = B_interp(t, M, vTime, vF)
    dq = np.dot(A,q)+B
    return dq
 
def odefun(t, dq):
    return dxdt(dq, t, A, vTime, vF)



# --------------------------------------------------------------------------------}
# --- Lineary system
# --------------------------------------------------------------------------------{
class LinearSystem():
    def __init__(self,A,B,C=None,D=None,q0=None):
        self.A=A
        self.B=B
        self.C=C
        self.D=C
        self.setStateInitialConditions(q0)

    def setStateInitialConditions(self,q0=None):
        self.q0 = np.zeros(self.nStates)
        if q0 is not None:
            if len(q0)!=self.nStates:
                raise Exception('Wrong dimension for q0 ({} instead of {} )'.format(len(q0),self.nStates))
            self.q0 = q0

    def setInputTimeSeries(self,vTime,vU):
        vTime = np.asarray(vTime)
        vU    = np.asarray(vU)
        if vU.shape[0]!=self.nInputs:
            raise Exception('Wrong first dimension for Inputs time series ({} instead of {} )'.format(vU.shape[0],self.nInputs))
        if vU.shape[1]!=len(vTime):
            raise Exception('Second dimension of Input time series does not match time dimension ({} instead of {} )'.format(vU.shape[1],len(vTime)))

        self._input_ts = vU
        self._time_ts  = vTime

    def setInputFunction(self,fn):
        self._input_fn = fn

    def Inputs(self,t,x=None):
        if hasattr(self,'_input_ts'):
            return vec_interp(t,self._time_ts,self._input_ts)
        elif hasattr(self,'_input_fn'):
            return self._input_fn(t,x)
        else:
            raise NotImplementedError('Please specify a time series of inputs first')

    def RHS(self,t,q=None):
        return np.dot(self.A, q) + np.dot(self.B, self.Inputs(t,q))
        #return np.dot(self.A, q) 

    def integrate(self, t_eval, method='RK4', y0=None, **options):
        #
        if y0 is not None:
            self.setStateInitialConditions(y0)
            
        odefun = lambda t, q : np.dot(self.A, q) + np.dot(self.B, self.Inputs(t,q))
        #odefun = lambda t, q : self.RHS(t,q)

        res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=self.q0, t_eval=t_eval, method=method, vectorized=False, **options)   

        if self.nOuputs>0:
            print('>>> Do something for outputs')


        return res

    @property
    def nStates(self):
        return self.A.shape[0]

    @property
    def nInputs(self):
        return self.B.shape[1]

    @property
    def nOuputs(self):
        if self.C is not None:
            return self.C.shape[1]
        else:
            return 0

    def __repr__(self):
        s='<MechSystem object>\n'
        s+='Read-only attributes:\n'
        s+=' - nState:{} \n'.format(self.nStates)
        if hasattr(self,'_force_ts'):
            if len(self._time_ts)>1:
                dt=self._time_ts[1]-self._time_ts[0]
            else:
                dt=np.nan
            s+='Force time series \n'
            s+=' - Time: [{} ... {}],  dt: {}, n: {} \n'.format(self._time_ts[0],self._time_ts[-1],dt,len(self._time_ts))
            s+=' - Force t0  : {} \n'.format(self._force_ts[:,0])
            s+=' - Force tend: {} \n'.format(self._force_ts[:,-1])
        s+='Attributes:\n'
        s+=' - A: State-State Matrix  \n'
        s+=str(self.A)+'\n'
        s+=' - B: State-Input Matrix  \n'
        s+=str(self.B)+'\n'
        s+=' - q0: Initial conditions (state) \n'
        s+=str(self.q0)+'\n'
        return s

# --------------------------------------------------------------------------------}
# --- Mech System Class
# --------------------------------------------------------------------------------{
class MechSystem():
    def __init__(self,M,C,K,x0=None,xdot0=None):
        self.M=M
        self.C=C
        self.K=K
        self.setInitialConditions(x0,xdot0)

    def setInitialConditions(self,x0,xdot0):
        self.q0 = np.zeros(2*self.nDOF)
        if x0 is not None:
            if len(x0)!=self.nDOF:
                raise Exception('Wrong dimension for x0 ({} instead of {} )'.format(len(x0),self.nDOF))
            self.q0[:self.nDOF] =  x0
        if xdot0 is not None:
            if len(xdot0)!=self.nDOF:
                raise Exception('Wrong dimension for xdot0 ({} instead of {} )'.format(len(xdot0),self.nDOF))
            self.q0[self.nDOF:] =  xdot0


    def setStateInitialConditions(self,q0):
        q0 = np.asarray(q0).ravel()
        if len(q0)!=self.nStates:
            raise Exception('Wrong dimension for q0 ({} instead of {} )'.format(len(q0),2*self.nStates))
        self.q0 = q0

    def setForceTimeSeries(self,vTime,vF):
        vTime = np.asarray(vTime)
        vF    = np.asarray(vF)
        if vF.shape[0]!=self.nDOF:
            raise Exception('Wrong first dimension for Force time series ({} instead of {} )'.format(vF.shape[0],self.nDOF))
        if vF.shape[1]!=len(vTime):
            raise Exception('Second dimension of Force time series does not match time dimension ({} instead of {} )'.format(vF.shape[1],len(vTime)))

        self._force_ts = vF
        self._time_ts  = vTime

    def setForceFunction(self,fn):
        self._force_fn = fn

    def Force(self,t,x=None):
        if hasattr(self,'_force_ts'):
            return vec_interp(t,self._time_ts,self._force_ts)
        elif hasattr(self,'_force_fn'):
            return self._force_fn(t,x)
        else:
            raise NotImplementedError('Please specify a time series of force first')

    def B(self,t,q):
        if hasattr(self,'_force_ts'):
            return  B_interp(t,self.M,self._time_ts,self._force_ts)
        elif hasattr(self,'_force_fn'):
            F = self._force_fn(t,q)
            return  B_reg(t,self.M, F)
        else:
            raise NotImplementedError('Please specify a time series of force first')


    def integrate(self,t_eval, method='RK45', y0=None, **options):
        #
        if y0 is not None:
            self.setStateInitialConditions(y0)
            
        A=self.A
        odefun = lambda t, q : np.dot(A,q)+self.B(t,q)

        res = solve_ivp(fun=odefun, t_span=[t_eval[0], t_eval[-1]], y0=self.q0, t_eval=t_eval, method=method, vectorized=True, **options)   
        return res

    @property
    def A(self):
        return StateMatrix(self.M,self.C,self.K)

    @property
    def nDOF(self):
        return len(self.M)

    @property
    def nStates(self):
        return 2*self.nDOF

    def __repr__(self):
        s='<MechSystem object>\n'
        s+='Read-only attributes:\n'
        s+=' - nDOF: {}, nState:{} \n'.format(self.nDOF,self.nStates)
        if hasattr(self,'_force_ts'):
            if len(self._time_ts)>1:
                dt=self._time_ts[1]-self._time_ts[0]
            else:
                dt=np.nan
            s+='Force time series \n'
            s+=' - Time: [{} ... {}],  dt: {}, n: {} \n'.format(self._time_ts[0],self._time_ts[-1],dt,len(self._time_ts))
            s+=' - Force t0  : {} \n'.format(self._force_ts[:,0])
            s+=' - Force tend: {} \n'.format(self._force_ts[:,-1])
        s+='Attributes:\n'
        s+=' - M: Mass Matrix  \n'
        s+=str(self.M)+'\n'
        s+=' - C: Damping Matrix  \n'
        s+=str(self.C)+'\n'
        s+=' - K: Stiffness Matrix  \n'
        s+=str(self.K)+'\n'
        s+=' - q0: Initial conditions (state) \n'
        s+=str(self.q0)+'\n'
        return s
