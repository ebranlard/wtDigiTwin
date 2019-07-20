import numpy as np
    

import scipy.optimize as sciopt
# res=sciopt.minimize_scalar(lambda k:np.abs(k*(1+k/(4*lambda_r**2))-Ct), bounds=[0,1.8], method='bounded')
# res=sciopt.minimize_scalar(lambda k:np.abs(k*(1+k/(4*lambda_r**2))-Ct), bounds=[0,1.8], method='bounded')

def BendingModesUniformBeam(Type,EI,rho,A,L,w=None,x=None,Mtop=0,norm='tip_norm',nModes=4):
    """
    returns Mode shapes and frequencies for a uniform beam

    References:
      Inman : Engineering variation
    
    Author: E. Branlard"""
    if x is None or len(x)==0:
        x = np.linspace(0,L,100)
    else:
        x = p.x
        if np.amax(x) != L:
            raise Exception('Max of x should be equal to L')

    freq = None
    ModesU = None
    ModesV = None
    ModesK = None
    # Dimensionless spanwise position
    x0 = x / L
    s = Type.split('-')
    if s[0].lower()=='transverse':
        if s[1].lower()=='unloaded':
            # --- "Theory" (clamped-free, vertical, no weight)
            # See Inman, p.335 or Nielsen1 p129
            if 'transverse-unloaded-clamped-free' == (Type.lower()):
                # NOTE: cosh(beta_n)cos(beta_n) =-1
                #    sigma_n = [ np.sinh(beta_n) - sin(beta_n) ]/[cosh(beta_n) + cos(beta_n)]
                #    for j>5, a good approx is B(j) = (2*j-1)np.pi/2  and S(j)=1;
                #B  = [1.87510407, 4.69409113, 7.85475744,10.99554073,14.13716839, (2*6-1)*np.pi/2];
                #S  = [0.734095514 1.018467319 0.999224497 1.000033553 0.999998550 1];
                B = np.zeros(nModes)
                for i in np.arange(len(B)):
                    B[i] = sciopt.fsolve(lambda x: 1 + np.cosh(x) * np.cos(x), (2*(i+1)-1)*np.pi/2)
            else:
                if 'transverse-unloaded-topmass-clamped-free' == (Type.lower()):
                    # The geometrical stiffning is not accounted for here
                    if Mtop is not None:
                        raise Exception('Please specify value for Mtop for %s',Type)
                    Mtop = p.Mtop
                    M = rho * A * L
                    B = np.zeros((1,nModes))
                    for i in np.arange(len(B)):
                        B[i] = sciopt.fsolve(lambda x: 1+np.cosh(x)*np.cos(x)-x*Mtop/M*(np.sin(x)*np.cosh(x)-np.cos(x)*np.sinh(x)),(2*(i+1)-1)*np.pi/2)
                else:
                    raise Exception('unknown type %s',Type)
            #S  = ( sinh(B)-sin(B) ) ./ ( cosh(B) + cos(B));  # Sigma
#C  = ( cosh(B)+cos(B) ) ./ ( sinh(B) + sin(B));  # Sigma
            SS = np.sinh(B) + np.sin(B)
            CC = np.cosh(B) + np.cos(B)
            # Frequency
            freq = (B / L) ** 2 / (2 * np.pi) * np.sqrt(EI / (rho * A))
            # --- Mode shapes
            print(freq)
            print(B.shape)
            print(x0.shape)
            ModesU = np.zeros((len(B),len(x0)))
            ModesV = np.zeros((len(B),len(x0)))
            ModesK = np.zeros((len(B),len(x0)))
            for i in np.arange(len(B)):
                ModesU[i,:] =            SS[i] * (np.cosh(B[i] * x0) - np.cos(B[i] * x0)) - CC[i] * (np.sinh(B[i] * x0) - np.sin(B[i] * x0))
                ModesV[i,:] = B[i]    * (SS[i] * (np.sinh(B[i] * x0) + np.sin(B[i] * x0)) - CC[i] * (np.cosh(B[i] * x0) - np.cos(B[i] * x0)))
                ModesK[i,:] = B[i]**2 * (SS[i] * (np.cosh(B[i] * x0) + np.cos(B[i] * x0)) - CC[i] * (np.sinh(B[i] * x0) + np.sin(B[i] * x0)))
                #  ModesU(i,:)  =        cosh(B[i]*x0)-cos(B[i]*x0) - S[i]*(sinh(B[i]*x0)-sin(B[i]*x0)) ;
                #  ModesV(i,:) = B[i]  *(sinh(B[i]*x0)+sin(B[i]*x0) - S[i]*(cosh(B[i]*x0)-cos(B[i]*x0)));
                #  ModesK(i,:) = B[i]^2*(cosh(B[i]*x0)+cos(B[i]*x0) - S[i]*(sinh(B[i]*x0)+sin(B[i]*x0)));
                #  ModesU(i,:)  =        cosh(B[i]*x0)-cos(B[i]*x0) - C[i]*(sinh(B[i]*x0)-sin(B[i]*x0)) ;
                #  ModesV(i,:) = B[i]  *(sinh(B[i]*x0)+sin(B[i]*x0) - C[i]*(cosh(B[i]*x0)-cos(B[i]*x0)));
                #  ModesK(i,:) = B[i]^2*(cosh(B[i]*x0)+cos(B[i]*x0) - C[i]*(sinh(B[i]*x0)+sin(B[i]*x0)));
        else:
            if s[1].lower()=='loaded':
                if 'transverse-loaded-clamped-free' == (Type.lower()):
                    if w is None:
                        w = A * rho
                    if L==0:
                        raise Exception('Please specify value for L for %s',Type)
                    B = np.array([1.875,4.694])
                    freq = (B / L) ** 2 / (2 * np.pi) * np.sqrt(EI / w)
                else:
                    raise Exception('unknown type %s',Type)
            else:
                raise Exception('Unknown %s',Type)
    ## Computation of derivatives if no analytical functions
    # # V=fgradient_regular(U(i,:),4,dx);
    # # K=fgradient_regular(V(i,:),4,dx);
    ## Going back to physical dimension
    x = x0 * L
    ModesV = ModesV/L
    ModesK = ModesK/L**2
    ## Normalization of modes
    if norm=='tip_norm':
        for i in np.arange(nModes):
            fact = 1 / ModesU[i,-1]
            ModesU[i,:] = ModesU[i,:] * fact
            ModesV[i,:] = ModesV[i,:] * fact
            ModesK[i,:] = ModesK[i,:] * fact
    else:
        raise Exception('Norm not implemented or incorrect: `%s`'%norm)

    return freq,x,ModesU,ModesV,ModesK


if __name__=='__main__':
    import matplotlib.pyplot as plt
    pass
    L = 100
    EI = 1868211939147.334
    m = 8828.201296825122
    _,x,U,_,_ = BendingModesUniformBeam('transverse-unloaded-clamped-free',EI,m,A=1,L=L)
    plt.figure
    plt.plot(x,U[0,:])
    plt.plot(x,U[1,:])
    plt.show()
