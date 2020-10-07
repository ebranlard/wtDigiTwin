import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# Local 
import weio


from mech_system.system import MechSystem

# --- Options
damping = True
damping = False

# --- Paramters
JzzA=3699.225
mp=443.9
kp=100
cp=1e4
l=5




df=weio.read('input_time_series_interface_joint.txt').toDataFrame()
print(df.columns)
accx=df['C13'].values
time=df['C0'].values
x_top = df['C1'].values


if damping:
    df_SD=weio.read('PendulumDamp.SD_ref.out').toDataFrame()
    df_OM=weio.read('PendulumDamp_OpenModelica.csv').toDataFrame()
    print(df_OM.keys())
    t_OM        = df_OM['Time [s]'].values
    x_OM_top    = df_OM['X1 (top) [m]'].values
    x_OM_bottom = df_OM['X3 (bottom) [m]'].values
    theta_OM    = df_OM['Rotated angle [deg]'].values
else:
    df_SD=weio.read('Pendulum.SD.out').toDataFrame()
t_SD = df_SD['Time_[]'].values
x_SD = df_SD['M2N2TDxss_[m]'].values
x_SD_top = df_SD['M1N1TDxss_[m]'].values


def Fx(t,x):
    theta=x[0,0]
    Force = np.zeros((1,1))
    accx_t = np.interp(t,time,accx)
    Force[0,0]= mp*l/2*accx_t* np.cos(theta)
    return Force

Fx_lin =  mp*l/2*accx
Fx_lin = Fx_lin.reshape((1,len(time)))


M=np.array([[JzzA]])
K=np.array([[kp]])
if damping:
    C=np.array([[cp]])
else:
    C=np.array([[0]])

sys=MechSystem(M,C,K)
# sys.setForceTimeSeries(time,Fx_lin)
sys.setForceFunction(Fx)
res=sys.integrate(time, method='RK45') # **options):

theta     = res.y[0,:]
theta_dot = res.y[1,:]

print(len(theta), len(x_SD_top), len(time))
x_end=-l*np.sin(theta) + x_top

print(res.keys())


print(sys)



import matplotlib.pyplot as plt
fig,ax = plt.subplots(1, 1, sharey=False, figsize=(6.4,4.8)) # (6.4,4.8)
fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.11, hspace=0.20, wspace=0.20)
# ax.plot(time,accx, label='acc')
# ax.plot(res.t,res.y[0,:], label='theta')
# ax.plot(res.t,res.y[1,:], label='theta dot')
if damping:
    ax.plot(res.t,x_end     , label='x')
    ax.plot(t_SD ,x_SD      , label='x_SD')
    ax.plot(t_OM ,x_OM_bottom, label='x_OM')

#     ax.plot(t_OM ,x_OM_top+np.sin(-theta_OM*np.pi/180)*l, label='x_OM')
#     ax.plot(t_SD ,x_SD_top  ,'--', label='x_SD_in')

#     ax.plot(res.t, res.y[0 ,:]*180/np.pi   ,label = 'theta')
#     ax.plot(t_OM ,theta_OM                ,label = 'theta_OM')

else:
    ax.plot(res.t,x_end     , label='x')
    ax.plot(t_SD ,x_SD      , label='x_SD')
#     ax.plot(t_SD ,x_SD_top  ,'--', label='x_SD_in')


ax.set_xlabel('')
ax.set_ylabel('')
ax.legend()
ax.tick_params(direction='in')
plt.show()


