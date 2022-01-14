import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

def analytic(p):
    taulwinf = 256
    gamma=0.04
    Fint = 1
    S = 1368/4
    
    taulw = taulwinf*p/p[-1]

    twosigt4 = Fint*(taulw + 1) + S*( 1 + 1/gamma + (gamma - 1/gamma)*np.exp(-gamma*taulw))

    return (twosigt4/2/5.67e-8)**0.25
#dat = np.genfromtxt('./PT_100000_300_H2O_CO2_Dry_100steps.dat')
#p = dat[:,0]
#T = dat[:,1]

#ds = xr.open_dataset('../fortran/output_new_spectrum_rayleigh.nc')
ds_ts =xr.open_dataset('../fortran/output_broad.nc')
ds_ts2 = xr.open_dataset('../fortran/output_broad2.nc')


#ds.Tf.plot(y='pfull')
ds_ts.Tf.plot(y='pfull')
ds_ts2.Tf.plot(y='pfull')
plt.plot(analytic(ds_ts.pfull), ds_ts.pfull)
#plt.plot(T,p)
#plt.plot(ds.Tf[-1]*(ds.pfull/ds.pfull[-1])**(2./7.), ds.pfull)
         
plt.gca().set_yscale('log')
#plt.gca().set_xscale('log')
plt.gca().invert_yaxis()

plt.figure()
ds_ts.q.plot(y='pedge')
plt.gca().invert_yaxis()
plt.gca().set_yscale('log')
plt.gca().set_xscale('log')

plt.show()

