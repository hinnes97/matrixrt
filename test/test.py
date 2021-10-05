import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

#dat = np.genfromtxt('./PT_100000_300_H2O_CO2_Dry_100steps.dat')
#p = dat[:,0]
#T = dat[:,1]

ds = xr.open_dataset('../fortran/output.nc')
ds_ts =xr.open_dataset('../fortran/output_timestep.nc')
ds.Tf.plot(y='pfull')
ds_ts.Tf.plot(y='pfull')
#plt.plot(T,p)
plt.gca().set_yscale('log')
plt.gca().invert_yaxis()
plt.show()

