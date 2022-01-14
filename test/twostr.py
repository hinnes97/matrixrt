import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

ds = xr.open_dataset('../fortran/output/output.nc')

ds.Tf.plot(y='pfull')
plt.gca().invert_yaxis()
plt.gca().set_yscale('log')

plt.savefig('TS.pdf')
