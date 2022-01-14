import xarray as xr
import matplotlib.pyplot as plt

root_dir = '../fortran/output/'

ds = xr.open_dataset(root_dir+'output.nc')

ds['q'].plot(y='pedge')

plt.gca().set_yscale('log')
plt.gca().set_xscale('log')
plt.gca().invert_yaxis()

plt.figure()
ds['Te'].plot(y='pedge')
plt.gca().set_yscale('log')
plt.gca().set_xscale('log')
plt.gca().invert_yaxis()

plt.show()
