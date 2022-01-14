import numpy as np
import xarray as xr
import plot_utils as pu
import matplotlib.pyplot as plt

# Constants
sb = 5.67e-8
tlwinf = 256
tswinf = 256*0.04
S0 = 1368/4
Fint = 0#sb*70**4

ds = xr.open_dataset('../fortran/output.nc')

# Calculate pressure and optical depths
pf = ds['pfull']
ph = ds['pedge']
tau_lw = tlwinf*(pf/ph[-1])
tau_sw = tswinf*(pf/ph[-1])
dtlw_dp = tlwinf/ph[-1]
dtsw_dp = tswinf/ph[-1]

t_an = pu.gg_analytic(tau_lw, tau_sw, dtlw_dp, dtsw_dp, S0, Fint)
input_t = xr.DataArray(t_an, coords={'pfull':pf}, name='Tf')
input_te = xr.DataArray(np.zeros_like(ph), coords={'pedge':ph}, name='Te')
ds_input = xr.Dataset(data_vars=dict(Tf=input_t, Te=input_te))
ds_input.to_netcdf('test_input.nc')

ds.Tf.plot(y='pfull', label='T computed')
plt.plot(t_an, pf, label='T analytic')

plt.gca().invert_yaxis()
plt.gca().set_yscale('log')
plt.legend()


fl_dn_test = pu.gg_flux_down(ds['Tf'].data, tau_lw.data)
fl_up_test = pu.gg_flux_up(ds['Tf'].data, tau_lw.data, sb*ds['Tf'].data[-1]**4)

fl_dn_test2 = pu.gg_flux_down(t_an, tau_lw.data)
fl_up_test2 = pu.gg_flux_up(t_an, tau_lw.data, sb*t_an[-1]**4)

fl_dn_test3 = pu.flux_down_shortchar(ds['Tf'].data, tau_lw.data)
fl_up_test3 = pu.flux_up_shortchar(ds['Tf'].data, tau_lw.data, bc=None)

fl_up, fl_dn = pu.flux_up_down(t_an, tau_lw, tau_sw, Fint, dtlw_dp, dtsw_dp, S0)
plt.figure()
ds['fdn'].plot(y='pedge', label='Fdn compute')
ds['fup'].plot(y='pedge', label='Fup compute')
plt.plot(fl_dn, pf, label='Fdn analytic')
plt.plot(fl_up, pf, label='Fup analytic')
plt.plot(fl_dn_test, pf, label='dn_test')
plt.plot(fl_up_test, pf, label='up_test')

plt.plot(fl_dn_test2, pf, label='dn_test2')
plt.plot(fl_up_test2, pf, label='up_test2')

plt.plot(fl_dn_test3, pf, label='dn_an_sc')
plt.plot(fl_up_test3, pf, label='up_an_sc')

plt.legend()
plt.gca().invert_yaxis()
plt.gca().set_yscale('log')


plt.show()


