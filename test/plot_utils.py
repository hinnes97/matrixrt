import numpy as np
import xarray as xr
from scipy.integrate import quad
import scipy.interpolate as spi

sb = 5.67e-8

def SW_down(tau_sw):
    return np.exp(-tau_sw)

def SW_down_diff(tau_sw, dtlw_dp, dtsw_dp):
    return -np.exp(-tau_sw)*dtsw_dp/dtlw_dp

def SW_down_integrand(tau_sw, dtlw_dp, dtsw_dp):
    return np.exp(-tau_sw)*dtlw_dp/dtsw_dp

def gg_analytic(tau_lw, tau_sw, dtlw_dp, dtsw_dp, S0, Fint):

    int_term = np.zeros_like(tau_lw)
    for i,tsw in enumerate(tau_sw):
        int_term[i] = S0*quad(SW_down_integrand, 0, tsw, args=(dtlw_dp, dtsw_dp))[0]

    diff_term = S0*SW_down_diff(tau_sw,dtlw_dp,dtsw_dp)

    twosigt4 = S0 - diff_term + int_term + Fint*(tau_lw + 1)

    return (twosigt4/2/sb)**0.25

def gg_flux_integrand(tau_lw_dummy, tau_lw, interp_fn):
    """ 
    Function can be used for both downwards and upwards integrands
    thanks to np.absolute() function
    """
    T = interp_fn(tau_lw_dummy)
    return sb*T**4*np.exp(-np.absolute((tau_lw - tau_lw_dummy)))

def gg_flux_down(T, tau_lw):
    interp_fn = spi.interp1d(tau_lw, T, fill_value='extrapolate')

    flux_down = np.zeros_like(tau_lw)

    for i,tlw in enumerate(tau_lw):
        flux_down[i] = quad(gg_flux_integrand, 0, tlw, args=(tlw, interp_fn))[0]

    return flux_down
    
def gg_flux_up(T, tau_lw, BC):
    interp_fn = spi.interp1d(tau_lw, T, fill_value='extrapolate')

    flux_up = np.zeros_like(tau_lw)

    # BC is upward boundary condition

    for i,tlw in enumerate(tau_lw):
        flux_up[i] = BC*np.exp(-(tau_lw[-1] - tlw)) + quad(gg_flux_integrand, tlw, tau_lw[-1], args=(tlw, interp_fn))[0]

    return flux_up

def flux_up_down(T, tau_lw, tau_sw, Fint, dtlw_dp, dtsw_dp, S0):
    flux_up =  0.5*(2*sb*T**4 + S0*SW_down(tau_sw) + Fint + S0*SW_down_diff(tau_sw, dtlw_dp, dtsw_dp))
    flux_down = flux_up - Fint - S0*SW_down(tau_sw)

    return flux_up, flux_down

def flux_down_shortchar(T, tau):
    # Testing this function
    #mus  = np.array([0.21132487, 0.78867513])
    #wus = np.array([0.5,0.5])

    mus = np.array([1])
    wus = np.array([1])
    wms = mus*mus

    S = sb*T**4
    #tau = ir_tau_water(atm)
    
    I_m_temp = np.zeros_like(T)
    I_m = np.zeros_like(T)

    print(I_m.shape, tau.shape, T.shape)
    for mu,wm in zip(mus,wms):
        
        dtau = np.diff(tau)/mu
        # Define constants for interpolation of source
        lc = linear_coeff(dtau)
        pc = parab_coeff(dtau)
    
        # Parabolic for every term but last (linear)
        dI_m = np.zeros((len(I_m) - 1))
        #print(S[:-1].shape, S[1:].shape, (S[:-1] + S[1:]).shape)
        for i in range(len(I_m)-1):
            dI_m[i] = lc.a_m[i]*S[i] + lc.b_m[i]*S[i+1]
        #dI_m = lc.a_m*S[:-1] + lc.b_m*S[1:]
        
        #dI_m[:-1] = pc.a_m*S[:-2] + pc.b_m*S[1:-1] + pc.g_m*S[2:]
        #dI_m[-1] = lc.a_m[-1]*S[-2] + lc.b_m[-1]*S[-1]
        
        for k in range(1,len(I_m)):
            I_m_temp[k] = I_m_temp[k-1]*np.exp(-dtau[k-1]) + dI_m[k-1]

        I_m = I_m + wm*I_m_temp
    return I_m

def flux_up_shortchar(T, tau, bc):
    mus = np.array([1])
    wus = np.array([1])
    
    wms = mus*mus

    I_p = np.zeros_like(T)
    I_p_temp = np.zeros_like(T)
    for mu,wm in zip(mus,wms):
        
        dtau = np.diff(tau)/mu

        lc = linear_coeff(dtau)
        pc = parab_coeff(dtau)

        S = sb*T**4
    
        if bc == 'flux':
            I_p[-1] = atm.Fint + atm.f_down[-1]
        else:
            I_p[-1] = S[-1]
        
        dI_p = np.zeros((len(I_p) - 1))

        for i in range(len(I_p) - 1):
            dI_p[i] = lc.b_p[i]*S[i] + lc.g_p[i]*S[i+1]
    
        #dI_p[1:] = pc.a_p*S[:-2]+ pc.b_p*S[1:-1] + pc.g_p*S[2:]
        dI_p[0] = lc.b_p[0]*S[0] + lc.g_p[0]*S[1]

        for k in range(len(I_p)-2,-1,-1):
            I_p_temp[k] = I_p_temp[k+1]*np.exp(-dtau[k]) + dI_p[k]

        I_p = I_p + wm*I_p_temp

    return I_p

class linear_coeff:
    def __init__(self, dtau):
        
        e0 = 1 - np.exp(-dtau)
        e1 = dtau - e0
        
        self.a_m = e0 - e1/dtau
        self.b_m = e1/dtau

        self.b_p = e1/dtau
        self.g_p = e0 - e1/dtau

class parab_coeff:
    def __init__(self, dtau):

        e0 = 1 - np.exp(-dtau)
        e1 = dtau - e0
        e2 = dtau**2 - 2*e1

        self.a_m = e0[:-1] + (e2[:-1] - (dtau[1:] + 2*dtau[:-1])*e1[:-1])/(dtau[:-1]*(dtau[1:] + dtau[:-1]))
        self.b_m = ((dtau[1:] + dtau[:-1])*e1[:-1] - e2[:-1])/(dtau[:-1]*dtau[1:])
        self.g_m = (e2[:-1] - dtau[:-1]*e1[:-1])/(dtau[1:]*(dtau[1:] + dtau[:-1]))

        self.a_p = (e2[1:] - dtau[1:]*e1[1:])/(dtau[:-1]*(dtau[1:] + dtau[:-1]))
        self.b_p = ((dtau[1:] + dtau[:-1])*e1[1:] - e2[1:])/(dtau[:-1]*dtau[1:])
        self.g_p = e0[1:] + (e2[1:] - (dtau[:-1] + 2*dtau[1:])*e1[1:])/(dtau[1:]*(dtau[1:] + dtau[:-1]))

