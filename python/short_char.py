import numpy as np
from moistad import satq
import phys

sig = 5.670374419e-8
taulwinf = 256
tauswinf = 256*0.04
ps = 1e6
Fint = 70
dry = phys.H2
wet = phys.H2O
q0 = 0.1

def ir_flux_down(Te, atm):
    # Testing this function
    #mus  = np.array([0.21132487, 0.78867513])
    #wus = np.array([0.5,0.5])

    mus = np.array([1])
    wus = np.array([1])
    wms = mus*mus

    pe = atm.pe
    
    tau = ir_tau(atm, taulwinf)
    S = sig*Te**4
    #tau = ir_tau_water(atm)
    
    I_m_temp = np.zeros_like(Te)
    I_m = np.zeros_like(Te)
    for mu,wm in zip(mus,wms):
        
        dtau = np.diff(tau)/mu

        # Define constants for interpolation of source
        lc = linear_coeff(dtau)
        pc = parab_coeff(dtau)
    
        # Parabolic for every term but last (linear)
        dI_m = np.zeros((len(I_m) - 1))
        dI_m = lc.a_m*S[:-1] + lc.b_m*S[1:]
    
        #dI_m[:-1] = pc.a_m*S[:-2] + pc.b_m*S[1:-1] + pc.g_m*S[2:]
        dI_m[-1] = lc.a_m[-1]*S[-2] + lc.b_m[-1]*S[-1]
    
        for k in range(1,len(I_m)):
            I_m_temp[k] = I_m_temp[k-1]*np.exp(-dtau[k-1]) + dI_m[k-1]

        I_m = I_m + wm*I_m_temp
    return I_m

def ir_flux_up(Te, atm, bc=None):

    #mus  = np.array([0.21132487, 0.78867513])
    #wus = np.array([0.5,0.5])
    mus = np.array([1])
    wus = np.array([1])
    
    wms = mus*mus

    pe = atm.pe
    tau = ir_tau(atm, taulwinf)
    #tau = ir_tau_water(atm)
    I_p = np.zeros_like(Te)
    I_p_temp = np.zeros_like(Te)
    for mu,wm in zip(mus,wms):
        
        dtau = np.diff(tau)/mu

        lc = linear_coeff(dtau)
        pc = parab_coeff(dtau)

        S = sig*Te**4
    
        if bc == 'flux':
            I_p[-1] = atm.Fint + atm.f_down[-1]
        else:
            I_p[-1] = S[-1]
        
        dI_p = np.zeros((len(I_p) - 1))
        dI_p = lc.b_p*S[:-1] + lc.g_p*S[1:]
    
        #dI_p[1:] = pc.a_p*S[:-2]+ pc.b_p*S[1:-1] + pc.g_p*S[2:]
        dI_p[0] = lc.b_p[0]*S[0] + lc.g_p[0]*S[1]

        for k in range(len(I_p)-2,-1,-1):
            I_p_temp[k] = I_p_temp[k+1]*np.exp(-dtau[k]) + dI_p[k]

        I_p = I_p + wm*I_p_temp

    return I_p

def sw_flux_down(atm):
    p = atm.pe
    S0 = atm.S0
    f = 1#0.8
    
    tau_sw = sw_tau(atm, tauswinf)
    tausw_2 = sw_tau(atm, 1)
    return S0*(f*np.exp(-tau_sw) + (1-f)*np.exp(-tausw_2))

def ir_tau(atm, tauinf):
    p = atm.pe
    ps = atm.ps
    tau_h2 = 180

    tau = np.zeros_like(p)
    dtau = np.diff(p)*atm.kappa/self.grav * (atm.q[1:]/2+atm.q[:-1]/2)
    tau = np.r_[0, np.cumsum(dtau)]

    return tau + atm.tau_0*p/ps
    #return tauinf*(p/ps)# + tau_h2*(p/ps)**2

def sw_tau(atm, tauinf):
    p = atm.pe
    ps = atm.ps
    kappa = 1.e-5
    tauinf = ps*atm.kappa/atm.grav
    return tauinf*(p/ps)[:]

def ir_tau_water(atm):
    p1 = 70
    #q = np.maximum(satq(T, p, dry, wet), q0)
    q = atm.q
    p = atm.pe
    dp = atm.dp
    tau0 = np.sum(q*dp)/p1

    return taulwinf*(p/ps) + tau0*(p/ps)**4
    
    
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
