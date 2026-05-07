import numpy as np
from . import params as p

def I_leak(V):
    return p.gl * (V - p.El)

def I_Na(aNa, V):
    return p.gNa * aNa * (V - p.ENa)

def aNa_inf(V):
    return 1.0 / (1.0 + np.exp(-p.SNa * (V - p.V0Na)))

def I_K(aK, V):
    return p.gK * aK * (V - p.EK)

def aK_inf(V):
    return 1.0 / (1.0 + np.exp(-p.SK * (V - p.V0K)))

def I_GLU(V, aGLU):
    return p.gGLU * aGLU * (V - p.EGLU)

def I_GLU_Trh(V, aGLU):
    return p.gGLU_Trh * aGLU * (V - p.EGLU)

def aGLU_inf(V):
    return 1.0 / (1.0 + np.exp(-p.Ssyn * (V - p.Vspike)))

def I_GABA(V, aGABA):
    return p.gGABA * aGABA * (V - p.EGABA)

def I_OREXIN(V, aOREXIN):
    return p.gOREXIN * aOREXIN * (V - p.EOREXIN)

def I_MC4R(V, aMC4R):
    return p.gMC4R * aMC4R * (V - p.EMC4R)

def I_INSR(V, aINSR, EINSR, gINSR):
    return gINSR * aINSR * (V - EINSR)

def I_GHSR(V, aGHSR, EGHSR):
    return p.gGHSR * aGHSR * (V - EGHSR)

def I_LEPR(V, aLEPR, ELEPR):
    return p.gLEPR * aLEPR * (V - ELEPR)

def I_circadian(t, wc=0.000073, Ac=0.045):
    phi = np.pi / 3.0
    return Ac * (
        0.97 * np.sin(wc * t - phi)
        + 0.22 * np.sin(2 * wc * t - phi)
        + 0.07 * np.sin(3 * wc * t - phi)
        + 0.03 * np.sin(4 * wc * t - phi)
        + 0.001 * np.sin(5 * wc * t - phi)
        + 1.0
    )

def I_SCN(t):
    return 0.045 * np.cos(2.44009422e-01 * (t / 3600.0) - 9.08555848e-01)

def I_tonic_SCN(t):
    CT = (t / 3600.0) % 24.0
    cfos = 497.67 * np.cos(0.2440 * CT - 0.9086) + 869.59
    norm = (cfos - 371.92) / (1367.26 - 371.92)
    return 0.3 + norm * (1.1 - 0.3)

def Ra_MC4R(aMSH, AgRP):
    numer = aMSH + p.Talpha * (1.0 + (AgRP / p.TA))
    denom = (
        p.Talpha * ((1.0 + p.m) + (1.0 / p.TA + p.m / p.ZA) * AgRP)
        + (1.0 + p.m * p.Talpha / p.Zalpha) * aMSH
            )
    return numer / denom

def C_release(V, maxval):
    return maxval / (1.0 + np.exp(-p.Ssyn * (V - p.Vspike)))

def Ra_INSR(Pin):
    return Pin / (p.EC50INSR + Pin)

def Ra_GHSR(Pghr):
    return Pghr / (p.EC50GHSR + Pghr)

def Ra_LEPR(Plep):
    return Plep / (p.EC50LEPR + Plep)

