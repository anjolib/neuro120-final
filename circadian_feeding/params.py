def gen(x):
    return float(x) * 3600.0

def break_(x):
    return float(x) * 3600.0 + 2340.0

def lun(x):
    return float(x) * 3600.0 + 1676.0

def din(x):
    return float(x) * 3600.0 + 1400.0

def num(x):
    return float(x) / 3600.0

def day(x):
    return float(x) * 3600.0

gl = 0.1
El = -60.0

gNa = 3.0
ENa = 50.0
SNa = 0.25
V0Na = -25.0

gK = 4.0
EK = -90.0
SK = 0.25
V0K = -25.0
tauK = 2.0

gGLU = 0.1
EGLU = 50.0
tauGLU = 30.0

gGABA = 0.1
EGABA = -65.0
tauGABA = 30.0

Ssyn = 1.0
Vspike = -20.0

gOREXIN = 0.2
EOREXIN = 50.0
tauOREXIN = 300.0
tauincrease = 27000.0
taudecrease = 3312.0
Mmax = 1.0

gMC4R = 0.1
EMC4R = -38.0
tauaco_aMSH = 400.0
taueli_aMSH = 300.0
tauaco_AgRP = 300.0
taueli_AgRP = 300.0
tauact_MC4R = 734.4
tauinac_MC4R = 2030.4

k12 = num(2.30)
k2x = num(0.71)
sigma = 0.20
khx = num(2.354)
beta = num(1591.0)
beta2 = num(1238.0)
gamma = 1.29

theta = num(528.63)
delta = 0.38
s = 4.0
kix = num(1.03)
rho = num(90.0)

epsilon = num(6231.39)
FM = 60000.0
BV = 5000.0

R = num(0.11)
eta = 24.0
l = 1.0
EC50LEP = 661.0

gINSRA = 0.1
gINSRP = 0.1
tauact_INSR = 300.0
tauinac_INSR = 300.0
EC50INSR = 661.0

gGHSR = 0.1
tauact_GHSR = 300.0
tauinac_GHSR = 300.0
EC50GHSR = 500.0

gLEPR = 0.045
tauact_LEPR = 300.0
tauinac_LEPR = 300.0
EC50LEPR = 16.0

basal = 0.5
m = 1.0

def _KK(kapp, maxval):
    return kapp * basal / maxval

def _kk(kapp, maxval):
    Lx = (1.0 - basal) / basal
    return _KK(kapp, maxval) * Lx * maxval / (1.0 - maxval)

TA = _KK(3.6e-9, 0.1)
ZA = _kk(3.6e-9, 0.1)
Talpha = _KK(0.16e-9, 0.9)
Zalpha = _kk(0.16e-9, 0.9)

I_tonic_S = 0.85
I_tonic_D = 0.82

gGLU_Trh = 0.15
