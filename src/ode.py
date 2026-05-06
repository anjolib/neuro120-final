from typing import Callable
from . import params as p
from .currents import (
    I_leak, I_Na, aNa_inf, I_K, aK_inf,
    aGLU_inf, I_GLU, I_GABA,
    I_OREXIN, I_MC4R,
    I_INSR, I_GHSR, I_LEPR,
    I_circadian,
    Ra_MC4R, C_release,
    Ra_INSR, Ra_GHSR, Ra_LEPR
)
from .hormones import (
    dGI1_dt, dGI2_dt,
    dGH_dt, dIN_dt, dLP_dt
)

Y0 = [
    -40.0, -60.0, -60.0, -60.0, -60.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0,
    0.0, 1.0,
    0.0, 0.0, 1.38,
    0.0, 0.01, 670.0, 79.0, 22.0,
    0.11, 0.6, 0.6,
]

STATE_NAMES = [
    "V", "A", "P", "L", "S",
    "aKV", "aKA", "aKP", "aKL", "aKS",
    "aGLUP", "aGLUV", "GABAA", "GABAL",
    "aOrexina", "Morexina",
    "alphaMSH", "AgRP", "AMC4R",
    "GI1", "GI2", "H1", "IN", "LP",
    "AINSR", "AGHSR", "ALEPR",
]

def make_ode(food_fn: Callable[[float], float]):
    def rhs(t: float, y: list) -> list:
        (V, A, P, L, S,
         aKV, aKA, aKP, aKL, aKS,
         aGLUP, aGLUV, GABAA, GABAL,
         aOrexina, Morexina,
         aMSH, AgRP, AMC4R,
         GI1, GI2, GH, IN, LP,
         AINSR, AGHSR, ALEPR) = y

        dV = (-I_leak(V)
              - I_Na(aNa_inf(V), V)
              - I_K(aKV, V)
              - I_GLU(V, aGLUP)
              - I_GABA(V, GABAA)
              - I_GABA(V, GABAL)
              - I_MC4R(V, AMC4R))

        dA = (-I_leak(A)
              - I_Na(aNa_inf(A), A)
              - I_K(aKA, A)
              - I_GLU(A, aGLUV)
              - I_OREXIN(A, aOrexina)
              - I_GHSR(A, AGHSR, -10.0)
              - I_LEPR(A, ALEPR, -90.0)
              - I_INSR(A, AINSR, -90.0, p.gINSRA))

        dP = (-I_leak(P)
              - I_Na(aNa_inf(P), P)
              - I_K(aKP, P)
              - I_GABA(P, GABAA)
              + I_OREXIN(P, aOrexina)
              - I_INSR(P, AINSR, -10.0, p.gINSRA)
              - I_LEPR(A, ALEPR, -20.0))

        dL = (-I_leak(L)
              - I_Na(aNa_inf(L), L)
              - I_K(aKL, L)
              - I_MC4R(L, AMC4R)
              + I_circadian(t))

        dS = (-I_leak(S)
              - I_Na(aNa_inf(S), S)
              - I_K(aKS, S)
              - I_circadian(t))

        daKV = (aK_inf(V) - aKV) / p.tauK
        daKA = (aK_inf(A) - aKA) / p.tauK
        daKP = (aK_inf(P) - aKP) / p.tauK
        daKL = (aK_inf(L) - aKL) / p.tauK
        daKS = (aK_inf(S) - aKS) / p.tauK

        daGLUP = (aGLU_inf(P) - aGLUP) / p.tauGLU
        daGLUV = (aGLU_inf(V) - aGLUV) / p.tauGLU
        dGABAA = (aGLU_inf(A) - GABAA) / p.tauGABA
        dGABAL = (aGLU_inf(L) - GABAL) / p.tauGABA

        s_L = aGLU_inf(L)
        daOrexina = (s_L * Morexina - aOrexina) / p.tauOREXIN
        dMorexina  = ((p.Mmax - Morexina) / p.tauincrease
                      - s_L * Morexina / p.taudecrease)

        daMSH  = C_release(P, 1e-9) / p.tauaco_aMSH - aMSH  / p.taueli_aMSH
        dAgRP  = C_release(A, 1e-8) / p.tauaco_AgRP - AgRP / p.taueli_AgRP
        dAMC4R = Ra_MC4R(aMSH, AgRP) / p.tauact_MC4R - AMC4R / p.tauinac_MC4R

        food  = food_fn(t)
        dGI1  = dGI1_dt(GI1, GI2) + food
        dGI2  = dGI2_dt(GI1, GI2)
        dGH   = dGH_dt(GH, GI2, t)
        dIN   = dIN_dt(GI2, IN)
        dLP   = dLP_dt(LP, IN)

        dAINSR = Ra_INSR(IN) / p.tauact_INSR - AINSR / p.tauinac_INSR
        dAGHSR = Ra_GHSR(GH)   / p.tauact_GHSR - AGHSR / p.tauinac_GHSR
        dALEPR = Ra_LEPR(LP) / p.tauact_LEPR - ALEPR / p.tauinac_LEPR

        return [
            dV, dA, dP, dL, dS,
            daKV, daKA, daKP, daKL, daKS,
            daGLUP, daGLUV, dGABAA, dGABAL,
            daOrexina, dMorexina,
            daMSH, dAgRP, dAMC4R,
            dGI1, dGI2, dGH, dIN, dLP,
            dAINSR, dAGHSR, dALEPR,
        ]

    return rhs

