import numpy as np
from dataclasses import dataclass
import inputs

@dataclass
class Params:
    k_A: float = 0.5 
    w_S: float = 1.0
    w_G: float = 1.0
    w_L: float = 1.0
    K_G: float = 1.0
    K_L: float = 4.0
    b:   float = 0.5

params = Params()

def dA_dt(t: float, A: float, p: Params = params, scn_on: bool = True) -> float:
    """
    AgRP rate-of-change at time t.

    Arguments:
        t (float): current time in hours
        A (float): current AgRP activity
        p (Params): model parameters
        scn_on (bool): if False, SCN ablated

    Returns:
        float: rate of change of AgRP activity
    """
    g = inputs.ghrelin(t, "h", False)
    l = inputs.leptin(t, False)
    scn_drive     = p.w_S * inputs.scn(t) if scn_on else 0.0
    ghrelin_drive = p.w_G * g / (p.K_G + g)
    leptin_drive  = p.w_L * l / (p.K_L + l)

    return -p.k_A * A + scn_drive + ghrelin_drive - leptin_drive + p.b

def euler_simulate(
    t0: float = 0.0,
    tf: float = 72.0,
    dt: float = 0.05,
    A0: float = 0.0,
    p: Params = params,
    scn_on: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Euler method for dA/dt.

    Arguments:
        t0 (float): initial time for simulation
        tf (float): final time for simulation
        dt (float): width of time step
        A0 (float): initial AgRP activity
        p (Params): model parameters
        scn_on (bool): if False, SCN ablated

    Returns:
        t (np.ndarray): time array
        A (np.ndarray): AgRP activity array
    """
    n = int(np.floor((tf - t0) / dt)) + 1
    t = np.linspace(t0, t0 + dt * (n - 1), n)
    A = np.empty(n)
    A[0] = A0
    for i in range(n - 1):
        A[i + 1] = A[i] + dt * dA_dt(t[i], A[i], p, scn_on)
    return t, A


