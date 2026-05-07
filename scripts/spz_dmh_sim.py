"""
Isolated simulation of the SPZ (S) and DMH (D) neurons.

State vector:  [S, D, aKS, aKD, GABAS]
Inputs:        I_circadian(t) drives S; GABAS (from S) inhibits D.
"""

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from scipy.integrate import solve_ivp

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from src import params as p
from src.currents import (
    I_leak, I_Na, aNa_inf,
    I_K, aK_inf,
    I_GABA, aGLU_inf,
    I_circadian,
)

# ---------------------------------------------------------------------------
# State layout
# ---------------------------------------------------------------------------
STATE_NAMES = ["S", "D", "aKS", "aKD", "GABAS"]

Y0 = [
    -60.0,  # S  (SPZ voltage, mV)
    -60.0,  # D  (DMH voltage, mV)
    0.0,    # aKS
    0.0,    # aKD
    0.0,    # GABAS
]


# ---------------------------------------------------------------------------
# ODE
# ---------------------------------------------------------------------------
def rhs(t: float, y: list) -> list:
    S, D, aKS, aKD, GABAS = y

    dS = (
        -I_leak(S)
        - I_Na(aNa_inf(S), S)
        - I_K(aKS, S)
        - I_circadian(t)
        + 1.0           # tonic excitatory bias (keeps SPZ near threshold)
    )

    dD = (
        -I_leak(D)
        - I_Na(aNa_inf(D), D)
        - I_K(aKD, D)
        - I_GABA(D, GABAS)
        + I_circadian(t)  # circadian excitation released from SPZ disinhibition
    )

    daKS = (aK_inf(S) - aKS) / p.tauK
    daKD = (aK_inf(D) - aKD) / p.tauK
    dGABAS = (aGLU_inf(S) - GABAS) / p.tauGABA

    return [dS, dD, daKS, daKD, dGABAS]


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
@dataclass
class SpzDmhResult:
    t_s: np.ndarray
    t_h: np.ndarray
    y:   np.ndarray

    S:     np.ndarray = field(init=False)
    D:     np.ndarray = field(init=False)
    aKS:   np.ndarray = field(init=False)
    aKD:   np.ndarray = field(init=False)
    GABAS: np.ndarray = field(init=False)

    def __post_init__(self):
        for i, name in enumerate(STATE_NAMES):
            object.__setattr__(self, name, self.y[i])


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
def run(
    duration_h: float = 72.0,
    dt: float = 0.1,
    y0: list | None = None,
    rtol: float = 1e-6,
    atol: float = 1e-9,
    max_step: float = 1.0,
) -> SpzDmhResult:
    if y0 is None:
        y0 = Y0

    t_end  = p.gen(duration_h)
    t_eval = np.arange(0.0, t_end, dt)

    sol = solve_ivp(
        rhs,
        [0.0, t_end],
        y0,
        method="LSODA",
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        max_step=max_step,
    )

    if not sol.success:
        raise RuntimeError(f"Solver failed: {sol.message}")

    return SpzDmhResult(t_s=sol.t, t_h=sol.t / 3600.0, y=sol.y)


# ---------------------------------------------------------------------------
# Quick plot (run as __main__)
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    res = run(duration_h=72.0)

    t = res.t_h
    circ = np.array([I_circadian(ti * 3600.0) for ti in t])

    fig, axes = plt.subplots(4, 1, figsize=(11, 10), sharex=True)

    axes[0].plot(t, circ, color="goldenrod")
    axes[0].set_ylabel("I_circ (nA)")
    axes[0].set_title("SCN circadian drive")
    axes[0].grid(alpha=0.3)

    axes[1].plot(t, res.S, color="steelblue")
    axes[1].set_ylabel("V (mV)")
    axes[1].set_title("SPZ voltage (S)")
    axes[1].grid(alpha=0.3)

    axes[2].plot(t, res.GABAS, color="mediumorchid")
    axes[2].set_ylabel("GABAS")
    axes[2].set_title("SPZ GABA release → DMH")
    axes[2].grid(alpha=0.3)

    axes[3].plot(t, res.D, color="tomato")
    axes[3].set_ylabel("V (mV)")
    axes[3].set_title("DMH voltage (D)")
    axes[3].grid(alpha=0.3)

    axes[3].set_xlabel("Time (hr)")
    plt.tight_layout()
    plt.show()
