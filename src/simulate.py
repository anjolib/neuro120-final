import numpy as np
from dataclasses import dataclass, field
from typing import Callable

from scipy.integrate import solve_ivp

from . import params as p
from .ode import make_ode, Y0, STATE_NAMES

@dataclass
class SimResult:
    t_s:  np.ndarray   # time in seconds
    t_h:  np.ndarray   # time in hours
    y:    np.ndarray   # raw (n_states × n_points) solution array

    # Voltages
    V: np.ndarray = field(init=False)
    A: np.ndarray = field(init=False)
    P: np.ndarray = field(init=False)
    L: np.ndarray = field(init=False)
    S: np.ndarray = field(init=False)
    D: np.ndarray = field(init=False)

    # K-channel gating
    aKV: np.ndarray = field(init=False)
    aKA: np.ndarray = field(init=False)
    aKP: np.ndarray = field(init=False)
    aKL: np.ndarray = field(init=False)
    aKS: np.ndarray = field(init=False)
    aKD: np.ndarray = field(init=False)

    # Ionotropic synaptic gating
    aGLUP: np.ndarray = field(init=False)
    aGLUV: np.ndarray = field(init=False)
    aGLUD: np.ndarray = field(init=False)
    GABAA: np.ndarray = field(init=False)
    GABAL: np.ndarray = field(init=False)
    GABAS: np.ndarray = field(init=False)

    # Orexin
    aOrexina: np.ndarray = field(init=False)
    Morexina: np.ndarray = field(init=False)

    # Neuropeptides / MC4R
    alphaMSH: np.ndarray = field(init=False)
    AgRP:     np.ndarray = field(init=False)
    AMC4R:    np.ndarray = field(init=False)

    # Hormone compartments
    GI1: np.ndarray = field(init=False)
    GI2: np.ndarray = field(init=False)
    H1:  np.ndarray = field(init=False)
    IN:  np.ndarray = field(init=False)
    LP:  np.ndarray = field(init=False)

    # Receptor gating
    AINSR: np.ndarray = field(init=False)
    AGHSR: np.ndarray = field(init=False)
    ALEPR: np.ndarray = field(init=False)

    def __post_init__(self):
        for i, name in enumerate(STATE_NAMES):
            object.__setattr__(self, name, self.y[i])

    def summary(self) -> str:
        lines = ["=" * 60, "Simulation summary", "=" * 60]
        pairs = [
            ("Plasma glucose GI2 (mmol/L)", self.GI2),
            ("Insulin IN (pM)",             self.IN),
            ("Ghrelin H1 (pg/mL)",          self.H1),
            ("Leptin LP (ng/mL)",           self.LP),
            ("α-MSH (nM)",                  self.alphaMSH * 1e9),
            ("AgRP (×1e-8 M)",              self.AgRP     * 1e8),
            ("MC4R activation",             self.AMC4R),
        ]
        fmt = "  {:<34s}  min={:8.2f}  mean={:8.2f}  max={:8.2f}"
        for name, arr in pairs:
            lines.append(fmt.format(name, arr.min(), arr.mean(), arr.max()))
        return "\n".join(lines)


def run(
    food_fn: Callable[[float], float],
    duration_h: float = 72.0,
    dt: float = 0.1,
    y0: list | None = None,
    rtol: float = 1e-6,
    atol: float = 1e-9,
    max_step: float = 1.0,
) -> SimResult:
    if y0 is None:
        y0 = Y0

    t_end  = p.gen(duration_h)
    t_eval = np.arange(0.0, t_end, dt)

    rhs = make_ode(food_fn)

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

    return SimResult(t_s=sol.t, t_h=sol.t / 3600.0, y=sol.y)

