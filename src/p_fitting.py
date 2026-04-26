import math
import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
import sys
sys.path.append("../src/")
from mylib import *
import importlib

import pandas as pd
df = pd.read_csv("../data/AgRP_control-Douglass_2G.csv")

t_data = df["time"].values
A_data = df[["m1", "m2", "m3", "m4", "m5"]].mean(axis=1).values
A_std = df[["m1", "m2", "m3", "m4", "m5"]].std(axis=1).values

import matplotlib.pyplot as plt

"""plt.plot(t_data, A_data, label="Mean")
plt.fill_between(t_data, A_data - A_std, A_data + A_std, alpha=0.3)
plt.legend()
plt.show()"""

A0 = A_data[0]

from scipy.integrate import solve_ivp

def simulate(params, t_span, t_eval, A0, scn_on=True):
    def ode(t, A):
        return dA_dt(t, A[0], params, scn_on)

    sol = solve_ivp(ode, t_span, [A0], t_eval=t_eval)
    return sol.y[0]

def loss(theta, t_data, A_data, A0):
    p = Params(*theta)
    
    A_model = simulate(p, (t_data[0], t_data[-1]), t_data, A0)
    
    return np.mean((A_model - A_data)**2)

from scipy.optimize import differential_evolution

bounds = [
    (0.01, 5.0),   # k_A
    (0.0, 5.0),    # w_S
    (0.0, 5.0),    # w_G
    (0.0, 5.0),    # w_L
    (0.1, 10.0),   # K_G
    (0.1, 10.0),   # K_L
    (0.0, 2.0)     # b
]

result = differential_evolution(
    lambda theta: loss(theta, t_data, A_data, A0),
    bounds
)

best_theta = result.x

np.save("../data/best_theta.npy", best_theta)

A_fit = simulate(Params(*best_theta), (t_data[0], t_data[-1]), t_data, A0)

import matplotlib.pyplot as plt
plt.plot(t_data, A_data, 'o', label="data")
plt.plot(t_data, A_fit, 'b:', label="model")
plt.legend()
plt.show()