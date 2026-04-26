import math
import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
import sys
sys.path.append("../src/")
from mylib import *
import importlib
import inputs

import pandas as pd
df_on = pd.read_csv("../data/AgRP_control-Douglass_2G.csv")
df_off = pd.read_csv("../data/AgRP_Caspase-Douglass_5H.csv")

t_data = df_on["time"].values
A_data_on = df_on[["m1", "m2", "m3", "m4", "m5"]].mean(axis=1).values
A_std = df_on[["m1", "m2", "m3", "m4", "m5"]].std(axis=1).values

t_data_off = df_off["time"].values
A_data_off = df_off[["m1", "m2", "m3", "m4", "m5"]].mean(axis=1).values
A_std_off = df_off[["m1", "m2", "m3", "m4", "m5"]].std(axis=1).values

import matplotlib.pyplot as plt

"""plt.plot(t_data, A_data, label="Mean")
plt.fill_between(t_data, A_data - A_std, A_data + A_std, alpha=0.3)
plt.legend()
plt.show()"""

A0_on = A_data_on[0]
A0_off = A_data_off[0]

from scipy.integrate import solve_ivp

def simulate(params, t_span, t_eval, A0, scn_on=True):
    def ode(t, A):
        return dA_dt(t, A[0], params, scn_on)

    sol = solve_ivp(ode, t_span, [A0], t_eval=t_eval)
    return sol.y[0]

def loss(theta, t_data, A_data_on, A_data_off, A0_on, A0_off):
    p = Params(*theta)

    A_on  = simulate(p, (t_data[0], t_data[-1]), t_data, A0_on, scn_on=True)
    A_off = simulate(p, (t_data[0], t_data[-1]), t_data, A0_off, scn_on=False)
    
    if np.any(np.isnan(A_on)) or np.any(np.isnan(A_off)):
        return 1e6
    
    loss_on  = np.mean((A_on - A_data_on)**2)
    loss_off = np.mean((A_off - A_data_off)**2)
    
    return loss_on + loss_off

from scipy.optimize import differential_evolution

bounds = [
    (0.01, 5.0),   # k_A
    (0.0, 5.0),    # w_S
    (0.0, 10.0),    # w_G
    (0.0, 10.0),    # w_L
    (0.1, 1.0),   # K_G
    (0.1, 1.0),   # K_L
    (0.0, 2.0)     # b
]

result = differential_evolution(
    lambda theta: loss(theta, t_data, A_data_on, A_data_off, A0_on, A0_off),
    bounds
)

best_theta = result.x

np.save("../data/best_theta.npy", best_theta)
print(best_theta)

A_fit_on = simulate(Params(*best_theta), (t_data[0], t_data[-1]), t_data, A0_on)
A_fit_off = simulate(Params(*best_theta), (t_data_off[0], t_data_off[-1]), t_data_off, A0_off, scn_on=False)

G_regular = []
L_regular = []
for i in t_data:
    G_regular.append(inputs.ghrelin_hypothalamus(t=i))
    L_regular.append(inputs.leptin(t=i))

import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(12, 5))



axes[0].plot(t_data, A_data_on, 'o', label="data mean ± STD")
axes[0].plot(t_data, A_fit_on, 'b:', label="model")
axes[0].fill_between(t_data, A_data_on - A_std, A_data_on + A_std, alpha=0.3)
axes[0].set_title(f"SCN ON vs Control")
axes[0].set_xlabel("Time (hr)")
axes[0].set_ylabel("AgRP activity")
#axes[0].plot(t_data, L_regular, 'r:', label="leptin")
#axes[0].plot(t_data, G_regular, 'y:', label="ghrelin")
axes[0].legend()

axes[1].plot(t_data, A_data_off, 'ro', label="data mean ± STD")
axes[1].plot(t_data, A_fit_off, 'r:', label="model")
axes[1].fill_between(t_data, A_data_off - A_std_off, A_data_off + A_std_off, color='red', alpha=0.3)
axes[1].set_title(f"SCN OFF vs Caspase")
axes[1].set_xlabel("Time (hr)")
axes[1].set_ylabel("AgRP activity")
#axes[1].plot(t_data, L_regular, 'r:', label="leptin")
#axes[1].plot(t_data, G_regular, 'y:', label="ghrelin")
axes[1].legend()

plt.suptitle("Model vs Experimental AgRP Activity (Douglass et al.)", fontsize=13)
plt.show()

