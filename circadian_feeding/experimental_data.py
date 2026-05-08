import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

def get_experimental_results(file: str):
    df = pd.read_csv(file)
    t = df["time"].to_numpy(float)
    data = df.drop("time", axis=1).to_numpy(float)
    return t, data
    
def cos_model(t, A, phi, C):
    return A * np.cos(2 * np.pi * t / 24 + phi) + C

def fit_cosine(t, y):
    p0 = [(y.max() - y.min()) / 2, 0.0, y.mean()]
    popt, _ = curve_fit(cos_model, t, y, p0=p0, maxfev=10000)
    return popt

def give_cosine(t, data, t0, tf):
    popt = fit_cosine(t, data)
    t_smooth = np.linspace(t0, tf, 1000)
    return cos_model(t_smooth, *popt)


#def rmse(pred, obs):
#    return np.sqrt(np.mean((pred - obs) ** 2))
#
#rmse_ctrl = rmse(A_ss_on_scaled,  cos_model(t_ss_on,  *ctrl_popt))
#rmse_casp = rmse(A_ss_off_scaled, cos_model(t_ss_off, *casp_popt))
#print(f"RMSE  model(SCN ON)  vs control cosine fit : {rmse_ctrl:.4f}")
#print(f"RMSE  model(SCN OFF) vs caspase cosine fit : {rmse_casp:.4f}")