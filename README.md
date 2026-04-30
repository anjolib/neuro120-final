# A Dynamical Systems Model of a Circadian Feeding Circuit: SCN-Driven Rhythmicity in the Arcuate Nucleus of the Hypothalamus

Final project for Neuro 120 (Spring 2026), by Antonino Libarnes, Soren James, Theo Tobel.

## Model

Population AgRP activity `A(t)` is governed by:

```
dA/dt = -k_A · A(t) + w_S · S(t) + u_G(t) - u_L(t) + b
```

where `S(t)` is normalized SCN drive and `u_G`, `u_L` are Michaelis-Menten functions for ghrelin and leptin. Parameters are fit jointly against control and SCN-ablated AgRP data (Douglass 2025), with hormone time series from Bodosi 2004 and SCN reference from Yamashita 2026.

## Layout

```
data/         CSV inputs + saved best-fit parameters (best_theta.npy)
src/
  mylib.py            Params, dA/dt, Euler integrator
  inputs.py           S(t), G(t), L(t)
  p_fitting.py        differential evolution fit, joint SCN-ON/OFF loss
  experimental_data.py  CSV loaders + cosinor helpers
scripts/      exploratory notebooks
flake.nix     Nix dev shell
```

## Running

```
nix develop                  # or: pip install numpy scipy pandas matplotlib jupyter
cd src
python p_fitting.py          # fits parameters, saves best_theta.npy
```

## Data

- **Douglass 2025** — AgRP fiber photometry, control + caspase-ablated SCN
- **Bodosi 2004** — circadian ghrelin and leptin (digitized via WebPlotDigitizer)
- **Yamashita 2026** — whole-brain c-Fos atlas, SCN reference (digitized via WebPlotDigitizer)
- **Sayar-Atasoy 2024** — held-out validation (data not in repo)

