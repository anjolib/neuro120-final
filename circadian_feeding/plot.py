import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from typing import Callable

from .simulate import SimResult
from .params import num, gen, break_, lun, din
from .currents import I_circadian
from .food import TRE

def _neuron_specs(neuron, res):
    match neuron:
        case "S":
            return (res.S, "SPZ", "blue")
        case "D":
            return (res.D, "DMH", "red")
        case "V":
            return (res.V, "PVH", "orange")
        case "A":
            return (res.A, "AgRP", "tomato")
        case "P":
            return (res.P, "Pomc", "seagreen")
        case "L":
            return (res.L, "LHA", "darkorange")
    return None

def _hormone_specs(hormone, res):
    match hormone:
        case "GI1":
            return (res.GI1, "GI1 (L)", "darkgreen")
        case "GI2":
            return (res.GI2, "GI2 (L)", "darkgreen")
        case "IN":
            return (res.IN, "Insulin (pM)", "steelblue")
        case "H1":
            return (res.H1, "Ghrelin (pg/mL)", "tomato")
        case "LP":
            return (res.LP, "Leptin plasma (ng/mL)", "purple")
    return None

def _peptide_specs(peptide, res):
    match peptide:
        case "AgRP":
            return (res.AgRP*1e8, "AgRP (×10⁻⁸ M)", "blue")
        case "aMSH":
            return (res.alphaMSH*1e9, "α-MSH (nM)", "green")
        case "AMC4R":
            return (res.AMC4R, "MC4R activation", "red")
    return None

_MEAL_SPANS_72H = [
    (gen(8) /3600, break_(8) /3600),  (gen(13)/3600, lun(13)/3600),
    (gen(17)/3600, din(17) /3600),
    (gen(32)/3600, break_(32)/3600),  (gen(37)/3600, lun(37)/3600),
    (gen(41)/3600, din(41) /3600),
    (gen(56)/3600, break_(56)/3600),  (gen(61)/3600, lun(61)/3600),
    (gen(65)/3600, din(65) /3600),
]

def _shade_meals(ax, meal_spans=_MEAL_SPANS_72H, alpha=0.13, color="gold"):
    for ms, me in meal_spans:
        ax.axvspan(ms, me, alpha=alpha, color=color)

def _day_lines(ax, days=(24, 48)):
    for d in days:
        ax.axvline(d, color="gray", lw=0.7, ls="--", alpha=0.5)

def plot_voltages(res: SimResult,
                  axes,
                  neurons: list = ["S", "D"],
                  meal_spans: list = _MEAL_SPANS_72H,
                  t_range: (float) = (0, 24.0)):
    """Plots voltages of neurons over time.

    Args:
        res: SimResult object containing results
        neurons: list of neurons to plot by variable name
        t_range: time range to plot (initial, final)

    Returns:
        pyplot figure containing voltage plots for each neuron
    """

    mask = (res.t_h >= t_range[0]) & (res.t_h <= t_range[1])
    t = res.t_h[mask]

    if len(neurons) == 1:
        data, label, color = _neuron_specs(neurons[0], res);
        data = data[mask]
        axes.plot(t, data, color=color, lw=0.5)
        axes.set_ylabel(f"{label}\n(mV)", fontsize=9)
        axes.set_ylim(-100, 70)
        axes.set_xlim(t_range)
        axes.axhline(0, color="gray", lw=0.4, ls="--")
        _shade_meals(axes, meal_spans)
        _day_lines(axes)
        axes.set_xlabel("Time (hours)")
        plt.tight_layout()
        return axes

    for ax, neuron in zip(axes, neurons):
        data, label, color = _neuron_specs(neuron, res)
        data = data[mask]
        ax.plot(t, data, color=color, lw=0.5)
        ax.set_ylabel(f"{label}\n(mV)", fontsize=9)
        ax.set_ylim(-100, 70)
        ax.set_xlim(t_range)
        ax.axhline(0, color="gray", lw=0.4, ls="--")
        _shade_meals(ax, meal_spans)
        _day_lines(ax)
        ax.set_xlabel("Time (hours)")

    return axes

def plot_hormones(res: SimResult,
                  axes,
                  hormones: list = ["GI2", "IN", "H1", "LP"],
                  meal_spans: list = _MEAL_SPANS_72H,
                  t_range: (float) = (0, 24.0)):
    """Glucose, insulin, ghrelin, and leptin over 72 h."""

    if len(hormones) == 1:
        data, label, color = _hormone_specs(hormones[0], res)
        axes.plot(res.t_h, data, color=color, lw=1.0)
        axes.set_ylabel(label, fontsize=9)
        _shade_meals(axes, meal_spans)
        _day_lines(axes)
        axes.set_xlabel("Time (hours)")
        axes.set_xlim(0, 72)
        axes.set_ylim(bottom=0)
        return axes

    for ax, hormone in zip(axes, hormones):
        data, label, color = _hormone_specs(hormone, res)
        ax.plot(res.t_h, data, color=color, lw=1.0)
        ax.set_ylabel(label, fontsize=9)
        _shade_meals(ax, meal_spans)
        _day_lines(ax)
        ax.set_xlabel("Time (hours)")
        ax.set_xlim(0, 72)
        ax.set_ylim(bottom=0)

    return axes

def plot_neuropeptides(res: SimResult,
                       axes,
                       peptides: list = ["aMSH", "AgRP", "AMC4R"],
                       meal_spans: list = _MEAL_SPANS_72H,
                       t_range: (float) = (0, 24.0)):
    """α-MSH, AgRP cleft concentrations and MC4R activation over 72 h."""

    if len(peptides) == 1:
        data, label, color = _peptide_specs(peptides[0], res)
        axes.plot(res.t_h, data, color=color, lw=1)
        axes.set_ylabel(label, fontsize=9)
        _shade_meals(axes, meal_spans)
        _day_lines(axes)
        axes.set_xlabel("Time (hours)")
        axes.set_xlim(t_range)
        axes.set_ylim(bottom=0)
        return axes

    for ax, peptide in zip(axes, peptides):
        data, label, color = _peptide_specs(peptide, res)
        ax.plot(res.t_h, data, color=color, lw=1)
        ax.set_ylabel(label, fontsize=9)
        _shade_meals(ax, meal_spans)
        _day_lines(ax)
        ax.set_xlabel("Time (hours)")
        ax.set_xlim(t_range)
        ax.set_ylim(bottom=0)

    return axes

def plot_spike_counts(
    res: SimResult,
    axes,
    neurons: list = ["S", "D"],
    meal_spans: list = _MEAL_SPANS_72H,
    t_range: (float) = (0, 24.0),
    window_h: float = 0.5,
    peak_height: float = 0.0,
    peak_distance: int = 10,
    zscore: bool = False
):

    edges   = np.arange(t_range[0], t_range[1], window_h)
    centers = (edges[:-1] + edges[1:]) / 2

    def _count(voltage):
        counts = []
        for i in range(len(edges) - 1):
            t0 = edges[i]     * 3600.0
            t1 = edges[i + 1] * 3600.0
            seg = voltage[(res.t_s >= t0) & (res.t_s < t1)]
            if len(seg):
                pks, _ = find_peaks(seg, height=peak_height, distance=peak_distance)
                counts.append(len(pks))
            else:
                counts.append(0)
        return counts

    def _zscore(x):
        x = np.array(x)
        return (x - np.mean(x)) / (np.std(x) + 1e-8)

    if len(neurons) == 1:
        data, label, color = _neuron_specs(neurons[0], res)
        counts = _count(data)
        axes.set_xlabel("Time (hours)")
        axes.set_xlim(t_range)
        _shade_meals(axes, meal_spans)
        _day_lines(axes)
        if zscore:
            counts = _zscore(counts)
            axes.plot(centers, counts, color=color)
            axes.set_ylabel(f"{label} (z-score)", fontsize=9)
        else:
            axes.bar(centers, counts, width=window_h * 0.9, color=color, alpha=0.8)
            axes.set_ylabel(f"{label} (spikes)", fontsize=9)
        return axes

    for ax, neurons in zip(axes, neurons):
        data, label, color = _neuron_specs(neurons, res)
        counts = _count(data)
        ax.set_xlabel("Time (hours)")
        ax.set_xlim(t_range)
        _shade_meals(ax, meal_spans)
        _day_lines(ax)
        if zscore:
            counts = _zscore(counts)
            ax.plot(centers, counts, color=color)
            ax.set_ylabel(f"{label} (z-score)", fontsize=9)
        else:
            ax.bar(centers, counts, width=window_h * 0.9, color=color, alpha=0.8)
            ax.set_ylabel(f"{label} (spikes)", fontsize=9)

    return axes

def plot_inputs(axes,
                food_fn: Callable[[float], float],
                food_rate: float = 2.1,
                t_range: (float) = (0, 24.0)):
    """Circadian signal and TRE food pattern over `hours` hours."""
    t_s = np.linspace(gen(t_range[0]), gen(t_range[1]), 10000)
    t_h = t_s / 3600.0

    from .params import num
    circ = I_circadian(t_s, wc=0.000073, Ac=1.0)
    food = np.array([food_fn(ti, food_rate) for ti in t_s])

    axes[0].fill_between(t_h, circ, color="goldenrod", lw=1.5, alpha=0.5)
    axes[0].axhline(1, color="gray", lw=0.5, ls="--")
    axes[0].set_ylabel("Circadian signal", fontsize=9)
    axes[0].set_ylim(bottom=0)

    axes[1].fill_between(t_h, food, color="sienna", alpha=0.7)
    axes[1].set_ylabel("Food intake (mmol/s)", fontsize=9)
    axes[1].set_xlabel("Time (hours)")
    axes[1].set_ylim(bottom=0)

    return axes

def plot_all(res: SimResult, save_dir: str | None = None) -> list[plt.Figure]:
    import os

    figs = {
        "neuron_voltages_24h":  plot_voltages(res),
        "hormones_72h":         plot_hormones(res),
        "neuropeptides_72h":    plot_neuropeptides(res),
        "synaptic_orexin_72h":  plot_synaptic(res),
        "spike_counts_24h":     plot_spike_counts(res),
        "inputs_24h":           plot_inputs(),
    }

    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        for name, fig in figs.items():
            path = os.path.join(save_dir, f"{name}.png")
            fig.savefig(path, dpi=150, bbox_inches="tight")
            print(f"Saved {path}")

    return list(figs.values())

