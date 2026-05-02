import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

from .simulate import SimResult
from .params import gen, break_, lun, din
from .currents import I_circadian
from .food import TRE


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

def plot_voltages(res: SimResult, hours: float = 24.0) -> plt.Figure:
    mask = res.t_h <= hours
    t = res.t_h[mask]

    meal_spans_24h = [
        (gen(8)/3600,  break_(8)/3600),
        (gen(13)/3600, lun(13)/3600),
        (gen(17)/3600, din(17)/3600),
    ]

    specs = [
        (res.V[mask], "PVH  (V)",  "steelblue"),
        (res.A[mask], "AgRP (A)",  "tomato"),
        (res.P[mask], "POMC (P)",  "seagreen"),
        (res.L[mask], "LHA  (L)",  "darkorange"),
    ]

    fig, axes = plt.subplots(4, 1, figsize=(14, 8), sharex=True)
    fig.suptitle(f"Membrane potentials – first {hours:.0f} h (TRE diet)", fontsize=12)

    for ax, (data, label, color) in zip(axes, specs):
        ax.plot(t, data, color=color, lw=0.5)
        ax.set_ylabel(f"{label}\n(mV)", fontsize=9)
        ax.set_ylim(-100, 70)
        ax.axhline(0, color="gray", lw=0.4, ls="--")
        _shade_meals(ax, meal_spans_24h)

    axes[-1].set_xlabel("Time (hours)")
    plt.tight_layout()
    return fig

def plot_hormones(res: SimResult) -> plt.Figure:
    """Glucose, insulin, ghrelin, and leptin over 72 h."""
    specs = [
        (res.GL2, "Plasma glucose (mmol/L)", "darkgreen"),
        (res.IN,  "Insulin (pM)",            "steelblue"),
        (res.H1,  "Ghrelin (pg/mL)",         "tomato"),
        (res.LP,  "Leptin plasma (ng/mL)",   "purple"),
    ]

    fig, axes = plt.subplots(4, 1, figsize=(14, 9), sharex=True)
    fig.suptitle("Hormone dynamics – 72 h  (gold = meal window)", fontsize=12)

    for ax, (data, label, color) in zip(axes, specs):
        ax.plot(res.t_h, data, color=color, lw=1.0)
        ax.set_ylabel(label, fontsize=9)
        _shade_meals(ax)
        _day_lines(ax)

    axes[-1].set_xlabel("Time (hours)")
    axes[-1].set_xlim(0, 72)
    plt.tight_layout()
    return fig

def plot_neuropeptides(res: SimResult) -> plt.Figure:
    """α-MSH, AgRP cleft concentrations and MC4R activation over 72 h."""
    fig, axes = plt.subplots(3, 1, figsize=(14, 7), sharex=True)
    fig.suptitle("Neuropeptide and MC4R dynamics – 72 h", fontsize=12)

    axes[0].plot(res.t_h, res.alphaMSH * 1e9, color="darkorchid", lw=1)
    axes[0].set_ylabel("α-MSH (nM)", fontsize=9)

    axes[1].plot(res.t_h, res.AgRP * 1e8, color="tomato", lw=1)
    axes[1].set_ylabel("AgRP (×10⁻⁸ M)", fontsize=9)

    axes[2].plot(res.t_h, res.AMC4R, color="navy", lw=1)
    axes[2].set_ylabel("MC4R activation", fontsize=9)

    for ax in axes:
        _shade_meals(ax)
        _day_lines(ax)

    axes[-1].set_xlabel("Time (hours)")
    axes[-1].set_xlim(0, 72)
    plt.tight_layout()
    return fig


def plot_synaptic(res: SimResult) -> plt.Figure:
    """Orexin, GABA, GLU gating and hormone receptor activations over 72 h."""
    fig, axes = plt.subplots(4, 1, figsize=(14, 9), sharex=True)
    fig.suptitle("Synaptic gating, orexin and receptor activations – 72 h", fontsize=12)

    axes[0].plot(res.t_h, res.aOrexina, label="aOrexin (gating)", color="darkorange")
    axes[0].plot(res.t_h, res.Morexina, label="Morexin (mediator)", color="orange", ls="--")
    axes[0].set_ylabel("Orexin", fontsize=9)
    axes[0].legend(fontsize=8, loc="upper right")

    axes[1].plot(res.t_h, res.GABAA, label="GABA from AgRP", color="tomato")
    axes[1].plot(res.t_h, res.GABAL, label="GABA from LHA",  color="darkred", ls="--")
    axes[1].set_ylabel("GABA gating", fontsize=9)
    axes[1].legend(fontsize=8, loc="upper right")

    axes[2].plot(res.t_h, res.aGLUP, label="GLU POMC→PVH", color="steelblue")
    axes[2].plot(res.t_h, res.aGLUV, label="GLU PVH→AgRP", color="blue", ls="--")
    axes[2].set_ylabel("GLU gating", fontsize=9)
    axes[2].legend(fontsize=8, loc="upper right")

    axes[3].plot(res.t_h, res.AINSR, label="AINSR (insulin)",  color="purple")
    axes[3].plot(res.t_h, res.AGHSR, label="AGHSR (ghrelin)", color="seagreen")
    axes[3].plot(res.t_h, res.ALEPR, label="ALEPR (leptin)",   color="teal")
    axes[3].set_ylabel("Receptor\nactivation", fontsize=9)
    axes[3].legend(fontsize=8, loc="upper right")

    for ax in axes:
        _day_lines(ax)

    axes[-1].set_xlabel("Time (hours)")
    axes[-1].set_xlim(0, 72)
    plt.tight_layout()
    return fig

def plot_spike_counts(
    res: SimResult,
    window_h: float = 0.5,
    hours: float = 24.0,
    peak_height: float = 0.0,
    peak_distance: int = 10,
) -> plt.Figure:
    edges   = np.arange(0.0, hours + window_h, window_h)
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

    specs = [
        (_count(res.V), "PVH",  "steelblue"),
        (_count(res.A), "AgRP", "tomato"),
        (_count(res.P), "POMC", "seagreen"),
        (_count(res.L), "LHA",  "darkorange"),
    ]

    fig, axes = plt.subplots(4, 1, figsize=(12, 7), sharex=True)
    fig.suptitle(f"Firing rate (spikes / {window_h*60:.0f}-min window) – first {hours:.0f} h",
                 fontsize=12)

    for ax, (counts, label, color) in zip(axes, specs):
        ax.bar(centers, counts, width=window_h * 0.9, color=color, alpha=0.8)
        ax.set_ylabel(f"{label}\n(spikes)", fontsize=9)

    axes[-1].set_xlabel("Time (hours)")
    axes[-1].set_xlim(0, hours)
    plt.tight_layout()
    return fig


def plot_inputs(hours: float = 24.0, food_rate: float = 2.1) -> plt.Figure:
    """Circadian signal and TRE food pattern over `hours` hours."""
    t_s = np.linspace(0, gen(hours), 10_000)
    t_h = t_s / 3600.0

    from .params import num
    circ = I_circadian(t_s, wc=0.000073, Ac=1.0)
    food = np.array([TRE(ti, food_rate) for ti in t_s])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 5), sharex=True)
    fig.suptitle(f"Input signals – {hours:.0f} h", fontsize=12)

    ax1.plot(t_h, circ, color="goldenrod", lw=1.5)
    ax1.axhline(1, color="gray", lw=0.5, ls="--")
    ax1.set_ylabel("Circadian signal\n(normalised)", fontsize=9)

    ax2.fill_between(t_h, food, color="sienna", alpha=0.7)
    ax2.set_ylabel("Food intake rate\n(TRE, mmol/s)", fontsize=9)
    ax2.set_xlabel("Time (hours)")

    plt.tight_layout()
    return fig


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

