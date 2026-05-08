import numpy as np
import matplotlib.pyplot as plt
from typing import Callable

from circadian_feeding.simulate import SimResult
from circadian_feeding.params import num, gen, break_, lun, din
from circadian_feeding.food import TRE
from circadian_feeding import plot

def main(res: SimResult):
    plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7})
    fig, axes = plt.subplots(10, 1, figsize=(10, 15), sharex=True)
    fig.tight_layout()
    twin = axes[0].twinx()
    plot.plot_inputs(axes=[axes[0], twin],
                food_fn=lambda t, food_rate: TRE(t, food_rate),
                t_range=(0.0,72.0))
    plot.plot_hormones(res=res,
                  axes=axes[1:6],
                  hormones=["GI1", "GI2", "H1", "IN", "LP"],
                  meal_spans=[],
                  t_range=(0.0,72.0))
    plot.plot_spike_counts(res=res,
                  axes=axes[6:],
                  neurons=["A", "P", "V", "L"],
                  meal_spans=[],
                  t_range=(0.0,72.0))
    
    for ax in list(axes) + [twin]:
        ax.set_xlabel("")
        if ax.get_lines() and ax != axes[0]:
            ax.get_lines()[0].set_color('black')
        for patch in ax.patches:
            patch.set_facecolor('steelblue')

    axes[3].set_ylim(top=1200)
    axes[4].set_ylim(top=600)
    axes[5].set_ylim(top=40)
        
    axes[-1].set_xlabel("Time (hours)")
    
    return fig, axes
