#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from itertools import combinations
from logging import getLogger, DEBUG, INFO, basicConfig

import numpy as np
import matplotlib.pyplot as plt

from vdwp.physconst import NkB, NA
from vdwp import drawLine


# basicConfig(level=DEBUG, format="%(levelname)s %(message)s")
basicConfig(level=INFO, format="%(levelname)s %(message)s")
logger = getLogger()
logger.debug("Debug mode.")


def plot_phase_diagram_Figure2(axes):
    # Karttunen, A. J., Fässler, T. F., Linnolahti, M. & Pakkanen, T. A. Structural principles of semiconducting Group 14 clathrate frameworks. Inorg. Chem. 50, 1733–1742 (2011)
    # Table 2
    elements = ["C", "Si", "Ge", "Sn"]
    karttunen = {element: dict() for element in elements}
    density = {element: dict() for element in elements}
    for line in """
    I       0.15 0.10 0.05 0.04   6.2 2.8 2.4 1.8 3.08 2.00 4.66 5.02
    II      0.11 0.08 0.04 0.03   6.0 2.8 1.9 1.4 3.06 1.98 4.61 4.96
    III     0.15 0.10 0.06 0.04   6.1 2.7 2.2 1.6 3.07 2.00 4.64 4.99
    IV      0.17 0.12 0.07 0.05   6.1 2.7 2.0 1.4 3.05 1.99 4.61 4.97
    V       0.12 0.09 0.04 0.03   6.0 2.8 1.9 1.4 3.05 1.98 4.61 4.96
    II-4H   0.12 0.09 0.04 0.03   6.0 2.8 1.9 1.4 3.05 1.98 4.61 4.96
    II+IV-a 0.20 0.13 0.08 0.06   5.9 2.6 1.5 1.0 3.04 1.98 4.60 4.96
    II+IV-b 0.15 0.10 0.06 0.04   6.0 2.7 2.0 1.5 3.05 1.98 4.61 4.96
    """.splitlines():
        cols = line.rstrip().split()
        if len(cols) > 0:
            phase = cols.pop(0)
            for i, element in enumerate(elements):
                karttunen[element][phase] = float(cols[i])
                density[element][phase] = float(cols[i + 8])

    print(karttunen)

    ratios = {
        "I": (1.0, 0.0, 0.0),
        "II": (0.0, 1.0, 0.0),
        "IV": (0.0, 0.0, 1.0),
        "III": (23 / 43, 0.0, 20 / 43),
        "V": (0.0, 1.0, 0.0),
        "II-4H": (0.0, 1.0, 0.0),
        "II+IV-a": (0, 0.739130435, 0.260869565),
        "II+IV-b": (0, 0.459459459, 0.540540541),
    }

    ((ax1, ax2), (ax3, ax4)) = axes

    xra = (-0.02, 0.06)
    yra = (0.0, 0.08)

    for element, panel in zip(elements, (ax1, ax2, ax3, ax4)):
        structures = karttunen[element].keys()
        colors = ("black", "red", "blue", "green", "pink", "orange", "purple", "brown")
        for x in np.linspace(*xra, 40):
            for y in np.linspace(*yra, 40):
                smin = -1
                emin = 1e10
                for i, s in enumerate(structures):
                    e = karttunen[element][s] + ratios[s][0] * x + ratios[s][1] * y
                    if e < emin:
                        emin = e
                        smin = i
                panel.plot(x, y, "o", color=colors[smin])

    for element, panel, label in zip(
        elements, (ax1, ax2, ax3, ax4), ("(a)", "(b)", "(c)", "(d)")
    ):
        structures = karttunen[element].keys()
        for s1, s2 in combinations(structures, 2):
            panel.set_xlim(*xra)
            panel.set_ylim(*yra)
            panel.set_xlabel(r"$\Delta\mu_c^{I} - \Delta\mu_c^{IV}$")
            panel.set_ylabel(r"$\Delta\mu_c^{II} - \Delta\mu_c^{IV}$")
            A = ratios[s1][0] - ratios[s2][0]
            B = ratios[s1][1] - ratios[s2][1]
            C = karttunen[element][s1] - karttunen[element][s2]
            drawLine(
                A,
                B,
                C,
                style="-",
                ax=panel,
                xtick=np.linspace(*xra, 100),
                ytick=np.linspace(*yra, 100),
            )
        panel.annotate(
            f"{label} {element}",  # this is the text
            (0.8, 0.2),  # these are the coordinates to position the label
            xycoords="axes fraction",  # how to position the text
            ha="center",
        )  # horizontal alignment can be left, right or center


if __name__ == "__main__":
    # plt.rcParams['text.usetex'] = True
    plt.rcParams["font.size"] = 14
    plt.rcParams["font.family"] = "sans-serif"

    fig = plt.figure(figsize=(7, 7))
    gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
    axes = gs.subplots(sharex="col", sharey="row")
    plot_phase_diagram_Figure2(axes)
    plt.show()
    fig.savefig("Figure2.pdf")
