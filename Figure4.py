#!/usr/bin/env python3
# -*- coding: utf-8 -*-


from logging import getLogger, DEBUG, INFO, basicConfig
from collections import namedtuple

import numpy as np
import matplotlib.pyplot as plt

from vdwp.physconst import NkB, NA
import vdwp.vdWP as vdWP
import vdwp.crystals as crystals
import vdwp.chempot as chempot
from ljd.ljd import fvalue
from LJparam import gases, inter, tip4pice


# basicConfig(level=DEBUG, format="%(levelname)s %(message)s")
basicConfig(level=INFO, format="%(levelname)s %(message)s")
logger = getLogger()
logger.debug("Debug mode.")


def determine_stable_phase(mu_e, Deltamu, structures):
    mumin = 0
    stmin = None
    for s in structures:
        if mu_e[s] + Deltamu[s] < mumin:
            mumin = mu_e[s] + Deltamu[s]
            stmin = s
    return stmin


"""
Nomenclature for variables
f, mu: free energy/chemical potential
h    : enthalpy

_w   : of water
_i   : of ice
_e   : of the empty lattice
_f   : of the filled lattice
_g   : of the guest
_c   : in the cage
"""


def calculate_mixed_clathrate_phase(
    g1,
    g2,
    beta,
    pressure,
    temperatures,
    mu_e,
    structures,
    composition_ratios=np.linspace(0.0, 1.0, 100),
    target_phase=None,
):
    f1 = dict()
    f2 = dict()
    for cage, radius in crystals.radii.items():
        # sigma and epsilon must be the intermolecular ones.
        f1[cage] = fvalue(
            {radius: crystals.nmemb[cage]}, g1.sig, g1.epsK * 8.314 / 1000, beta
        )
        f2[cage] = fvalue(
            {radius: crystals.nmemb[cage]}, g2.sig, g2.epsK * 8.314 / 1000, beta
        )

    phases = []
    for r in composition_ratios:
        p1 = (1 - r) * pressure
        p2 = r * pressure

        if p1 == 0:
            mu2 = chempot.ideal_gas_chemical_potential(
                temperatures, p2
            ) + chempot.phase_space_integration_correction(
                temperatures, degrees_of_freedom=0
            )
            Deltamu = vdWP.calculate_chemical_potential_by_occupation(
                temperatures, f2, mu2, structures
            )
        elif p2 == 0:
            mu1 = chempot.ideal_gas_chemical_potential(
                temperatures, p1
            ) + chempot.phase_space_integration_correction(
                temperatures, degrees_of_freedom=0
            )
            Deltamu = vdWP.calculate_chemical_potential_by_occupation(
                temperatures, f1, mu1, structures
            )
        else:
            mu1 = chempot.ideal_gas_chemical_potential(
                temperatures, p1
            ) + chempot.phase_space_integration_correction(
                temperatures, degrees_of_freedom=0
            )
            mu2 = chempot.ideal_gas_chemical_potential(
                temperatures, p2
            ) + chempot.phase_space_integration_correction(
                temperatures, degrees_of_freedom=0
            )
            Deltamu = vdWP.calculate_chemical_potential_by_occupation(
                temperatures, (f1, f2), (mu1, mu2), structures
            )

        phase = determine_stable_phase(mu_e, Deltamu, structures)
        phases.append(phase)
        if target_phase == phase:
            break
    return phases, phase


def plot_phase_marker(sig, epsK, phases, lastphase, ax):
    markers = {1: "o", 2: "+", 3: "^"}
    if lastphase == "CS1":
        color = "lightgreen"
    elif lastphase == "CS2":
        color = "#88f"
    else:
        color = "brown"
    # Distinguish I-III-I case (not actually found)
    if lastphase == "CS1" and "TS1" in phases and len(phases) == 2:
        ax.plot(sig, epsK, color="orange", marker=".")
    else:
        marker = markers[len(phases)]
        if len(phases) == 3 and lastphase == "CS1":
            color = "#0c0"  # slightly darker triangle
        ax.plot(sig, epsK, color=color, marker=marker)
        if len(phases) == 3:
            print(sig, epsK, phases, lastphase)


def plot_guest_molecules(ax, gases):
    """
    Put cities (i.e. locus of the molecules) on the map.
    """
    for name, gas in gases.items():
        sig = gas.sig
        eps = gas.epsK
        ax.plot(sig, eps, "ok")
        if name in ("CO2", "n-Butane"):
            ha = "right"
            xytext = (-5, 0)
        elif name == "N2O":
            ha = "left"
            xytext = (5, 0)
        else:
            ha = "center"
            xytext = (0, 5)
        ax.annotate(
            gas.TeX,  # this is the text
            (sig, eps),  # these are the coordinates to position the label
            textcoords="offset points",  # how to position the text
            xytext=xytext,  # distance from text to points (x,y)
            ha=ha,
        )  # horizontal alignment can be left, right or center


def plot_phase_contours(ax, X, Y, Z, Z2):
    # Phases of the 2nd component
    contours = ax.contour(X, Y, Z2, levels=[1.5, 2.5], colors="black", linewidths=2)
    # Number of phases by mixing
    contours = ax.contour(X, Y, Z, levels=[1.5, 2.5], colors="black", linewidths=1)
    cs = ax.contourf(
        X, Y, Z2, levels=[1.0, 2.0, 3.0], colors=["#8f8", "#ccf", "#fcc"], extend="min"
    )
    cs = ax.contourf(
        X,
        Y,
        Z,
        levels=[0.0, 1.0, 2.0],
        hatches=["", ".", "//"],
        colors="none",
        extend="max",
    )
    for collection in cs.get_children():
        collection.set_linewidth(0.0)
        collection.set_edgecolor("white")


def plot_phase_diagram_Figure4a(
    ax,
    temperatures,
    pressure,
    beta,
    mu_e,
    structures,
    inter,
    tip4pice,
    gases,
    composition_ratios,
    xyticks,
):
    # (a)
    # Calculate the minimum concentration that, when mixed with Q, produces II.
    guest = "Br2"
    ax.set_xlabel(r"$\sigma_g / \AA$")
    ax.set_ylabel(r"$\epsilon_g / K$")
    ax.set_xlim(3.45, 5.1)
    ax.set_ylim(120.0, 600)
    x = np.linspace(3.45, 5.1, xyticks * 2)
    y = np.linspace(120.0, 600, xyticks * 2)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    for ix, sig in enumerate(x):
        for iy, epsK in enumerate(y):
            inter2 = LJ(
                sig=(tip4pice.sig + sig) / 2, epsK=(tip4pice.epsK * epsK) ** 0.5
            )
            phases, lastphase = calculate_mixed_clathrate_phase(
                inter2,
                inter[guest],
                beta,
                pressure,
                temperatures,
                mu_e,
                structures,
                composition_ratios=composition_ratios,
                target_phase="CS2",
            )
            z = 1.0
            for i, phase in enumerate(phases):
                if phase == "CS2":
                    z = composition_ratios[i]
                    break
            Z[iy, ix] = z

    contours = ax.contour(
        X, Y, Z, levels=[0.0001, 0.001, 0.01, 0.1], colors="black", linewidths=1
    )
    ax.clabel(contours, inline=True, fontsize=8)
    contours = ax.contour(
        X,
        Y,
        Z,
        levels=[
            0.9999,
        ],
        colors="black",
        linestyles="dashed",
    )

    contours = ax.contour(
        X,
        Y,
        Z,
        levels=[
            0.0,
        ],
        colors="black",
        linewidths=2,
    )

    cs = ax.contourf(
        X,
        Y,
        Z,
        levels=[0.0, 0.0001, 0.001, 0.01, 0.1, 0.9999, 2.0],
        colors=["#ccf", "#dfd", "#cfc", "#bfb", "#afa", "#9f9", "#fcc"],
        extend="min",
    )

    plot_guest_molecules(ax=ax, gases=gases)

    ax.annotate(
        "(a) X + Q",  # this is the text
        xy=(0.8, 0.92),  # these are the coordinates to position the label
        xycoords="axes fraction",
        fontsize=20,
        ha="right",
    )


def plot_phase_diagram_Figure4b(
    ax,
    temperatures,
    pressure,
    beta,
    mu_e,
    structures,
    inter,
    tip4pice,
    gases,
    composition_ratios,
    xyticks,
):
    # (b)
    # Search for conditions for I-II-III-I transition when mixed with Me.
    guest = "Methane"
    ax.set_xlabel(r"$\sigma_g / \AA$")
    ax.set_xlim(3.45, 5.1)
    ax.set_ylim(120.0, 600)
    x = np.linspace(3.45, 5.1, xyticks * 2)
    y = np.linspace(120.0, 600, xyticks * 2)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    Z2 = np.zeros_like(X)
    for ix, sig in enumerate(x):
        for iy, epsK in enumerate(y):
            inter2 = LJ(
                sig=(tip4pice.sig + sig) / 2, epsK=(tip4pice.epsK * epsK) ** 0.5
            )
            phases, lastphase = calculate_mixed_clathrate_phase(
                inter[guest],
                inter2,
                beta,
                pressure,
                temperatures,
                mu_e,
                structures,
                composition_ratios=composition_ratios,
            )
            Z[iy, ix] = len(set(phases))
            if lastphase == "CS1":
                z = 1
            elif lastphase == "CS2":
                z = 2
            else:
                z = 3
            Z2[iy, ix] = z
    plot_guest_molecules(ax=ax, gases=gases)
    plot_phase_contours(ax, X, Y, Z, Z2)
    ax.annotate(
        "(b) Me + X",  # this is the text
        xy=(0.6, 0.92),  # these are the coordinates to position the label
        xycoords="axes fraction",
        fontsize=20,
        ha="right",
    )


def plot_phase_diagram_Figure4c(
    ax,
    temperatures,
    pressure,
    beta,
    mu_e,
    structures,
    inter,
    tip4pice,
    gases,
    composition_ratios,
    xyticks,
):
    # (c)
    # Search for conditions for I-II-III-I transition when mixed with Xe.
    guest = "Xe"
    ax.set_xlabel(r"$\sigma_g / \AA$")
    ax.set_xlim(4.55, 5.05)
    ax.set_ylim(120.0, 600)
    x = np.linspace(4.55, 5.05, xyticks * 2)
    y = np.linspace(120.0, 600, xyticks * 2)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)
    Z2 = np.zeros_like(X)
    for ix, sig in enumerate(x):
        for iy, epsK in enumerate(y):
            inter2 = LJ(
                sig=(tip4pice.sig + sig) / 2, epsK=(tip4pice.epsK * epsK) ** 0.5
            )
            phases, lastphase = calculate_mixed_clathrate_phase(
                inter[guest],
                inter2,
                beta,
                pressure,
                temperatures,
                mu_e,
                structures,
                composition_ratios=composition_ratios,
            )
            Z[iy, ix] = len(set(phases))
            if lastphase == "CS1":
                z = 1
            elif lastphase == "CS2":
                z = 2
            else:
                z = 3
            Z2[iy, ix] = z
    plot_guest_molecules(ax=ax, gases=gases)
    plot_phase_contours(ax, X, Y, Z, Z2)
    ax.annotate(
        "(c) Xe + X",  # this is the text
        xy=(0.4, 0.92),  # these are the coordinates to position the label
        xycoords="axes fraction",
        fontsize=20,
        ha="right",
    )


if __name__ == "__main__":
    # User variables
    pressure = 101325.00 * 50  # Pa
    temperatures = 273.15
    beta = 1.0 / (NkB * temperatures)

    LJ = namedtuple("LJ", ["sig", "epsK"])

    structures = ["CS2", "CS1", "HS1", "TS1"]

    ####### chemical potential of gas ######################################
    stericterm = chempot.molecular_chemical_potential_corrections(
        temperatures, mass=1.0, symmetry_number=1, moment_of_inertia=(0, 0, 0)
    )

    mu_g = (
        chempot.ideal_gas_chemical_potential(temperatures, pressure)
        + chempot.phase_space_integration_correction(temperatures, degrees_of_freedom=0)
        + stericterm
    )

    ####### structure-dependent terms ######################################

    mu_e = crystals.mu_e

    composition_ratios = np.concatenate(
        [np.array([0.0]), np.logspace(-4, 0.0, 300)]
    )  # 1e-4 .. 1e0

    plt.rcParams["font.size"] = 14
    plt.rcParams["font.family"] = "sans-serif"

    fig, axes = plt.subplots(
        nrows=1, ncols=3, figsize=(12, 5), sharey=True, gridspec_kw={"wspace": 0}
    )

    xyticks = 40

    plot_phase_diagram_Figure4a(
        axes[0],
        temperatures,
        pressure,
        beta,
        mu_e,
        structures,
        inter,
        tip4pice,
        gases,
        composition_ratios,
        xyticks,
    )
    plot_phase_diagram_Figure4b(
        axes[1],
        temperatures,
        pressure,
        beta,
        mu_e,
        structures,
        inter,
        tip4pice,
        gases,
        composition_ratios,
        xyticks,
    )
    plot_phase_diagram_Figure4c(
        axes[2],
        temperatures,
        pressure,
        beta,
        mu_e,
        structures,
        inter,
        tip4pice,
        gases,
        composition_ratios,
        xyticks,
    )
    plt.tight_layout()

    plt.show()
    fig.savefig("Figure4.pdf")
    fig.savefig("Figure4.png")
