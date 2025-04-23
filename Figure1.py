#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from itertools import combinations
from logging import getLogger, DEBUG, INFO, basicConfig

import matplotlib.pyplot as plt

import vdwp.vdWP as vdWP
import vdwp.crystals as crystals
from vdwp.physconst import NkB, NA
import vdwp.chempot as chempot
from vdwp import drawLine

# basicConfig(level=DEBUG, format="%(levelname)s %(message)s")
basicConfig(level=INFO, format="%(levelname)s %(message)s")
logger = getLogger()
logger.debug("Debug mode.")


# Nomenclature for variables
# f, mu: free energy/chemical potential
# h    : enthalpy
#
# _w   : of water
# _i   : of ice
# _e   : of the empty lattice
# _f   : of the filled lattice
# _g   : of the guest
# _c   : in the cage
#

plt.rcParams["font.size"] = 14
plt.rcParams["font.family"] = "sans-serif"

figure = plt.figure(figsize=(5, 5))

plt.xlim(-0.1, 0.5)
plt.ylim(-0.1, 0.5)
plt.xticks([0, 0.2, 0.4])
plt.xlabel(
    r"$(\Delta\mu_c^\mathrm{CS1} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$"
)
plt.ylabel(
    r"$(\Delta\mu_c^\mathrm{CS2} - \Delta\mu_c^\mathrm{HS1}) / \mathrm{kJ~mol}^{-1}$"
)
plt.axis("square")

temperatures = 273.15

####### chemical potential of gas ######################################
stericterm = chempot.molecular_chemical_potential_corrections(
    temperatures, mass=1.0, symmetry_number=1, moment_of_inertia=(0, 0, 0)
)

####### structure-dependent terms ######################################

mu_e = crystals.mu_e

for s1, s2 in combinations(crystals.names, 2):
    # The equation for coexistence lines
    # Ax + By + C = 0
    # A = xA - xB
    # B = yA - yB
    # C = mu_eA- mu_eB
    A = crystals.ratios[s1][0] - crystals.ratios[s2][0]
    B = crystals.ratios[s1][1] - crystals.ratios[s2][1]
    C = mu_e[s1] - mu_e[s2]
    drawLine(A, B, C, style="-k")

# manually labelled
plt.annotate(
    "I (CS-I)",  # this is the text
    xy=(0.2, 0.9),  # these are the coordinates to position the label
    xycoords="axes fraction",
    fontsize=24,
    ha="center",
)
plt.annotate(
    "II (CS-II)",  # this is the text
    xy=(0.8, 0.4),  # these are the coordinates to position the label
    xycoords="axes fraction",
    fontsize=24,
    ha="center",
)
plt.annotate(
    "III\n(TS-I)",  # this is the text
    xy=(0.57, 0.8),  # these are the coordinates to position the label
    xycoords="axes fraction",
    fontsize=24,
    ha="center",
)
plt.annotate(
    "IV\n(HS-I)",  # this is the text
    xy=(0.85, 0.8),  # these are the coordinates to position the label
    xycoords="axes fraction",
    fontsize=24,
    ha="center",
)

plt.tight_layout()
plt.show()
figure.savefig("Figure1.pdf")
