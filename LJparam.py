from vdwp.physconst import NkB, NA, kB
from collections import namedtuple

# LJ parameters

# Our paper MT2022
gastable = """
Methane  3.758 148.6 $\\mathrm{CH}_4$
Ethane   4.520 208.8
C2H4     4.232 205 $\\mathrm{C_2H_4}$
Ne       2.749  35.6
Ar       3.405 119.8
Kr       3.60  171.0
Xe       4.047 231.0
Br2      4.93   540 Q                        # manually set. Becomes sIII at 10 bar.
CO2      4.486 189.0 $\\mathrm{CO}_2$
CS2      4.438 488 $\\mathrm{CS_2}$
cC3H6    4.582733199595731 301.51400454201365 $c\\mathrm{C_3H_6}$ # from CP
""".splitlines()

LJ = namedtuple("LennardJones", ["sig", "epsK", "TeX"])
tip4pice = LJ(sig=3.1668, epsK=106.1, TeX=r"$\mathrm{tip4pice}$")

gases = dict()
inter = dict()
for gas in gastable:
    i = gas.find("#")
    if i >= 0:
        gas = gas[:i]
    cols = gas.rstrip().split(maxsplit=3)
    if len(cols) == 0:
        continue
    if len(cols) == 3:
        tex = r"$\mathrm{" + cols[0] + r"}$"
        cols.append(tex)
    # print(len(cols))
    name, sig, epsK, tex = cols
    gases[name] = LJ(sig=float(sig), epsK=float(epsK), TeX=tex)
    inter[name] = LJ(
        sig=(float(sig) + tip4pice.sig) / 2,
        epsK=(float(epsK) * tip4pice.epsK) ** 0.5,
        TeX=tex,
    )

if __name__ == "__main__":

    # MPa, K
    cp = """
    Methane 4.595 190.55
    Ethane 4.9 305.3
    Ne 2.76 44.4
    Ar 4.898 150.87
    Kr 5.5 209.41
    Xe 5.841 289.77
    Br2 10.34 588
    CO2 7.38 304.19
    cPropane 5.579700 398.3
    cPentane 4.44 511
    """.splitlines()
    # cPentane https://www.chemeo.com/cid/29-704-8/Cyclopentane

    # Estimate from the critical points
    for gas in cp:
        cols = gas.rstrip().split(maxsplit=3)
        if len(cols) == 0:
            continue
        if len(cols) == 1:
            print(cols)
        name, pc, tc = cols
        pc = float(pc) * 1e6
        tc = float(tc)
        Tc = 1.321  # T*, T*=kT/eps これで温度からepsを算出できる。
        eps = tc / Tc  # in kJ/mol
        Pc = 0.129  # P*, P* = P sig**3 / eps
        # epsをkJ/molからJ / moleculeに換算。これはそれらしい数値。
        sig = (Pc * eps * kB / pc) ** (1 / 3) * 1e10
        print(sig, eps)
        print(
            f"{name} sig: From CP {sig} ours {gases[name].sig} ratio {gases[name].sig/sig}"
        )
        print(
            f"{name} eps: From CP {eps} ours {gases[name].epsK} ratio {gases[name].epsK/eps}"
        )
