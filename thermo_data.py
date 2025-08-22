"""
thermo_data.py

Big, pragmatic table of standard thermochemical data for quick homework/app use.

Units
- ΔH°f (Hf): kJ/mol
- ΔG°f (Gf): kJ/mol
- S° (S): J/(mol·K)
- Temperature: 298.15 K unless noted; pressure 1 bar.

Notes
- Elements in their standard states have ΔH°f = ΔG°f = 0 by convention.
- Aqueous ions often omit S° (ionic convention); that's normal.
- Values vary a little between tables (± a few kJ/mol). If you need CRC-level precision,
  consider adding your own CSV (see load_additional_from_csv).
"""

from __future__ import annotations
from typing import Dict, Tuple, Optional, List, Iterable
import re
import csv
import math

# --------------------------------------------------------------------------------------
# CORE TABLE
# Key: (formula_str, phase_str) -> {"Hf": float, "Gf": float, "S": float}
# Phase: "g", "l", "s", "aq"
# --------------------------------------------------------------------------------------

THERMO: Dict[Tuple[str, str], Dict[str, float]] = {
    # ===== Elements (standard states) ==================================================
    ("H2","g"):  {"Hf":0.0,"Gf":0.0,"S":130.68},
    ("O2","g"):  {"Hf":0.0,"Gf":0.0,"S":205.15},
    ("N2","g"):  {"Hf":0.0,"Gf":0.0,"S":191.61},
    ("F2","g"):  {"Hf":0.0,"Gf":0.0,"S":202.78},
    ("Cl2","g"): {"Hf":0.0,"Gf":0.0,"S":223.08},
    ("Br2","l"): {"Hf":0.0,"Gf":0.0,"S":152.20},
    ("Br2","g"): {"Hf":30.91,"Gf":3.14,"S":245.46},  # Br2 gas form
    ("I2","s"):  {"Hf":0.0,"Gf":0.0,"S":116.10},
    ("I2","g"):  {"Hf":62.39,"Gf":19.38,"S":260.60}, # I2 gas form
    ("C","s"):   {"Hf":0.0,"Gf":0.0,"S":5.74},     # graphite
    ("S","s"):   {"Hf":0.0,"Gf":0.0,"S":31.80},
    ("Fe","s"):  {"Hf":0.0,"Gf":0.0,"S":27.30},
    ("Na","s"):  {"Hf":0.0,"Gf":0.0,"S":51.00},
    ("K","s"):   {"Hf":0.0,"Gf":0.0,"S":64.70},

    # ===== Monatomic gases (handy for atomization) ====================================
    ("H","g"):   {"Hf":218.00,"S":114.7},
    ("O","g"):   {"Hf":249.20,"Gf":231.70,"S":161.1},
    ("Na","g"):  {"Hf":107.50,"Gf":77.00,"S":153.7},
    ("K","g"):   {"Hf":89.99,"Gf":61.17,"S":160.2},
    ("Cl","g"):  {"Hf":121.70,"Gf":105.30,"S":165.2},

    # ===== Very common molecular species ==============================================
    ("H2O","l"): {"Hf":-285.83,"Gf":-237.13,"S":69.91},
    ("H2O","g"): {"Hf":-241.82,"Gf":-228.57,"S":188.83},
    ("CO2","g"): {"Hf":-393.51,"Gf":-394.36,"S":213.79},
    ("CO","g"):  {"Hf":-110.53,"Gf":-137.17,"S":197.66},
    ("O3","g"):  {"Hf":142.67,"Gf":163.20,"S":238.92},
    ("NO","g"):  {"Hf":90.25,"Gf":86.55,"S":210.76},
    ("NO2","g"): {"Hf":33.18,"Gf":51.84,"S":240.06},
    ("N2O","g"): {"Hf":82.05,"Gf":104.20,"S":219.90},
    ("NH3","g"): {"Hf":-46.11,"Gf":-16.45,"S":192.77},
    ("CH4","g"): {"Hf":-74.81,"Gf":-50.80,"S":186.26},
    ("C2H6","g"):{"Hf":-84.68,"Gf":-32.89,"S":229.49},
    ("C3H8","g"):{"Hf":-103.85,"Gf":-23.49,"S":269.91},
    ("C2H2","g"):{"Hf":226.73,"Gf":209.20,"S":200.94},
    ("H2S","g"): {"Hf":-20.60,"Gf":-33.56,"S":205.70},
    ("SO2","g"): {"Hf":-296.84,"Gf":-300.19,"S":248.22},
    ("SO3","g"): {"Hf":-395.72,"Gf":-371.06,"S":256.77},
    ("HCl","g"): {"Hf":-92.31,"Gf":-95.30,"S":186.90},
    ("HF","g"):  {"Hf":-271.10,"Gf":-273.30,"S":173.78},
    ("HBr","g"): {"Hf":-36.29,"Gf":-53.14,"S":198.70},
    ("HI","g"):  {"Hf":26.46,"Gf":1.30,"S":206.50},
    ("HCN","g"): {"Hf":135.10,"Gf":124.70,"S":201.60},

    # ===== Aromatics and related organics =============================================
    ("C6H6","l"):   {"Hf":49.04,"Gf":124.5,"S":173.3},    # benzene (liquid)
    ("C6H6","g"):   {"Hf":82.93,"Gf":124.5,"S":269.2},    # benzene (gas)
    ("C7H8","l"):   {"Hf":12.0,"Gf":124.7,"S":186.3},     # toluene (liquid)
    ("C7H8","g"):   {"Hf":50.1,"Gf":124.7,"S":282.7},     # toluene (gas)
    ("C8H10","l"):  {"Hf":-7.0,"Gf":117.0,"S":205.0},     # xylene (liquid, approx)
    ("C8H10","g"):  {"Hf":32.0,"Gf":117.0,"S":310.0},     # xylene (gas, approx)
    ("C6H5OH","l"): {"Hf":-165.0,"Gf":-44.0,"S":144.0},   # phenol (liquid)
    ("C6H5OH","g"): {"Hf":-110.0,"Gf":-44.0,"S":250.0},   # phenol (gas, approx)
    ("C6H5NH2","l"):{"Hf":49.0,"Gf":73.0,"S":151.0},      # aniline (liquid)
    ("C6H5NH2","g"):{"Hf":97.0,"Gf":73.0,"S":260.0},      # aniline (gas, approx)
    ("C6H5CHO","l"):{"Hf":-24.0,"Gf":-8.0,"S":170.0},     # benzaldehyde (liquid)
    ("C6H5CHO","g"):{"Hf":20.0,"Gf":-8.0,"S":280.0},      # benzaldehyde (gas, approx)
    ("C6H5COOH","s"):{"Hf":-385.0,"Gf":-314.0,"S":197.0}, # benzoic acid (solid)
    ("C6H5COOH","g"):{"Hf":-320.0,"Gf":-314.0,"S":320.0}, # benzoic acid (gas, approx)
    ("C6H5NO2","s"):{"Hf":-9.0,"Gf":-4.0,"S":180.0},      # nitrobenzene (solid)
    ("C6H5NO2","g"):{"Hf":40.0,"Gf":-4.0,"S":290.0},      # nitrobenzene (gas, approx)
    ("C10H8","s"):  {"Hf":78.7,"Gf":101.0,"S":170.0},     # naphthalene (solid)
    ("C10H8","g"):  {"Hf":126.0,"Gf":101.0,"S":320.0},    # naphthalene (gas, approx)
    ("C9H8","l"):   {"Hf":-12.0,"Gf":90.0,"S":220.0},     # styrene (liquid, approx)
    ("C9H8","g"):   {"Hf":30.0,"Gf":90.0,"S":340.0},      # styrene (gas, approx)

    # ---- Nitrosyl chloride (explicit request) ----
    ("NOCl","g"):{"Hf":51.71,"Gf":66.08,"S":261.68},  # NIST/Janaf + UWisc

    # ===== Liquids / extra solvents & organics ========================================
    ("CH3OH","l"):{"Hf":-238.7,"Gf":-166.2,"S":126.8},
    ("CH3OH","g"):{"Hf":-201.0,"Gf":-162.0,"S":239.8},
    ("C2H5OH","l"):{"Hf":-277.7,"Gf":-174.8,"S":160.7},
    ("C2H5OH","g"):{"Hf":-235.1,"Gf":-168.6,"S":282.6},

    # ===== Common solids (salts, oxides, carbonates, etc.) =============================
    ("NaCl","s"):   {"Hf":-411.12,"Gf":-384.14,"S":72.11},
    ("KCl","s"):    {"Hf":-436.7,"Gf":-408.5,"S":82.6},   # Gf from EngToolbox
    ("KBr","s"):    {"Hf":-392.0,"Gf":-380.7},           # typical table value (Gf from EngToolbox)
    ("NaBr","s"):   {"Hf":-361.0,"Gf":-354.7},
    ("NaI","s"):    {"Hf":-287.1,"Gf":-284.6},
    ("KI","s"):     {"Hf":-328.0,"Gf":-325.5},
    ("AgCl","s"):   {"Hf":-127.0,"Gf":-109.8,"S":96.3},  # EngToolbox
    ("AgBr","s"):   {"Hf":-100.4,"Gf":-96.9,"S":107.1},
    ("AgI","s"):    {"Hf":-61.8,"Gf":-66.2,"S":115.5},

    ("CaCO3","s"):  {"Hf":-1206.9,"Gf":-1128.8,"S":92.90}, # calcite ~
    ("CaO","s"):    {"Hf":-634.9,"Gf":-603.3,"S":38.1},    # EngToolbox
    ("MgO","s"):    {"Hf":-601.6,"Gf":-596.3,"S":26.9},
    ("Al2O3","s"):  {"Hf":-1675.7,"Gf":-1582.3,"S":50.9},
    ("Fe2O3","s"):  {"Hf":-824.2,"Gf":-742.2,"S":87.4},    # representative; small table variations exist
    ("SiO2","s"):   {"Hf":-910.9,"Gf":-856.4,"S":41.5},    # quartz
    ("CaSO4","s"):  {"Hf":-1434.5,"Gf":-1322.0,"S":106.5}, # EngToolbox
    ("Na2SO4","s"): {"Hf":-1387.1,"Gf":-1270.2,"S":149.6},
    ("NaNO3","s"):  {"Hf":-467.9,"Gf":-367.1,"S":116.5},
    ("KNO3","s"):   {"Hf":-494.5,"Gf":-394.7,"S":132.5},
    ("CaF2","s"):   {"Hf":-1228.0,"Gf":-1175.6,"S":68.5},

    # ===== Aqueous ions (298 K) =======================================================
    ("H+","aq"):     {"Hf":0.0,"Gf":0.0,"S":0.0},           # convention
    ("OH-","aq"):    {"Hf":-230.0,"Gf":-157.24,"S":-10.75},
    ("F-","aq"):     {"Hf":-329.1},
    ("Cl-","aq"):    {"Hf":-167.16,"Gf":-131.23},
    ("Br-","aq"):    {"Hf":-121.4,"Gf":-104.0},
    ("I-","aq"):     {"Hf":-55.9},
    ("NO3-","aq"):   {"Hf":-205.0,"Gf":-111.3},
    ("NO2-","aq"):   {"Hf":-104.6},
    ("SO4^2-","aq"): {"Hf":-909.3,"Gf":-744.6},
    ("HSO4-","aq"):  {"Hf":-885.8,"Gf":-738.3},
    ("SO3^2-","aq"): {"Hf":-635.1,"Gf":-557.1},
    ("HSO3-","aq"):  {"Hf":-587.0,"Gf":-530.0},
    ("HCO3-","aq"):  {"Hf":-691.99,"Gf":-586.85},
    ("CO3^2-","aq"): {"Hf":-677.1,"Gf":-527.9},
    ("CH3COO-","aq"):{"Hf":-486.0,"Gf":-369.4},            # acetate (representative)
    ("C2H5COO-","aq"):{"Hf":-546.0,"Gf":-425.0},           # propionate (representative)
    ("CN-","aq"):    {"Hf":150.0,"Gf":111.4},

    ("NH4+","aq"):   {"Hf":-132.5,"Gf":-79.31},
    ("NH3","aq"):    {"Hf":-80.29,"Gf":-26.57},
    ("Na+","aq"):    {"Hf":-240.12,"Gf":-261.9},
    ("K+","aq"):     {"Hf":-251.2,"Gf":-283.3},
    ("Ca^2+","aq"):  {"Hf":-542.8,"Gf":-553.6},
    ("Mg^2+","aq"):  {"Hf":-466.9,"Gf":-454.8},
    ("Zn^2+","aq"):  {"Hf":-153.9},
    ("Cu^2+","aq"):  {"Hf":  64.4},
    ("Al^3+","aq"):  {"Hf":-524.7},
    ("Fe^2+","aq"):  {"Hf": -89.1,"Gf": -78.9},
    ("Fe^3+","aq"):  {"Hf": -48.5,"Gf":  -4.1},

    # ===== Acids/bases in solution ====================================================
    ("HF","aq"):     {"Hf":-320.1,"Gf":-278.9},
    ("HNO3","aq"):   {"Hf":-207.0},
    ("H2SO4","l"):   {"Hf":-814.0},                         # conc. acid (l)
    ("HNO2","aq"):   {"Hf": -104.0,"Gf": -45.3},           # representative values
    ("HCl","aq"):    {"Hf":-167.2,"Gf":-131.2},

    # ===== Solids/liquids from your original table ====================================
    ("C6H12O6","s"): {"Hf":-1273.3,"Gf":-910.5,"S":212.10}, # glucose (approx)
}

# --------------------------------------------------------------------------------------
# OPTIONAL: more entries that are commonly handy (safely merged into THERMO)
# --------------------------------------------------------------------------------------
THERMO.update({
    # Nitrates / chlorates / perchlorates
    ("NaNO3","aq"): {"Hf":-468.8,"Gf":-373.2},
    ("KNO3","aq"):  {"Hf":-503.0,"Gf":-407.9},
    ("NaClO","s"):  {"Hf":-365.0,"Gf":-314.0},
    ("NaClO3","s"): {"Hf":-364.0,"Gf":-296.3},
    ("KClO3","s"):  {"Hf":-397.7,"Gf":-296.3},            # Gf per EngToolbox
    ("NaClO4","s"): {"Hf":-382.9,"Gf":-307.0},
    ("KClO4","s"):  {"Hf":-432.8,"Gf":-303.1},

    # Alkali/alkaline earth hydroxides (solid)
    ("NaOH","s"):   {"Hf":-425.6,"Gf":-379.4,"S":64.5},
    ("KOH","s"):    {"Hf":-424.7,"Gf":-379.1,"S":78.7},
    ("Ca(OH)2","s"):{"Hf":-985.2,"Gf":-897.5,"S":83.4},

    # Carbonates/bicarbonates (solid)
    ("Na2CO3","s"): {"Hf":-1130.7,"Gf":-1044.4,"S":136.0},
    ("NaHCO3","s"): {"Hf":-947.7,"Gf":-851.9,"S":101.7},
    ("K2CO3","s"):  {"Hf":-1150.0,"Gf":-1092.3,"S":154.0},

    # Some additional gases/liquids you may hit
    ("N2O4","g"):   {"Hf": 9.16,"Gf": 97.9,"S":304.3},
    ("H2O2","l"):   {"Hf":-187.8,"Gf":-120.4},
    ("H2O2","g"):   {"Hf":-136.1},  # often only Hf tabulated

    # Interhalogen compounds
    ("BrCl","g"):   {"Hf":8.6,"Gf":13.6,"S":220.0},
    ("ICl","g"):    {"Hf":107.2,"Gf":99.0,"S":256.0},
    ("ICl","s"):    {"Hf":106.6,"Gf":98.4,"S":151.0},
    ("IBr","g"):    {"Hf":66.2,"Gf":57.6,"S":282.0},
    ("IBr","s"):    {"Hf":65.7,"Gf":57.1,"S":164.0},
    ("ClF","g"):    {"Hf":-56.3,"Gf":-13.6,"S":222.0},
    ("BrF","g"):    {"Hf":13.6,"Gf":22.0,"S":247.0},
    ("BrF3","s"):   {"Hf":-189.0,"Gf":-176.0,"S":149.0},
    ("BrF5","s"):   {"Hf":-277.0,"Gf":-261.0,"S":176.0},
    ("IF","g"):     {"Hf":-37.6,"Gf":-29.0,"S":260.0},
    ("IF5","s"):    {"Hf":-361.0,"Gf":-340.0,"S":207.0},
    ("IF7","g"):    {"Hf":-437.0,"Gf":-413.0,"S":298.0},

    # More oxides
    ("NO","g"):     {"Hf":90.25,"Gf":86.55,"S":210.76},
    ("NO2","g"):    {"Hf":33.18,"Gf":51.84,"S":240.06},
    ("N2O","g"):    {"Hf":82.05,"Gf":104.20,"S":219.90},
    ("N2O3","g"):   {"Hf":83.7,"Gf":111.0,"S":281.0},
    ("N2O4","g"):   {"Hf":9.16,"Gf":97.9,"S":304.3},
    ("N2O5","s"):   {"Hf":11.3,"Gf":56.0,"S":140.0},
    ("SO2","g"):    {"Hf":-296.84,"Gf":-300.19,"S":248.22},
    ("SO3","g"):    {"Hf":-395.72,"Gf":-371.06,"S":256.77},
    ("P4O10","s"):  {"Hf":-2984.0,"Gf":-2940.0,"S":163.0},
    ("Cl2O","g"):   {"Hf":80.3,"Gf":101.0,"S":223.0},
    ("ClO2","g"):   {"Hf":101.0,"Gf":108.0,"S":261.0},
    ("Cl2O7","l"):  {"Hf":65.0,"Gf":80.0,"S":234.0},

    # More acids
    ("H2SO4","aq"): {"Hf":-909.27,"Gf":-744.6,"S":20.1},
    ("H3PO4","aq"): {"Hf":-1281.0,"Gf":-1122.0,"S":110.0},
    ("HClO4","aq"): {"Hf":-467.9,"Gf":-287.0,"S":110.0},
    ("H2CO3","aq"): {"Hf":-699.7,"Gf":-623.1,"S":93.0},
    ("HNO3","l"):   {"Hf":-207.4,"Gf":-79.9,"S":146.4},

    # More bases
    ("NH4OH","aq"): {"Hf":-80.29,"Gf":-26.57,"S":63.0},
    ("Ba(OH)2","s"):{"Hf":-1044.4,"Gf":-945.0,"S":110.0},
    ("Mg(OH)2","s"):{"Hf":-924.7,"Gf":-833.6,"S":62.0},

    # More salts
    ("BaSO4","s"):  {"Hf":-1465.0,"Gf":-1361.0,"S":148.0},
    ("PbSO4","s"):  {"Hf":-1017.0,"Gf":-813.0,"S":148.0},
    ("CuSO4","s"):  {"Hf":-769.0,"Gf":-661.0,"S":110.0},
    ("ZnSO4","s"):  {"Hf":-980.0,"Gf":-879.0,"S":120.0},
    ("FeSO4","s"):  {"Hf":-928.0,"Gf":-740.0,"S":110.0},
    ("Al2(SO4)3","s"):{"Hf":-3445.0,"Gf":-3160.0,"S":343.0},

    # More organics
    ("C2H4","g"):   {"Hf":52.47,"Gf":68.11,"S":219.56},
    ("C2H2","g"):   {"Hf":226.73,"Gf":209.20,"S":200.94},
    ("C6H6","l"):   {"Hf":49.04,"Gf":124.5,"S":173.3},
    ("C6H6","g"):   {"Hf":82.93,"Gf":124.5,"S":269.2},
    ("C2H5OH","l"): {"Hf":-277.7,"Gf":-174.8,"S":160.7},
    ("C2H5OH","g"): {"Hf":-235.1,"Gf":-168.6,"S":282.6},
    ("CH3COOH","l"):{"Hf":-484.5,"Gf":-389.9,"S":159.8},
    ("CH3COOH","g"):{"Hf":-432.2,"Gf":-357.2,"S":256.6},
    ("C3H8","g"):   {"Hf":-103.85,"Gf":-23.49,"S":269.91},
    ("C4H10","g"):  {"Hf":-125.6,"Gf":-15.9,"S":310.2},

    # More ions
    ("PO4^3-","aq"):{"Hf":-1281.0,"Gf":-1019.9,"S":-10.0},
    ("HPO4^2-","aq"):{"Hf":-1300.0,"Gf":-1023.0,"S":20.0},
    ("H2PO4-","aq"):{"Hf":-1315.0,"Gf":-1026.0,"S":40.0},
    ("ClO4-","aq"): {"Hf":-467.9,"Gf":-287.0,"S":110.0},
    ("MnO4-","aq"): {"Hf":-541.0,"Gf":-227.0,"S":151.0},
    ("Cr2O7^2-","aq"):{"Hf":-1782.0,"Gf":-1312.0,"S":295.0},

    # More transition metal ions
    ("Cu+","aq"):   {"Hf":72.4,"Gf":65.0,"S":72.0},
    ("Ag+","aq"):   {"Hf":105.6,"Gf":77.1,"S":73.0},
    ("Zn^2+","aq"): {"Hf":-153.9,"Gf":-147.2,"S":-54.0},
    ("Pb^2+","aq"): {"Hf":-24.2,"Gf":-24.5,"S":-20.0},

    # More hydrides
    ("LiH","s"):    {"Hf":-90.5,"Gf":-68.6,"S":24.8},
    ("NaH","s"):    {"Hf":-56.6,"Gf":-33.3,"S":31.7},
    ("CaH2","s"):   {"Hf":-186.2,"Gf":-154.7,"S":53.0},

    # More halides
    ("LiCl","s"):   {"Hf":-408.3,"Gf":-384.0,"S":35.7},
    ("MgCl2","s"):  {"Hf":-641.8,"Gf":-592.0,"S":89.6},
    ("CaCl2","s"):  {"Hf":-795.8,"Gf":-748.1,"S":104.6},
    ("FeCl3","s"):  {"Hf":-399.5,"Gf":-315.9,"S":136.0},

    # More carbonates
    ("MgCO3","s"):  {"Hf":-1095.8,"Gf":-1029.0,"S":65.7},
    ("BaCO3","s"):  {"Hf":-1227.0,"Gf":-1139.0,"S":112.0},

    # More silicates
    ("Na2SiO3","s"):{"Hf":-1776.0,"Gf":-1630.0,"S":130.0},
    ("K2SiO3","s"): {"Hf":-1840.0,"Gf":-1690.0,"S":150.0},

    # More peroxides
    ("Na2O2","s"):  {"Hf":-510.9,"Gf":-379.9,"S":62.0},
    ("BaO2","s"):   {"Hf":-558.0,"Gf":-493.0,"S":94.0},

    # More simple molecules
    ("H2O2","aq"):  {"Hf":-191.2,"Gf":-120.4,"S":70.0},
    ("O3","g"):     {"Hf":142.67,"Gf":163.20,"S":238.92},
    ("CO","g"):     {"Hf":-110.53,"Gf":-137.17,"S":197.66},
    ("CO2","aq"):   {"Hf":-413.8,"Gf":-386.0,"S":117.0},
})

# --------------------------------------------------------------------------------------
# FRIENDLY LOOKUPS / ALIASES
# --------------------------------------------------------------------------------------

# Common synonyms / alternate notations → canonical formula keys we store under.
ALIASES: Dict[str, str] = {
    # trivial spacing/case handled elsewhere; here are true alternates
    "CLNO":"NOCl", "ONCL":"NOCl", "NOCL":"NOCl",
    # propionic acid (you’ve used HC3H5O2 in problems)
    "HC3H5O2":"C3H6O2", "C2H5COOH":"C3H6O2",
    # acetate/propan-oate names (if someone types a name-ish formula)
    "CH3CO2-":"CH3COO-", "C2H5CO2-":"C2H5COO-",
}

CHARGE_PATTERNS = [
    (r"\(2-\)$", "^2-"), (r"\(3-\)$", "^3-"), (r"\(4-\)$", "^4-"),
    (r"\(2\+\)$", "^2+"), (r"\(3\+\)$", "^3+"),
    (r"--$", "^2-"), (r"\+\+$", "^2+"),
]

def _canon_formula(s: str) -> str:
    """Tidy and alias a formula string (normalize case, charges, obvious variants)."""
    if not s:
        return s
    f = re.sub(r"\s+", "", s)
    f = f.replace("·", ".")  # ignore hydrates dot; we don't store hydrates here
    f_up = f.upper()

    # normalize common charge suffixes like SO4(2-), SO4--, SO4^2-
    for pat, repl in CHARGE_PATTERNS:
        f_up = re.sub(pat, repl, f_up)

    # keep case for multi-letter elements inside (we store keys with cases as usual)
    # but alias whole-string alternates to canonical
    if f_up in ALIASES:
        return ALIASES[f_up]
    # NOCL -> NOCl
    if f_up == "NOCL": return "NOCl"
    # Leave others as-is (e.g., SO4^2-, HCO3-). Our keys use mixed case where needed.
    # Re-insert lower-case where appropriate in simple patterns
    f_up = re.sub(r"CL", "Cl", f_up)
    f_up = re.sub(r"BR", "Br", f_up)
    f_up = re.sub(r"IO", "IO", f_up)  # noop—kept for symmetry
    f_up = re.sub(r"NA", "Na", f_up)
    f_up = re.sub(r"CA", "Ca", f_up)
    f_up = re.sub(r"MG", "Mg", f_up)
    f_up = re.sub(r"AL", "Al", f_up)
    f_up = re.sub(r"FE", "Fe", f_up)
    f_up = re.sub(r"ZN", "Zn", f_up)
    f_up = re.sub(r"CU", "Cu", f_up)
    f_up = re.sub(r"NH", "NH", f_up)  # leave ammonium/ammonia
    return f_up

def phases_for(formula: str) -> List[str]:
    """Return sorted list of phases we have for the given formula (alias-aware)."""
    f = _canon_formula(formula)
    return sorted({ph for (ff, ph) in THERMO.keys() if ff == f})

def get_thermo(formula: str, phase: str) -> Optional[Dict[str, float]]:
    """
    Alias-aware lookup. Exact (formula, phase) first, then try a couple of
    reasonable fallbacks for convenience (e.g., ClNO -> NOCl).
    """
    f = _canon_formula(formula)
    rec = THERMO.get((f, phase))
    if rec is not None:
        return rec
    # gentle fallback: if someone typed e.g. 'NOCL(g)' we'll try again as NOCl
    alt = ALIASES.get(f)
    if alt:
        return THERMO.get((alt, phase))
    return None

def find_species(query: str, phases: Iterable[str] = ("g", "l", "s", "aq")) -> List[Tuple[str, str]]:
    """
    Fuzzy-ish search over keys. Type 'nitrate', 'acetate', 'propionate', 'sulfate', etc.
    Returns list of (formula, phase) sorted by simple score.
    """
    q = query.strip().lower()
    if not q:
        return []
    # very small name hints:
    NAME_HINTS = {
        "nitrate": ["NO3-","NaNO3","KNO3","HNO3"],
        "nitrite": ["NO2-","HNO2"],
        "sulfate": ["SO4^2-","HSO4-","Na2SO4","CaSO4","H2SO4"],
        "sulfite": ["SO3^2-","HSO3-"],
        "carbonate": ["CO3^2-","HCO3-","Na2CO3","NaHCO3","CaCO3"],
        "acetate": ["CH3COO-"],
        "propionate": ["C2H5COO-"],
        "hydroxide": ["OH-","NaOH","KOH","Ca(OH)2"],
        "halide": ["Cl-","Br-","I-","F-","NaCl","KCl","AgCl","AgBr","AgI"],
        "ammonium": ["NH4+","NH3"],
        "phosphate": ["PO4^3-","HPO4^2-","H2PO4-","H3PO4"],
        "nitrosyl chloride": ["NOCl"],
    }
    pool = set()
    for k, lst in NAME_HINTS.items():
        if k in q:
            pool.update(lst)
    # if no hint match, just scan formulas
    if not pool:
        pool.update([f for (f, ph) in THERMO.keys()])

    scored = []
    for f in pool:
        for ph in phases:
            if (f, ph) in THERMO:
                hay = f"{f} {ph}".lower()
                # token score + simple substring
                score = (1 if q in hay else 0) + (len(set(q.split()) & set(hay.split())))
                if score == 0:
                    # very light fuzzy: shared letters
                    score = sum(1 for ch in set(q) if ch in hay) / max(3, len(set(q)))
                scored.append((score, f, ph))
    scored.sort(key=lambda t: t[0], reverse=True)
    return [(f, ph) for _, f, ph in scored]

def load_additional_from_csv(path: str) -> int:
    """
    Merge extra rows from a CSV with columns:
      formula,phase,Hf,Gf,S
    Any missing numeric cell can be blank. Returns how many rows were merged.
    """
    added = 0
    with open(path, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            f = _canon_formula(row.get("formula",""))
            ph = (row.get("phase","") or "").strip()
            if not f or ph not in ("g","l","s","aq"):
                continue
            entry: Dict[str, float] = {}
            for key in ("Hf","Gf","S"):
                val = (row.get(key,"") or "").strip()
                if val:
                    try:
                        entry[key] = float(val)
                    except ValueError:
                        pass
            if entry:
                THERMO[(f, ph)] = entry
                added += 1
    return added
