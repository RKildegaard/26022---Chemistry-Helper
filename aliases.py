# aliases.py
from typing import Dict, List, Tuple, Optional
import difflib
import math

# Variable registry with human-friendly names, descriptions, aliases, and units.
# Each var:
#  "units": {
#     "base": "<base unit>",
#     "choices": {
#        "<UNIT>": {"factor": <float>, "offset": <float>}  # base = value*factor + offset
#     }
#  }
VARIABLES: Dict[str, Dict] = {
    "Q": {
        "name": "heat (energy)",
        "desc": "Thermal energy transferred.",
        "aliases": ["q", "heat", "energy"],
        "units": {
            "base": "J",
            "choices": {
                "J": {"factor": 1.0, "offset": 0.0},
                "kJ": {"factor": 1e3, "offset": 0.0},
                "cal": {"factor": 4.184, "offset": 0.0},
                "kcal": {"factor": 4184.0, "offset": 0.0},
            },
        },
    },
    "m": {
        "name": "mass",
        "desc": "Amount of matter in a sample.",
        "aliases": ["m", "mass"],
        "units": {
            "base": "kg",
            "choices": {
                "kg": {"factor": 1.0, "offset": 0.0},
                "g": {"factor": 1e-3, "offset": 0.0},
            },
        },
    },
    "c": {
        "name": "specific heat capacity",
        "desc": "Energy needed to raise 1 kg of substance by 1 K.",
        "aliases": ["c", "specific heat", "specificheat", "shc"],
        "units": {
            "base": "J/(kg·K)",
            "choices": {
                "J/(kg·K)": {"factor": 1.0, "offset": 0.0},
                "J/(g·K)": {"factor": 1000.0, "offset": 0.0},
                "cal/(g·°C)": {"factor": 4184.0, "offset": 0.0},
            },
        },
    },
    "ΔT": {
        "name": "temperature change",
        "desc": "Change in temperature (final - initial).",
        "aliases": [
            "deltat",
            "delta_t",
            "dt",
            "dtemp",
            "delta temp",
            "Δt",
            "deltatemp",
            "deltatemperature",
            "deltaT",
            "delta temperature",
        ],
        "units": {
            "base": "K",
            "choices": {
                "K": {"factor": 1.0, "offset": 0.0},
                "°C": {"factor": 1.0, "offset": 0.0},  # ΔK == Δ°C
            },
        },
    },
    "p": {
        "name": "pressure",
        "desc": "Force per unit area exerted by particles.",
        "aliases": ["p", "pressure", "pres", "press"],
        "units": {
            "base": "Pa",
            "choices": {
                "Pa": {"factor": 1.0, "offset": 0.0},
                "kPa": {"factor": 1e3, "offset": 0.0},
                "bar": {"factor": 1e5, "offset": 0.0},
                "atm": {"factor": 101325.0, "offset": 0.0},
                "psi": {"factor": 6894.757293, "offset": 0.0},
            },
        },
    },
    "V": {
        "name": "volume",
        "desc": "Space occupied by the sample.",
        "aliases": ["v", "vol", "volume"],
        "units": {
            "base": "m³",
            "choices": {
                "m³": {"factor": 1.0, "offset": 0.0},
                "L": {"factor": 1e-3, "offset": 0.0},
                "mL": {"factor": 1e-6, "offset": 0.0},
                "cm³": {"factor": 1e-6, "offset": 0.0},
            },
        },
    },
    "n": {
        "name": "amount (moles)",
        "desc": "Amount of substance in moles.",
        "aliases": ["n", "moles", "amount", "mol"],
        "units": {
            "base": "mol",
            "choices": {
                "mol": {"factor": 1.0, "offset": 0.0},
                "mmol": {"factor": 1e-3, "offset": 0.0},
            },
        },
    },
    "R": {
        "name": "gas constant",
        "desc": "Proportionality constant in ideal gas law.",
        "aliases": ["r", "gas constant", "gasconstant", "rgas"],
        "units": {
            "base": "J/(mol·K)",
            "choices": {
                "J/(mol·K)": {"factor": 1.0, "offset": 0.0},
                "L·atm/(mol·K)": {
                    "factor": 101325.0,
                    "offset": 0.0,
                },  # 1 L·atm = 101,325 J
            },
        },
    },
    "T": {
        "name": "temperature",
        "desc": "Absolute temperature.",
        "aliases": ["t", "temp", "temperature"],
        "units": {
            "base": "K",
            "choices": {
                "K": {"factor": 1.0, "offset": 0.0},
                "°C": {"factor": 1.0, "offset": 273.15},  # K = °C + 273.15
            },
        },
    },
    "ρ": {
        "name": "density",
        "desc": "Mass per unit volume.",
        "aliases": ["rho", "density", "ϱ", "dens"],
        "units": {
            "base": "kg/m³",
            "choices": {
                "kg/m³": {"factor": 1.0, "offset": 0.0},
                "g/mL": {"factor": 1000.0, "offset": 0.0},  # 1 g/mL = 1000 kg/m³
                "g/cm³": {"factor": 1000.0, "offset": 0.0},
            },
        },
    },
    "c_m": {
        "name": "molarity (concentration)",
        "desc": "Solute amount per liter of solution.",
        "aliases": ["molarity", "concentration", "cm", "c_m", "conc"],
        "units": {
            "base": "mol/L",
            "choices": {
                "mol/L": {"factor": 1.0, "offset": 0.0},
                "mmol/L": {"factor": 1e-3, "offset": 0.0},
            },
        },
    },
    # ---- Boyle variables (explicit units) ----
    "p1": {
        "name": "initial pressure",
        "desc": "Starting pressure.",
        "aliases": ["p1", "p_1", "p initial", "p init"],
        "units": {
            "base": "Pa",
            "choices": {
                "Pa": {"factor": 1.0, "offset": 0.0},
                "kPa": {"factor": 1e3, "offset": 0.0},
                "bar": {"factor": 1e5, "offset": 0.0},
                "atm": {"factor": 101325.0, "offset": 0.0},
                "psi": {"factor": 6894.757293, "offset": 0.0},
            },
        },
    },
    "v1": {
        "name": "initial volume",
        "desc": "Starting volume.",
        "aliases": ["v1", "v_1", "v initial", "v init"],
        "units": {
            "base": "m³",
            "choices": {
                "m³": {"factor": 1.0, "offset": 0.0},
                "L": {"factor": 1e-3, "offset": 0.0},
                "mL": {"factor": 1e-6, "offset": 0.0},
                "cm³": {"factor": 1e-6, "offset": 0.0},
            },
        },
    },
    "p2": {
        "name": "final pressure",
        "desc": "Ending pressure.",
        "aliases": ["p2", "p_2", "p final"],
        "units": {
            "base": "Pa",
            "choices": {
                "Pa": {"factor": 1.0, "offset": 0.0},
                "kPa": {"factor": 1e3, "offset": 0.0},
                "bar": {"factor": 1e5, "offset": 0.0},
                "atm": {"factor": 101325.0, "offset": 0.0},
                "psi": {"factor": 6894.757293, "offset": 0.0},
            },
        },
    },
    "v2": {
        "name": "final volume",
        "desc": "Ending volume.",
        "aliases": ["v2", "v_2", "v final"],
        "units": {
            "base": "m³",
            "choices": {
                "m³": {"factor": 1.0, "offset": 0.0},
                "L": {"factor": 1e-3, "offset": 0.0},
                "mL": {"factor": 1e-6, "offset": 0.0},
                "cm³": {"factor": 1e-6, "offset": 0.0},
            },
        },
    },
    # ---- Dalton's law variables ----
    "p_total": {
        "name": "total pressure",
        "desc": "Sum of all component partial pressures.",
        "aliases": ["p_total", "ptotal", "p tot", "total pressure", "Ptot"],
        "units": {
            "base": "Pa",
            "choices": {
                "Pa": {"factor": 1.0, "offset": 0.0},
                "kPa": {"factor": 1e3, "offset": 0.0},
                "bar": {"factor": 1e5, "offset": 0.0},
                "atm": {"factor": 101325.0, "offset": 0.0},
                "psi": {"factor": 6894.757293, "offset": 0.0},
            },
        },
    },
    "pA": {
        "name": "partial pressure (A)",
        "desc": "Partial pressure of component A.",
        "aliases": ["pa", "pA", "partial pressure a", "p_A"],
        "units": {
            "base": "Pa",
            "choices": {
                "Pa": {"factor": 1.0, "offset": 0.0},
                "kPa": {"factor": 1e3, "offset": 0.0},
                "bar": {"factor": 1e5, "offset": 0.0},
                "atm": {"factor": 101325.0, "offset": 0.0},
                "psi": {"factor": 6894.757293, "offset": 0.0},
            },
        },
    },
    "pB": {
        "name": "partial pressure (B)",
        "desc": "Partial pressure of component B.",
        "aliases": ["pb", "pB", "partial pressure b", "p_B"],
        "units": {
            "base": "Pa",
            "choices": {
                "Pa": {"factor": 1.0, "offset": 0.0},
                "kPa": {"factor": 1e3, "offset": 0.0},
                "bar": {"factor": 1e5, "offset": 0.0},
                "atm": {"factor": 101325.0, "offset": 0.0},
                "psi": {"factor": 6894.757293, "offset": 0.0},
            },
        },
    },
    "pC": {
        "name": "partial pressure (C)",
        "desc": "Partial pressure of component C.",
        "aliases": ["pc", "pC", "partial pressure c", "p_C"],
        "units": {
            "base": "Pa",
            "choices": {
                "Pa": {"factor": 1.0, "offset": 0.0},
                "kPa": {"factor": 1e3, "offset": 0.0},
                "bar": {"factor": 1e5, "offset": 0.0},
                "atm": {"factor": 101325.0, "offset": 0.0},
                "psi": {"factor": 6894.757293, "offset": 0.0},
            },
        },
    },
    "y_i": {
        "name": "mole fraction (i)",
        "desc": "Component i fraction of total moles; unitless (0–1).",
        "aliases": ["y", "yi", "y_i", "mole fraction", "molefraction"],
        "units": {
            "base": "",
            "choices": {
                "": {"factor": 1.0, "offset": 0.0}  # unitless
            },
        },
    },
    "n_i": {
        "name": "component moles (i)",
        "desc": "Moles of component i.",
        "aliases": ["ni", "n_i", "component moles", "moles i"],
        "units": {
            "base": "mol",
            "choices": {
                "mol": {"factor": 1.0, "offset": 0.0},
                "mmol": {"factor": 1e-3, "offset": 0.0},
            },
        },
    },
    "n_total": {
        "name": "total moles",
        "desc": "Sum of component moles.",
        "aliases": ["ntotal", "n_total", "total moles", "sum moles"],
        "units": {
            "base": "mol",
            "choices": {
                "mol": {"factor": 1.0, "offset": 0.0},
                "mmol": {"factor": 1e-3, "offset": 0.0},
            },
        },
    },
}


# --- Thermodynamics: core state-function variables ---
VARIABLES.update(
    {
        "ΔH°rxn": {
            "name": "standard reaction enthalpy",
            "desc": "ΔH° for the balanced reaction (per mole of reaction).",
            "aliases": [
                "deltaHrxn",
                "ΔHrxn",
                "ΔH°",
                "reaction enthalpy",
                "dHrxn",
                "del H rxn",
            ],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "ΔS°rxn": {
            "name": "standard reaction entropy",
            "desc": "ΔS° for the balanced reaction (per mole of reaction).",
            "aliases": ["deltaSrxn", "ΔSrxn", "reaction entropy", "dSrxn"],
            "units": {
                "base": "J/(mol·K)",
                "choices": {
                    "J/(mol·K)": {"factor": 1.0, "offset": 0.0},
                    "kJ/(mol·K)": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "ΔG°rxn": {
            "name": "standard reaction Gibbs energy",
            "desc": "ΔG° for the balanced reaction (per mole of reaction).",
            "aliases": ["deltaGrxn", "ΔGrxn", "reaction gibbs", "dGrxn", "ΔG°"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
    }
)

# --- Formation data & standard molar entropies (species properties, summed with stoichiometry) ---
VARIABLES.update(
    {
        "ΔH°f": {
            "name": "standard enthalpy of formation",
            "desc": "Per species; use stoichiometric sum for a full reaction.",
            "aliases": ["deltaHf", "ΔHf", "formation enthalpy", "dHf", "ΔH°f"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "ΔG°f": {
            "name": "standard Gibbs energy of formation",
            "desc": "Per species; use stoichiometric sum for a full reaction.",
            "aliases": ["deltaGf", "ΔGf", "formation gibbs", "dGf", "ΔG°f"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "S°": {
            "name": "standard molar entropy",
            "desc": "Per species; use stoichiometric sum for a full reaction.",
            "aliases": ["Sstd", "Sdeg", "standard entropy", "molar entropy"],
            "units": {
                "base": "J/(mol·K)",
                "choices": {
                    "J/(mol·K)": {"factor": 1.0, "offset": 0.0},
                    "kJ/(mol·K)": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        # Pre-summed helpers: user can type the two sums directly (fewer inputs)
        "sum_Hf_prod": {
            "name": "ΣνΔH°f(products)",
            "desc": "Stoichiometric sum of formation enthalpies for all products.",
            "aliases": ["sum hf products", "Σ hf prod", "Hf sum prod"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "sum_Hf_react": {
            "name": "ΣνΔH°f(reactants)",
            "desc": "Stoichiometric sum of formation enthalpies for all reactants.",
            "aliases": ["sum hf reactants", "Σ hf react", "Hf sum react"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "sum_Gf_prod": {
            "name": "ΣνΔG°f(products)",
            "desc": "Stoichiometric sum of formation Gibbs energies for all products.",
            "aliases": ["sum gf products", "Σ gf prod", "Gf sum prod"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "sum_Gf_react": {
            "name": "ΣνΔG°f(reactants)",
            "desc": "Stoichiometric sum of formation Gibbs energies for all reactants.",
            "aliases": ["sum gf reactants", "Σ gf react", "Gf sum react"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "sum_S_prod": {
            "name": "ΣνS°(products)",
            "desc": "Stoichiometric sum of standard molar entropies for products.",
            "aliases": ["sum S products", "Σ S prod", "entropy sum prod"],
            "units": {
                "base": "J/(mol·K)",
                "choices": {
                    "J/(mol·K)": {"factor": 1.0, "offset": 0.0},
                    "kJ/(mol·K)": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "sum_S_react": {
            "name": "ΣνS°(reactants)",
            "desc": "Stoichiometric sum of standard molar entropies for reactants.",
            "aliases": ["sum S reactants", "Σ S react", "entropy sum react"],
            "units": {
                "base": "J/(mol·K)",
                "choices": {
                    "J/(mol·K)": {"factor": 1.0, "offset": 0.0},
                    "kJ/(mol·K)": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
    }
)

# --- Equilibrium & reaction extent / heat bookkeeping ---
VARIABLES.update(
    {
        "K": {
            "name": "equilibrium constant",
            "desc": "Dimensionless equilibrium constant (activities).",
            "aliases": ["equilibrium constant", "K_eq", "Keq"],
            "units": {"base": "", "choices": {"": {"factor": 1.0, "offset": 0.0}}},
        },
        "lnK": {
            "name": "natural log of K",
            "desc": "ln(K). Often used with ΔG° = −RT lnK.",
            "aliases": ["ln K", "log K", "log_e K"],
            "units": {"base": "", "choices": {"": {"factor": 1.0, "offset": 0.0}}},
        },
        "ξ": {
            "name": "extent of reaction",
            "desc": "Moles of reaction advanced (per balanced equation).",
            "aliases": ["xi", "extent", "n_rxn", "n reaction"],
            "units": {
                "base": "mol",
                "choices": {
                    "mol": {"factor": 1.0, "offset": 0.0},
                    "mmol": {"factor": 1e-3, "offset": 0.0},
                },
            },
        },
        "ν_i": {
            "name": "stoichiometric coefficient (i)",
            "desc": "Positive for products, positive input here; use separate sign bookkeeping if needed.",
            "aliases": ["nu_i", "nu", "stoich i", "stoichiometric coefficient"],
            "units": {"base": "", "choices": {"": {"factor": 1.0, "offset": 0.0}}},
        },
        "q_rxn": {
            "name": "heat of reaction (system)",
            "desc": "Heat when the reaction advances by ξ at constant p.",
            "aliases": [
                "q_rxn",
                "q(reaction)",
                "reaction heat",
                "q reaction",
                "q",
            ],
            "units": {
                "base": "J",
                "choices": {
                    "J": {"factor": 1, "offset": 0},
                    "kJ": {"factor": 1e3, "offset": 0},
                },
            },
        },
    }
)

# --- Convenient phase-change enthalpies (optional but handy for H2O(l) → H2O(g) etc.) ---
VARIABLES.update(
    {
        "ΔH_vap": {
            "name": "enthalpy of vaporization",
            "desc": "Per mole; equals ΔH°f(g) − ΔH°f(l) for the same substance.",
            "aliases": ["Hvap", "deltaHvap", "ΔHvap", "enthalpy vaporization"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "ΔH_fus": {
            "name": "enthalpy of fusion",
            "desc": "Per mole; equals ΔH°f(solid) − ΔH°f(liquid) with sign convention.",
            "aliases": ["Hfus", "deltaHfus", "ΔHfus", "enthalpy fusion"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
    }
)

VARIABLES.update(
    {
        "ΔH°f_gas": {
            "name": "standard enthalpy of formation (gas)",
            "desc": "ΔH°f of a substance in the gas phase (per mol). Use with H2O(g).",
            "aliases": ["Hf gas", "deltaHf gas", "ΔHf(g)", "ΔH°f(g)", "Hf_gas"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
        "ΔH°f_liq": {
            "name": "standard enthalpy of formation (liquid)",
            "desc": "ΔH°f of a substance in the liquid phase (per mol). Use with H2O(l).",
            "aliases": ["Hf liq", "deltaHf liq", "ΔHf(l)", "ΔH°f(l)", "Hf_liq"],
            "units": {
                "base": "J/mol",
                "choices": {
                    "J/mol": {"factor": 1.0, "offset": 0.0},
                    "kJ/mol": {"factor": 1e3, "offset": 0.0},
                },
            },
        },
    }
)

# --- Photons / waves ---
VARIABLES.update(
    {
        "ν_freq": {
            "name": "frequency",
            "desc": "Wave frequency.",
            "aliases": ["nu", "v_freq", "freq", "frequency", "f"],
            "units": {
                "base": "Hz",
                "choices": {
                    "Hz": {"factor": 1.0, "offset": 0.0},
                    "kHz": {"factor": 1e3, "offset": 0.0},
                    "MHz": {"factor": 1e6, "offset": 0.0},
                    "GHz": {"factor": 1e9, "offset": 0.0},
                    "THz": {"factor": 1e12, "offset": 0.0},
                },
            },
        },
        "λ": {
            "name": "wavelength",
            "desc": "Wave/photonic wavelength.",
            "aliases": ["lambda", "wavelength", "lam"],
            "units": {
                "base": "m",
                "choices": {
                    "m": {"factor": 1.0, "offset": 0.0},
                    "nm": {"factor": 1e-9, "offset": 0.0},
                    "Å": {"factor": 1e-10, "offset": 0.0},
                    "pm": {"factor": 1e-12, "offset": 0.0},
                },
            },
        },
        "h_planck": {  # was "h"
            "name": "Planck constant",
            "desc": "Planck's constant.",
            "aliases": ["h", "planck", "planck constant"],
            "units": {
                "base": "J·s",
                "choices": {"J·s": {"factor": 1.0, "offset": 0.0}},
            },
        },
        "c0": {  # was "c" (speed of light)
            "name": "speed of light",
            "desc": "Speed of light in vacuum.",
            "aliases": ["c", "lightspeed", "speed of light"],
            "units": {
                "base": "m/s",
                "choices": {"m/s": {"factor": 1.0, "offset": 0.0}},
            },
        },
        "E_ph": {  # more explicit photon energy
            "name": "photon energy",
            "desc": "Energy of a single photon.",
            "aliases": ["E", "E_ph", "photon energy", "energy photon"],
            "units": {
                "base": "J",
                "choices": {
                    "J": {"factor": 1, "offset": 0},
                    "kJ": {"factor": 1e3, "offset": 0},
                    "eV": {"factor": 1.602_176_634e-19, "offset": 0},
                    "keV": {"factor": 1.602_176_634e-16, "offset": 0},
                },
            },
        },
    }
)

# --- Hydrogen-like transitions / Bohr/Rydberg ---
VARIABLES.update(
    {
        "Z": {
            "name": "atomic number (Z)",
            "desc": "Nuclear charge for hydrogen-like ion.",
            "aliases": ["Z", "atomic number", "charge number"],
            "units": {"base": "", "choices": {"": {"factor": 1.0, "offset": 0.0}}},
        },
        "n1": {
            "name": "lower principal quantum number",
            "desc": "Lower energy level index (n₁).",
            "aliases": ["n1", "n lower", "n initial"],
            "units": {"base": "", "choices": {"": {"factor": 1.0, "offset": 0.0}}},
        },
        "n2": {
            "name": "upper principal quantum number",
            "desc": "Upper energy level index (n₂).",
            "aliases": ["n2", "n upper", "n final"],
            "units": {"base": "", "choices": {"": {"factor": 1.0, "offset": 0.0}}},
        },
        "ΔE": {
            "name": "transition energy",
            "desc": "Energy difference between two levels.",
            "aliases": ["dE", "deltaE", "transition energy"],
            "units": {
                "base": "J",
                "choices": {
                    "J": {"factor": 1.0, "offset": 0.0},
                    "eV": {"factor": 1.602_176_634e-19, "offset": 0.0},
                    "keV": {"factor": 1.602_176_634e-16, "offset": 0.0},
                },
            },
        },
        "R∞": {
            "name": "Rydberg constant (∞ mass)",
            "desc": "R∞ ≈ 1.097×10⁷ m⁻¹.",
            "aliases": ["Rinf", "Rydberg", "Rydberg constant"],
            "units": {
                "base": "1/m",
                "choices": {"1/m": {"factor": 1.0, "offset": 0.0}},
            },
        },
    }
)

# --- Bragg / diffraction: unique canonical keys and indices ---
VARIABLES.update(
    {
        "theta": {  # was "θ"
            "name": "Bragg angle θ",
            "desc": "Incidence angle (half of 2θ).",
            "aliases": ["θ", "theta", "bragg angle", "th"],
            "units": {
                "base": "rad",
                "choices": {
                    "rad": {"factor": 1, "offset": 0},
                    "deg": {"factor": math.pi / 180, "offset": 0},
                },
            },
        },
        "two_theta": {  # was "2θ"
            "name": "diffraction angle 2θ",
            "desc": "Twice the Bragg angle; powder XRD reports this.",
            "aliases": ["2θ", "2theta", "two theta"],
            "units": {
                "base": "rad",
                "choices": {
                    "rad": {"factor": 1, "offset": 0},
                    "deg": {"factor": math.pi / 180, "offset": 0},
                },
            },
        },
        "d_spacing": {  # was "d"
            "name": "interplanar spacing d",
            "desc": "Spacing between lattice planes (hkl).",
            "aliases": ["d", "d-spacing", "interplanar spacing"],
            "units": {
                "base": "m",
                "choices": {
                    "m": {"factor": 1, "offset": 0},
                    "nm": {"factor": 1e-9, "offset": 0},
                    "Å": {"factor": 1e-10, "offset": 0},
                    "pm": {"factor": 1e-12, "offset": 0},
                },
            },
        },
        "a_cubic": {  # was "a"
            "name": "cubic lattice parameter a",
            "desc": "Cubic unit-cell edge length.",
            "aliases": ["a", "lattice parameter", "lattice constant"],
            "units": {
                "base": "m",
                "choices": {
                    "m": {"factor": 1, "offset": 0},
                    "nm": {"factor": 1e-9, "offset": 0},
                    "Å": {"factor": 1e-10, "offset": 0},
                    "pm": {"factor": 1e-12, "offset": 0},
                },
            },
        },
        "n_bragg": {  # Bragg order; do NOT alias plain "n"
            "name": "Bragg order",
            "desc": "Diffraction order (1,2, …).",
            "aliases": ["order", "bragg order", "n_bragg"],
            "units": {"base": "", "choices": {"": {"factor": 1, "offset": 0}}},
        },
        "h_mi": {
            "name": "Miller index h",
            "desc": "First Miller index.",
            "aliases": ["h_mi", "miller h", "h(index)"],
            "units": {"base": "", "choices": {"": {"factor": 1, "offset": 0}}},
        },
        "k_mi": {
            "name": "Miller index k",
            "desc": "Second Miller index.",
            "aliases": ["k_mi", "miller k", "k(index)"],
            "units": {"base": "", "choices": {"": {"factor": 1, "offset": 0}}},
        },
        "l_mi": {
            "name": "Miller index l",
            "desc": "Third Miller index.",
            "aliases": ["l_mi", "miller l", "l(index)"],
            "units": {"base": "", "choices": {"": {"factor": 1, "offset": 0}}},
        },
        "q_scat": {  # was "q" (scattering vector)
            "name": "scattering vector magnitude q",
            "desc": "q = 4π sinθ / λ = 2π/d.",
            "aliases": ["q", "q_scat", "scattering vector"],
            "units": {
                "base": "1/m",
                "choices": {
                    "1/m": {"factor": 1, "offset": 0},
                    "Å⁻¹": {"factor": 1e10, "offset": 0},
                    "1/Å": {"factor": 1e10, "offset": 0},
                },
            },
        },
        "s_recip": {  # was "s"
            "name": "reciprocal spacing s = 1/d",
            "desc": "s = 2 sinθ / λ.",
            "aliases": ["s", "s_recip", "1/d", "reciprocal spacing"],
            "units": {
                "base": "1/m",
                "choices": {
                    "1/m": {"factor": 1, "offset": 0},
                    "Å⁻¹": {"factor": 1e10, "offset": 0},
                    "1/Å": {"factor": 1e10, "offset": 0},
                },
            },
        },
    }
)

# ---------- Solutions & colligative properties ----------
VARIABLES.update({
    # liquid mole fraction of component i
    "x_i_liq": {
        "name": "mole fraction (liquid, i)",
        "desc": "x_i in the liquid phase (unitless).",
        "aliases": ["x_i", "x(liq)", "mole fraction liquid", "x"],
        "units": {"base":"", "choices":{"":{"factor":1.0,"offset":0.0}}},
    },
    # gas-phase mole fraction of component i
    "y_i_gas": {
        "name": "mole fraction (gas, i)",
        "desc": "y_i in the gas phase (unitless).",
        "aliases": ["y_i", "y(gas)", "mole fraction gas", "y"],
        "units": {"base":"", "choices":{"":{"factor":1.0,"offset":0.0}}},
    },
    # pure-component saturation vapor pressure of i
    "P_sat_i": {
        "name": "pure-component vapor pressure (i)",
        "desc": "P_i* at the given T for component i.",
        "aliases": ["P*_i", "Psat", "Psat_i", "vapour pressure i"],
        "units": {
            "base":"Pa",
            "choices":{"Pa":{"factor":1.0,"offset":0.0},"kPa":{"factor":1e3,"offset":0.0},"bar":{"factor":1e5,"offset":0.0},"atm":{"factor":101325.0,"offset":0.0}}
        },
    },
    # component partial pressure
    "p_i": {
        "name": "partial pressure (i)",
        "desc": "Partial pressure of component i.",
        "aliases": ["p_i","partial pressure i"],
        "units": {
            "base":"Pa",
            "choices":{"Pa":{"factor":1.0,"offset":0.0},"kPa":{"factor":1e3,"offset":0.0},"bar":{"factor":1e5,"offset":0.0},"atm":{"factor":101325.0,"offset":0.0}}
        },
    },
    # total vapor pressure of mixture
    "P_vap_total": {
        "name": "total vapor pressure (mixture)",
        "desc": "Sum of vapor-phase partial pressures.",
        "aliases": ["P_total_vap","Pt,vap","vapour pressure total"],
        "units": {
            "base":"Pa",
            "choices":{"Pa":{"factor":1.0,"offset":0.0},"kPa":{"factor":1e3,"offset":0.0},"bar":{"factor":1e5,"offset":0.0},"atm":{"factor":101325.0,"offset":0.0}}
        },
    },
    # Henry’s constants (two common conventions)
    "kH_Px": {
        "name": "Henry constant (p = kH·x)",
        "desc": "Pressure-based Henry constant: p_gas = kH·x_solute.",
        "aliases": ["kH", "k_H", "henry constant", "kH(px)"],
        "units": {"base":"Pa", "choices":{"Pa":{"factor":1.0,"offset":0.0},"kPa":{"factor":1e3,"offset":0.0},"bar":{"factor":1e5,"offset":0.0},"atm":{"factor":101325.0,"offset":0.0}}},
    },
    "kH_cP": {
        "name": "Henry constant (c = kH·p)",
        "desc": "Concentration-based Henry constant: c = kH·p.",
        "aliases": ["kH(cp)","henry c=p*k","solubility constant"],
        "units": {"base":"mol/(m³·Pa)", "choices":{"mol/(m³·Pa)":{"factor":1.0,"offset":0.0},"mol/(L·atm)":{"factor":1.0/(1e3*101325.0),"offset":0.0}}},
    },
    # van ’t Hoff factor
    "i_vH": {
        "name": "van ’t Hoff factor",
        "desc": "Effective number of dissolved particles per formula unit.",
        "aliases": ["i", "vant hoff factor", "van't hoff factor"],
        "units": {"base":"", "choices":{"":{"factor":1.0,"offset":0.0}}},
    },
    # molality (to avoid collision with your c_m molarity)
    "m_molal": {
        "name": "molality",
        "desc": "mol solute / kg solvent.",
        "aliases": ["molality","m (molal)","b"],
        "units": {"base":"mol/kg","choices":{"mol/kg":{"factor":1.0,"offset":0.0}}},
    },
    "Kb": {
        "name": "ebullioscopic constant",
        "desc": "Boiling-point elevation constant for solvent.",
        "aliases": ["Kb","K_b","ebullioscopic"],
        "units": {"base":"K·kg/mol","choices":{"K·kg/mol":{"factor":1.0,"offset":0.0}}},
    },
    "Kf": {
        "name": "cryoscopic constant",
        "desc": "Freezing-point depression constant for solvent.",
        "aliases": ["Kf","K_f","cryoscopic"],
        "units": {"base":"K·kg/mol","choices":{"K·kg/mol":{"factor":1.0,"offset":0.0}}},
    },
    "ΔTb": {
        "name": "boiling-point elevation",
        "desc": "ΔT_b = i Kb m.",
        "aliases": ["delta Tb","ΔT_b"],
        "units": {"base":"K","choices":{"K":{"factor":1.0,"offset":0.0},"°C":{"factor":1.0,"offset":0.0}}},
    },
    "ΔTf": {
        "name": "freezing-point depression",
        "desc": "ΔT_f = i Kf m.",
        "aliases": ["delta Tf","ΔT_f"],
        "units": {"base":"K","choices":{"K":{"factor":1.0,"offset":0.0},"°C":{"factor":1.0,"offset":0.0}}},
    },
    "π_osm": {
        "name": "osmotic pressure",
        "desc": "π = i M R T (van ’t Hoff).",
        "aliases": ["pi", "osmotic pressure", "π"],
        "units": {
            "base":"Pa",
            "choices":{"Pa":{"factor":1.0,"offset":0.0},"kPa":{"factor":1e3,"offset":0.0},"bar":{"factor":1e5,"offset":0.0},"atm":{"factor":101325.0,"offset":0.0}}
        },
    },
    # optional: total ion concentration (for quick display)
    "c_total_ions": {
        "name": "total ion concentration",
        "desc": "Sum of molar concentrations of all ions in solution.",
        "aliases": ["ion total conc","total ions"],
        "units": {"base":"mol/L","choices":{"mol/L":{"factor":1.0,"offset":0.0},"mmol/L":{"factor":1e-3,"offset":0.0}}},
    },
})


VARIABLES.update(
    {
        "r_part": {
            "name": "particle radius",
            "desc": "Radius of a (spherical) particle.",
            "aliases": ["radius", "particle radius", "r_part", "r_sphere"],
            "units": {
                "base": "m",
                "choices": {
                    "m": {"factor": 1.0, "offset": 0.0},
                    "nm": {"factor": 1e-9, "offset": 0.0},
                    "Å": {"factor": 1e-10, "offset": 0.0},
                    "pm": {"factor": 1e-12, "offset": 0.0},
                },
            },
        },
        "M_molar": {
            "name": "molar mass",
            "desc": "Mass per mole of species.",
            "aliases": ["M", "molar mass", "MW", "Mr", "molecular weight"],
            "units": {
                "base": "kg/mol",
                "choices": {
                    "kg/mol": {"factor": 1.0, "offset": 0.0},
                    "g/mol": {"factor": 1e-3, "offset": 0.0},
                },
            },
        },
        "N_A": {
            "name": "Avogadro constant",
            "desc": "Avogadro's constant.",
            "aliases": ["NA", "N_A", "avogadro", "avogadro constant"],
            "units": {"base": "1/mol", "choices": {"1/mol": {"factor": 1.0, "offset": 0.0}}},
        },
        "N": {
            "name": "number of particles",
            "desc": "Count of atoms/molecules/particles; unitless.",
            "aliases": ["N", "atoms", "molecules", "count", "# particles"],
            "units": {"base": "", "choices": {"": {"factor": 1.0, "offset": 0.0}}},
        },
    }
)

_ALIAS_TO_CANON: Dict[str, str] = {}
for canon, meta in VARIABLES.items():
    for a in meta["aliases"]:
        _ALIAS_TO_CANON[a.lower()] = canon
    _ALIAS_TO_CANON[meta["name"].lower()] = canon
    _ALIAS_TO_CANON[canon.lower()] = canon


def _startswith_rank(query: str, target: str) -> int:
    return 0 if target.startswith(query) else -1


# ---------- normalization & suggestions ----------
def normalize_var(name: str) -> str:
    s = name.strip()
    if not s:
        return s
    key = s.lower()

    if key in _ALIAS_TO_CANON:
        return _ALIAS_TO_CANON[key]

    starts: List[Tuple[int, str]] = []
    for canon, meta in VARIABLES.items():
        r = _startswith_rank(key, meta["name"].lower())
        if r != -1:
            starts.append((r, canon))
            continue
        for a in meta["aliases"]:
            r = _startswith_rank(key, a.lower())
            if r != -1:
                starts.append((r, canon))
                break
    if starts:
        starts.sort(key=lambda x: (x[0], x[1]))
        return starts[0][1]

    candidates = list(_ALIAS_TO_CANON.keys())
    best = difflib.get_close_matches(key, candidates, n=1, cutoff=0.75)
    if best:
        return _ALIAS_TO_CANON[best[0]]

    return s


def suggestions(query: str, limit: int = 7) -> List[Tuple[str, str]]:
    q = query.strip().lower()
    if not q:
        base = ["p", "V", "n", "T", "R", "m", "c", "ΔT", "ρ", "c_m"]
        out = []
        for c in base:
            if c in VARIABLES:
                out.append((c, f"{c} — {VARIABLES[c]['name']}"))
        return out[:limit]

    starts, contains = [], []
    for canon, meta in VARIABLES.items():
        name = meta["name"].lower()
        aliases = [a.lower() for a in meta["aliases"]]
        label = f"{canon} — {meta['name']}"
        if name.startswith(q) or any(a.startswith(q) for a in aliases):
            starts.append((canon, label))
            continue
        if q in name or any(q in a for a in aliases):
            contains.append((canon, label))

    fuzzy = []
    if len(starts) + len(contains) < limit:
        search_space = []
        canon_for_token = {}
        for canon, meta in VARIABLES.items():
            tokens = [canon, meta["name"], *meta["aliases"]]
            for t in tokens:
                tok = t.lower()
                search_space.append(tok)
                canon_for_token[tok] = canon
        best = difflib.get_close_matches(q, search_space, n=limit, cutoff=0.75)
        seen = set(c for c, _ in starts + contains)
        for b in best:
            c = canon_for_token[b]
            if c in seen:
                continue
            fuzzy.append((c, f"{c} — {VARIABLES[c]['name']}"))
            seen.add(c)

    seen = set()
    result = []
    for group in (starts, contains, fuzzy):
        for c, label in group:
            if c in seen:
                continue
            result.append((c, label))
            seen.add(c)
    return result[:limit]


# ---------- metadata helpers ----------
def var_name(canonical: str) -> str:
    return VARIABLES.get(canonical, {}).get("name", canonical)


def var_desc(canonical: str) -> str:
    return VARIABLES.get(canonical, {}).get("desc", "")


def units_for(canonical: str) -> List[str]:
    meta = VARIABLES.get(canonical, {})
    units = meta.get("units", {})
    choices = units.get("choices", {})
    return list(choices.keys())


def base_unit(canonical: str) -> Optional[str]:
    return VARIABLES.get(canonical, {}).get("units", {}).get("base")


def preferred_unit_for_display(canonical: str) -> Optional[str]:
    """
    Choose a friendly display unit (different from base when useful).
    """
    # Defaults; tweak to your taste
    preferred = {
        "V": "L",  # show liters instead of m³
        "p": "kPa",  # friendlier than Pa
        "m": "kg",
        "Q": "J",
        "n": "mol",
        "T": "K",
        "ΔT": "K",
        "c": "J/(kg·K)",
        "ρ": "kg/m³",
        "c_m": "mol/L",
    }
    preferred.update(
        {
            "λ": "nm",
            "ν": "Hz",
            "E": "eV",
        }
    )
    u = preferred.get(canonical)
    if not u:
        return base_unit(canonical)
    # ensure it's a defined unit
    return u if u in units_for(canonical) else base_unit(canonical)


# ---------- conversions ----------
def convert_to_base(canonical: str, value: float, unit: Optional[str]) -> float:
    meta = VARIABLES.get(canonical, {})
    units = meta.get("units", {})
    choices = units.get("choices", {})
    if not unit or unit not in choices:
        return value
    f = choices[unit].get("factor", 1.0)
    off = choices[unit].get("offset", 0.0)
    return value * f + off


def convert_from_base(canonical: str, base_value: float, unit: Optional[str]) -> float:
    meta = VARIABLES.get(canonical, {})
    units = meta.get("units", {})
    choices = units.get("choices", {})
    if not unit or unit not in choices:
        return base_value
    f = choices[unit].get("factor", 1.0)
    off = choices[unit].get("offset", 0.0)
    # base = value * f + off  => value = (base - off) / f
    return (base_value - off) / f


def pretty_label(canonical: str) -> str:
    return f"{canonical} — {var_name(canonical)}"
