# equations.py

from typing import Dict, Any, List
import math


# Each equation has:
# - key: short identifier
# - name: human-friendly name
# - variables: list of canonical variable keys used by this equation
# - formula: a pretty string
# - solve: dict mapping target variable -> function(values_dict) -> float
#   (values_dict uses canonical variable keys; constants will be injected at solve time)
def _asin_clamped(x: float) -> float:
    return math.asin(max(-1.0, min(1.0, x)))


EQUATIONS: List[Dict[str, Any]] = [
    {
        "key": "specific_heat",
        "name": "Specific Heat Capacity",
        "variables": ["Q", "m", "c", "ΔT"],
        "formula": "Q = m · c · ΔT",
        "solve": {
            "Q": lambda v: v["m"] * v["c"] * v["ΔT"],
            "m": lambda v: v["Q"] / (v["c"] * v["ΔT"]),
            "c": lambda v: v["Q"] / (v["m"] * v["ΔT"]),
            "ΔT": lambda v: v["Q"] / (v["m"] * v["c"]),
        },
        "notes": "Typical units: Q in J, m in kg (or g), c in J/(kg·K) (or J/(g·K)), ΔT in K or °C.",
    },
    {
        "key": "ideal_gas_law",
        "name": "Ideal Gas Law",
        "variables": ["p", "V", "n", "R", "T"],
        "formula": "p · V = n · R · T",
        "solve": {
            "p": lambda v: v["n"] * v["R"] * v["T"] / v["V"],
            "V": lambda v: v["n"] * v["R"] * v["T"] / v["p"],
            "n": lambda v: v["p"] * v["V"] / (v["R"] * v["T"]),
            "R": lambda v: v["p"] * v["V"] / (v["n"] * v["T"]),
            "T": lambda v: v["p"] * v["V"] / (v["n"] * v["R"]),
        },
        "notes": "Use consistent units. Default R is SI (J/(mol·K)); if you prefer L·atm, supply your own R.",
    },
    {
        "key": "density",
        "name": "Density",
        "variables": ["ρ", "m", "V"],
        "formula": "ρ = m / V",
        "solve": {
            "ρ": lambda v: v["m"] / v["V"],
            "m": lambda v: v["ρ"] * v["V"],
            "V": lambda v: v["m"] / v["ρ"],
        },
        "notes": "ρ in kg/m³ (or g/mL), m in kg (or g), V in m³ (or mL).",
    },
    # Sphere volume from radius
    {
        "key": "sphere_volume",
        "name": "Sphere volume",
        "variables": ["V", "r_part"],
        "formula": "V = (4/3) π r_part^3",
        "solve": {
            "V": lambda v: (4.0 / 3.0) * math.pi * v["r_part"] ** 3,
            "r_part": lambda v: ((3.0 * v["V"]) / (4.0 * math.pi)) ** (1.0 / 3.0),
        },
        "notes": "r_part in m, V in m³.",
    },
    # Mass–moles–molar mass
    {
        "key": "mass_moles_molar_mass",
        "name": "Mass–moles–molar mass",
        "variables": ["m", "n", "M_molar"],
        "formula": "m = n · M_molar",
        "solve": {
            "m": lambda v: v["n"] * v["M_molar"],
            "n": lambda v: v["m"] / v["M_molar"],
            "M_molar": lambda v: v["m"] / v["n"],
        },
        "notes": "m in kg, M_molar in kg/mol, n in mol.",
    },
    # Particle count from moles (Avogadro)
    {
        "key": "particle_count",
        "name": "Particle count",
        "variables": ["N", "n", "N_A"],
        "formula": "N = n · N_A",
        "solve": {
            "N": lambda v: v["n"] * v["N_A"],
            "n": lambda v: v["N"] / v["N_A"],
            "N_A": lambda v: v["N"] / v["n"],
        },
        "notes": "N is unitless, n in mol, N_A in 1/mol.",
    },
    # One-step: atoms in a spherical particle from ρ, r, M, N_A
    {
        "key": "particle_count_from_radius",
        "name": "Atoms in spherical particle",
        "variables": ["N", "ρ", "r_part", "M_molar", "N_A"],
        "formula": "N = (ρ · (4/3)π r_part^3 / M_molar) · N_A",
        "solve": {
            "N": lambda v: (
                v["ρ"] * (4.0 / 3.0) * math.pi * v["r_part"] ** 3 / v["M_molar"]
            )
            * v["N_A"],
            "r_part": lambda v: (
                ((v["N"] * v["M_molar"]) / (v["ρ"] * v["N_A"]))
                * (3.0 / (4.0 * math.pi))
            )
            ** (1.0 / 3.0),
            "ρ": lambda v: (v["N"] * v["M_molar"])
            / (((4.0 / 3.0) * math.pi * v["r_part"] ** 3) * v["N_A"]),
            "M_molar": lambda v: (
                v["ρ"] * (4.0 / 3.0) * math.pi * v["r_part"] ** 3 * v["N_A"]
            )
            / v["N"],
            "N_A": lambda v: (v["N"] * v["M_molar"])
            / (v["ρ"] * (4.0 / 3.0) * math.pi * v["r_part"] ** 3),
        },
        "notes": "Use consistent units: ρ in kg/m³ (or g/cm³), r_part in m (nm ok), M_molar in kg/mol (or g/mol).",
    },
    {
        "key": "molarity",
        "name": "Molarity",
        "variables": ["c_m", "n", "V"],
        "formula": "c = n / V",
        "solve": {
            "c_m": lambda v: v["n"] / v["V"],
            "n": lambda v: v["c_m"] * v["V"],
            "V": lambda v: v["n"] / v["c_m"],
        },
        "notes": "c (molarity) in mol/L, n in mol, V in L.",
    },
    {
        "key": "boyle",
        "name": "Boyle's Law (isothermal)",
        "variables": ["p1", "v1", "p2", "v2"],
        "formula": "p₁ · V₁ = p₂ · V₂",
        "solve": {
            "p1": lambda v: v["p2"] * v["v2"] / v["v1"],
            "v1": lambda v: v["p2"] * v["v2"] / v["p1"],
            "p2": lambda v: v["p1"] * v["v1"] / v["v2"],
            "v2": lambda v: v["p1"] * v["v1"] / v["p2"],
        },
        "notes": "Keep pressure and volume units consistent.",
    },
    {
        "key": "dalton_partial",
        "name": "Dalton's Law (partial pressure)",
        "variables": ["p_i", "y_i", "p_total"],
        "formula": "pᵢ = yᵢ · p_total",
        "solve": {
            "p_i": lambda v: v["y_i"] * v["p_total"],
            "y_i": lambda v: v["p_i"] / v["p_total"],
            "p_total": lambda v: v["p_i"] / v["y_i"],
        },
        "notes": "Mole fraction yᵢ = nᵢ / n_total (unitless). Keep pressure units consistent.",
    },
    {
        "key": "mole_fraction",
        "name": "Mole Fraction",
        "variables": ["y_i", "n_i", "n_total"],
        "formula": "yᵢ = nᵢ / n_total",
        "solve": {
            "y_i": lambda v: v["n_i"] / v["n_total"],
            "n_i": lambda v: v["y_i"] * v["n_total"],
            "n_total": lambda v: v["n_i"] / v["y_i"],
        },
        "notes": "All amounts in mol; yᵢ is unitless and between 0 and 1.",
    },
    {
        "key": "dalton_sum2",
        "name": "Dalton's Law (sum of partials, 2 components)",
        "variables": ["p_total", "pA", "pB"],
        "formula": "p_total = pA + pB",
        "solve": {
            "p_total": lambda v: v["pA"] + v["pB"],
            "pA": lambda v: v["p_total"] - (v["pB"]),
            "pB": lambda v: v["p_total"] - (v["pA"]),
        },
        "notes": "Total pressure equals sum of component partial pressures. Use any consistent pressure unit.",
    },
    {
        "key": "dalton_sum3",
        "name": "Dalton's Law (sum of partials, 3 components)",
        "variables": ["p_total", "pA", "pB", "pC"],
        "formula": "p_total = pA + pB + pC",
        "solve": {
            "p_total": lambda v: v["pA"] + v["pB"] + v["pC"],
            "pA": lambda v: v["p_total"] - (v["pB"] + v["pC"]),
            "pB": lambda v: v["p_total"] - (v["pA"] + v["pC"]),
            "pC": lambda v: v["p_total"] - (v["pA"] + v["pB"]),
        },
        "notes": "Total pressure equals sum of component partial pressures. Use any consistent pressure unit.",
    },
    # Gibbs relation (covers spontaneity tests & "at what T?" by solving for T when ΔG°=0)
    {
        "key": "gibbs_standard",
        "name": "Gibbs relation (standard)",
        "variables": ["ΔG°rxn", "ΔH°rxn", "T", "ΔS°rxn"],
        "formula": "ΔG° = ΔH° − T·ΔS°",
        "solve": {
            "ΔG°rxn": lambda v: v["ΔH°rxn"] - v["T"] * v["ΔS°rxn"],
            "ΔH°rxn": lambda v: v["ΔG°rxn"] + v["T"] * v["ΔS°rxn"],
            "T": lambda v: (v["ΔH°rxn"] - v["ΔG°rxn"]) / v["ΔS°rxn"],
            "ΔS°rxn": lambda v: (v["ΔH°rxn"] - v["ΔG°rxn"]) / v["T"],
        },
        "notes": "T in K; ΔH° and ΔG° in J/mol; ΔS° in J/(mol·K).",
    },
    # ΔG° ↔ K (equilibrium link)
    {
        "key": "equilibrium_dg",
        "name": "Equilibrium ↔ Gibbs (standard)",
        "variables": ["ΔG°rxn", "K", "R", "T"],
        "formula": "ΔG° = − R·T·ln K",
        "solve": {
            "ΔG°rxn": lambda v: -v["R"] * v["T"] * math.log(v["K"]),
            "K": lambda v: math.exp(-v["ΔG°rxn"] / (v["R"] * v["T"])),
        },
        "notes": "Use activities; K is unitless. R default from constants.py.",
    },
    # Reaction enthalpy from formation enthalpies
    {
        "key": "rxn_enthalpy_from_formation",
        "name": "Reaction enthalpy from formation data",
        "variables": ["ΔH°rxn", "sum_Hf_prod", "sum_Hf_react"],
        "formula": "ΔH°rxn = ΣνΔH°f(products) − ΣνΔH°f(reactants)",
        "solve": {
            "ΔH°rxn": lambda v: v["sum_Hf_prod"] - v["sum_Hf_react"],
            "sum_Hf_prod": lambda v: v["ΔH°rxn"] + v["sum_Hf_react"],
            "sum_Hf_react": lambda v: v["sum_Hf_prod"] - v["ΔH°rxn"],
        },
        "notes": "Enter stoichiometric sums (ν·ΔH°f) for each side, then solve ΔH°rxn.",
    },
    # Reaction entropy from standard molar entropies
    {
        "key": "rxn_entropy_from_S",
        "name": "Reaction entropy from standard molar entropies",
        "variables": ["ΔS°rxn", "sum_S_prod", "sum_S_react"],
        "formula": "ΔS°rxn = ΣνS°(products) − ΣνS°(reactants)",
        "solve": {
            "ΔS°rxn": lambda v: v["sum_S_prod"] - v["sum_S_react"],
            "sum_S_prod": lambda v: v["ΔS°rxn"] + v["sum_S_react"],
            "sum_S_react": lambda v: v["sum_S_prod"] - v["ΔS°rxn"],
        },
        "notes": "Use S° values at 298 K unless specified otherwise.",
    },
    # Reaction Gibbs from formation Gibbs
    {
        "key": "rxn_gibbs_from_formation",
        "name": "Reaction Gibbs from formation data",
        "variables": ["ΔG°rxn", "sum_Gf_prod", "sum_Gf_react"],
        "formula": "ΔG°rxn = ΣνΔG°f(products) − ΣνΔG°f(reactants)",
        "solve": {
            "ΔG°rxn": lambda v: v["sum_Gf_prod"] - v["sum_Gf_react"],
            "sum_Gf_prod": lambda v: v["ΔG°rxn"] + v["sum_Gf_react"],
            "sum_Gf_react": lambda v: v["sum_Gf_prod"] - v["ΔG°rxn"],
        },
        "notes": "Combine with ΔG° ↔ K to get equilibrium constants.",
    },
    # Heat released/absorbed from reaction extent
    {
        "key": "heat_from_extent",
        "name": "Heat from reaction extent",
        "variables": ["q", "ξ", "ΔH°rxn"],
        "formula": "q = ξ · ΔH°rxn",
        "solve": {
            "q": lambda v: v["ξ"] * v["ΔH°rxn"],
            "ξ": lambda v: v["q"] / v["ΔH°rxn"],
            "ΔH°rxn": lambda v: v["q"] / v["ξ"],
        },
        "notes": "Exothermic reactions have negative ΔH°rxn (system convention).",
    },
    # Extent from a component amount and its stoichiometric coefficient
    {
        "key": "extent_from_component",
        "name": "Extent from component moles",
        "variables": ["ξ", "n_i", "ν_i"],
        "formula": "ξ = n_i / ν_i",
        "solve": {
            "ξ": lambda v: v["n_i"] / v["ν_i"],
            "n_i": lambda v: v["ξ"] * v["ν_i"],
            "ν_i": lambda v: v["n_i"] / v["ξ"],
        },
        "notes": "Use the coefficient of the species in the balanced equation.",
    },
    # (Optional) Vaporization enthalpy from formation enthalpies of same substance
    {
        "key": "vap_from_formation",
        "name": "Enthalpy of vaporization from formation data",
        "variables": ["ΔH_vap", "ΔH°f_gas", "ΔH°f_liq"],
        "formula": "ΔH_vap = ΔH°f(g) − ΔH°f(l)",
        "solve": {
            "ΔH_vap": lambda v: v["ΔH°f_gas"] - v["ΔH°f_liq"],
            "ΔH°f_gas": lambda v: v["ΔH_vap"] + v["ΔH°f_liq"],
            "ΔH°f_liq": lambda v: v["ΔH°f_gas"] - v["ΔH_vap"],
        },
        "notes": "Handy shortcut for e.g. H₂O(l) → H₂O(g).",
    },
    # (Optional) Scaling/reversal helpers for enthalpy bookkeeping
    {
        "key": "enthalpy_scaling",
        "name": "Scale reaction enthalpy by factor",
        "variables": ["ΔH°rxn_scaled", "ΔH°rxn", "scale_n"],
        "formula": "ΔH°(scaled) = n · ΔH°(base)",
        "solve": {
            "ΔH°rxn_scaled": lambda v: v["scale_n"] * v["ΔH°rxn"],
            "ΔH°rxn": lambda v: v["ΔH°rxn_scaled"] / v["scale_n"],
            "scale_n": lambda v: v["ΔH°rxn_scaled"] / v["ΔH°rxn"],
        },
        "notes": "If you double the reaction, you double ΔH°. Flip sign when reversing.",
    },
    {
        "key": "enthalpy_reverse",
        "name": "Reverse reaction enthalpy",
        "variables": ["ΔH°rxn_rev", "ΔH°rxn"],
        "formula": "ΔH°(reverse) = − ΔH°(forward)",
        "solve": {
            "ΔH°rxn_rev": lambda v: -v["ΔH°rxn"],
            "ΔH°rxn": lambda v: -v["ΔH°rxn_rev"],
        },
        "notes": "Use with scaling to handle arbitrary reaction manipulations.",
    },
    {
        "key": "combustion_energy",
        "name": "Energy released by combustion",
        "variables": ["q", "n", "ΔG°rxn"],
        "formula": "q = n · (−ΔG°rxn)",
        "solve": {
            "q": lambda v: v["n"] * (-v["ΔG°rxn"]),
            "n": lambda v: v["q"] / (-v["ΔG°rxn"]),
            "ΔG°rxn": lambda v: -v["q"] / v["n"],
        },
        "notes": "q is energy released (J); n is mol of substance combusted; ΔG°rxn is standard Gibbs energy of reaction (J/mol, negative for spontaneous combustion).",
    },
    # Planck / wave
    {
        "key": "planck_hnu",
        "name": "Photon energy from frequency",
        "variables": ["E_ph", "h_planck", "ν_freq"],
        "formula": "E = h · ν",
        "solve": {
            "E_ph": lambda v: v["h_planck"] * v["ν_freq"],
            "ν_freq": lambda v: v["E_ph"] / v["h_planck"],
            "h_planck": lambda v: v["E_ph"] / v["ν_freq"],
        },
        "notes": "E in J (or eV), ν in Hz.",
    },
    {
        "key": "wave_speed",
        "name": "Wave speed relation",
        "variables": ["c0", "λ", "ν"],
        "formula": "c = λ · ν",
        "solve": {
            "c0": lambda v: v["λ"] * v["ν"],
            "λ": lambda v: v["c0"] / v["ν"],
            "ν": lambda v: v["c0"] / v["λ"],
        },
    },
    {
        "key": "planck_hc_over_lambda",
        "name": "Photon energy from wavelength",
        "variables": ["E_ph", "h_planck", "c0", "λ"],
        "formula": "E = h · c / λ",
        "solve": {
            "E_ph": lambda v: v["h_planck"] * v["c0"] / v["λ"],
            "λ": lambda v: v["h_planck"] * v["c0"] / v["E_ph"],
        },
    },
    # --- Hydrogen-like lines (Rydberg/Bohr) ---
    {
        "key": "rydberg_lambda",
        "name": "Hydrogen-like wavelength (Rydberg)",
        "variables": ["λ", "R∞", "Z", "n1", "n2"],
        "formula": "1/λ = R∞ · Z² · (1/n₁² − 1/n₂²), n₂>n₁",
        "solve": {
            "λ": lambda v: 1.0
            / (v["R∞"] * (v["Z"] ** 2) * (1.0 / (v["n1"] ** 2) - 1.0 / (v["n2"] ** 2))),
        },
        "notes": "λ in m. For H, Z=1. Use with E = h·c/λ to get energy.",
    },
    {
        "key": "bohr_transition",
        "name": "Hydrogen-like transition energy (Bohr)",
        "variables": ["ΔE", "Z", "n1", "n2", "E_H"],
        "formula": "ΔE = E_H · Z² · (1/n₁² − 1/n₂²)  (in eV)",
        "solve": {
            # return J so units system is consistent; E_H is in eV in constants
            "ΔE": lambda v: (
                v["E_H"] * (v["Z"] ** 2) * (1.0 / (v["n1"] ** 2) - 1.0 / (v["n2"] ** 2))
            )
            * 1.602_176_634e-19,
        },
        "notes": "Result returned in J; display in eV by choosing the eV unit.",
    },
    # Bragg family
    {
        "key": "bragg",
        "name": "Bragg’s law",
        "variables": ["n_bragg", "λ", "d_spacing", "theta"],
        "formula": "n·λ = 2·d·sinθ",
        "solve": {
            "n_bragg": lambda v: 2.0 * v["d_spacing"] * math.sin(v["theta"]) / v["λ"],
            "λ": lambda v: 2.0 * v["d_spacing"] * math.sin(v["theta"]) / v["n_bragg"],
            "d_spacing": lambda v: (v["n_bragg"] * v["λ"])
            / (2.0 * math.sin(v["theta"])),
            "theta": lambda v: _asin_clamped(
                v["n_bragg"] * v["λ"] / (2.0 * v["d_spacing"])
            ),
        },
    },
    {
        "key": "bragg_2theta",
        "name": "Bragg with 2θ",
        "variables": ["n_bragg", "λ", "d_spacing", "two_theta"],
        "formula": "n·λ = 2·d·sin(2θ/2)",
        "solve": {
            "d_spacing": lambda v: (v["n_bragg"] * v["λ"])
            / (2.0 * math.sin(0.5 * v["two_theta"])),
            "λ": lambda v: 2.0
            * v["d_spacing"]
            * math.sin(0.5 * v["two_theta"])
            / v["n_bragg"],
            "n_bragg": lambda v: 2.0
            * v["d_spacing"]
            * math.sin(0.5 * v["two_theta"])
            / v["λ"],
            "two_theta": lambda v: 2.0
            * _asin_clamped(v["n_bragg"] * v["λ"] / (2.0 * v["d_spacing"])),
        },
    },
    {
        "key": "cubic_d_hkl",
        "name": "Cubic: d from a and (hkl)",
        "variables": ["d_spacing", "a_cubic", "h_mi", "k_mi", "l_mi"],
        "formula": "d = a / √(h²+k²+l²)",
        "solve": {
            "d_spacing": lambda v: v["a_cubic"]
            / math.sqrt(v["h_mi"] ** 2 + v["k_mi"] ** 2 + v["l_mi"] ** 2),
            "a_cubic": lambda v: v["d_spacing"]
            * math.sqrt(v["h_mi"] ** 2 + v["k_mi"] ** 2 + v["l_mi"] ** 2),
        },
    },
    {
        "key": "s_from_bragg",
        "name": "Reciprocal spacing s from Bragg",
        "variables": ["s_recip", "λ", "theta"],
        "formula": "s = 2·sinθ / λ",
        "solve": {
            "s_recip": lambda v: 2.0 * math.sin(v["theta"]) / v["λ"],
            "theta": lambda v: _asin_clamped(0.5 * v["s_recip"] * v["λ"]),
            "λ": lambda v: 2.0 * math.sin(v["theta"]) / v["s_recip"],
        },
    },
    {
        "key": "q_from_bragg",
        "name": "Scattering vector q",
        "variables": ["q_scat", "λ", "theta"],
        "formula": "q = 4π·sinθ / λ = 2π / d",
        "solve": {
            "q_scat": lambda v: 4.0 * math.pi * math.sin(v["theta"]) / v["λ"],
            "theta": lambda v: _asin_clamped(0.25 * v["q_scat"] * v["λ"] / math.pi),
            "λ": lambda v: 4.0 * math.pi * math.sin(v["theta"]) / v["q_scat"],
        },
    },
    {
        "key": "q_s_d_links",
        "name": "Links: q, s, d",
        "variables": ["q_scat", "s_recip", "d_spacing"],
        "formula": "s = 1/d ; q = 2π·s = 2π/d",
        "solve": {
            "s_recip": lambda v: 1.0 / v["d_spacing"]
            if "d_spacing" in v
            else v["q_scat"] / (2.0 * math.pi),
            "q_scat": lambda v: 2.0 * math.pi / v["d_spacing"]
            if "d_spacing" in v
            else 2.0 * math.pi * v["s_recip"],
            "d_spacing": lambda v: 1.0 / v["s_recip"]
            if "s_recip" in v
            else (2.0 * math.pi) / v["q_scat"],
        },
    },
    # ---------- Raoult’s law (two-component mix, ideal) ----------
    {
        "key": "raoult_2comp",
        "name": "Raoult’s law (2 components)",
        "variables": ["p_i", "x_i_liq", "P_sat_i", "P_vap_total"],
        "formula": "p_i = x_i · P_i* ;  P_total = Σ p_i",
        "solve": {
            "p_i": lambda v: v["x_i_liq"] * v["P_sat_i"],
            "x_i_liq": lambda v: v["p_i"] / v["P_sat_i"],
            "P_sat_i": lambda v: v["p_i"] / v["x_i_liq"],
            # If both components’ p_i are provided in values dict, compute sum:
            "P_vap_total": lambda v: sum(val for k, val in v.items() if k == "p_i")
            if isinstance(v.get("p_i"), (list, tuple))
            else v["p_i"],
        },
        "notes": "Use per-component p_i, x_i, P_i*. For totals, sum p_i of all components.",
    },
    # ---------- Henry’s law (two conventions) ----------
    {
        "key": "henry_px",
        "name": "Henry: p = kH·x",
        "variables": ["p", "kH_Px", "x_i_liq"],
        "formula": "p_gas = kH · x_liq",
        "solve": {
            "p": lambda v: v["kH_Px"] * v["x_i_liq"],
            "kH_Px": lambda v: v["p"] / v["x_i_liq"],
            "x_i_liq": lambda v: v["p"] / v["kH_Px"],
        },
        "notes": "kH in Pa (or bar/kPa).",
    },
    {
        "key": "henry_cp",
        "name": "Henry: c = kH·p",
        "variables": ["c_m", "kH_cP", "p"],
        "formula": "c = kH · p",
        "solve": {
            "c_m": lambda v: v["kH_cP"]
            * v["p"],  # treat c_m as mol/L; kH_cP supports mol/(L·atm)
            "kH_cP": lambda v: v["c_m"] / v["p"],
            "p": lambda v: v["c_m"] / v["kH_cP"],
        },
        "notes": "Mind units for kH; here c uses mol/L if kH is in mol/(L·atm).",
    },
    # ---------- Osmotic pressure ----------
    {
        "key": "osmotic_pressure",
        "name": "Osmotic pressure (van ’t Hoff)",
        "variables": ["π_osm", "i_vH", "c_m", "R", "T"],
        "formula": "π = i·M·R·T",
        "solve": {
            "π_osm": lambda v: v["i_vH"] * v["c_m"] * v["R"] * v["T"],
            "i_vH": lambda v: v["π_osm"] / (v["c_m"] * v["R"] * v["T"]),
            "c_m": lambda v: v["π_osm"] / (v["i_vH"] * v["R"] * v["T"]),
        },
        "notes": "M in mol/L is fine if R uses L·atm; otherwise convert. Default R in J/(mol·K) → use Pa for π.",
    },
    # ---------- Colligative properties ----------
    {
        "key": "boiling_elevation",
        "name": "Boiling-point elevation",
        "variables": ["ΔTb", "i_vH", "Kb", "m_molal"],
        "formula": "ΔT_b = i · K_b · m",
        "solve": {
            "ΔTb": lambda v: v["i_vH"] * v["Kb"] * v["m_molal"],
            "i_vH": lambda v: v["ΔTb"] / (v["Kb"] * v["m_molal"]),
            "Kb": lambda v: v["ΔTb"] / (v["i_vH"] * v["m_molal"]),
            "m_molal": lambda v: v["ΔTb"] / (v["i_vH"] * v["Kb"]),
        },
        "notes": "Use solvent-specific Kb (K·kg/mol).",
    },
    {
        "key": "freezing_depression",
        "name": "Freezing-point depression",
        "variables": ["ΔTf", "i_vH", "Kf", "m_molal"],
        "formula": "ΔT_f = i · K_f · m",
        "solve": {
            "ΔTf": lambda v: v["i_vH"] * v["Kf"] * v["m_molal"],
            "i_vH": lambda v: v["ΔTf"] / (v["Kf"] * v["m_molal"]),
            "Kf": lambda v: v["ΔTf"] / (v["i_vH"] * v["m_molal"]),
            "m_molal": lambda v: v["ΔTf"] / (v["i_vH"] * v["Kf"]),
        },
        "notes": "Use solvent-specific Kf (K·kg/mol).",
    },
    # ---------- Total ion concentration (quick relation) ----------
    {
        "key": "total_ion_conc",
        "name": "Total ion concentration (via i)",
        "variables": ["c_total_ions", "i_vH", "c_m"],
        "formula": "c_total_ions = i · M",
        "solve": {
            "c_total_ions": lambda v: v["i_vH"] * v["c_m"],
            "i_vH": lambda v: v["c_total_ions"] / v["c_m"],
            "c_m": lambda v: v["c_total_ions"] / v["i_vH"],
        },
        "notes": "For ideal dissociation i = number of ions per formula unit.",
    },
]
