import json, os
from typing import Dict, Tuple

# Units:
#  - heat capacities: J/(kg·K)
#  - latent heats:    J/kg
#  - transition temps: °C (at ~1 atm; edit if you need different p)
DEFAULT_PHASE_CATALOG: Dict[str, Dict[str, float]] = {
    "water (H2O)": {        # H2O
        "formula": "H2O",
        "c_solid": 2090.0,
        "c_liquid": 4184.0,
        "c_gas": 1996.0,
        "H_fus": 333_550.0,
        "H_vap": 2_256_000.0,
        "T_melt_C": 0.0,
        "T_boil_C": 100.0,
    },
    "ethanol (CH3OH)": {      # C2H5OH (typical values at 1 atm, ~room T)
        "formula": "C2H5OH",
        "c_solid": 1600.0,
        "c_liquid": 2440.0,
        "c_gas": 1430.0,
        "H_fus": 108_000.0,
        "H_vap": 840_000.0,
        "T_melt_C": -114.1,
        "T_boil_C": 78.37,
    },
    "methanol (CH3OH)": {     # CH3OH
        "formula": "CH3OH",
        "c_solid": 1500.0,
        "c_liquid": 2510.0,
        "c_gas": 1500.0,
        "H_fus": 100_000.0,
        "H_vap": 1_100_000.0,
        "T_melt_C": -97.6,
        "T_boil_C": 64.7,
    },
    "acetone (C3H6O)": {      # C3H6O
        "formula": "C3H6O",
        "c_solid": 1300.0,
        "c_liquid": 2180.0,
        "c_gas": 1200.0,
        "H_fus": 98_000.0,
        "H_vap": 500_000.0,
        "T_melt_C": -94.7,
        "T_boil_C": 56.05,
    },
    "benzene (C6H6)": {      # C6H6
        "formula": "C6H6",
        "c_solid": 1200.0,
        "c_liquid": 1740.0,
        "c_gas": 1100.0,
        "H_fus": 126_000.0,
        "H_vap": 394_000.0,
        "T_melt_C": 5.5,
        "T_boil_C": 80.1,
    },
    "ammonia (NH3)": {      # NH3
        "formula": "NH3",
        "c_solid": 2300.0,
        "c_liquid": 4700.0,
        "c_gas": 2080.0,
        "H_fus": 332_000.0,
        "H_vap": 1_370_000.0,
        "T_melt_C": -77.7,
        "T_boil_C": -33.34,
    },
    "Sodium chloride (NaCl)": {         # Sodium chloride
        "formula": "NaCl",
        "c_solid": 864.0,
        "c_liquid": 850.0,
        "c_gas": 820.0,
        "H_fus": 28_160.0,
        "H_vap": 502_000.0,
        "T_melt_C": 801.0,
        "T_boil_C": 1_413.0,
    },
    "Potassium chloride (KCl)": {          # Potassium chloride
        "formula": "KCl",
        "c_solid": 860.0,
        "c_liquid": 900.0,
        "c_gas": 800.0,
        "H_fus": 25_800.0,
        "H_vap": 437_000.0,
        "T_melt_C": 770.0,
        "T_boil_C": 1_420.0,
    },
    "Magnesium chloride (MgCl2)": {        # Magnesium chloride
        "formula": "MgCl2",
        "c_solid": 850.0,
        "c_liquid": 1_100.0,
        "c_gas": 900.0,
        "H_fus": 35_000.0,
        "H_vap": 641_000.0,
        "T_melt_C": 714.0,
        "T_boil_C": 1_412.0,
    },
    "Calcium oxide (CaO)": {          # Calcium oxide
        "formula": "CaO",
        "c_solid": 750.0,
        "c_liquid": 1_100.0,
        "c_gas": 1_000.0,
        "H_fus": 63_700.0,
        "H_vap": 515_000.0,
        "T_melt_C": 2_572.0,
        "T_boil_C": 2_850.0,
    },
    "Aluminum oxide (Al2O3)": {        # Aluminum oxide
        "formula": "Al2O3",
        "c_solid": 880.0,
        "c_liquid": 1_200.0,
        "c_gas": 1_100.0,
        "H_fus": 1_093_000.0,
        "H_vap": 4_800_000.0,
        "T_melt_C": 2_072.0,
        "T_boil_C": 2_977.0,
    },
}

def _catalog_path(path: str | None = None) -> str:
    return path or os.path.join(os.path.dirname(__file__), "phase_constants.json")

def load_phase_catalog(path: str | None = None) -> Tuple[Dict[str, Dict[str, float]], str]:
    p = _catalog_path(path)
    if os.path.exists(p):
        with open(p, "r", encoding="utf-8") as f:
            return json.load(f), p
    # first run: write defaults
    with open(p, "w", encoding="utf-8") as f:
        json.dump(DEFAULT_PHASE_CATALOG, f, indent=2, ensure_ascii=False)
    return json.loads(json.dumps(DEFAULT_PHASE_CATALOG)), p  # deep copy

def save_phase_catalog(cat: Dict[str, Dict[str, float]], path: str | None = None) -> str:
    p = _catalog_path(path)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(cat, f, indent=2, ensure_ascii=False)
    return p
