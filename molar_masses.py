from typing import Dict
import re
from typing import List, Tuple

# Atomic weights (g/mol). For synthetic/no-stable-isotope elements we use the most common
# mass number (rounded) as a conventional placeholder for calculations.
ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.00794,   "He": 4.002602,
    "Li": 6.941,    "Be": 9.012182, "B": 10.811,   "C": 12.0107,  "N": 14.0067,
    "O": 15.9994,   "F": 18.9984032,"Ne": 20.1797, "Na": 22.98976928, "Mg": 24.3050,
    "Al": 26.9815386,"Si": 28.0855, "P": 30.973762,"S": 32.065,   "Cl": 35.453,
    "Ar": 39.948,   "K": 39.0983,   "Ca": 40.078,  "Sc": 44.955912,"Ti": 47.867,
    "V": 50.9415,   "Cr": 51.9961,  "Mn": 54.938045,"Fe": 55.845, "Co": 58.933195,
    "Ni": 58.6934,  "Cu": 63.546,   "Zn": 65.38,   "Ga": 69.723,  "Ge": 72.64,
    "As": 74.92160, "Se": 78.96,    "Br": 79.904,  "Kr": 83.798,  "Rb": 85.4678,
    "Sr": 87.62,    "Y": 88.90585,  "Zr": 91.224,  "Nb": 92.90638,"Mo": 95.96,
    "Tc": 98.0,     "Ru": 101.07,   "Rh": 102.90550,"Pd": 106.42, "Ag": 107.8682,
    "Cd": 112.411,  "In": 114.818,  "Sn": 118.710, "Sb": 121.760, "Te": 127.60,
    "I": 126.90447, "Xe": 131.293,  "Cs": 132.9054519,"Ba": 137.327,"La": 138.90547,
    "Ce": 140.116,  "Pr": 140.90765,"Nd": 144.242, "Pm": 145.0,   "Sm": 150.36,
    "Eu": 151.964,  "Gd": 157.25,   "Tb": 158.92535,"Dy": 162.500,"Ho": 164.93032,
    "Er": 167.259,  "Tm": 168.93421,"Yb": 173.054, "Lu": 174.9668,"Hf": 178.49,
    "Ta": 180.94788,"W": 183.84,    "Re": 186.207, "Os": 190.23,  "Ir": 192.217,
    "Pt": 195.084,  "Au": 196.966569,"Hg": 200.59, "Tl": 204.3833,"Pb": 207.2,
    "Bi": 208.98040,"Po": 209.0,    "At": 210.0,   "Rn": 222.0,   "Fr": 223.0,
    "Ra": 226.0,    "Ac": 227.0,    "Th": 232.03806,"Pa": 231.03588,"U": 238.02891,
    "Np": 237.0,    "Pu": 244.0,    "Am": 243.0,   "Cm": 247.0,   "Bk": 247.0,
    "Cf": 251.0,    "Es": 252.0,    "Fm": 257.0,   "Md": 258.0,   "No": 259.0,
    "Lr": 262.0,    "Rf": 267.0,    "Db": 270.0,   "Sg": 271.0,   "Bh": 270.0,
    "Hs": 277.0,    "Mt": 276.0,    "Ds": 281.0,   "Rg": 280.0,   "Cn": 285.0,
    "Nh": 284.0,    "Fl": 289.0,    "Mc": 288.0,   "Lv": 293.0,   "Ts": 294.0,
    "Og": 294.0,
}

# Common and alternate names → symbol (lowercased keys)
NAMES_TO_SYMBOL = {
    # Main-group/early
    "hydrogen": "H", "helium": "He", "lithium": "Li", "beryllium": "Be", "boron": "B",
    "carbon": "C", "nitrogen": "N", "oxygen": "O", "fluorine": "F", "neon": "Ne",
    "sodium": "Na", "magnesium": "Mg", "aluminum": "Al", "aluminium": "Al",
    "silicon": "Si", "phosphorus": "P", "phosphorous": "P", "sulfur": "S", "sulphur": "S",
    "chlorine": "Cl", "argon": "Ar", "potassium": "K", "calcium": "Ca",
    "scandium": "Sc", "titanium": "Ti", "vanadium": "V", "chromium": "Cr",
    "manganese": "Mn", "iron": "Fe", "cobalt": "Co", "nickel": "Ni",
    "copper": "Cu", "zinc": "Zn", "gallium": "Ga", "germanium": "Ge",
    "arsenic": "As", "selenium": "Se", "bromine": "Br", "krypton": "Kr",
    "rubidium": "Rb", "strontium": "Sr", "yttrium": "Y", "zirconium": "Zr",
    "niobium": "Nb", "columbium": "Nb", "molybdenum": "Mo", "technetium": "Tc",
    "ruthenium": "Ru", "rhodium": "Rh", "palladium": "Pd", "silver": "Ag",
    "cadmium": "Cd", "indium": "In", "tin": "Sn", "stannum": "Sn",
    "antimony": "Sb", "stibium": "Sb", "tellurium": "Te", "iodine": "I", "xenon": "Xe",
    "cesium": "Cs", "caesium": "Cs", "barium": "Ba", "lanthanum": "La",
    "cerium": "Ce", "praseodymium": "Pr", "neodymium": "Nd", "promethium": "Pm",
    "samarium": "Sm", "europium": "Eu", "gadolinium": "Gd", "terbium": "Tb",
    "dysprosium": "Dy", "holmium": "Ho", "erbium": "Er", "thulium": "Tm",
    "ytterbium": "Yb", "lutetium": "Lu",
    "hafnium": "Hf", "tantalum": "Ta", "tungsten": "W", "wolfram": "W",
    "rhenium": "Re", "osmium": "Os", "iridium": "Ir", "platinum": "Pt",
    "gold": "Au", "mercury": "Hg", "quicksilver": "Hg", "thallium": "Tl",
    "lead": "Pb", "bismuth": "Bi", "polonium": "Po", "astatine": "At",
    "radon": "Rn", "francium": "Fr", "radium": "Ra", "actinium": "Ac",
    "thorium": "Th", "protactinium": "Pa", "uranium": "U", "neptunium": "Np",
    "plutonium": "Pu", "americium": "Am", "curium": "Cm", "berkelium": "Bk",
    "californium": "Cf", "einsteinium": "Es", "fermium": "Fm", "mendelevium": "Md",
    "nobelium": "No", "lawrencium": "Lr", "rutherfordium": "Rf", "dubnium": "Db",
    "seaborgium": "Sg", "bohrium": "Bh", "hassium": "Hs", "meitnerium": "Mt",
    "darmstadtium": "Ds", "roentgenium": "Rg", "copernicium": "Cn",
    "nihonium": "Nh", "flerovium": "Fl", "moscovium": "Mc", "livermorium": "Lv",
    "tennessine": "Ts", "oganesson": "Og",
}




# Common substances and their chemical formulas (lowercased keys)
SPECIES: Dict[str, str] = {
    # Elements (monatomic & diatomic where common)
    "hydrogen": "H", "dihydrogen": "H2",
    "nitrogen": "N", "dinitrogen": "N2",
    "oxygen": "O", "dioxygen": "O2", "ozone": "O3",
    "fluorine": "F2", "chlorine": "Cl2", "bromine": "Br2", "iodine": "I2",

    # Simple inorganics / gases
    "water": "H2O", "heavy water": "D2O",  # parser will choke on D (deuterium) unless added to table; keep H2O preferred
    "ammonia": "NH3", "methane": "CH4", "ethane": "C2H6", "propane": "C3H8", "butane": "C4H10",
    "carbon monoxide": "CO", "carbon dioxide": "CO2", "nitric oxide": "NO",
    "nitrosyl chloride": "NOCl", "nitrosylchlorid": "NOCl",
    "phosphine": "PH3", "hydrazine": "N2H4", "carbon disulfide": "CS2",
    "hydrogen sulfide": "H2S", "sulfur hexafluoride": "SF6",
    "tetrafluoromethane": "CF4", "silicon dioxide": "SiO2",
    "aluminum oxide": "Al2O3", "magnesium oxide": "MgO",
    "nitrogen dioxide": "NO2", "dinitrogen monoxide": "N2O", "nitrous oxide": "N2O",
    "sulfur dioxide": "SO2", "sulfur trioxide": "SO3", "hydrogen sulfide": "H2S",
    "hydrogen chloride": "HCl", "hydrochloric acid": "HCl",
    "hydrogen fluoride": "HF", "hydrogen bromide": "HBr", "hydrogen iodide": "HI",
    "hydrogen peroxide": "H2O2",

    # Acids / bases / salts
    "sulfuric acid": "H2SO4", "sulphuric acid": "H2SO4",
    "nitric acid": "HNO3", "phosphoric acid": "H3PO4",
    "acetic acid": "CH3COOH", "formic acid": "HCOOH", "citric acid": "C6H8O7",
    "sodium hydroxide": "NaOH", "potassium hydroxide": "KOH", "ammonium hydroxide": "NH4OH",
    "sodium chloride": "NaCl", "potassium chloride": "KCl",
    "calcium carbonate": "CaCO3", "sodium bicarbonate": "NaHCO3", "baking soda": "NaHCO3",
    "sodium carbonate": "Na2CO3", "magnesium sulfate": "MgSO4",
    "calcium sulfate": "CaSO4", "aluminum sulfate": "Al2(SO4)3",

    # Organics / common lab solvents etc.
    "ethanol": "C2H5OH", "methanol": "CH3OH", "isopropanol": "C3H8O", "isopropyl alcohol": "C3H8O",
    "acetone": "C3H6O", "acetaldehyde": "C2H4O", "formaldehyde": "CH2O",
    "glucose": "C6H12O6", "sucrose": "C12H22O11", "fructose": "C6H12O6",
    "benzene": "C6H6", "toluene": "C7H8", "acetonitrile": "C2H3N",
    "acrylonitrile": "C3H3N", "acrylic acid": "C3H4O2", "lactic acid": "C3H6O3",
    "urea": "CH4N2O", "aniline": "C6H7N", "phenol": "C6H6O", "chloroform": "CHCl3",
    "dimethyl sulfoxide": "C2H6OS", "dmso": "C2H6OS",

    # Aqueous ions (as formula units; users often want these)
    "ammonium nitrate": "NH4NO3", "ammonium chloride": "NH4Cl",
    "sodium sulfate": "Na2SO4", "potassium nitrate": "KNO3",
    "calcium hydroxide": "Ca(OH)2", "magnesium hydroxide": "Mg(OH)2",

    # Gases misc
    "phosgene": "COCl2", "hydrogen cyanide": "HCN",
    "silane": "SiH4", "diborane": "B2H6",

    # Misc favorites
    "sodium acetate": "CH3COONa", "sodium lactate": "C3H5NaO3",
    "sodium citrate": "Na3C6H5O7",
    "acetylene": "C2H2",
}

SPECIES.update({
    "hydrogen ion": "H+",
    "proton": "H+",
    "hydroxide": "OH-",
    "oxide ion": "O2-",      # not tabulated in water; see note
    "chloride": "Cl-",
    "nitrate": "NO3-",
    "sulfate": "SO4^2-",
    "bisulfate": "HSO4-",
    "bicarbonate": "HCO3-",
    "carbonate": "CO3^2-",
    "ammonium": "NH4+",
    "sodium ion": "Na+",
    "potassium ion": "K+",
    "calcium ion": "Ca^2+",
    "magnesium ion": "Mg^2+",
    "iron(ii)": "Fe^2+",
    "iron(iii)": "Fe^3+",
})

# More salts / alcohols / commons
SPECIES.update({
    # salts
    "sodium sulfate": "Na2SO4",
    "sodium sulphate": "Na2SO4",
    "na2so4": "Na2SO4",   # handy alias so typing 'na2...' also works

    # alcohols
    "methanol": "CH3OH",
    "methyl alcohol": "CH3OH",
    "wood alcohol": "CH3OH",

    "ethanol": "C2H5OH",
    "ethyl alcohol": "C2H5OH",
})

SPECIES.update({
    # ions (phase picked in the UI; names help the search)
    "sodium ion": "Na+",
    "sodium ion gas": "Na+",
    "potassium ion": "K+",
    "chloride ion": "Cl-",
    "chloride ion gas": "Cl-",
    "fluoride ion": "F-",
    "iodide ion": "I-",
    "nitrite": "NO2-",
    "nitrate": "NO3-",
    "chlorate": "ClO3-",
    "perchlorate": "ClO4-",
    "ammonium": "NH4+",
    "copper(ii) ion": "Cu2+",
    "zinc(ii) ion": "Zn2+",
    "aluminum(iii) ion": "Al3+",
    "barium(ii) ion": "Ba2+",
    "lithium ion": "Li+",

    # more common salts & molecules
    "potassium chloride": "KCl",
    "potassium nitrate": "KNO3",
    "sodium nitrate": "NaNO3",
    "calcium chloride": "CaCl2",
    "magnesium chloride": "MgCl2",
    "hydrogen peroxide": "H2O2",
    "nitric acid": "HNO3",
    "sulfuric acid": "H2SO4",
})

SPECIES.update({
    # Elements (standard states)
    "dihydrogen": "H2", "dioxygen": "O2", "dinitrogen": "N2", "difluorine": "F2", "dichlorine": "Cl2",
    "dibromine": "Br2", "diiodine": "I2", "graphite": "C", "sulfur": "S", "iron": "Fe", "sodium": "Na", "potassium": "K",
    # Monatomic gases
    "atomic hydrogen": "H", "atomic oxygen": "O", "atomic sodium": "Na", "atomic potassium": "K", "atomic chlorine": "Cl",
    # Common molecular species
    "water (liquid)": "H2O", "water (gas)": "H2O", "carbon dioxide": "CO2", "carbon monoxide": "CO", "ozone": "O3",
    "nitric oxide": "NO", "nitrogen dioxide": "NO2", "dinitrogen monoxide": "N2O", "ammonia": "NH3", "methane": "CH4",
    "ethane": "C2H6", "propane": "C3H8", "acetylene": "C2H2", "hydrogen sulfide": "H2S", "sulfur dioxide": "SO2",
    "sulfur trioxide": "SO3", "hydrogen chloride": "HCl", "hydrogen fluoride": "HF", "hydrogen bromide": "HBr", "hydrogen iodide": "HI",
    "hydrogen cyanide": "HCN", "nitrosyl chloride": "NOCl",
    # Liquids/solvents
    "methanol (liquid)": "CH3OH", "methanol (gas)": "CH3OH", "ethanol (liquid)": "C2H5OH", "ethanol (gas)": "C2H5OH",
    # Common solids
    "sodium chloride": "NaCl", "potassium chloride": "KCl", "potassium bromide": "KBr", "sodium bromide": "NaBr", "sodium iodide": "NaI",
    "potassium iodide": "KI", "silver chloride": "AgCl", "silver bromide": "AgBr", "silver iodide": "AgI", "calcium carbonate": "CaCO3",
    "calcium oxide": "CaO", "magnesium oxide": "MgO", "aluminum oxide": "Al2O3", "iron(iii) oxide": "Fe2O3", "silicon dioxide": "SiO2",
    "calcium sulfate": "CaSO4", "sodium sulfate": "Na2SO4", "sodium nitrate": "NaNO3", "potassium nitrate": "KNO3",
    "calcium fluoride": "CaF2", "glucose": "C6H12O6",
    # Aqueous ions
    "hydrogen ion": "H+", "hydroxide ion": "OH-", "fluoride ion": "F-", "chloride ion": "Cl-", "bromide ion": "Br-",
    "iodide ion": "I-", "nitrate ion": "NO3-", "nitrite ion": "NO2-", "sulfate ion": "SO4^2-",
    "bisulfate ion": "HSO4-", "sulfite ion": "SO3^2-", "bisulfite ion": "HSO3-", "bicarbonate ion": "HCO3-",
    "carbonate ion": "CO3^2-", "acetate ion": "CH3COO-", "propionate ion": "C2H5COO-", "cyanide ion": "CN-",
    "ammonium ion": "NH4+", "ammonia (aq)": "NH3", "sodium ion": "Na+", "potassium ion": "K+",
    "calcium ion": "Ca^2+", "magnesium ion": "Mg^2+", "zinc ion": "Zn^2+", "copper(ii) ion": "Cu^2+",
    "aluminum ion": "Al^3+", "iron(ii) ion": "Fe^2+", "iron(iii) ion": "Fe^3+",
    # Acids/bases in solution
    "hydrofluoric acid (aq)": "HF", "nitric acid (aq)": "HNO3", "sulfuric acid (liquid)": "H2SO4", "nitrous acid (aq)": "HNO2", "hydrochloric acid (aq)": "HCl",
    # More entries from update
    "sodium nitrate (aq)": "NaNO3", "potassium nitrate (aq)": "KNO3", "sodium hypochlorite": "NaClO", "sodium chlorate": "NaClO3",
    "potassium chlorate": "KClO3", "sodium perchlorate": "NaClO4", "potassium perchlorate": "KClO4", "sodium hydroxide (solid)": "NaOH",
    "potassium hydroxide (solid)": "KOH", "calcium hydroxide (solid)": "Ca(OH)2", "sodium carbonate (solid)": "Na2CO3", "sodium bicarbonate (solid)": "NaHCO3",
    "potassium carbonate (solid)": "K2CO3", "dinitrogen tetroxide": "N2O4", "hydrogen peroxide (liquid)": "H2O2", "hydrogen peroxide (gas)": "H2O2",
    "bromine monochloride": "BrCl", "iodine monochloride (gas)": "ICl", "iodine monochloride (solid)": "ICl", "iodine monobromide (gas)": "IBr", "iodine monobromide (solid)": "IBr",
    "chlorine monofluoride": "ClF", "bromine monofluoride": "BrF", "bromine trifluoride": "BrF3", "bromine pentafluoride": "BrF5", "iodine monofluoride": "IF",
    "iodine pentafluoride": "IF5", "iodine heptafluoride": "IF7", "dinitrogen trioxide": "N2O3", "dinitrogen pentoxide": "N2O5",
    "tetraphosphorus decaoxide": "P4O10", "dichlorine monoxide": "Cl2O", "chlorine dioxide": "ClO2", "dichlorine heptoxide": "Cl2O7",
    "sulfuric acid (aq)": "H2SO4", "phosphoric acid (aq)": "H3PO4", "perchloric acid (aq)": "HClO4", "carbonic acid (aq)": "H2CO3",
    "nitric acid (liquid)": "HNO3", "ammonium hydroxide (aq)": "NH4OH", "barium hydroxide (solid)": "Ba(OH)2", "magnesium hydroxide (solid)": "Mg(OH)2",
    "barium sulfate": "BaSO4", "lead(ii) sulfate": "PbSO4", "copper(ii) sulfate": "CuSO4", "zinc sulfate": "ZnSO4",
    "iron(ii) sulfate": "FeSO4", "aluminum sulfate": "Al2(SO4)3", "ethylene": "C2H4", "acetylene": "C2H2",
    "benzene (liquid)": "C6H6", "benzene (gas)": "C6H6", "acetic acid (liquid)": "CH3COOH", "acetic acid (gas)": "CH3COOH", "butane": "C4H10",
    "phosphate ion": "PO4^3-", "hydrogen phosphate ion": "HPO4^2-", "dihydrogen phosphate ion": "H2PO4-", "perchlorate ion": "ClO4-",
    "permanganate ion": "MnO4-", "dichromate ion": "Cr2O7^2-", "copper(i) ion": "Cu+", "silver ion": "Ag+",
    "lead(ii) ion": "Pb^2+", "lithium hydride": "LiH", "sodium hydride": "NaH",
    "calcium hydride": "CaH2", "lithium chloride": "LiCl", "magnesium chloride": "MgCl2", "calcium chloride": "CaCl2",
    "iron(iii) chloride": "FeCl3", "magnesium carbonate": "MgCO3", "barium carbonate": "BaCO3", "sodium metasilicate": "Na2SiO3",
    "potassium metasilicate": "K2SiO3", "sodium peroxide": "Na2O2", "barium peroxide": "BaO2", "hydrogen peroxide (aq)": "H2O2",
    "carbon dioxide (aq)": "CO2",
})

_token = re.compile(r"([A-Z][a-z]?|\(|\)|\d+)")


def _parse(tokens):
    """Recursive descent parse for parentheses and counts."""
    total = {}
    while tokens:
        t = tokens[0]
        if t == ")":
            break
        if t == "(":
            tokens.pop(0)  # consume '('
            inner = _parse(tokens)
            if not tokens or tokens[0] != ")":
                raise ValueError("Unmatched '(' in formula")
            tokens.pop(0)  # consume ')'
            # optional count
            count = 1
            if tokens and tokens[0].isdigit():
                count = int(tokens.pop(0))
            for sym, n in inner.items():
                total[sym] = total.get(sym, 0) + n * count
        else:
            # element symbol
            if not re.match(r"[A-Z][a-z]?$", t):
                raise ValueError(f"Unexpected token: {t}")
            sym = t
            tokens.pop(0)
            # optional count
            count = 1
            if tokens and tokens[0].isdigit():
                count = int(tokens.pop(0))
            total[sym] = total.get(sym, 0) + count
    return total


def molar_mass(formula_or_name_or_symbol: str) -> float:
    """
    Return molar mass (g/mol) for:
      - element name: "oxygen"
      - element symbol: "O"
      - simple formula: "O2", "CO2", "H2O", "Ca(OH)2", etc.
    """
    s = formula_or_name_or_symbol.strip()
    if not s:
        raise ValueError("Empty formula")
    # Name?
    key = s.lower()
    if key in NAMES_TO_SYMBOL:
        sym = NAMES_TO_SYMBOL[key]
        if sym not in ATOMIC_WEIGHTS:
            raise ValueError(f"No atomic weight for {sym}")
        return ATOMIC_WEIGHTS[sym]
    # Symbol?
    if s in ATOMIC_WEIGHTS:
        return ATOMIC_WEIGHTS[s]
    # Formula parsing
    toks = _token.findall(s)
    composition = _parse(toks)
    if toks:
        raise ValueError("Unexpected trailing tokens in formula")
    mm = 0.0
    for sym, count in composition.items():
        if sym not in ATOMIC_WEIGHTS:
            raise ValueError(f"Unknown element {sym}")
        mm += ATOMIC_WEIGHTS[sym] * count
    return mm


# Add element names (from NAMES_TO_SYMBOL) as suggestions too
for _name, _sym in NAMES_TO_SYMBOL.items():
    SPECIES.setdefault(_name, _sym)

# Helper: resolve a user entry to a chemical formula if it's a known name; else assume it's a formula
def name_to_formula(text: str) -> str:
    s = text.strip()
    if not s:
        return s
    key = s.lower()
    if key in SPECIES:
        return SPECIES[key]
    # treat symbols like "O" or "Na" directly
    if s in ATOMIC_WEIGHTS:
        return s
    # otherwise assume it's already a formula (parser will validate elsewhere)
    return s

def suggest_substances(query: str, limit: int = 50) -> List[Tuple[str, str]]:
    """
    Return a ranked list of (label, formula) suggestions.
    Label examples: "oxygen — O", "dioxygen — O2", "carbon dioxide — CO2"
    Ranking: name startswith > formula startswith > name contains > formula contains.
    """
    q = (query or "").strip().lower()
    items: List[Tuple[str, str]] = []
    for name, formula in SPECIES.items():
        label = f"{name} — {formula}"
        items.append((name, formula, label))

    if not q:
        # some common starters
        common = ["oxygen", "dioxygen", "water", "carbon dioxide", "nitrogen", "ammonia", "methane",
                  "hydrogen", "sulfuric acid", "sodium chloride", "ethanol", "glucose"]
        out = []
        for n in common:
            if n in SPECIES:
                out.append((f"{n} — {SPECIES[n]}", SPECIES[n]))
        # pad with first N alphabetically if needed
        if len(out) < 15:
            rest = sorted(set(f"{n} — {f}" for n, f, _ in items if n not in common))
            out.extend(rest[: (15 - len(out))])
        return out[:limit]

    starts_name = []
    starts_formula = []
    contains_name = []
    contains_formula = []

    for name, formula, label in items:
        n = name.lower()
        f = formula.lower()
        if n.startswith(q):
            starts_name.append((label, formula))
        elif f.startswith(q):
            starts_formula.append((label, formula))
        elif q in n:
            contains_name.append((label, formula))
        elif q in f:
            contains_formula.append((label, formula))

    ranked = starts_name + starts_formula + contains_name + contains_formula
    # Deduplicate by label
    seen = set()
    out: List[Tuple[str, str]] = []
    for label, formula in ranked:
        if label in seen:
            continue
        out.append((label, formula))
        seen.add(label)
        if len(out) >= limit:
            break
    return out
