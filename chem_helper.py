# chem_helper.py

from typing import Dict, List, Optional, Tuple, Any

from equations import EQUATIONS
from constants import CONSTANTS, CONSTANT_NOTES


def merged_values(user_values: Dict[str, float]) -> Dict[str, float]:
    """
    Merge constants into user-provided values.
    User values override constants if both are present.
    """
    merged = dict(CONSTANTS)
    merged.update(user_values)
    return merged


def available_keys(user_values: Dict[str, float]) -> set:
    """Set of keys considered 'known' (user + constants)."""
    return set(merged_values(user_values).keys())


def find_applicable_equations(known: dict, min_overlap: int = 1):
    """
    Return a list of (equation_dict, missing_vars) for equations that have at least
    `min_overlap` variables known. Constants (e.g., R) count as known automatically.
    Sorted: solvable first, then by overlap (desc), then by name.
    """
    matches = []
    known_keys = set(known.keys())

    for eq in EQUATIONS:
        vars_set = set(eq["variables"])

        # variables covered by user-known values
        overlap_user = len(vars_set & known_keys)

        # constants that this equation needs and we can auto-fill
        overlap_consts = sum(1 for v in vars_set if v in CONSTANTS and v not in known_keys)

        overlap_total = overlap_user + overlap_consts

        # variables still missing (ignore constants, since we can supply them)
        missing = [v for v in eq["variables"] if v not in known_keys and v not in CONSTANTS]

        if overlap_total >= min_overlap:
            matches.append((eq, missing, overlap_total))

    def is_solvable(eq, missing):
        # solvable if exactly one unknown AND we have a solver for it
        return 1 if (len(missing) == 1 and missing[0] in eq.get("solve", {})) else 0

    # sort: solvable first, then by overlap desc, then alpha name
    matches.sort(key=lambda item: (-is_solvable(item[0], item[1]), -item[2], item[0]["name"].lower()))

    # return original shape (eq, missing)
    return [(eq, missing) for (eq, missing, _overlap) in matches]


def can_solve(eq: Dict[str, Any], known_vars: Dict[str, float]) -> Optional[str]:
    """
    If exactly one variable in eq is unknown (after constants) AND we have a solver for it, return that variable.
    """
    needed = set(eq["variables"])
    unknowns = list(needed - available_keys(known_vars))
    if len(unknowns) == 1 and unknowns[0] in eq["solve"]:
        return unknowns[0]
    return None


def solve_for(eq: Dict[str, Any], target: str, values: Dict[str, float]) -> float:
    """
    Solve equation 'eq' for 'target' using numeric values (with constants merged).
    """
    solver = eq["solve"].get(target)
    if not solver:
        raise ValueError(f"No solver available for {target} in {eq['name']}")
    v = merged_values(values)
    # Restrict values to required keys to avoid surprises
    required = {k: v[k] for k in eq["variables"] if k in v}
    return float(solver(required))


def constant_notes_for_eq(eq: Dict[str, Any]) -> List[str]:
    """Return any relevant constant notes for this equation."""
    notes = []
    for var in eq["variables"]:
        if var in CONSTANT_NOTES:
            notes.append(CONSTANT_NOTES[var])
    return notes
