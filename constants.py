# Default constants (user values override these if provided)
CONSTANTS = {
    # Ideal gas constant (SI): Pa·m³/(mol·K) == J/(mol·K)
    "R": 8.314462618,
    "h": 6.626_070_15e-34,     # Planck constant, J·s (exact)
    "c": 299_792_458.0,        # speed of light in vacuum, m/s (exact)
    "e": 1.602_176_634e-19,    # elementary charge, J per eV (exact)
    "R∞": 1.097_373_156_816e7, # Rydberg constant for infinite mass, 1/m
    "E_H": 13.605_693_122994,  # Hydrogen ground-state |energy|, eV (for Bohr/ΔE)
    "λ_CuKα": 1.54060e-10,  # m  (Cu Kα average)
    "λ_CuKα1": 1.54059e-10, # m
    "λ_CuKα2": 1.54443e-10, # m
    "λ_MoKα": 7.093e-11,    # m  (Mo Kα ≈ 0.7093 Å)
}

CONSTANT_NOTES = {
    "R": "Ideal gas constant in SI units: 8.314462618 J/(mol·K). "
         "If you use other units (e.g., L·atm), provide R yourself (e.g., 0.082057 L·atm/(mol·K))."
}
