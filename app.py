# app.py

from aliases import normalize_var
from chem_helper import (
    find_applicable_equations,
    can_solve,
    solve_for,
    constant_notes_for_eq,
)

BANNER = r"""
Chem Formula Finder — MVP
Enter known variables; I'll show equations that fit and what you can solve.
Type 'help' for tips, or press Enter on variable name to finish input.
"""

HELP = """
Tips:
- You can type aliases like: deltaT, dt, temp, vol, mass, pressure, molarity, rho...
- Examples:
    Variable name: m
    Value: 2.5
    Variable name: deltaT      (or dt, delta temperature)
    Value: 30
    Variable name: c
    Value: 4.18
- Constants like R are built-in. If you want a different unit system (e.g., L·atm),
  you can enter your own R and it will override the default.
- Units are YOUR responsibility—keep them consistent within an equation.
"""

def prompt_known_vars():
    print(BANNER)
    known = {}
    while True:
        var = input("Variable name (Enter to finish): ").strip()
        if var == "":
            break
        if var.lower() in ("help", "?"):
            print(HELP)
            continue
        canon = normalize_var(var)
        try:
            val_str = input(f"Value for {canon}: ").strip()
            if val_str == "":
                print("No value entered; skipping.")
                continue
            val = float(val_str.replace(",", "."))  # allow commas
        except ValueError:
            print("Please enter a numeric value.")
            continue
        known[canon] = val
    return known


def main():
    while True:
        known = prompt_known_vars()
        if not known:
            print("No variables entered. Exiting.")
            return

        matches = find_applicable_equations(known)
        if not matches:
            print("\nNo matching equations found yet. Try adding more variables.\n")
        else:
            print("\nWith that you can find the following:\n")
            numbered = []
            for i, (eq, missing) in enumerate(matches, start=1):
                unknown = can_solve(eq, known)
                status = "Solvable" if unknown else f"Missing: {', '.join(missing) if missing else '—'}"
                print(f"{i}) {eq['name']}: {eq['formula']}  [{status}]")
                # Show equation-specific notes and related constant notes
                if 'notes' in eq and eq['notes']:
                    print(f"   Notes: {eq['notes']}")
                const_notes = constant_notes_for_eq(eq)
                for note in const_notes:
                    print(f"   Const: {note}")
                numbered.append((eq, unknown))

            # Offer to compute if any solvable
            solvable_options = [(idx+1, eq, unk) for idx, (eq, unk) in enumerate(numbered) if unk]
            if solvable_options:
                try:
                    choice = input("\nEnter the number to compute the unknown (or press Enter to skip): ").strip()
                    if choice:
                        num = int(choice)
                        # Validate choice
                        if 1 <= num <= len(numbered) and numbered[num-1][1]:
                            eq, target = numbered[num-1]
                            result = solve_for(eq, target, known)
                            print(f"\nResult: {target} = {result}\n")
                        else:
                            print("That selection is not solvable or out of range.")
                except ValueError:
                    print("Not a valid selection; skipping calculation.")

        again = input("Start over? (y/N): ").strip().lower()
        if again != "y":
            break

    print("Goodbye!")


if __name__ == "__main__":
    main()
