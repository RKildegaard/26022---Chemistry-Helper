# gui.py
import tkinter as tk
from tkinter import ttk, messagebox
from molar_masses import molar_mass, suggest_substances, name_to_formula
from thermo_data import get_thermo, phases_for
from constants import CONSTANTS
import math
from phase_catalog import load_phase_catalog, save_phase_catalog, DEFAULT_PHASE_CATALOG


from aliases import (
    normalize_var, suggestions, var_name, var_desc, units_for,
    convert_to_base, convert_from_base, pretty_label, base_unit, preferred_unit_for_display
)
from chem_helper import (
    find_applicable_equations,
    can_solve,
    solve_for,
    constant_notes_for_eq,
)
from molar_masses import molar_mass

APP_TITLE = "Chem Formula Finder — GUI"

class ChemGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("1200x680")
        self.minsize(980, 620)

        # State (values in base units)
        self.known_base = {}
        # UI memory (unit + display_value)
        self.known_ui = {}

        self._build_layout()
        self._refresh_equation_list()
        self.phase_catalog, self.phase_catalog_path = load_phase_catalog()


    # ---------------- UI BUILD ----------------
    def _build_layout(self):
        # Menu bar
        menubar = tk.Menu(self)
        convert_menu = tk.Menu(menubar, tearoff=0)
        convert_menu.add_command(label="Mass ↔ Moles (element/formula)…", command=self._open_mol_converter)
        convert_menu.add_command(label="Pressure units…", command=self._open_pressure_converter)
        convert_menu.add_command(label="Thermo lookup (ΔH°f, ΔG°f, S°)…", command=self._open_thermo_lookup)
        convert_menu.add_command(label="Reaction builder (ΔH°, ΔS°, ΔG°)…", command=self._open_reaction_builder)
        convert_menu.add_command(label="Atoms in nanoparticle…", command=self._open_nanoparticle_atoms_tool)


        explain_menu = tk.Menu(menubar, tearoff=0)
        explain_menu.add_command(label="ΔG° ↔ K (show steps)…", command=self._open_equilibrium_dg_steps)
        explain_menu.add_command(label="Unpaired d-electrons (octahedral)…",command=self._open_unpaired_d_electrons_tool)


        convert_menu.add_command(label="Photon color from E/λ…", command=self._open_photon_color_tool)
        convert_menu.add_command(label="Heating/Cooling (phase changes)…",command=self._open_phase_change_tool)
        
        menu_solutions = tk.Menu(menubar, tearoff=0)
        menu_solutions.add_command(label="Percent & Mixing…", command=self._open_percent_mixer_tool)
        menu_solutions.add_command(label="Henry’s Law…", command=self._open_henry_tool)
        menu_solutions.add_command(label="Raoult’s Law (2 comp)…", command=self._open_raoult_tool)
        menu_solutions.add_command(label="Ion Concentrations / van ’t Hoff…", command=self._open_ion_vant_hoff_tool)
        menu_solutions.add_command(label="Colligative ΔT (Kb/Kf)…", command=self._open_colligative_tool)
        menu_solutions.add_command(label="Osmotic Pressure…", command=self._open_osmotic_tool)

        # after you set up menubar / other cascades:
        lookup_menu = tk.Menu(menubar, tearoff=0)
        menubar.add_cascade(label="Lookup", menu=lookup_menu)
        lookup_menu.add_command(label="Kb / Kf (common solvents)…", command=self._open_kb_kf_lookup)

        # --- Kinetics menu ----------------------------------------------------------
        kin_menu = tk.Menu(menubar, tearoff=0)
        kin_menu.add_command(label="Rate laws (0/1/2 order)…", command=self._open_rate_law_tool)
        kin_menu.add_command(label="Arrhenius (k, Ea, A, T)…", command=self._open_arrhenius_tool)
        kin_menu.add_command(label="Initial Rates (find m, n)…", command=self._open_initial_rates_tool)

        menubar.add_cascade(label="Kinetics", menu=kin_menu)

        # --- Equilibrium menu ---
        menu_eq = tk.Menu(menubar, tearoff=False)
        menubar.add_cascade(label="Equilibrium", menu=menu_eq)
        menu_eq.add_command(label="ICE Solver (Kc/Kp, Q, shift)…", command=self._open_equilibrium_ice_tool)
        menu_eq.add_command(label="Kp ↔ Kc Converter…",        command=self._open_kp_kc_converter)
        menu_eq.add_command(label="Combine Equilibria (overall K)…", command=self._open_equilibrium_combine_tool)
        menu_eq.add_command(label="Q & Le Châtelier Helper…",  command=self._open_q_direction_tool)
        # under your Equilibrium menu
        menu_eq.add_command(label="K-expression (builder)…", command=self._open_k_expression_builder)
        menu_eq.add_command(label="Ksp & Heterogeneous Helper…", command=self._open_ksp_tool)
        menu_eq.add_command(label="van ’t Hoff (ΔH° from K vs T)…",command=self._open_vant_hoff_tool)



        acid_base_menu = tk.Menu(menubar, tearoff=False)
        menubar.add_cascade(label="Acid & Base", menu=acid_base_menu)
        acid_base_menu.add_command(label="Ph (strong)", command=self._open_ph_modal)
        acid_base_menu.add_command(label="Ph (weak)", command=self._open_weak_acid_base_modal)

        # --- Oxidation menu ---
        m_ox = tk.Menu(menubar, tearoff=0)
        m_ox.add_command(label="Oxidation number (single species)…", command=self._open_oxidation_number_tool)
        m_ox.add_command(label="Redox analysis (reaction)…", command=self._open_redox_analysis_tool)  # optional
        menubar.add_cascade(label="Oxidation", menu=m_ox)

        m_rxn = tk.Menu(menubar, tearoff=0)
        m_rxn.add_command(label="Net-ionic equation…", command=self._open_net_ionic_tool)
        m_rxn.add_command(label="Organic Structures", command=self._open_functional_group_finder)
        m_rxn.add_command(label="Draw Organic Structures", command=self._open_structure_drawer)
        menubar.add_cascade(label="Reactions", menu=m_rxn)





        menubar.add_cascade(label="Convert", menu=convert_menu)
        menubar.add_cascade(label="Solutions", menu=menu_solutions)
        menubar.add_cascade(label="Explain", menu=explain_menu)
        self.config(menu=menubar)

                # --- Search / Command Palette -----------------------------------------
        search_menu = tk.Menu(menubar, tearoff=0)
        search_menu.add_command(label="Search tools… (Ctrl+K / ⌘K)",
                                command=self._open_command_palette)
        menubar.add_cascade(label="Search", menu=search_menu)

        # Keyboard shortcuts for quick open
        self.bind_all("<Control-k>", lambda e: self._open_command_palette())
        self.bind_all("<Command-k>", lambda e: self._open_command_palette())  # macOS


        title = ttk.Label(self, text="Chem Formula Finder", font=("Segoe UI", 18, "bold"))
        title.pack(pady=(10, 0))

        subtitle = ttk.Label(
            self,
            text="Add known variables with descriptions and units. Type aliases (e.g., 'pres' → pressure). "
                 "Equations use base SI internally. Use Convert → Mass ↔ Moles when needed.",
            wraplength=980, justify="center"
        )
        subtitle.pack(pady=(2, 10))

        paned = ttk.Panedwindow(self, orient=tk.HORIZONTAL)
        paned.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

        left = ttk.Frame(paned, padding=8)
        paned.add(left, weight=2)
        self._build_known_vars_panel(left)

        right = ttk.Frame(paned, padding=8)
        paned.add(right, weight=3)
        self._build_equations_panel(right)

        help_text = (
            "Tips:\n"
            "• We convert inputs to base SI (e.g., kPa→Pa, L→m³, °C→K) automatically.\n"
            "• Results show base units and a friendly unit when helpful (e.g., L for volume).\n"
            "• For molecules like O2 or H2O, use the converter and type the chemical formula."
        )
        ttk.Label(self, text=help_text, foreground="#444", wraplength=980, justify="left").pack(
            padx=12, pady=(0, 10), anchor="w"
        )


    def _add_compute_bar(self, win, compute):
        bar = ttk.Frame(win)
        bar.pack(side=tk.BOTTOM, fill=tk.X, padx=12, pady=8)  # PACK on win
        ttk.Button(bar, text="Compute", command=compute).pack(side=tk.LEFT)
        win.bind("<Return>", lambda e: compute())



    def _build_known_vars_panel(self, parent):
        header = ttk.Label(parent, text="Known Variables", font=("Segoe UI", 12, "bold"))
        header.pack(anchor="w", pady=(0, 6))

        form = ttk.Frame(parent)
        form.pack(fill=tk.X, pady=(0, 6))

        ttk.Label(form, text="Variable:").grid(row=0, column=0, padx=(0, 6), pady=2, sticky="w")
        self.var_name = ttk.Entry(form, width=22)
        self.var_name.grid(row=0, column=1, padx=(0, 8), pady=2, sticky="w")
        self.var_name.bind("<KeyRelease>", self._on_var_type)
        self.var_name.bind("<Down>", self._focus_suggest_list)
        self.var_name.bind("<Return>", self._accept_from_entry)

        ttk.Label(form, text="Unit:").grid(row=0, column=2, padx=(6, 6), pady=2, sticky="w")
        self.var_unit = ttk.Combobox(form, width=16, state="readonly", values=[])
        self.var_unit.grid(row=0, column=3, padx=(0, 8), pady=2, sticky="w")

        ttk.Label(form, text="Value:").grid(row=0, column=4, padx=(6, 6), pady=2, sticky="w")
        self.var_value = ttk.Entry(form, width=16)
        self.var_value.grid(row=0, column=5, padx=(0, 10), pady=2, sticky="w")
        self.var_value.bind("<Return>", lambda e: self._on_add_known())

        add_btn = ttk.Button(form, text="Add / Update", command=self._on_add_known)
        add_btn.grid(row=0, column=6, padx=(0, 10), pady=2, sticky="w")

        clear_btn = ttk.Button(form, text="Clear All", command=self._on_clear_all)
        clear_btn.grid(row=0, column=7, padx=(0, 0), pady=2, sticky="w")

        # Suggestion list
        self.suggest_list = tk.Listbox(parent, height=0, width=60)
        self.suggest_list.pack(fill=tk.X, padx=2, pady=(0, 6))
        self.suggest_scrollx = ttk.Scrollbar(parent, orient="horizontal", command=self.suggest_list.xview)
        self.suggest_list.configure(xscrollcommand=self.suggest_scrollx.set)
        self.suggest_scrollx.pack(fill=tk.X, padx=2, pady=(0, 6))
        self.suggest_list.bind("<Double-Button-1>", self._accept_suggestion)
        self.suggest_list.bind("<Return>", self._accept_suggestion)
        self.suggest_list.bind("<Escape>", lambda e: self._hide_suggestions())
        self.suggest_list.bind("<Up>", self._listbox_up)
        self.suggest_list.bind("<Down>", self._listbox_down)
        self._hide_suggestions()

        # Table
        self.known_tree = ttk.Treeview(parent, columns=("var", "name", "val", "unit"), show="headings", height=11)
        self.known_tree.heading("var", text="Symbol")
        self.known_tree.heading("name", text="Name")
        self.known_tree.heading("val", text="Value")
        self.known_tree.heading("unit", text="Unit")
        self.known_tree.column("var", width=90, anchor="w")
        self.known_tree.column("name", width=210, anchor="w")
        self.known_tree.column("val", width=120, anchor="e")
        self.known_tree.column("unit", width=100, anchor="w")
        self.known_tree.pack(fill=tk.BOTH, expand=True, pady=(6, 4))

        del_btn = ttk.Button(parent, text="Delete Selected", command=self._on_delete_selected)
        del_btn.pack(anchor="e", pady=(0, 6))

    def _build_equations_panel(self, parent):
        header = ttk.Label(parent, text="Applicable Equations", font=("Segoe UI", 12, "bold"))
        header.pack(anchor="w", pady=(0, 6))

        cols = ("name", "formula", "status")
        self.eq_tree = ttk.Treeview(parent, columns=cols, show="headings", height=15)
        for c, w, anchor in zip(cols, (260, 420, 180), ("w", "w", "center")):
            self.eq_tree.heading(c, text=c.capitalize())
            self.eq_tree.column(c, width=w, anchor=anchor)
        self.eq_tree.pack(fill=tk.BOTH, expand=True)

        self.eq_info = tk.Text(parent, height=7, wrap="word")
        self.eq_info.pack(fill=tk.BOTH, expand=False, pady=(6, 6))
        self.eq_info.configure(state="disabled")

        btns = ttk.Frame(parent)
        btns.pack(fill=tk.X)

        show_btn = ttk.Button(btns, text="Show Details", command=self._on_show_details)
        show_btn.pack(side=tk.LEFT, padx=(0, 8))

        compute_btn = ttk.Button(btns, text="Compute Unknown", command=self._on_compute)
        compute_btn.pack(side=tk.LEFT, padx=(0, 8))

        refresh_btn = ttk.Button(btns, text="Refresh", command=self._refresh_equation_list)
        refresh_btn.pack(side=tk.LEFT, padx=(0, 8))

        self.result_lbl = ttk.Label(parent, text="", font=("Segoe UI", 11, "bold"), foreground="#0a5")
        self.result_lbl.pack(anchor="w", pady=(8, 0))

        self.eq_tree.bind("<<TreeviewSelect>>", lambda e: self._on_tree_select())

    # ---------------- Suggestion logic ----------------
    def _on_var_type(self, event=None):
        text = self.var_name.get()
        items = suggestions(text, limit=7)
        self.suggest_list.delete(0, tk.END)
        for canon, label in items:
            self.suggest_list.insert(tk.END, label)
        if items:
            self.suggest_list.configure(height=min(7, len(items)))
            self.suggest_list.pack_configure()
            self.suggest_list.selection_clear(0, tk.END)
            self.suggest_list.selection_set(0)
            self.suggest_list.activate(0)
        else:
            self._hide_suggestions()
        self._update_units_for_current_entry()

    def _update_units_for_current_entry(self):
        text = self.var_name.get().strip()
        if not text:
            self.var_unit["values"] = []
            self.var_unit.set("")
            return
        canon = normalize_var(text)
        choices = units_for(canon)
        if choices:
            self.var_unit["values"] = choices
            # Prefer base if present, else first
            base = base_unit(canon)
            default = base if base in choices else choices[0]
            if not self.var_unit.get() or self.var_unit.get() not in choices:
                self.var_unit.set(default)
        else:
            self.var_unit["values"] = []
            self.var_unit.set("")

    def _hide_suggestions(self):
        self.suggest_list.forget()

    def _focus_suggest_list(self, event=None):
        if str(self.suggest_list.winfo_manager()) != "pack":
            return
        self.suggest_list.focus_set()
        if self.suggest_list.size() > 0:
            self.suggest_list.selection_clear(0, tk.END)
            self.suggest_list.selection_set(0)
            self.suggest_list.activate(0)

    def _accept_suggestion(self, event=None):
        sel = self.suggest_list.curselection()
        if not sel:
            return
        label = self.suggest_list.get(sel[0])
        canon = label.split(" — ", 1)[0].strip()
        self.var_name.delete(0, tk.END)
        self.var_name.insert(0, canon)
        self._update_units_for_current_entry()
        self.var_value.focus_set()
        self._hide_suggestions()

    def _listbox_up(self, event=None):
        idxs = self.suggest_list.curselection()
        if not idxs:
            return "break"
        i = idxs[0]
        if i > 0:
            self.suggest_list.selection_clear(0, tk.END)
            self.suggest_list.selection_set(i - 1)
            self.suggest_list.activate(i - 1)
        return "break"

    def _listbox_down(self, event=None):
        idxs = self.suggest_list.curselection()
        if not idxs:
            return "break"
        i = idxs[0]
        if i < self.suggest_list.size() - 1:
            self.suggest_list.selection_clear(0, tk.END)
            self.suggest_list.selection_set(i + 1)
            self.suggest_list.activate(i + 1)
        return "break"

    def _accept_from_entry(self, event=None):
        text = self.var_name.get().strip()
        if not text:
            return
        canon = normalize_var(text)
        self.var_name.delete(0, tk.END)
        self.var_name.insert(0, canon)
        self._update_units_for_current_entry()
        self.var_value.focus_set()
        self._hide_suggestions()

    # ---------------- Known variable handlers ----------------
    def _on_add_known(self):
        name = self.var_name.get().strip()
        if not name:
            messagebox.showinfo("Input needed", "Please enter a variable name.")
            return
        canon = normalize_var(name)

        val_str = self.var_value.get().strip().replace(",", ".")
        try:
            display_val = float(val_str)
        except ValueError:
            messagebox.showerror("Invalid value", "Please enter a numeric value.")
            return

        unit = self.var_unit.get().strip() or None
        base_val = convert_to_base(canon, display_val, unit)

        self.known_base[canon] = base_val
        self.known_ui[canon] = {"unit": unit or "", "display_value": display_val}

        self._refresh_known_table()
        self._refresh_equation_list()
        self.result_lbl.config(text="")
        self.var_name.delete(0, tk.END)
        self.var_value.delete(0, tk.END)
        self.var_unit.set("")
        self.var_name.focus_set()

    def _on_clear_all(self):
        if messagebox.askyesno("Clear all", "Remove all known variables?"):
            self.known_base.clear()
            self.known_ui.clear()
            self._refresh_known_table()
            self._refresh_equation_list()
            self.result_lbl.config(text="")
            self._hide_suggestions()

    def _on_delete_selected(self):
        sel = self.known_tree.selection()
        if not sel:
            return
        for item in sel:
            var = self.known_tree.item(item, "values")[0]
            self.known_base.pop(var, None)
            self.known_ui.pop(var, None)
        self._refresh_known_table()
        self._refresh_equation_list()
        self.result_lbl.config(text="")

    # ---------------- Equation panel handlers ----------------
    def _on_tree_select(self):
        self._on_show_details()

    def _on_show_details(self):
        eq = self._get_selected_equation()
        self._set_info("")
        if not eq:
            return
        unknown = can_solve(eq, self.known_base)
        missing = [v for v in eq["variables"] if v not in set(self.known_base.keys())]
        txt = []
        txt.append(f"Name: {eq['name']}")
        txt.append(f"Formula: {eq['formula']}")
        txt.append(f"Variables: {', '.join(eq['variables'])}")
        if unknown:
            txt.append(f"Status: Solvable — can compute {unknown}")
        else:
            txt.append(f"Status: Missing -> {', '.join(missing) if missing else '—'}")
        if eq.get("notes"):
            txt.append(f"Notes: {eq['notes']}")
        txt.append("")
        txt.append("Variable help:")
        for v in eq["variables"]:
            desc = var_desc(v)
            if desc:
                txt.append(f"  • {v} — {var_name(v)}: {desc}")
        const_notes = constant_notes_for_eq(eq)
        for note in const_notes:
            txt.append(f"Constant: {note}")
        self._set_info("\n".join(txt))

    def _on_compute(self):
        eq = self._get_selected_equation()
        self.result_lbl.config(text="")
        if not eq:
            messagebox.showinfo("No selection", "Please select an equation first.")
            return
        target = can_solve(eq, self.known_base)
        if not target:
            messagebox.showwarning("Not solvable", "This equation is not solvable with current inputs.")
            return
        try:
            base_result = solve_for(eq, target, self.known_base)
        except KeyError as e:
            messagebox.showerror("Missing value", f"Missing required variable: {e}")
            return
        except ZeroDivisionError:
            messagebox.showerror("Math error", "Division by zero encountered. Check your inputs.")
            return
        except Exception as ex:
            messagebox.showerror("Error", f"Could not compute: {ex}")
            return

        # Units on result (base + friendly)
        bu = base_unit(target) or ""
        friendly_u = preferred_unit_for_display(target)
        text = f"Result: {target} = {base_result:g} {bu}".strip()
        if friendly_u and friendly_u != bu:
            friendly_val = convert_from_base(target, base_result, friendly_u)
            text += f"  (≈ {friendly_val:g} {friendly_u})"
        self.result_lbl.config(text=text)

    # ---------------- Helpers ----------------
    def _refresh_known_table(self):
        for row in self.known_tree.get_children():
            self.known_tree.delete(row)
        for k in sorted(self.known_base.keys()):
            ui = self.known_ui.get(k, {"unit": "", "display_value": self.known_base[k]})
            self.known_tree.insert("", tk.END, values=(k, var_name(k), f"{ui['display_value']:g}", ui["unit"]))

    def _refresh_equation_list(self):
        for row in self.eq_tree.get_children():
            self.eq_tree.delete(row)

        matches = find_applicable_equations(self.known_base)
        self._row_to_eq = {}
        for eq, missing in matches:
            unknown = can_solve(eq, self.known_base)
            status = "Solvable" if unknown else (f"Missing: {', '.join(missing)}" if missing else "—")
            iid = self.eq_tree.insert("", tk.END, values=(eq["name"], eq["formula"], status))
            self._row_to_eq[iid] = eq

        self._set_info("")
        self.result_lbl.config(text="")

    def _get_selected_equation(self):
        sel = self.eq_tree.selection()
        if not sel:
            return None
        iid = sel[0]
        return self._row_to_eq.get(iid)

    def _set_info(self, text):
        self.eq_info.configure(state="normal")
        self.eq_info.delete("1.0", tk.END)
        self.eq_info.insert("1.0", text)
        self.eq_info.configure(state="disabled")

    # ---------------- Converter modal ----------------
    def _open_mol_converter(self):
        win = tk.Toplevel(self)
        win.title("Mass ↔ Moles Converter")
        win.transient(self)
        win.grab_set()
        win.geometry("620x420")

        frm = ttk.Frame(win, padding=12)
        frm.pack(fill=tk.BOTH, expand=True)

        # Substance entry + suggestions
        ttk.Label(frm, text="Substance (name/symbol/formula):").grid(row=0, column=0, sticky="w")
        sub_entry = ttk.Entry(frm, width=36)
        sub_entry.grid(row=0, column=1, columnspan=3, sticky="we", padx=8, pady=4)
        sub_entry.insert(0, "oxy")  # helpful default for demo

        # Suggestion list
        sugg = tk.Listbox(frm, height=0)
        sugg.grid(row=1, column=0, columnspan=4, sticky="we", padx=(0,0))
        sugg_scrollx = ttk.Scrollbar(frm, orient="horizontal", command=sugg.xview)
        sugg.configure(xscrollcommand=sugg_scrollx.set)
        sugg_scrollx.grid(row=2, column=0, columnspan=4, sticky="we", pady=(0,6))
        sugg.grid_remove()
        sugg_scrollx.grid_remove()

        # Formula preview
        formula_var = tk.StringVar(value="")
        mm_var = tk.StringVar(value="")
        preview = ttk.Label(frm, textvariable=formula_var, font=("Segoe UI", 10, "bold"))
        preview.grid(row=3, column=0, columnspan=4, sticky="w", pady=(2,8))
        mm_lbl = ttk.Label(frm, textvariable=mm_var, foreground="#444")
        mm_lbl.grid(row=4, column=0, columnspan=4, sticky="w", pady=(0,8))

        # Mode & inputs
        mode_var = tk.StringVar(value="m2n")  # m2n or n2m
        mode_frame = ttk.Frame(frm)
        mode_frame.grid(row=5, column=0, columnspan=4, sticky="w", pady=(2, 8))
        ttk.Radiobutton(mode_frame, text="Mass → Moles", variable=mode_var, value="m2n").pack(side=tk.LEFT, padx=(0, 12))
        ttk.Radiobutton(mode_frame, text="Moles → Mass", variable=mode_var, value="n2m").pack(side=tk.LEFT)

        ttk.Label(frm, text="Amount:").grid(row=6, column=0, sticky="w")
        amt_entry = ttk.Entry(frm, width=14)
        amt_entry.grid(row=6, column=1, sticky="w", padx=8, pady=4)

        unit_combo = ttk.Combobox(frm, width=12, state="readonly")
        unit_combo.grid(row=6, column=2, sticky="w", pady=4)
        # default units based on mode
        def _set_units_for_mode(*_):
            if mode_var.get() == "m2n":
                unit_combo.config(values=["g", "kg", "mg"])
                if unit_combo.get() not in ("g", "kg", "mg"):
                    unit_combo.set("g")
            else:
                unit_combo.config(values=["mol"])
                unit_combo.set("mol")
        _set_units_for_mode()
        mode_var.trace_add("write", _set_units_for_mode)

        out_var = tk.StringVar(value="")
        out_lbl = ttk.Label(frm, textvariable=out_var, font=("Segoe UI", 10, "bold"))
        out_lbl.grid(row=7, column=0, columnspan=4, sticky="w", pady=(10,6))

        def _update_preview():
            user_text = sub_entry.get().strip()
            if not user_text:
                formula_var.set("")
                mm_var.set("")
                return
            f = name_to_formula(user_text)
            formula_var.set(f"Formula: {f}")
            try:
                mm = molar_mass(f)
                mm_var.set(f"Molar mass = {mm:g} g/mol")
            except Exception:
                mm_var.set("Molar mass = (unknown — check formula/elements)")

        def _refresh_suggestions(_ev=None):
            q = sub_entry.get()
            items = suggest_substances(q, limit=60)
            sugg.delete(0, tk.END)
            for label, _formula in items:
                sugg.insert(tk.END, label)
            if items:
                sugg.configure(height=min(8, len(items)))
                sugg.grid()
                sugg_scrollx.grid()
            else:
                sugg.grid_remove()
                sugg_scrollx.grid_remove()
            _update_preview()

        def _accept_suggestion(_ev=None):
            sel = sugg.curselection()
            if not sel:
                return
            label = sugg.get(sel[0])  # "name — FORMULA"
            # put the name part back into entry, but we can also directly insert the formula
            # we'll insert the name so preview shows both name & parsed formula consistently
            name_part = label.split(" — ", 1)[0]
            sub_entry.delete(0, tk.END)
            sub_entry.insert(0, name_part)
            _update_preview()
            amt_entry.focus_set()
            sugg.grid_remove()
            sugg_scrollx.grid_remove()

        sub_entry.bind("<KeyRelease>", _refresh_suggestions)
        sub_entry.bind("<Down>", lambda e: (sugg.focus_set(), sugg.selection_clear(0, tk.END), sugg.selection_set(0), sugg.activate(0)) if sugg.winfo_ismapped() else None)
        sugg.bind("<Return>", _accept_suggestion)
        sugg.bind("<Double-Button-1>", _accept_suggestion)
        sugg.bind("<Escape>", lambda e: (sugg.grid_remove(), sugg_scrollx.grid_remove()))

        def compute():
            substance_text = sub_entry.get().strip()
            if not substance_text:
                messagebox.showerror("Input needed", "Enter a substance (e.g., oxygen, O2, water, CO2).", parent=win)
                return
            f = name_to_formula(substance_text)
            try:
                mm = molar_mass(f)  # g/mol
            except Exception as ex:
                messagebox.showerror("Invalid substance", str(ex), parent=win)
                return

            try:
                amt = float(amt_entry.get().strip().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid amount", "Enter a numeric amount.", parent=win)
                return

            u = unit_combo.get().strip()
            if mode_var.get() == "m2n":
                # mass → mol
                if u == "kg":
                    mass_g = amt * 1000.0
                elif u == "mg":
                    mass_g = amt * 0.001
                else:
                    mass_g = amt
                n = mass_g / mm
                out_var.set(f"{amt:g} {u} {substance_text} → {n:g} mol   (MM = {mm:g} g/mol; formula {f})")
                return n, mm, f
            else:
                # mol → mass
                n = amt
                mass_g = n * mm
                out_var.set(f"{n:g} mol {substance_text} → {mass_g:g} g   (MM = {mm:g} g/mol; formula {f})")
                return n, mm, f

        self._add_compute_bar(win, compute)


        def add_to_known():
            res = compute()
            if not res:
                return
            n, mm, f = res
            if mode_var.get() == "m2n":
                # add n (mol)
                self.known_base["n"] = float(n)
                self.known_ui["n"] = {"unit": "mol", "display_value": float(n)}
            else:
                # add m in kg (our base is kg)
                mass_g = n * mm
                mass_kg = mass_g / 1000.0
                self.known_base["m"] = mass_kg
                self.known_ui["m"] = {"unit": "kg", "display_value": mass_kg}
            self._refresh_known_table()
            self._refresh_equation_list()
            messagebox.showinfo("Added", "Value added to known variables.", parent=win)

        btns = ttk.Frame(frm)
        btns.grid(row=8, column=0, columnspan=4, sticky="w", pady=(10, 0))
        ttk.Button(btns, text="Compute", command=compute).pack(side=tk.LEFT, padx=(0, 8))
        ttk.Button(btns, text="Add to Known Variables", command=add_to_known).pack(side=tk.LEFT)

        # initial populate
        _refresh_suggestions()

    # ---------------- end modal ----------------

    # ---------------- start of pressure converter modal ----------------
    def _open_pressure_converter(self):
        win = tk.Toplevel(self)
        win.title("Pressure Converter")
        win.transient(self)
        win.grab_set()
        win.geometry("520x260")

        frm = ttk.Frame(win, padding=12)
        frm.pack(fill=tk.BOTH, expand=True)

        # Input
        ttk.Label(frm, text="Value:").grid(row=0, column=0, sticky="w", pady=(0,4))
        val_entry = ttk.Entry(frm, width=12)
        val_entry.grid(row=0, column=1, sticky="w", padx=(6,12), pady=(0,4))
        val_entry.insert(0, "101.325")  # friendly default

        units = units_for("p") or ["Pa", "kPa", "bar", "atm", "psi"]

        ttk.Label(frm, text="From:").grid(row=1, column=0, sticky="w")
        from_combo = ttk.Combobox(frm, width=10, state="readonly", values=units)
        from_combo.grid(row=1, column=1, sticky="w", padx=(6,12))
        from_combo.set("kPa" if "kPa" in units else units[0])

        ttk.Label(frm, text="To:").grid(row=1, column=2, sticky="w")
        to_combo = ttk.Combobox(frm, width=10, state="readonly", values=units)
        to_combo.grid(row=1, column=3, sticky="w", padx=(6,0))
        to_combo.set("bar" if "bar" in units else units[-1])

        # Output
        out_var = tk.StringVar(value="")
        out_lbl = ttk.Label(frm, textvariable=out_var, font=("Segoe UI", 10, "bold"))
        out_lbl.grid(row=2, column=0, columnspan=4, sticky="w", pady=(10,4))

        def compute():
            try:
                x = float(val_entry.get().strip().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid value", "Please enter a numeric value.", parent=win)
                return
            u_from = from_combo.get().strip()
            u_to = to_combo.get().strip()
            # Convert via base (Pa)
            base = convert_to_base("p", x, u_from)
            y = convert_from_base("p", base, u_to)
            out_var.set(f"{x:g} {u_from} = {y:g} {u_to}")
            return y, u_to, base  # also return base value for adding to known vars
        
        self._add_compute_bar(win, compute)


        # Action buttons
        btns = ttk.Frame(frm)
        btns.grid(row=3, column=0, columnspan=4, sticky="w", pady=(10,0))
        ttk.Button(btns, text="Convert", command=compute).pack(side=tk.LEFT, padx=(0,8))

        # Add-to-known controls
        ttk.Label(frm, text="Add as:").grid(row=4, column=0, sticky="w", pady=(12,0))
        add_as_var = tk.StringVar(value="p")
        add_row = ttk.Frame(frm)
        add_row.grid(row=4, column=1, columnspan=3, sticky="w", pady=(8,0))
        for label, key in [("p (pressure)", "p"), ("p1 (initial pressure)", "p1"), ("p2 (final pressure)", "p2")]:
            ttk.Radiobutton(add_row, text=label, variable=add_as_var, value=key).pack(side=tk.LEFT, padx=(0,12))

        def add_to_known():
            res = compute()
            if not res:
                return
            y, u_to, base = res
            canon = add_as_var.get()  # 'p', 'p1', or 'p2'
            # Store base (Pa) internally, but reflect user's chosen unit/value in UI
            self.known_base[canon] = float(base)
            self.known_ui[canon] = {"unit": u_to, "display_value": float(y)}
            self._refresh_known_table()
            self._refresh_equation_list()
            messagebox.showinfo("Added", f"Added {canon} = {y:g} {u_to}", parent=win)

        ttk.Button(btns, text="Add to Known Variables", command=add_to_known).pack(side=tk.LEFT)
    # ---------------- end modal ----------------
    """ def _open_hf_water_modal(self):
        win = tk.Toplevel(self)
        win.title("Formation Enthalpy ΔH°f — H₂O(g)")
        win.transient(self)
        win.grab_set()
        win.geometry("560x420")

        frm = ttk.Frame(win, padding=12)
        frm.pack(fill=tk.BOTH, expand=True)

        # Mode select
        mode = tk.StringVar(value="liq_plus_vap")  # or "from_rxn"
        mrow = ttk.Frame(frm)
        mrow.pack(anchor="w", pady=(0,8))
        ttk.Radiobutton(mrow, text="From H₂O(l) and ΔHᵥₐₚ", variable=mode, value="liq_plus_vap").pack(side=tk.LEFT, padx=(0,14))
        ttk.Radiobutton(mrow, text="From formation reaction (elements → H₂O(g))", variable=mode, value="from_rxn").pack(side=tk.LEFT)

        # --- Inputs for method A: Hf(liq) + ΔHvap ---
        boxA = ttk.LabelFrame(frm, text="Method A")
        boxA.pack(fill=tk.X, pady=(6,6))

        a1 = ttk.Frame(boxA); a1.pack(fill=tk.X, pady=4)
        ttk.Label(a1, text="ΔH°f[H₂O(l)]").grid(row=0, column=0, sticky="w")
        hf_liq_val = ttk.Entry(a1, width=14); hf_liq_val.grid(row=0, column=1, padx=6)
        hf_liq_unit = ttk.Combobox(a1, width=10, state="readonly", values=["kJ/mol","J/mol"])
        hf_liq_unit.grid(row=0, column=2); hf_liq_unit.set("kJ/mol")

        a2 = ttk.Frame(boxA); a2.pack(fill=tk.X, pady=4)
        ttk.Label(a2, text="ΔHᵥₐₚ (same substance)").grid(row=0, column=0, sticky="w")
        hvap_val = ttk.Entry(a2, width=14); hvap_val.grid(row=0, column=1, padx=6)
        hvap_unit = ttk.Combobox(a2, width=10, state="readonly", values=["kJ/mol","J/mol"])
        hvap_unit.grid(row=0, column=2); hvap_unit.set("kJ/mol")

        # --- Inputs for method B: reaction enthalpy equals ΔHf° ---
        boxB = ttk.LabelFrame(frm, text="Method B")
        boxB.pack(fill=tk.X, pady=(6,6))

        b1 = ttk.Frame(boxB); b1.pack(fill=tk.X, pady=4)
        ttk.Label(b1, text="ΔH°rxn for  H₂(g) + ½ O₂(g) → H₂O(g)").grid(row=0, column=0, sticky="w")
        rxn_val = ttk.Entry(b1, width=14); rxn_val.grid(row=0, column=1, padx=6)
        rxn_unit = ttk.Combobox(b1, width=10, state="readonly", values=["kJ/mol","J/mol"])
        rxn_unit.grid(row=0, column=2); rxn_unit.set("kJ/mol")

        # Output
        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(10,4))

        def compute():
            try:
                if mode.get() == "liq_plus_vap":
                    # Convert both to base J/mol
                    hfl = float(hf_liq_val.get().strip().replace(",", "."))
                    hv  = float(hvap_val.get().strip().replace(",", "."))
                    hfl_base = convert_to_base("ΔH°f_liq", hfl, hf_liq_unit.get())
                    hv_base  = convert_to_base("ΔH_vap",    hv,  hvap_unit.get())
                    hfg_base = hfl_base + hv_base
                else:
                    # ΔH°f[H2O(g)] = ΔH°rxn for formation from elements (per 1 mol product)
                    dhr = float(rxn_val.get().strip().replace(",", "."))
                    hfg_base = convert_to_base("ΔH°rxn", dhr, rxn_unit.get())
            except ValueError:
                messagebox.showerror("Invalid input", "Please enter numeric values.", parent=win)
                return

            # Display both J/mol and kJ/mol
            hfg_kJmol = convert_from_base("ΔH°f_gas", hfg_base, "kJ/mol")
            out.set(f"ΔH°f[H₂O(g)] = {hfg_kJmol:g} kJ/mol   (= {hfg_base:g} J/mol)")
            return hfg_base

        # Buttons
        btns = ttk.Frame(frm); btns.pack(anchor="w", pady=(8,0))
        ttk.Button(btns, text="Compute", command=compute).pack(side=tk.LEFT, padx=(0,8))

        def add_known():
            hfg_base = compute()
            if hfg_base is None:
                return
            # Save as ΔH°f_gas
            self.known_base["ΔH°f_gas"] = float(hfg_base)
            disp = convert_from_base("ΔH°f_gas", hfg_base, "kJ/mol")
            self.known_ui["ΔH°f_gas"] = {"unit": "kJ/mol", "display_value": float(disp)}
            self._refresh_known_table()
            self._refresh_equation_list()
            messagebox.showinfo("Added", "Added ΔH°f[H₂O(g)] to Known Variables.", parent=win)

        ttk.Button(btns, text="Add to Known Variables", command=add_known).pack(side=tk.LEFT) """
    
    def _open_thermo_lookup(self):
        win = tk.Toplevel(self)
        win.title("Thermo lookup — ΔH°f, ΔG°f, S° (298 K)")
        win.transient(self); win.grab_set(); win.geometry("640x420")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # Substance entry + suggestions (reuse the same list UX as mol converter)
        ttk.Label(frm, text="Substance (name/symbol/formula):").grid(row=0, column=0, sticky="w")
        ent = ttk.Entry(frm, width=36); ent.grid(row=0, column=1, columnspan=3, sticky="we", padx=8, pady=4)
        ent.insert(0, "water")

        sugg = tk.Listbox(frm, height=0); sugg.grid(row=1, column=0, columnspan=4, sticky="we")
        scx = ttk.Scrollbar(frm, orient="horizontal", command=sugg.xview)
        sugg.configure(xscrollcommand=scx.set); scx.grid(row=2, column=0, columnspan=4, sticky="we", pady=(0,6))
        sugg.grid_remove(); scx.grid_remove()

        # Phase dropdown
        ttk.Label(frm, text="Phase:").grid(row=3, column=0, sticky="w", pady=(4,0))
        phase = ttk.Combobox(frm, width=8, state="readonly", values=["g","l","s","aq"])
        phase.grid(row=3, column=1, sticky="w", padx=8, pady=(4,0))
        phase.set("l")

        # Output labels
        val_Hf = tk.StringVar(value="ΔH°f: —")
        val_Gf = tk.StringVar(value="ΔG°f: —")
        val_S  = tk.StringVar(value="S°: —")
        ttk.Label(frm, textvariable=val_Hf, font=("Segoe UI", 10, "bold")).grid(row=4, column=0, columnspan=4, sticky="w", pady=(12,2))
        ttk.Label(frm, textvariable=val_Gf).grid(row=5, column=0, columnspan=4, sticky="w", pady=2)
        ttk.Label(frm, textvariable=val_S).grid(row=6, column=0, columnspan=4, sticky="w", pady=2)

        # Helper: current formula & refresh
        current_formula = {"f": ""}

        def _refresh_suggestions(_=None):
            items = suggest_substances(ent.get(), limit=60)
            sugg.delete(0, tk.END)
            for label, f in items:
                sugg.insert(tk.END, label)
            if items:
                sugg.configure(height=min(8, len(items))); sugg.grid(); scx.grid()
            else:
                sugg.grid_remove(); scx.grid_remove()
            _refresh_formula_and_phases()

        def _accept_suggestion(_=None):
            sel = sugg.curselection()
            if not sel: return
            name = sugg.get(sel[0]).split(" — ", 1)[0]
            ent.delete(0, tk.END); ent.insert(0, name)
            ent.focus_set(); sugg.grid_remove(); scx.grid_remove()
            _refresh_formula_and_phases()

        def _refresh_formula_and_phases():
            f = name_to_formula(ent.get().strip())
            current_formula["f"] = f
            phs = phases_for(f) or ["g","l","s","aq"]
            phase.config(values=phs)
            if phase.get() not in phs:
                phase.set(phs[0])
            _show_values()

        def _show_values(*_):
            f = current_formula["f"]; ph = phase.get()
            data = get_thermo(f, ph)
            if data:
                val_Hf.set(f"ΔH°f: {data['Hf']:g} kJ/mol")
                val_Gf.set(f"ΔG°f: {data['Gf']:g} kJ/mol")
                val_S.set(f"S°: {data['S']:g} J/(mol·K)")
            else:
                val_Hf.set("ΔH°f: (no data)")
                val_Gf.set("ΔG°f: (no data)")
                val_S.set("S°: (no data)")

        ent.bind("<KeyRelease>", _refresh_suggestions)
        ent.bind("<Down>", lambda e: (sugg.focus_set(), sugg.selection_clear(0, tk.END), sugg.selection_set(0), sugg.activate(0)) if sugg.winfo_ismapped() else None)
        sugg.bind("<Return>", _accept_suggestion)
        sugg.bind("<Double-Button-1>", _accept_suggestion)
        phase.bind("<<ComboboxSelected>>", _show_values)

        # Buttons: copy or add to known (optional)
        btns = ttk.Frame(frm); btns.grid(row=7, column=0, columnspan=4, sticky="w", pady=(12,0))
        def _add_known(kind: str):
            f = current_formula["f"]; ph = phase.get(); data = get_thermo(f, ph)
            if not data:
                tk.messagebox.showerror("No data", "No table value for that species/phase.", parent=win); return
            if kind == "Hf":
                base = data["Hf"] * 1000.0  # -> J/mol
                self.known_base["ΔH°f"] = base
                self.known_ui["ΔH°f"] = {"unit": "kJ/mol", "display_value": data["Hf"]}
            elif kind == "Gf":
                base = data["Gf"] * 1000.0
                self.known_base["ΔG°f"] = base
                self.known_ui["ΔG°f"] = {"unit": "kJ/mol", "display_value": data["Gf"]}
            elif kind == "S":
                base = data["S"]  # already J/(mol·K)
                self.known_base["S°"] = base
                self.known_ui["S°"] = {"unit": "J/(mol·K)", "display_value": data["S"]}
            self._refresh_known_table(); self._refresh_equation_list()
            tk.messagebox.showinfo("Added", "Value added to Known Variables.", parent=win)

        ttk.Button(btns, text="Add ΔH°f", command=lambda: _add_known("Hf")).pack(side=tk.LEFT, padx=(0,6))
        ttk.Button(btns, text="Add ΔG°f", command=lambda: _add_known("Gf")).pack(side=tk.LEFT, padx=(0,6))
        ttk.Button(btns, text="Add S°",    command=lambda: _add_known("S")).pack(side=tk.LEFT)

        _refresh_suggestions()


    def _open_reaction_builder(self):
        import tkinter as tk
        from tkinter import ttk, messagebox
        import math

        win = tk.Toplevel(self)
        win.title("Reaction builder — ΔH°, ΔS°, ΔG° (298 K)")
        win.transient(self); win.grab_set()
        win.geometry("1620x960")  # wider to fit per-species thermo columns

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # ---------- Picker row (choose species, phase, coeff) ----------
        ttk.Label(frm, text="Substance:").grid(row=0, column=0, sticky="w")
        ent = ttk.Entry(frm, width=32); ent.grid(row=0, column=1, sticky="w", padx=(6,8))
        ent.insert(0, "oxygen")

        sugg = tk.Listbox(frm, height=0); sugg.grid(row=1, column=0, columnspan=8, sticky="we")
        scx = ttk.Scrollbar(frm, orient="horizontal", command=sugg.xview)
        sugg.configure(xscrollcommand=scx.set); scx.grid(row=2, column=0, columnspan=8, sticky="we", pady=(0,6))
        sugg.grid_remove(); scx.grid_remove()

        ttk.Label(frm, text="Phase:").grid(row=0, column=2, sticky="w")
        phase = ttk.Combobox(frm, width=8, state="readonly", values=["g","l","s","aq"])
        phase.grid(row=0, column=3, sticky="w"); phase.set("g")

        ttk.Label(frm, text="Coeff (ν):").grid(row=0, column=4, sticky="w")
        coeff = ttk.Entry(frm, width=10); coeff.grid(row=0, column=5, sticky="w"); coeff.insert(0, "1")

        # live preview of table values for the picker
        pv_Hf = tk.StringVar(value="ΔH°f: —")
        pv_Gf = tk.StringVar(value="ΔG°f: —")
        pv_S  = tk.StringVar(value="S°: —")
        ttk.Label(frm, textvariable=pv_Hf).grid(row=0, column=6, sticky="w", padx=(14,0))
        ttk.Label(frm, textvariable=pv_Gf).grid(row=0, column=7, sticky="w", padx=(10,0))
        ttk.Label(frm, textvariable=pv_S).grid(row=0, column=8, sticky="w", padx=(10,0))

        # ---------- Two panels with Treeviews (show ν, species, ΔH°f, ΔG°f, S°) ----------
        left = ttk.LabelFrame(frm, text="Reactants (left)")
        right = ttk.LabelFrame(frm, text="Products (right)")
        left.grid(row=3, column=0, columnspan=4, sticky="nsew", padx=(0,6), pady=(6,6))
        right.grid(row=3, column=4, columnspan=5, sticky="nsew", padx=(6,0), pady=(6,6))
        frm.grid_rowconfigure(3, weight=1)
        for c in (0,1,2,3,4,5,6,7,8):
            frm.grid_columnconfigure(c, weight=1)

        cols = ("nu","species","Hf","Gf","S")
        head = {
            "nu": "ν",
            "species": "Species (phase)",
            "Hf": "ΔH°f (kJ/mol)",
            "Gf": "ΔG°f (kJ/mol)",
            "S": "S° (J/mol·K)"
        }

        def _make_tree(parent):
            tv = ttk.Treeview(parent, columns=cols, show="headings", height=10)
            tv.pack(fill=tk.BOTH, expand=True)
            # columns sizing
            tv.column("nu", width=60, anchor="e")
            tv.column("species", width=240, anchor="w")
            tv.column("Hf", width=120, anchor="e")
            tv.column("Gf", width=120, anchor="e")
            tv.column("S", width=120, anchor="e")
            for c in cols:
                tv.heading(c, text=head[c])
            # vertical scrollbar
            vs = ttk.Scrollbar(parent, orient="vertical", command=tv.yview)
            tv.configure(yscrollcommand=vs.set)
            vs.place(relx=1.0, rely=0, relheight=1.0, anchor="ne")
            return tv

        treeL = _make_tree(left)
        treeR = _make_tree(right)

        # Internal storage: list of dicts {side, formula, phase, nu, Hf, Gf, S}
        items = []

        # ---------- Suggestion helpers ----------
        def _refresh_suggestions(_=None):
            from molar_masses import suggest_substances, name_to_formula
            data = suggest_substances(ent.get(), limit=60)
            sugg.delete(0, tk.END)
            for label, f in data:
                sugg.insert(tk.END, label)
            if data:
                sugg.configure(height=min(8, len(data))); sugg.grid(); scx.grid()
            else:
                sugg.grid_remove(); scx.grid_remove()
            _update_phase_and_preview()

        def _accept_suggestion(_=None):
            label = sugg.get(sugg.curselection()[0]) if sugg.curselection() else ""
            name = label.split(" — ", 1)[0] if label else ent.get()
            ent.delete(0, tk.END); ent.insert(0, name)
            sugg.grid_remove(); scx.grid_remove()
            _update_phase_and_preview()

        def _update_phase_and_preview():
            from molar_masses import name_to_formula
            f = name_to_formula(ent.get().strip())
            phs = phases_for(f) or ["g","l","s","aq"]
            phase.config(values=phs)
            if phase.get() not in phs: phase.set(phs[0])
            _preview_values(f, phase.get())

        def _preview_values(f, ph):
            data = get_thermo(f, ph)
            if data:
                pv_Hf.set(f"ΔH°f: {data.get('Hf','—'):g} kJ/mol" if "Hf" in data else "ΔH°f: —")
                pv_Gf.set(f"ΔG°f: {data.get('Gf','—'):g} kJ/mol" if "Gf" in data else "ΔG°f: —")
                pv_S.set( f"S°: {data.get('S','—'):g} J/(mol·K)"  if "S"  in data else "S°: —")
            else:
                pv_Hf.set("ΔH°f: —"); pv_Gf.set("ΔG°f: —"); pv_S.set("S°: —")

        phase.bind("<<ComboboxSelected>>", lambda e: _update_phase_and_preview())

        # ---------- Add/Remove rows ----------
        def _add(side: str):
            from molar_masses import name_to_formula
            f = name_to_formula(ent.get().strip())
            ph = phase.get()
            try:
                nu = float(coeff.get().strip().replace(",", "."))
            except ValueError:
                messagebox.showerror("Bad ν", "Coefficient must be a number.", parent=win); return
            data = get_thermo(f, ph) or {}
            Hf = data.get("Hf"); Gf = data.get("Gf"); S = data.get("S")

            row = {"side": side, "formula": f, "phase": ph, "nu": nu, "Hf": Hf, "Gf": Gf, "S": S}
            items.append(row)

            vals = (
                f"{nu:g}",
                f"{f}({ph})",
                f"{Hf:g}" if Hf is not None else "—",
                f"{Gf:g}" if Gf is not None else "—",
                f"{S:g}" if S  is not None else "—",
            )
            (treeL if side == "L" else treeR).insert("", tk.END, values=vals)

        def _remove(side: str):
            tv = treeL if side == "L" else treeR
            sel = tv.selection()
            if not sel:
                return
            idx = tv.index(sel[0])
            # remove idx-th item of this side from items
            count = -1
            for i, it in enumerate(items):
                if it["side"] == side:
                    count += 1
                    if count == idx:
                        items.pop(i); break
            tv.delete(sel[0])

        btnrow = ttk.Frame(frm); btnrow.grid(row=4, column=0, columnspan=9, sticky="w")
        ttk.Button(btnrow, text="Add → Reactant", command=lambda: _add("L")).pack(side=tk.LEFT, padx=(0,8), pady=(6,0))
        ttk.Button(btnrow, text="Add → Product", command=lambda: _add("R")).pack(side=tk.LEFT, padx=(0,8), pady=(6,0))
        ttk.Button(btnrow, text="Remove Reactant", command=lambda: _remove("L")).pack(side=tk.LEFT, padx=(16,8))
        ttk.Button(btnrow, text="Remove Product", command=lambda: _remove("R")).pack(side=tk.LEFT, padx=(8,8))

        # ---------- Outputs (numbers) ----------
        outH = tk.StringVar(value="ΔH°rxn: —")
        outS = tk.StringVar(value="ΔS°rxn: —")
        outG = tk.StringVar(value="ΔG°rxn: —")
        ttk.Label(frm, textvariable=outH, font=("Segoe UI", 10, "bold")).grid(row=5, column=0, columnspan=9, sticky="w", pady=(10,2))
        ttk.Label(frm, textvariable=outS).grid(row=6, column=0, columnspan=9, sticky="w", pady=2)
        ttk.Label(frm, textvariable=outG).grid(row=7, column=0, columnspan=9, sticky="w", pady=2)

        # --- Nonstandard ΔG (ΔG = ΔG° + RT ln Q) ------------------------------------
        ns_box = ttk.LabelFrame(frm, text="ΔG at non-standard conditions (298 K)")
        ns_box.grid(row=8, column=0, columnspan=9, sticky="we", pady=(8,0))

        mode = tk.StringVar(value="P")  # P for partial pressure (atm), C for concentration (M)
        ttk.Radiobutton(ns_box, text="Use partial pressures (atm) → Qp", variable=mode, value="P").grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(ns_box, text="Use concentrations (M) → Qc",       variable=mode, value="C").grid(row=0, column=1, sticky="w", padx=(8,0))

        tbl_ns = ttk.Frame(ns_box); tbl_ns.grid(row=1, column=0, columnspan=9, sticky="we", pady=(6,2))
        ns_entries = []  # list of dicts: {"sp": "SO2(g)", "nu": +/−power, "entry": Entry}

        def _build_ns_table():
            for w in tbl_ns.winfo_children(): w.destroy()
            ns_entries.clear()
            ttk.Label(tbl_ns, text="Species").grid(row=0, column=0, sticky="w")
            ttk.Label(tbl_ns, text="Value (atm or M)").grid(row=0, column=1, sticky="w")
            r = 1
            # collect species with sign: products positive, reactants negative
            for side in ("L","R"):
                for it in items:
                    if it["side"] != side: continue
                    sp = f"{it['formula']}({it['phase']})"
                    power = it["nu"] if side=="R" else -it["nu"]
                    ttk.Label(tbl_ns, text=f"{sp}   (power {power:g})").grid(row=r, column=0, sticky="w")
                    e = ttk.Entry(tbl_ns, width=12); e.grid(row=r, column=1, sticky="w"); e.insert(0, "1.0")
                    ns_entries.append({"sp": sp, "nu": power, "entry": e})
                    r += 1
            if r == 1:
                ttk.Label(tbl_ns, text="(Add reactants/products above, then click ‘Use current reaction’.)").grid(row=1, column=0, columnspan=2, sticky="w")

        ttk.Button(ns_box, text="Use current reaction", command=_build_ns_table).grid(row=0, column=2, padx=(12,0))

        out_nonstd = tk.StringVar(value="")
        ttk.Label(ns_box, textvariable=out_nonstd, font=("Segoe UI", 10, "bold")).grid(row=2, column=0, columnspan=9, sticky="w", pady=(4,2))

        def _compute_nonstandard():
            # Prefer manual ΔG°rxn if provided; otherwise compute from species
            txt = (man_dG0.get() or "").strip()
            if txt:
                try:
                    dG0 = float(txt.replace(",", "."))  # kJ/mol
                except ValueError:
                    messagebox.showerror("Bad ΔG°rxn", "ΔG°rxn must be numeric (kJ/mol).", parent=win); return
                source_msg = f"manual ΔG°rxn = {dG0:.6g} kJ/mol"
            else:
                res = _compute()  # updates ΔG°rxn from table if available
                if not res or res.get("dG_kJ") is None:
                    messagebox.showerror("Need ΔG°rxn",
                                        "Enter ΔG°rxn above or add ΔG°f values for all species.",
                                        parent=win); return
                dG0 = res["dG_kJ"]  # kJ/mol
                source_msg = f"ΔG°rxn (from table) = {dG0:.6g} kJ/mol"

            if not ns_entries:
                messagebox.showerror("No species",
                                    "Click ‘Use current reaction’ and enter P or C values.",
                                    parent=win); return

            # Build Q (products^ν / reactants^ν). Entries already store +ν for products, −ν for reactants
            try:
                Q = 1.0
                for row in ns_entries:
                    val = float(row["entry"].get().replace(",", "."))
                    if val <= 0: raise ValueError
                    Q *= val ** row["nu"]
            except ValueError:
                messagebox.showerror("Bad input", "All P/C values must be positive numbers.", parent=win); return

            # Thermo math
            R_J = 8.314462618      # J mol^-1 K^-1
            RkJ = 0.008314462618   # kJ mol^-1 K^-1
            T   = 298.15
            dG  = dG0 + RkJ*T*math.log(Q)              # kJ/mol
            K   = math.exp(-(dG0*1000.0)/(R_J*T))      # dimensionless

            # Direction from Q vs K
            rel = abs(math.log(Q/K))                   # tolerance in log space
            if rel < 1e-6:
                dir_txt = "Q ≈ K — at equilibrium (no net shift)"
                arrow = "↔"
            elif Q < K:
                dir_txt = "Q < K — proceeds to the right (toward products)"
                arrow = "→"
            else:
                dir_txt = "Q > K — proceeds to the left (toward reactants)"
                arrow = "←"

            kind = "Qp" if mode.get() == "P" else "Qc"
            verdict = "spontaneous" if dG < 0 else ("non-spontaneous" if dG > 0 else "at equilibrium tendency")

            out_nonstd.set(
                f"{source_msg}\n"
                f"{kind} = {Q:.6g};  K (298 K) = {K:.3e}  →  {dir_txt} {arrow}\n"
                f"ΔG = ΔG° + RT ln Q = {dG0:.6g} + {(RkJ*T):.6g}·ln({Q:.6g}) = {dG:.6g} kJ/mol  — {verdict} at 298 K"
            )



        ttk.Button(ns_box, text="Compute ΔG (non-standard)", command=_compute_nonstandard).grid(row=0, column=3, padx=(10,0))

        # Manual ΔG°rxn override (useful when a problem gives ΔG° directly)
        man_dG0 = tk.StringVar()
        ttk.Label(ns_box, text="Given ΔG°rxn (kJ/mol) [optional]:").grid(row=0, column=4, sticky="e", padx=(12,4))
        ttk.Entry(ns_box, width=12, textvariable=man_dG0).grid(row=0, column=5, sticky="w")

        def _compute_K_from_dG0():
            val = (man_dG0.get() or "").strip()
            if not val:
                messagebox.showerror("Need ΔG°rxn", "Enter ΔG°rxn (kJ/mol) or compute it above.", parent=win); return
            try:
                dG0 = float(val.replace(",", "."))
            except ValueError:
                messagebox.showerror("Bad ΔG°rxn", "ΔG°rxn must be numeric (kJ/mol).", parent=win); return
            R = 8.314462618  # J/(mol·K)
            T = 298.15
            K = math.exp(-(dG0*1000.0)/(R*T))
            out_nonstd.set(f"K (298 K) from ΔG° = {dG0:.6g} kJ/mol  →  K ≈ {K:.3e}  (≈ {K:.6g})")

        def _compute_Q_only():
            # Build Q from table entries
            Q_num = 1.0; Q_den = 1.0
            try:
                for row in ns_entries:
                    val = float(row["entry"].get().replace(",", "."))
                    if val <= 0: raise ValueError
                    power = row["nu"]
                    if power > 0: Q_num *= val**power
                    elif power < 0: Q_den *= val**(-power)
            except ValueError:
                messagebox.showerror("Bad input", "All P/C values must be positive numbers.", parent=win); return
            Q = Q_num / Q_den
            kind = "Qp" if mode.get()=="P" else "Qc"
            out_nonstd.set(f"{kind} = {Q:.6g}")
        
        ttk.Button(ns_box, text="Compute K from ΔG°", command=_compute_K_from_dG0).grid(row=0, column=6, padx=(10,0))
        ttk.Button(ns_box, text="Compute Q only", command=_compute_Q_only).grid(row=0, column=7, padx=(10,0))





        # ---------- Explanation panel ----------
        explain = tk.Text(frm, height=8, wrap="word")
        explain.grid(row=10, column=0, columnspan=9, sticky="nsew", pady=(8,0))
        frm.grid_rowconfigure(9, weight=1)
        sv = ttk.Scrollbar(frm, orient="vertical", command=explain.yview)
        explain.configure(yscrollcommand=sv.set); sv.grid(row=9, column=9, sticky="ns")
        explain.configure(state="disabled")

        def _write_explain(lines):
            explain.configure(state="normal")
            explain.delete("1.0", "end")
            explain.insert("1.0", "\n".join(lines))
            explain.configure(state="disabled")

        # ---------- Compute & Add to Known ----------
        def _compute():
            import math

            # sums (kJ/mol for H,G; J/(mol·K) for S) and how many valid terms we had
            sum_H_L = sum_G_L = 0.0
            sum_S_L = 0.0
            sum_H_R = sum_G_R = 0.0
            sum_S_R = 0.0
            cnt_H = cnt_S = cnt_G = 0
            missing = []

            # pretty per-side contribution strings (uses items already in scope)
            def contrib_str(side):
                lines_H, lines_S, lines_G = [], [], []
                for it in items:
                    if it["side"] != side: 
                        continue
                    sp = f"{it['formula']}({it['phase']})"
                    nu = it["nu"]
                    Hf = it.get("Hf"); Gf = it.get("Gf"); S = it.get("S")
                    lines_H.append(f"  {nu:g}×ΔH°f[{sp}] = {nu*Hf:.6g} kJ" if Hf is not None else f"  {nu:g}×ΔH°f[{sp}] = —")
                    lines_S.append(f"  {nu:g}×S°[{sp}] = {nu*S:.6g} J/K"   if S  is not None else f"  {nu:g}×S°[{sp}] = —")
                    lines_G.append(f"  {nu:g}×ΔG°f[{sp}] = {nu*Gf:.6g} kJ" if Gf is not None else f"  {nu:g}×ΔG°f[{sp}] = —")
                return lines_H, lines_S, lines_G

            # accumulate sums from the table (use get_thermo to re-read authoritative values)
            for it in items:
                data = get_thermo(it["formula"], it["phase"]) or {}
                nu = it["nu"]
                isL = (it["side"] == "L")

                if "Hf" in data:
                    (sum_H_L if isL else sum_H_R)  # just to choose var
                    if isL: sum_H_L += nu * data["Hf"]
                    else:   sum_H_R += nu * data["Hf"]
                    cnt_H += 1
                else:
                    missing.append(f"ΔH°f {it['formula']}({it['phase']})")

                if "Gf" in data:
                    if isL: sum_G_L += nu * data["Gf"]
                    else:   sum_G_R += nu * data["Gf"]
                    cnt_G += 1
                else:
                    missing.append(f"ΔG°f {it['formula']}({it['phase']})")

                if "S" in data:
                    if isL: sum_S_L += nu * data["S"]
                    else:   sum_S_R += nu * data["S"]
                    cnt_S += 1
                else:
                    missing.append(f"S° {it['formula']}({it['phase']})")

            if missing:
                tk.messagebox.showwarning(
                    "Missing data",
                    "No/partial table data for: " + ", ".join(sorted(set(missing))) +
                    "\n(You can add them in thermo_data.py.)",
                    parent=win
                )

            # decide availability and compute deltas
            haveH = cnt_H > 0
            haveS = cnt_S > 0
            haveG = cnt_G > 0
            dH = (sum_H_R - sum_H_L) if haveH else None   # kJ/mol
            dS = (sum_S_R - sum_S_L) if haveS else None   # J/(mol·K)
            dG = (sum_G_R - sum_G_L) if haveG else None   # kJ/mol

            # headline numbers + quick meanings
            if dH is not None:
                heat = "exothermic (releases heat)" if dH < 0 else ("endothermic (absorbs heat)" if dH > 0 else "thermally neutral")
                outH.set(f"ΔH°rxn = {dH:.6g} kJ/mol   — {heat}")
            else:
                outH.set("ΔH°rxn: (not enough ΔH°f data)")

            if dS is not None:
                entro = "entropy increases (more disorder)" if dS > 0 else ("entropy decreases (more order)" if dS < 0 else "no net entropy change")
                outS.set(f"ΔS°rxn = {dS:.6g} J/(mol·K)   — {entro}")
            else:
                outS.set("ΔS°rxn: (not enough S° data)")

            if dG is not None:
                R = 8.314462618  # J/(mol·K)
                T = 298.15
                K = math.exp(-(dG*1000.0)/(R*T))
                spont = "spontaneous at 298 K" if dG < 0 else ("non-spontaneous at 298 K" if dG > 0 else "at equilibrium tendency at 298 K")
                outG.set(f"ΔG°rxn = {dG:.6g} kJ/mol   — {spont}   (K≈{K:.3e})")
            else:
                outG.set("ΔG°rxn: (not enough ΔG°f data)")

            # algebra + narrative in the explanation box
            lines = []
            if missing:
                lines.append("Missing table entries (ignored in sums): " + ", ".join(sorted(set(missing))) + "\n")

            LH, LS, LG = contrib_str("L")
            RH, RS, RG = contrib_str("R")

            if dH is not None:
                lines.append("ΔH°rxn = ΣνΔH°f(products) − ΣνΔH°f(reactants)")
                lines += ["  Reactants:"] + LH + [f"  Σ_L = {sum_H_L:.6g} kJ"]
                lines += ["  Products:"] + RH + [f"  Σ_R = {sum_H_R:.6g} kJ"]
                lines.append(f"  ⇒ ΔH°rxn = {sum_H_R:.6g} − {sum_H_L:.6g} = {dH:.6g} kJ\n")
            if dS is not None:
                lines.append("ΔS°rxn = ΣνS°(products) − ΣνS°(reactants)")
                lines += ["  Reactants:"] + LS + [f"  Σ_L = {sum_S_L:.6g} J/K"]
                lines += ["  Products:"] + RS + [f"  Σ_R = {sum_S_R:.6g} J/K"]
                lines.append(f"  ⇒ ΔS°rxn = {sum_S_R:.6g} − {sum_S_L:.6g} = {dS:.6g} J/K\n")
            if dG is not None:
                lines.append("ΔG°rxn = ΣνΔG°f(products) − ΣνΔG°f(reactants)")
                lines += ["  Reactants:"] + LG + [f"  Σ_L = {sum_G_L:.6g} kJ"]
                lines += ["  Products:"] + RG + [f"  Σ_R = {sum_G_R:.6g} kJ"]
                lines.append(f"  ⇒ ΔG°rxn = {sum_G_R:.6g} − {sum_G_L:.6g} = {dG:.6g} kJ\n")

            lines.append("What these mean (at ~298 K):")
            if dH is not None:
                lines.append(f"  • ΔH°rxn {dH:+.6g} kJ/mol → " +
                            ("exothermic (releases heat)" if dH < 0 else
                            "endothermic (absorbs heat)" if dH > 0 else "thermally neutral"))
            if dS is not None:
                lines.append(f"  • ΔS°rxn {dS:+.6g} J/(mol·K) → " +
                            ("entropy increases (more disorder)" if dS > 0 else
                            "entropy decreases (more order)" if dS < 0 else "no net entropy change"))
            if dG is not None:
                R = 8.314462618; T = 298.15
                K = math.exp(-(dG*1000.0)/(R*T))
                lines.append(f"  • ΔG°rxn {dG:+.6g} kJ/mol → " +
                            ("spontaneous at 298 K (thermodynamically favorable)" if dG < 0 else
                            "non-spontaneous at 298 K (requires driving force)" if dG > 0 else
                            "at equilibrium tendency at 298 K"))
                lines.append(f"  • Estimated K (298 K): K ≈ exp(−ΔG°/RT) = {K:.3e}  (≈ {K:.6g})")
                lines.append("    Rule-of-thumb: |ΔG°| ≈ 5.7 kJ/mol ↔ 10× change in K at 298 K.")
            else:
                lines.append("  • ΔG°rxn not available → can’t comment on spontaneity/K.")

            _write_explain(lines)

            return {"dH_kJ": dH, "dS_J": dS, "dG_kJ": dG}


        def _add_known():
            res = _compute()
            if res is None: return
            if res.get("dH_kJ") is not None:
                base = res["dH_kJ"] * 1000.0
                self.known_base["ΔH°rxn"] = base
                self.known_ui["ΔH°rxn"] = {"unit": "kJ/mol", "display_value": res["dH_kJ"]}
            if res.get("dS_J") is not None:
                self.known_base["ΔS°rxn"] = res["dS_J"]
                self.known_ui["ΔS°rxn"] = {"unit": "J/(mol·K)", "display_value": res["dS_J"]}
            if res.get("dG_kJ") is not None:
                base = res["dG_kJ"] * 1000.0
                self.known_base["ΔG°rxn"] = base
                self.known_ui["ΔG°rxn"] = {"unit": "kJ/mol", "display_value": res["dG_kJ"]}
            self._refresh_known_table(); self._refresh_equation_list()
            messagebox.showinfo("Added", "Reaction values added to Known Variables.", parent=win)

        ctl = ttk.Frame(frm); ctl.grid(row=9, column=0, columnspan=9, sticky="w", pady=(8,0))
        ttk.Button(ctl, text="Compute", command=_compute).pack(side=tk.LEFT, padx=(0,8))
        ttk.Button(ctl, text="Add to Known Variables", command=_add_known).pack(side=tk.LEFT)

        # bindings
        ent.bind("<KeyRelease>", _refresh_suggestions)
        ent.bind("<Down>", lambda e: (sugg.focus_set(), sugg.selection_clear(0, tk.END), sugg.selection_set(0), sugg.activate(0)) if sugg.winfo_ismapped() else None)
        sugg.bind("<Return>", _accept_suggestion)
        sugg.bind("<Double-Button-1>", _accept_suggestion)

        _refresh_suggestions()


    def _open_equilibrium_dg_steps(self):
        win = tk.Toplevel(self)
        win.title("ΔG° ↔ K — step-by-step")
        win.transient(self)
        win.grab_set()
        win.geometry("720x520")

        frm = ttk.Frame(win, padding=12)
        frm.pack(fill=tk.BOTH, expand=True)

        # Pull inputs from Known Variables (base units)
        if "ΔG°rxn" not in self.known_base:
            messagebox.showerror("Missing ΔG°rxn", "Enter ΔG°rxn first (J/mol or kJ/mol).", parent=win)
            return
        if "T" not in self.known_base:
            messagebox.showerror("Missing T", "Enter temperature T in K.", parent=win)
            return

        dG_base = float(self.known_base["ΔG°rxn"])          # J/mol (base)
        dG_disp = self.known_ui.get("ΔG°rxn", {"unit": "kJ/mol", "display_value": dG_base/1000.0})
        dG_kJ = dG_base / 1000.0

        T = float(self.known_base["T"])                      # K
        T_disp = self.known_ui.get("T", {"unit": "K", "display_value": T})

        # Gas constant: user-provided R or default
        R = float(self.known_base.get("R", CONSTANTS.get("R", 8.314462618)))
        R_disp = self.known_ui.get("R", {"unit": "J/(mol·K)", "display_value": R})

        # Compute
        x = - dG_base / (R * T)          # exponent
        K = math.exp(x)

        # Handy reference deltas at this T
        RT = R * T                        # J/mol
        rtln10_kJ = (RT * math.log(10.0)) / 1000.0
        rtln2_kJ  = (RT * math.log(2.0))  / 1000.0

        # UI: big read-only text area
        txt = tk.Text(frm, height=20, wrap="word")
        txt.pack(fill=tk.BOTH, expand=True)
        vs = ttk.Scrollbar(frm, orient="vertical", command=txt.yview)
        txt.configure(yscrollcommand=vs.set)
        vs.place(relx=1.0, rely=0, relheight=1.0, anchor="ne")

        def fmt(x):
            # compact scientific (up to 6 sig figs)
            return f"{x:.6g}"

        steps = []
        steps.append("Given (using Known Variables):\n")
        steps.append(f"  ΔG°rxn = {dG_disp['display_value']:.6g} {dG_disp['unit']}  (= {dG_base:.6g} J/mol)\n")
        steps.append(f"  T      = {T_disp['display_value']:.6g} {T_disp['unit']}\n")
        steps.append(f"  R      = {R_disp['display_value']:.6g} {R_disp['unit']}\n\n")
        steps.append("Compute exponent x = −ΔG°rxn / (R·T):\n")
        steps.append(f"  x = −({dG_base:.6g} J/mol) / ({R:.6g} J/(mol·K) × {T:.6g} K)\n")
        steps.append(f"    = {x:.6g}\n\n")
        steps.append("Then K = exp(x):\n")
        steps.append(f"  K = exp({x:.6g}) = {K:.6g}\n")
        steps.append(f"  → scientific: {K:.3e}\n\n")
        steps.append("Useful rule-of-thumb at this T:\n")
        steps.append(f"  RT·ln(10) = {rtln10_kJ:.3f} kJ/mol  (ΔG° change per 10× change in K)\n")
        steps.append(f"  RT·ln(2)  = {rtln2_kJ:.3f} kJ/mol   (ΔG° change per 2× change in K)\n")
        txt.insert("1.0", "".join(steps))
        txt.configure(state="disabled")

        # Buttons
        btns = ttk.Frame(frm)
        btns.pack(anchor="w", pady=(10,0))
        def add_K():
            # Save K as known, unitless
            self.known_base["K"] = float(K)
            self.known_ui["K"] = {"unit": "", "display_value": float(K)}
            self._refresh_known_table()
            self._refresh_equation_list()
            messagebox.showinfo("Added", "Added K to Known Variables.", parent=win)

        ttk.Button(btns, text="Add K to Known Variables", command=add_K).pack(side=tk.LEFT)

    def _open_photon_color_tool(self):
        win = tk.Toplevel(self)
        win.title("Photon color from energy / wavelength")
        win.transient(self)
        win.grab_set()
        win.geometry("560x560")  # taller to fit the step-by-step log

        frm = ttk.Frame(win, padding=12)
        frm.pack(fill=tk.BOTH, expand=True)

        ttk.Label(frm, text="Enter E or λ (leave the other empty):").grid(row=0, column=0, columnspan=4, sticky="w")

        ttk.Label(frm, text="E:").grid(row=1, column=0, sticky="e")
        e_val = ttk.Entry(frm, width=14)
        e_val.grid(row=1, column=1, sticky="w", padx=6)
        e_unit = ttk.Combobox(frm, width=8, state="readonly", values=["J", "eV"])
        e_unit.grid(row=1, column=2, sticky="w")
        e_unit.set("J")

        ttk.Label(frm, text="λ:").grid(row=2, column=0, sticky="e")
        lam_val = ttk.Entry(frm, width=14)
        lam_val.grid(row=2, column=1, sticky="w", padx=6)
        lam_unit = ttk.Combobox(frm, width=8, state="readonly", values=["m", "nm", "Å"])
        lam_unit.grid(row=2, column=2, sticky="w")
        lam_unit.set("nm")

        # Short formulas info (still useful)
        info_text = (
            "Formulas:\n"
            "• E = h·c / λ    • λ = h·c / E\n"
            "Photoelectric effect: KE = E_photon − φ (if E_photon > φ)."
        )
        ttk.Label(frm, text=info_text, justify="left", wraplength=520, foreground="gray25").grid(
            row=3, column=0, columnspan=4, sticky="w", pady=(8, 10)
        )

        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 10, "bold")).grid(
            row=4, column=0, columnspan=4, sticky="w", pady=(0, 6)
        )

        # Work function inputs
        ttk.Label(frm, text="Work function φ (kJ/mol):").grid(row=5, column=0, sticky="e")
        wf_entry = ttk.Entry(frm, width=14)
        wf_entry.grid(row=5, column=1, sticky="w", padx=6)
        wf_entry.insert(0, "221")  # Potassium default
        ttk.Label(frm, text="(photoelectric: E_photon > φ)").grid(row=5, column=2, columnspan=2, sticky="w")

        ke_out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=ke_out, font=("Segoe UI", 10, "bold")).grid(
            row=6, column=0, columnspan=4, sticky="w", pady=(4, 6)
        )

        # Step-by-step log (monospace label; easy and no extra imports)
        ttk.Label(frm, text="Calculation details:", font=("Segoe UI", 9, "bold")).grid(
            row=7, column=0, columnspan=4, sticky="w", pady=(6, 2)
        )
        steps_var = tk.StringVar(value="")
        steps_lbl = ttk.Label(
            frm, textvariable=steps_var, justify="left", wraplength=520, font=("Courier New", 9)
        )
        steps_lbl.grid(row=8, column=0, columnspan=4, sticky="w")

        def compute_ke(E_J, steps):
            """Compute KE from photon energy and work function field. Writes to ke_out and appends steps."""
            txt = wf_entry.get().strip()
            if not txt:
                ke_out.set("")  # no work function provided → clear
                return steps
            try:
                wf_kJmol = float(txt.replace(",", "."))
            except ValueError:
                ke_out.set("Invalid work function.")
                return steps + "\n[φ] Invalid work function; cannot compute KE."

            NA = 6.022_140_76e23  # Avogadro's number
            wf_J = wf_kJmol * 1000.0 / NA  # J per electron
            steps += (
                f"\n[Photoelectric]\n"
                f"  φ_given = {wf_kJmol:g} kJ/mol\n"
                f"  Convert φ → J/electron: φ_J = (φ_kJmol × 1000) / N_A\n"
                f"    = ({wf_kJmol:g} × 1000) / {NA:.6e} = {wf_J:.6e} J\n"
                f"  KE = E_photon − φ_J = {E_J:.6e} − {wf_J:.6e} = {E_J - wf_J:.6e} J\n"
            )

            KE = E_J - wf_J
            if KE <= 0:
                ke_out.set("Photon energy < work function: no electrons ejected.")
                steps += "  Result: E_photon ≤ φ → no electrons ejected.\n"
                return steps
            ke_out.set(f"Kinetic energy of ejected electron: {KE:.3e} J  ({KE/1.602_176_634e-19:.3g} eV)")
            steps += f"  Result: KE = {KE:.6e} J = {KE/1.602_176_634e-19:.6g} eV\n"
            return steps

        def compute():
            lam_txt = lam_val.get().strip()
            e_txt = e_val.get().strip()

            # constants
            h = float(self.known_base.get("h", 6.626_070_15e-34))
            c = float(self.known_base.get("c", 299_792_458.0))
            q_e = 1.602_176_634e-19

            eJ_entered = None
            lam_m_entered = None
            steps = ""

            # Parse inputs and log
            try:
                if e_txt:
                    raw_E = float(e_txt.replace(",", "."))
                    if e_unit.get() == "eV":
                        steps += (
                            "[Input E]\n"
                            f"  E_entered = {raw_E:g} eV\n"
                            f"  Convert eV → J: E_J = E_eV × e = {raw_E:g} × {q_e:.6e} = {raw_E*q_e:.6e} J\n"
                        )
                        eJ_entered = raw_E * q_e
                    else:
                        steps += f"[Input E]\n  E_entered = {raw_E:g} J\n"
                        eJ_entered = raw_E
                    if eJ_entered <= 0:
                        raise ValueError

                if lam_txt:
                    raw_lam = float(lam_txt.replace(",", "."))
                    u = lam_unit.get()
                    if u == "nm":
                        lam_m_entered = raw_lam * 1e-9
                        steps += (
                            "[Input λ]\n"
                            f"  λ_entered = {raw_lam:g} nm\n"
                            f"  Convert nm → m: λ_m = {raw_lam:g} × 1e-9 = {lam_m_entered:.6e} m\n"
                        )
                    elif u == "Å":
                        lam_m_entered = raw_lam * 1e-10
                        steps += (
                            "[Input λ]\n"
                            f"  λ_entered = {raw_lam:g} Å\n"
                            f"  Convert Å → m: λ_m = {raw_lam:g} × 1e-10 = {lam_m_entered:.6e} m\n"
                        )
                    else:
                        lam_m_entered = raw_lam
                        steps += f"[Input λ]\n  λ_entered = {raw_lam:g} m\n"
                    if lam_m_entered <= 0:
                        raise ValueError
            except ValueError:
                messagebox.showerror("Invalid", "Enter positive numeric values for E or λ.", parent=win)
                return

            steps += (
                "\n[Constants]\n"
                f"  h = {h:.8e} J·s\n"
                f"  c = {c:.6e} m/s\n"
            )

            # Compute missing quantity (and cross-check if both given)
            if lam_m_entered is not None:
                lam_m = lam_m_entered
                eJ = h * c / lam_m
                steps += (
                    "\n[Compute E from λ]\n"
                    f"  E = h·c / λ = ({h:.8e} × {c:.6e}) / {lam_m:.6e} = {eJ:.6e} J\n"
                    f"  E_eV = E / e = {eJ:.6e} / {q_e:.6e} = {eJ/q_e:.6g} eV\n"
                )
                if eJ_entered is not None:
                    lam_from_E = h * c / eJ_entered
                    rel_err = abs(lam_from_E - lam_m) / lam_m if lam_m else 0.0
                    steps += (
                        "\n[Consistency check]\n"
                        f"  λ_from_E = h·c / E_entered = ({h:.8e} × {c:.6e}) / {eJ_entered:.6e} = {lam_from_E:.6e} m\n"
                        f"  Relative mismatch = |λ_from_E − λ_entered| / λ_entered = {rel_err:.3%}\n"
                    )
                    if rel_err > 0.02:
                        messagebox.showwarning(
                            "Inconsistent inputs",
                            f"The entered E and λ disagree by {rel_err*100:.1f}%. Using λ to compute E.",
                            parent=win,
                        )
            elif eJ_entered is not None:
                eJ = eJ_entered
                lam_m = h * c / eJ
                steps += (
                    "\n[Compute λ from E]\n"
                    f"  λ = h·c / E = ({h:.8e} × {c:.6e}) / {eJ:.6e} = {lam_m:.6e} m\n"
                    f"  λ_nm = λ × 1e9 = {lam_m*1e9:.6g} nm\n"
                )
            else:
                messagebox.showerror("Missing input", "Enter E or λ.", parent=win)
                return

            # Convert m → nm (for display)
            lam_nm = lam_m * 1e9

            # Color label
            def nm_to_color(nm: float) -> str:
                if nm < 380 or nm > 780: return "outside visible range"
                if nm < 450: return "violet"
                if nm < 495: return "blue"
                if nm < 570: return "green"
                if nm < 590: return "yellow"
                if nm < 620: return "orange"
                return "red"

            color = nm_to_color(lam_nm)
            out.set(
                f"λ = {lam_nm:.1f} nm → {color};  "
                f"E = {eJ/q_e:.3g} eV  ({eJ:.3e} J)"
            )

            # Photoelectric effect (optional, with steps)
            steps = compute_ke(eJ, steps)

            # Final color note
            steps += f"\n[Color]\n  λ_nm = {lam_nm:.3f} → {color}\n"

            # Push the step-by-step log to UI
            steps_var.set(steps)

            return lam_m, eJ

        self._add_compute_bar(win, compute)


        #----------------------------------------------------------------------------------------------#

    def _open_phase_change_tool(self):
        win = tk.Toplevel(self)
        win.title("Heating/Cooling with Phase Changes (catalog-backed)")
        win.transient(self); win.grab_set(); win.geometry("820x640")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # ---------- Substance selector ----------
        ttk.Label(frm, text="Substance:").grid(row=0, column=0, sticky="e")
        sub = tk.StringVar(value="water")
        sub_box = ttk.Combobox(frm, textvariable=sub, width=28, values=sorted(self.phase_catalog.keys()))
        sub_box.grid(row=0, column=1, columnspan=2, sticky="w", padx=6)
        ttk.Label(frm, text="(You can type a new name and Save)").grid(row=0, column=3, sticky="w")

        # Mass / temps
        ttk.Label(frm, text="Mass:").grid(row=1, column=0, sticky="e")
        mass_val = ttk.Entry(frm, width=10); mass_val.grid(row=1, column=1, sticky="w", padx=6); mass_val.insert(0, "1.00")
        mass_unit = ttk.Combobox(frm, width=6, state="readonly", values=["kg","g"]); mass_unit.grid(row=1, column=2, sticky="w"); mass_unit.set("kg")

        ttk.Label(frm, text="Initial T:").grid(row=2, column=0, sticky="e")
        Ti_val = ttk.Entry(frm, width=10); Ti_val.grid(row=2, column=1, sticky="w", padx=6); Ti_val.insert(0, "-20")
        Ti_unit = ttk.Combobox(frm, width=6, state="readonly", values=["°C","K"]); Ti_unit.grid(row=2, column=2, sticky="w"); Ti_unit.set("°C")

        ttk.Label(frm, text="Final T:").grid(row=3, column=0, sticky="e")
        Tf_val = ttk.Entry(frm, width=10); Tf_val.grid(row=3, column=1, sticky="w", padx=6); Tf_val.insert(0, "80")
        Tf_unit = ttk.Combobox(frm, width=6, state="readonly", values=["°C","K"]); Tf_unit.grid(row=3, column=2, sticky="w"); Tf_unit.set("°C")

        # ---------- Constants panel ----------
        box = ttk.LabelFrame(frm, text="Constants (editable; saved per substance)")
        box.grid(row=4, column=0, columnspan=4, sticky="nsew", pady=(8,8))
        for i in range(6): box.grid_columnconfigure(i, weight=1)

        # ui helpers
        def mk_field(r, label, unit, width=12):
            ttk.Label(box, text=label).grid(row=r, column=0, sticky="w")
            e = ttk.Entry(box, width=width); e.grid(row=r, column=1, sticky="w", padx=6)
            ttk.Label(box, text=unit).grid(row=r, column=2, sticky="w")
            return e

        f_formula = mk_field(0, "Formula (optional)", "")
        f_c_solid = mk_field(1, "c_solid", "J/(kg·K)")
        f_c_liq   = mk_field(2, "c_liquid", "J/(kg·K)")
        f_c_gas   = mk_field(3, "c_gas", "J/(kg·K)")
        f_H_fus   = mk_field(4, "ΔH_fus", "J/kg")
        f_H_vap   = mk_field(5, "ΔH_vap", "J/kg")
        f_Tm      = mk_field(6, "T_melt", "°C")
        f_Tb      = mk_field(7, "T_boil", "°C")

        # Output
        txt = tk.Text(frm, height=14, wrap="word")
        txt.grid(row=5, column=0, columnspan=4, sticky="nsew", pady=(6,4))
        frm.grid_rowconfigure(5, weight=1)
        scy = ttk.Scrollbar(frm, orient="vertical", command=txt.yview); scy.grid(row=5, column=4, sticky="ns")
        txt.configure(yscrollcommand=scy.set)
        total = tk.StringVar(value="")
        ttk.Label(frm, textvariable=total, font=("Segoe UI",10,"bold")).grid(row=6, column=0, columnspan=4, sticky="w")

        # ----- helpers -----
        def load_current():
            name = sub.get().strip()
            rec = self.phase_catalog.get(name) or DEFAULT_PHASE_CATALOG.get(name) or {}
            f_formula.delete(0,"end"); f_formula.insert(0, rec.get("formula",""))
            for e, key, default in [
                (f_c_solid, "c_solid", 2090.0),
                (f_c_liq,   "c_liquid",4184.0),
                (f_c_gas,   "c_gas",   1996.0),
                (f_H_fus,   "H_fus",   333_550.0),
                (f_H_vap,   "H_vap",   2_256_000.0),
                (f_Tm,      "T_melt_C",0.0),
                (f_Tb,      "T_boil_C",100.0),
            ]:
                e.delete(0,"end"); e.insert(0, str(rec.get(key, default)))

        def current_record_from_fields():
            def num(e):
                return float(e.get().strip().replace(",", "."))
            return {
                "formula": f_formula.get().strip(),
                "c_solid": num(f_c_solid),
                "c_liquid": num(f_c_liq),
                "c_gas": num(f_c_gas),
                "H_fus": num(f_H_fus),
                "H_vap": num(f_H_vap),
                "T_melt_C": num(f_Tm),
                "T_boil_C": num(f_Tb),
            }

        def to_C(txt_val, unit_combo):
            t = float(txt_val.replace(",", "."))
            return t if unit_combo.get() == "°C" else (t - 273.15)

        def add_line(s):
            txt.insert("end", s + "\n")

        # ----- compute (same physics you already had) -----
        def compute():
            txt.configure(state="normal"); txt.delete("1.0","end")
            try:
                m_in = float(mass_val.get().replace(",", "."))
                m_kg = convert_to_base("m", m_in, mass_unit.get())
                Ti_C = to_C(Ti_val.get(), Ti_unit)
                Tf_C = to_C(Tf_val.get(), Tf_unit)
                rec  = current_record_from_fields()
            except ValueError:
                messagebox.showerror("Invalid input", "Enter numeric values.", parent=win); return

            if m_kg <= 0:
                messagebox.showerror("Mass error","Mass must be > 0.", parent=win); return

            c_s, c_l, c_g = rec["c_solid"], rec["c_liquid"], rec["c_gas"]
            Hfus, Hvap    = rec["H_fus"], rec["H_vap"]
            Tm, Tb        = rec["T_melt_C"], rec["T_boil_C"]

            def phase(Tc):
                return "solid" if Tc < Tm else ("liquid" if Tc < Tb else "gas")

            Q = 0.0
            T = Ti_C
            add_line(f"{sub.get()}  (formula: {rec.get('formula') or '—'})")
            add_line(f"m = {m_kg:.6g} kg;  path: {Ti_C:.2f} °C → {Tf_C:.2f} °C")
            add_line("Steps:")

            def heat_to(target, c, label):
                nonlocal T, Q
                if abs(target - T) < 1e-12: return
                dT = target - T
                q = m_kg * c * dT
                Q += q
                add_line(f"  • {label}: ΔT = {dT:+.2f} K → q = {q/1000:.3f} kJ")
                T = target

            def phase_change(latent, label, sign=+1):
                nonlocal Q
                q = sign * m_kg * latent
                Q += q
                add_line(f"  • {label}: q = {q/1000:.3f} kJ")

            # march upward or downward
            direction = 1 if Tf_C >= Ti_C else -1
            while True:
                ph = phase(T)
                if direction > 0:
                    if ph == "solid":
                        heat_to(min(Tm, Tf_C), c_s, "heat solid")
                        if T < Tf_C: phase_change(Hfus, "melt @ T_melt")
                    elif ph == "liquid":
                        heat_to(min(Tb, Tf_C), c_l, "heat liquid")
                        if T < Tf_C: phase_change(Hvap, "boil @ T_boil")
                    else:
                        heat_to(Tf_C, c_g, "heat gas")
                else:
                    if ph == "gas":
                        heat_to(max(Tb, Tf_C), c_g, "cool gas")
                        if T > Tf_C: phase_change(Hvap, "condense @ T_boil", sign=-1)
                    elif ph == "liquid":
                        heat_to(max(Tm, Tf_C), c_l, "cool liquid")
                        if T > Tf_C: phase_change(Hfus, "freeze @ T_melt", sign=-1)
                    else:
                        heat_to(Tf_C, c_s, "cool solid")
                if abs(T - Tf_C) < 1e-12: break

            txt.configure(state="disabled")
            total.set(f"Total Q = {Q/1000:.3f} kJ  (= {Q:.3e} J)")
            return Q

        self._add_compute_bar(win, compute)

        # ----- Save / Reset / Compute buttons -----
        btns = ttk.Frame(frm); btns.grid(row=7, column=0, columnspan=4, sticky="w", pady=(6,0))
        def save_consts():
            name = sub.get().strip()
            if not name:
                messagebox.showerror("Name required","Type a substance name to save.", parent=win); return
            try:
                rec = current_record_from_fields()
            except ValueError:
                messagebox.showerror("Bad number","Check constants fields.", parent=win); return
            self.phase_catalog[name] = rec
            save_phase_catalog(self.phase_catalog, self.phase_catalog_path)
            sub_box.configure(values=sorted(self.phase_catalog.keys()))
            messagebox.showinfo("Saved", f"Saved constants for '{name}' to phase_constants.json", parent=win)

        def reset_to_file():
            load_current()
            messagebox.showinfo("Reloaded", "Reloaded values for this substance.", parent=win)

        ttk.Button(btns, text="Compute", command=compute).pack(side=tk.LEFT, padx=(0,8))
        ttk.Button(btns, text="Save constants to catalog", command=save_consts).pack(side=tk.LEFT, padx=(0,8))
        ttk.Button(btns, text="Reset fields from catalog", command=reset_to_file).pack(side=tk.LEFT)

        # initial populate + change handler
        def on_sub_change(_=None): load_current()
        sub_box.bind("<<ComboboxSelected>>", on_sub_change)
        load_current()


    def _open_percent_mixer_tool(self):
        """
        Concentration & % Converter
        - Pick % type (w/w, w/v, v/v) and enter %.
        - Give solute formula (e.g., HCl) and (optionally) solvent formula (default H2O).
        - Optionally give densities: solution rho (g/mL) and solute rho (g/mL).
        - Computes: mass/volume basis, m (molality), M (molarity, when possible), mass %,
        mole fraction (if solvent formula given), and explanatory details.
        """
        import re
        try:
            from constants import ATOMIC_WEIGHTS  # your project file with atom masses
        except Exception:
            # Minimal fallback (extend if you don't import constants)
            ATOMIC_WEIGHTS = {"H":1.00794,"C":12.0107,"N":14.0067,"O":15.9994,"Cl":35.453,"Na":22.98976928}

        def _parse_num(s, i):
            m = re.match(r"\d+", s[i:])
            return (int(m.group(0)), i + len(m.group(0))) if m else (1, i)

        def _parse_formula(s, i=0):
            comp = {}
            while i < len(s):
                ch = s[i]
                if ch == "(":
                    sub, i = _parse_formula(s, i + 1)
                    mult, i = _parse_num(s, i)
                    for el, n in sub.items():
                        comp[el] = comp.get(el, 0) + n * mult
                elif ch == ")":
                    return comp, i + 1
                else:
                    m = re.match(r"[A-Z][a-z]?", s[i:])
                    if not m:
                        break
                    el = m.group(0)
                    i += len(el)
                    mult, i = _parse_num(s, i)
                    comp[el] = comp.get(el, 0) + mult
            return comp, i

        def molar_mass(formula: str) -> float:
            comp, _ = _parse_formula(formula.strip())
            return sum(ATOMIC_WEIGHTS[el] * n for el, n in comp.items())

        win = tk.Toplevel(self); win.title("Concentration & % Converter")
        win.transient(self); win.grab_set(); win.geometry("680x520")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # --- Inputs ---------------------------------------------------------------
        kind = tk.StringVar(value="w/w")
        ttk.Label(frm, text="Given as:").grid(row=0, column=0, sticky="e")
        ttk.Combobox(frm, textvariable=kind, state="readonly", width=8,
                    values=["w/w","w/v","v/v"]).grid(row=0, column=1, sticky="w")

        ttk.Label(frm, text="% value:").grid(row=1, column=0, sticky="e")
        pct_e = ttk.Entry(frm, width=10); pct_e.grid(row=1, column=1, sticky="w"); pct_e.insert(0, "30")

        ttk.Label(frm, text="Solute formula:").grid(row=2, column=0, sticky="e")
        solute_e = ttk.Entry(frm, width=12); solute_e.grid(row=2, column=1, sticky="w"); solute_e.insert(0, "HCl")

        ttk.Label(frm, text="Solvent formula:").grid(row=3, column=0, sticky="e")
        solvent_e = ttk.Entry(frm, width=12); solvent_e.grid(row=3, column=1, sticky="w"); solvent_e.insert(0, "H2O")

        ttk.Label(frm, text="ρ_solution (g/mL) [optional]:").grid(row=0, column=3, sticky="e")
        rho_soln_e = ttk.Entry(frm, width=10); rho_soln_e.grid(row=0, column=4, sticky="w")  # e.g., 1.15 for strong HCl

        ttk.Label(frm, text="ρ_solute (g/mL) [optional]:").grid(row=1, column=3, sticky="e")
        rho_solute_e = ttk.Entry(frm, width=10); rho_solute_e.grid(row=1, column=4, sticky="w")

        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 10, "bold")).grid(
            row=4, column=0, columnspan=5, sticky="w", pady=(8, 4)
        )
        txt = tk.Text(frm, height=18, wrap="word")
        txt.grid(row=5, column=0, columnspan=5, sticky="nsew"); frm.grid_rowconfigure(5, weight=1)

        # --- Compute --------------------------------------------------------------
        def compute():
            txt.delete("1.0", "end")
            try:
                pct = float(pct_e.get().replace(",", ".")) / 100.0
            except ValueError:
                messagebox.showerror("Invalid", "Enter a numeric percent.", parent=win); return

            solute = solute_e.get().strip()
            solvent = (solvent_e.get().strip() or "H2O")
            try:
                MM_solute = molar_mass(solute)
            except Exception:
                messagebox.showerror("Unknown element", f"Cannot parse formula: {solute}", parent=win); return

            MM_solvent = None
            try:
                MM_solvent = molar_mass(solvent)
            except Exception:
                pass  # fine; only needed for mole fraction

            rho_soln = float(rho_soln_e.get().replace(",", ".")) if rho_soln_e.get() else None
            rho_solute = float(rho_solute_e.get().replace(",", ".")) if rho_solute_e.get() else None

            # Canonical bases:
            mass_soln_g = None
            vol_soln_mL = None
            mass_solute_g = None

            # Results we can compute
            m_molality = None
            M_molarity = None
            x_solute = None

            if kind.get() == "w/w":
                # 100 g solution basis
                mass_soln_g = 100.0
                mass_solute_g = pct * mass_soln_g
                mass_solvent_g = mass_soln_g - mass_solute_g
                n_solute = mass_solute_g / MM_solute
                m_molality = n_solute / (mass_solvent_g / 1000.0)  # mol/kg
                if rho_soln:
                    vol_soln_mL = mass_soln_g / rho_soln
                    M_molarity = n_solute / (vol_soln_mL / 1000.0)  # mol/L
                if MM_solvent is not None:
                    n_solvent = mass_solvent_g / MM_solvent
                    x_solute = n_solute / (n_solute + n_solvent)

            elif kind.get() == "w/v":
                # By definition: (% g per 100 mL solution)
                vol_soln_mL = 100.0
                mass_solute_g = pct * 100.0
                n_solute = mass_solute_g / MM_solute
                M_molarity = n_solute / 0.1  # mol/L
                if rho_soln:
                    mass_soln_g = rho_soln * vol_soln_mL
                    mass_solvent_g = mass_soln_g - mass_solute_g
                    m_molality = n_solute / (mass_solvent_g / 1000.0)
                    if MM_solvent is not None and mass_solvent_g > 0:
                        n_solvent = mass_solvent_g / MM_solvent
                        x_solute = n_solute / (n_solute + n_solvent)

            else:  # v/v
                # (% mL solute per 100 mL solution). Need ρ_solute to get mass-based things.
                vol_soln_mL = 100.0
                vol_solute_mL = pct * 100.0
                if rho_solute:
                    mass_solute_g = rho_solute * vol_solute_mL
                    n_solute = mass_solute_g / MM_solute
                    M_molarity = n_solute / 0.1
                    if rho_soln:
                        mass_soln_g = rho_soln * vol_soln_mL
                        mass_solvent_g = mass_soln_g - mass_solute_g
                        if mass_solvent_g > 0:
                            m_molality = n_solute / (mass_solvent_g / 1000.0)
                            if MM_solvent is not None:
                                n_solvent = mass_solvent_g / MM_solvent
                                x_solute = n_solute / (n_solute + n_solvent)

            # --- Output -----------------------------------------------------------
            lines = []
            if mass_soln_g is not None:
                lines.append(f"Basis: solution mass = {mass_soln_g:.3g} g")
            if vol_soln_mL is not None:
                lines.append(f"Basis: solution volume = {vol_soln_mL:.3g} mL")
            if mass_solute_g is not None:
                lines.append(f"Solute mass = {mass_solute_g:.3g} g  (Mᵣ={MM_solute:.3f} g/mol)")

            if m_molality is not None:
                lines.append(f"Molality m = {m_molality:.4g} mol/kg")
            else:
                lines.append("Molality m = (needs mass of solvent ⇒ either w/w, or w/v with ρ_solution, or v/v with both densities)")

            if M_molarity is not None:
                lines.append(f"Molarity M = {M_molarity:.4g} mol/L")
            else:
                lines.append("Molarity M = (needs solution volume ⇒ w/v, or w/w with ρ_solution, or v/v plus densities)")

            # Mass percent is always whatever you entered, but re-state clearly:
            if kind.get() == "w/w":
                lines.append(f"Mass percent (w/w) = {pct*100:.4g}%")
            elif kind.get() == "w/v":
                lines.append(f"Mass/volume percent (w/v) = {pct*100:.4g}%  (={pct*100:.4g} g per 100 mL)")
            else:
                lines.append(f"Volume percent (v/v) = {pct*100:.4g}%  ({pct*100:.4g} mL per 100 mL)")

            if x_solute is not None:
                lines.append(f"Mole fraction x_solute = {x_solute:.4f}  (solvent={solvent})")

            out.set("Computed:")
            txt.insert("end", "\n".join(lines) + "\n")

        self._add_compute_bar(win, compute)




    def _open_henry_tool(self):
        win = tk.Toplevel(self); win.title("Henry’s Law"); win.transient(self); win.grab_set(); win.geometry("560x360")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        mode = tk.StringVar(value="px")
        ttk.Radiobutton(frm, text="p = kH·x", variable=mode, value="px").grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(frm, text="c = kH·p", variable=mode, value="cp").grid(row=0, column=1, sticky="w")

        p_e   = ttk.Entry(frm, width=12); x_e   = ttk.Entry(frm, width=12); k_e = ttk.Entry(frm, width=12)
        ttk.Label(frm, text="p (Pa / kPa / atm):").grid(row=1,column=0,sticky="e"); p_e.grid(row=1,column=1,sticky="w")
        ttk.Label(frm, text="x (mole fraction):").grid(row=2,column=0,sticky="e"); x_e.grid(row=2,column=1,sticky="w")
        ttk.Label(frm, text="kH:").grid(row=3,column=0,sticky="e");             k_e.grid(row=3,column=1,sticky="w")
        out = tk.StringVar(value=""); ttk.Label(frm, textvariable=out, font=("Segoe UI",10,"bold")).grid(row=4,column=0,columnspan=3,sticky="w",pady=(8,0))

        def compute():
            try:
                p  = float(p_e.get().replace(",", ".")) if p_e.get() else None
                x  = float(x_e.get().replace(",", ".")) if x_e.get() else None
                kH = float(k_e.get().replace(",", ".")) if k_e.get() else None
            except ValueError:
                messagebox.showerror("Invalid","Numbers only.", parent=win); return

            if mode.get()=="px":
                if p is None: p = kH * x
                elif x is None: x = p / kH
                elif kH is None: kH = p / x
                out.set(f"p = {p:.4g} ;  x = {x:.4g} ;  kH = {kH:.4g} (Pa)")
            else:
                # c = kH * p  (c as mol/L if kH entered in mol/(L·atm) and p in atm)
                if p is None or kH is None:
                    if p is None and kH is not None and x is not None:
                        # allow quick convert: if x given, approximate c ≈ x·(n/V) via ideal gas? keep simple: require c or p
                        pass
                c = float(x_e.get().replace(",", ".")) if (mode.get()=="cp" and x_e.get()) else None  # repurpose x_e as c field
                # if user put c in that box:
                if c is None and p is not None and kH is not None: c = kH * p
                elif p is None and c is not None and kH is not None: p = c / kH
                elif kH is None and p is not None and c is not None: kH = c / p
                out.set(f"c = {c:.4g} ;  p = {p:.4g} ;  kH = {kH:.4g} (units per your entry)")
        
        self._add_compute_bar(win, compute)


    def _open_raoult_tool(self):
        win = tk.Toplevel(self); win.title("Raoult’s Law — 2 Components"); win.transient(self); win.grab_set(); win.geometry("640x430")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # Component A
        ttk.Label(frm, text="x_A:").grid(row=0,column=0,sticky="e"); xA=ttk.Entry(frm,width=10); xA.grid(row=0,column=1,sticky="w")
        ttk.Label(frm, text="P*_A:").grid(row=0,column=2,sticky="e"); PsA=ttk.Entry(frm,width=10); PsA.grid(row=0,column=3,sticky="w"); ttk.Label(frm,text="(same units for A/B)").grid(row=0,column=4,sticky="w")
        # Component B
        ttk.Label(frm, text="x_B:").grid(row=1,column=0,sticky="e"); xB=ttk.Entry(frm,width=10); xB.grid(row=1,column=1,sticky="w")
        ttk.Label(frm, text="P*_B:").grid(row=1,column=2,sticky="e"); PsB=ttk.Entry(frm,width=10); PsB.grid(row=1,column=3,sticky="w")

        out = tk.StringVar(value=""); ttk.Label(frm, textvariable=out, font=("Segoe UI",10,"bold")).grid(row=2,column=0,columnspan=5,sticky="w",pady=(8,0))
        txt = tk.Text(frm, height=12, wrap="word"); txt.grid(row=3,column=0,columnspan=5,sticky="nsew"); frm.grid_rowconfigure(3, weight=1)

        def compute():
            txt.delete("1.0","end")
            try:
                x_A = float(xA.get().replace(",", ".")) if xA.get() else None
                x_B = float(xB.get().replace(",", ".")) if xB.get() else None
                P_A = float(PsA.get().replace(",", "."))
                P_B = float(PsB.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid","Enter numeric values.", parent=win); return
            if x_A is None and x_B is None: messagebox.showerror("Need a composition","Provide x_A or x_B.", parent=win); return
            if x_A is None: x_A = 1.0 - x_B
            if x_B is None: x_B = 1.0 - x_A
            pA = x_A * P_A; pB = x_B * P_B; Ptot = pA + pB
            yA = pA / Ptot; yB = pB / Ptot
            out.set(f"P_total = {Ptot:.4g}   (p_A={pA:.4g}, p_B={pB:.4g})")
            txt.insert("end", f"Gas composition: y_A = {yA:.4f}, y_B = {yB:.4f}\n")
            return
        
        self._add_compute_bar(win, compute)


    def _open_ion_vant_hoff_tool(self):
        win = tk.Toplevel(self); win.title("Ion Concentrations & van ’t Hoff"); win.transient(self); win.grab_set(); win.geometry("720x520")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # minimal dissociation dictionary
        DISS = {
            "NaCl": {"Na+":1, "Cl-":1},
            "KNO3":{"K+":1, "NO3-":1},
            "K2SO4":{"K+":2, "SO4^2-":1},
            "CaCl2":{"Ca^2+":1, "Cl-":2},
            "Al2(SO4)3":{"Al^3+":2, "SO4^2-":3},
            "K3PO4":{"K+":3, "PO4^3-":1},
            "MgCl2":{"Mg^2+":1, "Cl-":2},
        }

        ttk.Label(frm, text="Solute formula:").grid(row=0,column=0,sticky="e"); fE=ttk.Entry(frm,width=18); fE.grid(row=0,column=1,sticky="w"); fE.insert(0,"K2SO4")
        ttk.Label(frm, text="Formal molarity M:").grid(row=1,column=0,sticky="e"); ME=ttk.Entry(frm,width=10); ME.grid(row=1,column=1,sticky="w"); ME.insert(0,"0.15")
        ttk.Label(frm, text="van ’t Hoff factor i (optional):").grid(row=2,column=0,sticky="e"); iE=ttk.Entry(frm,width=10); iE.grid(row=2,column=1,sticky="w")
        out = tk.StringVar(value=""); ttk.Label(frm,textvariable=out, font=("Segoe UI",10,"bold")).grid(row=3,column=0,columnspan=3,sticky="w",pady=(8,0))
        txt = tk.Text(frm, height=16, wrap="word"); txt.grid(row=4,column=0,columnspan=3,sticky="nsew"); frm.grid_rowconfigure(4, weight=1)

        def compute():
            txt.delete("1.0","end")
            formula = fE.get().strip()
            try: M = float(ME.get().replace(",", ".")); 
            except ValueError: messagebox.showerror("Invalid","Give a numeric molarity.", parent=win); return
            st = DISS.get(formula)
            if not st:
                txt.insert("end","Not in internal dissociation map. Add it there if you use it often.\n"); return
            ideal_i = sum(st.values())
            concs = {ion: coeff*M for ion,coeff in st.items()}
            total_ideal = sum(concs.values())
            out.set(f"Ideal i = {ideal_i} ;  total ideal ion conc = {total_ideal:.3g} mol/L")
            txt.insert("end","Per-ion concentrations (ideal):\n")
            for ion, c in concs.items():
                txt.insert("end", f"  [{ion}] = {c:.4g} M\n")
            if iE.get():
                try: iobs = float(iE.get().replace(",", ".")); 
                except ValueError: iobs = None
                if iobs:
                    total_obs = iobs * M
                    txt.insert("end", f"\nWith observed i = {iobs}: total ion conc = {total_obs:.4g} M\n")
            return
        
        self._add_compute_bar(win, compute)


    def _open_colligative_tool(self):
        import tkinter as tk
        from tkinter import ttk, messagebox
        from molar_masses import molar_mass, name_to_formula  # you already import these above

        win = tk.Toplevel(self)
        win.title("Colligative Properties (ΔTb, ΔTf)")
        win.transient(self); win.grab_set(); win.geometry("780x700")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # ---- Solvent presets (auto-fill Kb/Kf and pure Tb/Tf) -------------------
        KBKF = {
            # Aqueous standard
            "water (H2O)":            {"Kb": 0.512, "Kf": 1.86, "Tb_C": 100.0,  "Tf_C": 0.0},

            # Aromatics
            "benzene (C6H6)":         {"Kb": 2.53,  "Kf": 5.12, "Tb_C": 80.1,   "Tf_C": 5.5},
            "toluene (C7H8)":         {"Kb": 3.40,  "Kf": 7.60, "Tb_C": 110.6,  "Tf_C": -95.0},
            "nitrobenzene (C6H5NO2)": {"Kb": 5.57,  "Kf": 7.00, "Tb_C": 210.9,  "Tf_C": 5.7},

            # Alcohols / glycols
            "ethanol (C2H5OH)":       {"Kb": 1.22,  "Kf": 1.99, "Tb_C": 78.37,  "Tf_C": -114.1},
            "methanol (CH3OH)":       {"Kb": 0.78,  "Kf": 1.90, "Tb_C": 64.7,   "Tf_C": -97.6},
            "ethylene glycol (C2H6O2)":{"Kb": 2.56, "Kf": 3.90, "Tb_C": 197.3,  "Tf_C": -12.9},  # Kf commonly tabulated ~3.9
            "glycerol (C3H8O3)":      {"Kb": 3.73,  "Kf": 5.10, "Tb_C": 290.0,  "Tf_C": 18.2},

            # Halogenated solvents
            "chloroform (CHCl3)":     {"Kb": 3.63,  "Kf": 4.68, "Tb_C": 61.2,   "Tf_C": -63.5},


            "carbon tetrachloride (CCl4)":{"Kb": 5.03,"Kf": 29.8,"Tb_C": 76.7,  "Tf_C": -22.9},
            "dichloromethane (CH2Cl2)":{"Kb": 3.38, "Kf": 4.90, "Tb_C": 39.6,   "Tf_C": -95.0},

            # Ketones / esters / nitriles / amides
            "acetone (C3H6O)":        {"Kb": 1.71,  "Kf": 2.39, "Tb_C": 56.1,   "Tf_C": -94.7},
            "ethyl acetate (C4H8O2)": {"Kb": 2.77,  "Kf": 4.50, "Tb_C": 77.1,   "Tf_C": -83.6},
            "acetonitrile (C2H3N)":   {"Kb": 0.74,  "Kf": 3.90, "Tb_C": 81.6,   "Tf_C": -45.7},
            "DMF (C3H7NO)":           {"Kb": 2.37,  "Kf": 3.70, "Tb_C": 153.0,  "Tf_C": -60.0},
            "DMSO (C2H6OS)":          {"Kb": 2.85,  "Kf": 4.20, "Tb_C": 189.0,  "Tf_C": 18.5},

            # Ethers / cycloalkanes
            "diethyl ether (C4H10O)": {"Kb": 2.02,  "Kf": 1.79, "Tb_C": 34.6,   "Tf_C": -116.3},
            "cyclohexane (C6H12)":    {"Kb": 2.79,  "Kf": 20.0, "Tb_C": 80.7,   "Tf_C": 6.5},

            # Carboxylic acids / others
            "acetic acid (CH3COOH)":  {"Kb": 2.93,  "Kf": 3.90, "Tb_C": 118.1,  "Tf_C": 16.6},

            # Solids used for cryoscopy (very large Kf; Tb often irrelevant)
            "camphor (C10H16O)":      {"Kb": None,  "Kf": 40.0, "Tb_C": 204.0,  "Tf_C": 179.8},  # classic cryoscopic solvent
            "naphthalene (C10H8)":    {"Kb": None,  "Kf": 6.9,  "Tb_C": 218.0,  "Tf_C": 80.2},
        }


        ttk.Label(frm, text="Solvent preset:").grid(row=0, column=0, sticky="e")
        solvent = tk.StringVar(value="water (H2O)")
        solvent_box = ttk.Combobox(frm, width=28, state="readonly",
                                values=sorted(KBKF.keys()), textvariable=solvent)
        solvent_box.grid(row=0, column=1, sticky="w", padx=(6,12))

        # ---- Main inputs --------------------------------------------------------
        ttk.Label(frm, text="Molality m (mol/kg):").grid(row=1, column=0, sticky="e")
        mE = ttk.Entry(frm, width=12); mE.grid(row=1, column=1, sticky="w")
        ttk.Label(frm, text="van ’t Hoff i:").grid(row=2, column=0, sticky="e")
        iE = ttk.Entry(frm, width=12); iE.grid(row=2, column=1, sticky="w"); iE.insert(0, "1")

        # A small hint label under i
        hint_i = tk.StringVar(value="Tip: choose a solute to auto-fill i")
        ttk.Label(frm, textvariable=hint_i, foreground="#555").grid(row=3, column=0, columnspan=2, sticky="w", pady=(2,6))

        # --- Solute type (large searchable picker) --------------------------------
        ttk.Label(frm, text="Solute type:").grid(row=2, column=2, sticky="e")
        choose_i_btn = ttk.Button(frm, text="Choose… (Ctrl+L)")
        choose_i_btn.grid(row=2, column=3, sticky="w", padx=(6,0))

        SOLUTE_I_PRESETS = [
            # --- Non-electrolytes (i ≈ 1) ------------------------------------------
            ("Glucose (non-electrolyte)", "C6H12O6", 1.0, "no dissociation; i≈1"),
            ("Sucrose (non-electrolyte)", "C12H22O11", 1.0, "no dissociation; i≈1"),
            ("Fructose (non-electrolyte)", "C6H12O6", 1.0, "no dissociation; i≈1"),
            ("Urea (non-electrolyte)", "CH4N2O", 1.0, "no dissociation; i≈1"),
            ("Ethanol (non-electrolyte)", "C2H5OH", 1.0, "no dissociation; i≈1"),
            ("Methanol (non-electrolyte)", "CH3OH", 1.0, "no dissociation; i≈1"),
            ("Glycerol (non-electrolyte)", "C3H8O3", 1.0, "no dissociation; i≈1"),
            ("Ethylene glycol (non-electrolyte)", "C2H6O2", 1.0, "no dissociation; i≈1"),
            ("Propylene glycol (non-electrolyte)", "C3H8O2", 1.0, "no dissociation; i≈1"),

            # --- Strong acids/bases (nominally i≈2; activity < 2 in practice) -------
            ("HCl (strong acid)",  "HCl",   2.0, ""),
            ("HBr (strong acid)",  "HBr",   2.0, ""),
            ("HI (strong acid)",   "HI",    2.0, ""),
            ("HNO3 (strong acid)", "HNO3",  2.0, ""),
            ("HClO4 (strong acid)","HClO4", 2.0, ""),
            ("NaOH (strong base)", "NaOH",  2.0, ""),
            ("KOH (strong base)",  "KOH",   2.0, ""),
            ("LiOH (strong base)", "LiOH",  2.0, ""),

            # --- Weak acids/bases (i a little > 1; depends on α) --------------------
            ("HF (weak acid)",        "HF",      1.0, "i slightly >1; strongly solvent-dependent"),
            ("CH3COOH (acetic acid)", "CH3COOH", 1.0, "i slightly >1; depends on dissociation"),
            ("H2CO3 (carbonic acid)", "H2CO3",   1.0, "polyprotic weak acid; i>1 but small in dilute soln"),
            ("H3PO4 (phosphoric acid)","H3PO4",  1.0, "triprotic weak; i between 1 and 4 depending on pH"),
            ("NH3(aq) (weak base)",   "NH3",     1.0, "forms NH4+ + OH−; i>1 but small in dilute soln"),
            ("NH4OH (aq ammonia)",    "NH4OH",   1.0, "convenience formula; weak base"),

            # --- 1:1 salts (→ 2 ions; i≈2) ------------------------------------------
            ("NaCl", "NaCl", 2.0, ""), ("KCl", "KCl", 2.0, ""), ("LiCl", "LiCl", 2.0, ""),
            ("NaBr", "NaBr", 2.0, ""), ("KBr", "KBr", 2.0, ""), ("NaI", "NaI", 2.0, ""), ("KI", "KI", 2.0, ""),
            ("NaF",  "NaF",  2.0, ""), ("KF",  "KF",  2.0, ""),
            ("NaNO3","NaNO3",2.0, ""), ("KNO3","KNO3",2.0, ""), ("NH4NO3","NH4NO3",2.0, ""),
            ("NaClO4","NaClO4",2.0,""), ("KClO4","KClO4",2.0,""),
            ("NaOAc (sodium acetate)","CH3COONa", 2.0, "salt of weak acid; still dissociates to 2 ions"),
            ("NH4Cl", "NH4Cl", 2.0, ""), ("NaHCO3","NaHCO3",2.0,""), ("NaHSO3","NaHSO3",2.0,""),

            # --- 2:1 or 1:2 salts (→ 3 ions; i≈3) -----------------------------------
            ("MgCl2", "MgCl2", 3.0, ""), ("CaCl2", "CaCl2", 3.0, ""), ("SrCl2", "SrCl2", 3.0, ""), ("BaCl2","BaCl2",3.0,""),
            ("ZnCl2","ZnCl2",3.0,""), ("FeCl2","FeCl2",3.0,""), ("CuCl2","CuCl2",3.0,""), ("Pb(NO3)2","Pb(NO3)2",3.0,""),
            ("Na2SO4","Na2SO4",3.0,""), ("K2SO4","K2SO4",3.0,""), ("(NH4)2SO4","(NH4)2SO4",3.0,""),
            ("Na2CO3","Na2CO3",3.0,""), ("K2CO3","K2CO3",3.0,""), ("(NH4)2CO3","(NH4)2CO3",3.0,""),
            ("Ca(NO3)2","Ca(NO3)2",3.0,""), ("Mg(NO3)2","Mg(NO3)2",3.0,""), ("Ba(NO3)2","Ba(NO3)2",3.0,""),
            ("Ba(OH)2","Ba(OH)2",3.0,"strong base; limited solubility"),
            ("Ca(OH)2","Ca(OH)2",3.0,"strong base; limited solubility"),

            # --- 3:1 or 1:3 salts (→ 4 ions; i≈4) -----------------------------------
            ("AlCl3","AlCl3",4.0,"hydrolyzes somewhat in water (real i < 4)"),
            ("FeCl3","FeCl3",4.0,"hydrolysis lowers effective i"),
            ("CrCl3","CrCl3",4.0,""),
            ("Al(NO3)3","Al(NO3)3",4.0,""),
            ("Fe(NO3)3","Fe(NO3)3",4.0,""),
            ("Na3PO4","Na3PO4",4.0,"basic solution; phosphate equilibria in real systems"),

            # --- 3:2 salts (→ 5 ions; i≈5) ------------------------------------------
            ("Al2(SO4)3","Al2(SO4)3",5.0,""),
            ("Fe2(SO4)3","Fe2(SO4)3",5.0,""),
            ("Cr2(SO4)3","Cr2(SO4)3",5.0,""),

            # --- Complex ion salts (fully dissociate to counterions + complex) -------
            ("K3[Fe(CN)6]","K3Fe(CN)6",4.0,"3 K+ + [Fe(CN)6]3− → 4 particles"),
            ("K4[Fe(CN)6]","K4Fe(CN)6",5.0,"4 K+ + [Fe(CN)6]4− → 5 particles"),

            # --- Polyprotic strong (upper bounds) ------------------------------------
            ("H2SO4 (upper-bound ideal)","H2SO4",3.0,"dilute realistic 2<i<3 (2nd dissociation incomplete)"),

            # --- “Real world” caveats ------------------------------------------------
            ("CuSO4 (association/hydration)","CuSO4",2.0,"hydrates; ion pairing can lower i"),
            ("AgNO3 (complexation in some media)","AgNO3",2.0,"can form complexes; i<2 in some conditions"),
            ("CaSO4 (very low solubility)","CaSO4",2.0,"sparingly soluble; colligative effects negligible"),
        ]


        # ---- Kb/Kf and pure Tb/Tf fields ----------------------------------------
        ttk.Label(frm, text="K_b (K·kg/mol):").grid(row=4, column=0, sticky="e")
        KbE = ttk.Entry(frm, width=12); KbE.grid(row=4, column=1, sticky="w")
        ttk.Label(frm, text="K_f (K·kg/mol):").grid(row=4, column=2, sticky="e")
        KfE = ttk.Entry(frm, width=12); KfE.grid(row=4, column=3, sticky="w")

        ttk.Label(frm, text="Pure Tb (°C):").grid(row=5, column=0, sticky="e")
        Tb0E = ttk.Entry(frm, width=12); Tb0E.grid(row=5, column=1, sticky="w")
        ttk.Label(frm, text="Pure Tf (°C):").grid(row=5, column=2, sticky="e")
        Tf0E = ttk.Entry(frm, width=12); Tf0E.grid(row=5, column=3, sticky="w")

        # ---- Quick molality calculator -----------------------------------------
        box = ttk.LabelFrame(frm, text="Quick molality calculator (optional)")
        box.grid(row=6, column=0, columnspan=4, sticky="we", pady=(10,6))
        for c in range(4): box.grid_columnconfigure(c, weight=1)

        ttk.Label(box, text="Solute (name/formula):").grid(row=0, column=0, sticky="e")
        solute_txt = ttk.Entry(box, width=20); solute_txt.grid(row=0, column=1, sticky="w", padx=(6,8))
        ttk.Label(box, text="Mr (g/mol):").grid(row=0, column=2, sticky="e")
        MrE = ttk.Entry(box, width=12); MrE.grid(row=0, column=3, sticky="w")

        mm_info = tk.StringVar(value="")
        ttk.Label(box, textvariable=mm_info, foreground="#555").grid(row=1, column=0, columnspan=4, sticky="w", pady=(0,4))

        ttk.Label(box, text="Solute mass (g):").grid(row=2, column=0, sticky="e")
        msE = ttk.Entry(box, width=12); msE.grid(row=2, column=1, sticky="w", padx=(6,8))
        ttk.Label(box, text="Solvent mass (g):").grid(row=2, column=2, sticky="e")
        msolvE = ttk.Entry(box, width=12); msolvE.grid(row=2, column=3, sticky="w")

        ttk.Label(box, text="or Solvent volume (mL) × ρ (g/mL):").grid(row=3, column=0, sticky="e")
        vsolvE = ttk.Entry(box, width=12); vsolvE.grid(row=3, column=1, sticky="w", padx=(6,8))
        rhoE = ttk.Entry(box, width=12); rhoE.grid(row=3, column=2, sticky="w")
        rhoE.insert(0, "1.0")  # water default

        def refresh_Mr_from_solute(_=None):
            txt = solute_txt.get().strip()
            if not txt:
                mm_info.set(""); return
            try:
                f = name_to_formula(txt)
                Mr = molar_mass(f)
                MrE.delete(0, "end"); MrE.insert(0, f"{Mr:g}")
                mm_info.set(f"Formula: {f} — Mr = {Mr:g} g/mol")
            except Exception:
                mm_info.set("Unrecognized name/formula. Enter Mr manually.")

        solute_txt.bind("<KeyRelease>", refresh_Mr_from_solute)

        def compute_molality():
            # Mr from auto or manual
            Mr = None
            if MrE.get():
                try: Mr = float(MrE.get().replace(",", "."))
                except ValueError:
                    messagebox.showerror("Invalid Mr", "Enter a numeric molar mass.", parent=win); return
            else:
                refresh_Mr_from_solute()
                if MrE.get():
                    Mr = float(MrE.get().replace(",", "."))
            if Mr is None:
                messagebox.showerror("Need molar mass", "Type a formula/name I can parse or Mr (g/mol).", parent=win); return

            try:
                ms = float(msE.get().replace(",", "."))  # g
            except ValueError:
                messagebox.showerror("Invalid mass", "Enter solute mass (g).", parent=win); return

            # solvent mass
            if msolvE.get():
                try:
                    msolv_g = float(msolvE.get().replace(",", "."))
                except ValueError:
                    messagebox.showerror("Invalid solvent mass", "Enter solvent mass (g).", parent=win); return
            elif vsolvE.get():
                try:
                    V = float(vsolvE.get().replace(",", "."))
                    rho = float(rhoE.get().replace(",", "."))
                    msolv_g = V * rho
                except ValueError:
                    messagebox.showerror("Invalid volume/density", "Enter numeric volume and density.", parent=win); return
            else:
                messagebox.showerror("Need solvent amount","Give solvent mass or volume×density.", parent=win); return

            if msolv_g <= 0:
                messagebox.showerror("Invalid solvent", "Solvent mass must be > 0.", parent=win); return

            n = ms / Mr                     # mol
            m = n / (msolv_g / 1000.0)      # mol/kg
            mE.delete(0, "end"); mE.insert(0, f"{m:g}")
            messagebox.showinfo("Molality set", f"m = {m:g} mol/kg", parent=win)

        ttk.Button(box, text="Compute m →", command=compute_molality).grid(row=4, column=0, columnspan=4, sticky="w", pady=(6,0))

        # ---- Prefill from Known Variables (if saved earlier) --------------------
        if not KbE.get():
            for key in ("Kb","K_b","K_boil"):
                if key in self.known_ui:
                    try: KbE.insert(0, f"{float(self.known_ui[key]['display_value']):g}"); break
                    except: pass
        if not KfE.get():
            for key in ("Kf","K_f","K_freeze"):
                if key in self.known_ui:
                    try: KfE.insert(0, f"{float(self.known_ui[key]['display_value']):g}"); break
                    except: pass

        # ---- Fill solvent fields (Kb/Kf/Tb/Tf) ----------------------------------
        def _fill_from_solvent(force=False, *_):
            d = KBKF.get(solvent.get())
            if not d: return
            def put(entry, value):
                if force or not entry.get():
                    entry.delete(0, "end"); entry.insert(0, str(value))
            put(KbE,  d["Kb"])
            put(KfE,  d["Kf"])
            put(Tb0E, d["Tb_C"])
            put(Tf0E, d["Tf_C"])

        solvent_box.bind("<<ComboboxSelected>>", lambda e: _fill_from_solvent(True))
        _fill_from_solvent(False)  # initial prefill only where empty

        # ---- Solute i picker dialog ---------------------------------------------
        def _open_i_picker():
            pop = tk.Toplevel(win)
            pop.title("Choose solute type (van ’t Hoff i)")
            pop.transient(win); pop.grab_set(); pop.geometry("720x560")

            top = ttk.Frame(pop, padding=10); top.pack(fill=tk.BOTH, expand=True)
            ttk.Label(top, text="Search (name or formula):").pack(anchor="w")
            q = ttk.Entry(top, width=36); q.pack(fill=tk.X, pady=(2,8))

            lst = tk.Listbox(top, height=18); lst.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
            vs = ttk.Scrollbar(top, orient="vertical", command=lst.yview)
            lst.configure(yscrollcommand=vs.set); vs.pack(side=tk.RIGHT, fill=tk.Y, padx=(6,0))

            info = tk.StringVar(value="Double-click or press Enter to choose. Esc to cancel.")
            ttk.Label(pop, textvariable=info, foreground="#555").pack(anchor="w", padx=10, pady=6)

            state = {"rows": SOLUTE_I_PRESETS[:]}

            def refresh(_=None):
                needle = q.get().strip().lower()
                rows = [r for r in SOLUTE_I_PRESETS
                        if needle in r[0].lower() or needle in r[1].lower()]
                state["rows"] = rows
                lst.delete(0, tk.END)
                for name, formula, i_val, note in rows:
                    tail = f" — {note}" if note else ""
                    lst.insert(tk.END, f"{name}  [{formula}]  →  i={i_val:g}{tail}")
                if rows:
                    lst.selection_clear(0, tk.END); lst.selection_set(0); lst.activate(0)

            def use_selected(_=None):
                if not state["rows"]: return
                sel = lst.curselection()
                if not sel: return
                name, formula, i_val, note = state["rows"][sel[0]]
                iE.delete(0, "end"); iE.insert(0, f"{i_val:g}")
                hint_i.set(f"{name}: i={i_val:g}" + (f" — {note}" if note else ""))
                # also drop formula → Mr for Quick-m
                solute_txt.delete(0, "end"); solute_txt.insert(0, formula)
                refresh_Mr_from_solute()
                pop.destroy()

            def cancel(_=None): pop.destroy()

            q.bind("<KeyRelease>", refresh)
            lst.bind("<Double-Button-1>", use_selected)
            lst.bind("<Return>", use_selected)
            pop.bind("<Return>", use_selected)
            pop.bind("<Escape>", cancel)

            refresh(); q.focus_set()

        choose_i_btn.configure(command=_open_i_picker)
        win.bind("<Control-l>", lambda e: _open_i_picker())

        # ---- Output / explanation ----------------------------------------------
        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 10, "bold")).grid(
            row=7, column=0, columnspan=4, sticky="w", pady=(10,4)
        )

        expl = tk.Text(frm, height=12, wrap="word")
        expl.grid(row=8, column=0, columnspan=4, sticky="nsew")
        frm.grid_rowconfigure(8, weight=1)

        def compute():
            expl.configure(state="normal"); expl.delete("1.0", "end")
            try:
                m  = float(mE.get().replace(",", "."))
                i  = float(iE.get().replace(",", "."))
                Kb = float(KbE.get().replace(",", ".")) if KbE.get() else 0.0
                Kf = float(KfE.get().replace(",", ".")) if KfE.get() else 0.0
            except ValueError:
                messagebox.showerror("Invalid", "Numbers only in m, i, Kb, Kf.", parent=win); return

            Tb0 = float(Tb0E.get().replace(",", ".")) if Tb0E.get() else None
            Tf0 = float(Tf0E.get().replace(",", ".")) if Tf0E.get() else None

            dTb = i*Kb*m if Kb>0 else None
            dTf = i*Kf*m if Kf>0 else None

            lines = []
            lines.append("Using ΔTb = i·Kb·m and ΔTf = i·Kf·m\n")
            lines.append(f"Inputs:  m={m:g} mol/kg,  i={i:g},  Kb={Kb:g} K·kg/mol,  Kf={Kf:g} K·kg/mol\n")
            if dTb is not None:
                lines.append(f"ΔTb = {i:g} × {Kb:g} × {m:g} = {dTb:.6g} K\n")
                if Tb0 is not None:
                    lines.append(f"Predicted boiling point: Tb(soln) = {Tb0:g} + {dTb:.6g} = {Tb0 + dTb:.6g} °C\n")
            if dTf is not None:
                lines.append(f"ΔTf = {i:g} × {Kf:g} × {m:g} = {dTf:.6g} K\n")
                if Tf0 is not None:
                    lines.append(f"Predicted freezing point: Tf(soln) = {Tf0:g} − {dTf:.6g} = {Tf0 - dTf:.6g} °C\n")

            if not (dTb or dTf):
                out.set("Enter Kb and/or Kf.")
            else:
                msg = []
                if dTb is not None:
                    msg.append(f"ΔT_b = {dTb:.4g} K")
                    if Tb0 is not None: msg.append(f"Tb ≈ {Tb0 + dTb:.4g} °C")
                if dTf is not None:
                    msg.append(f"ΔT_f = {dTf:.4g} K")
                    if Tf0 is not None: msg.append(f"Tf ≈ {Tf0 - dTf:.4g} °C")
                out.set(" ;  ".join(msg))

            expl.insert("end", "".join(lines))
            expl.configure(state="disabled")

        # ✅ Put the compute bar on the FRAME (not the toplevel) so it always shows
        self._add_compute_bar(win, compute)




    def _open_osmotic_tool(self):
        import tkinter as tk
        from tkinter import ttk, messagebox
        from molar_masses import molar_mass, name_to_formula, suggest_substances

        win = tk.Toplevel(self)
        win.title("Osmotic Pressure — π = i·M·R·T")
        win.transient(self); win.grab_set(); win.geometry("760x560")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # ---------------- Core inputs ----------------
        ttk.Label(frm, text="Molarity M (mol/L):").grid(row=0, column=0, sticky="e")
        ME = ttk.Entry(frm, width=12); ME.grid(row=0, column=1, sticky="w", padx=(6,14))

        ttk.Label(frm, text="van ’t Hoff i:").grid(row=0, column=2, sticky="e")
        iE = ttk.Entry(frm, width=10); iE.grid(row=0, column=3, sticky="w", padx=(6,14)); iE.insert(0, "1")

        ttk.Label(frm, text="Temperature T:").grid(row=1, column=0, sticky="e", pady=(4,0))
        TE = ttk.Entry(frm, width=12); TE.grid(row=1, column=1, sticky="w", padx=(6,6), pady=(4,0)); TE.insert(0, "298.15")
        Tunit = ttk.Combobox(frm, width=6, state="readonly", values=["K","°C"])
        Tunit.grid(row=1, column=2, sticky="w", pady=(4,0)); Tunit.set("K")

        # ---------------- Quick M calculator (optional) ----------------
        box = ttk.LabelFrame(frm, text="Quick M from mass & volume (optional)")
        box.grid(row=2, column=0, columnspan=4, sticky="nsew", pady=(10,8))
        for c in range(6): box.grid_columnconfigure(c, weight=1)

        ttk.Label(box, text="Solute (name/symbol/formula):").grid(row=0, column=0, sticky="w")
        f_entry = ttk.Entry(box, width=28); f_entry.grid(row=0, column=1, columnspan=3, sticky="we", padx=6)
        f_entry.insert(0, "sucrose")

        # suggestions
        sugg = tk.Listbox(box, height=0); sugg.grid(row=1, column=0, columnspan=4, sticky="we")
        scx = ttk.Scrollbar(box, orient="horizontal", command=sugg.xview)
        sugg.configure(xscrollcommand=scx.set); scx.grid(row=2, column=0, columnspan=4, sticky="we", pady=(0,6))
        sugg.grid_remove(); scx.grid_remove()

        # mass, volume
        ttk.Label(box, text="Solute mass:").grid(row=3, column=0, sticky="e")
        massE = ttk.Entry(box, width=10); massE.grid(row=3, column=1, sticky="w", padx=6); massE.insert(0, "15")
        massU = ttk.Combobox(box, width=6, state="readonly", values=["g","mg","kg"]); massU.grid(row=3, column=2, sticky="w"); massU.set("g")

        ttk.Label(box, text="Solution volume:").grid(row=3, column=3, sticky="e")
        volE = ttk.Entry(box, width=10); volE.grid(row=3, column=4, sticky="w", padx=6); volE.insert(0, "100")
        volU = ttk.Combobox(box, width=6, state="readonly", values=["mL","L"]); volU.grid(row=3, column=5, sticky="w"); volU.set("mL")

        # preview
        p_formula = tk.StringVar(value="")
        p_mm = tk.StringVar(value="")
        ttk.Label(box, textvariable=p_formula, font=("Segoe UI", 9, "bold")).grid(row=4, column=0, columnspan=6, sticky="w", pady=(4,0))
        ttk.Label(box, textvariable=p_mm, foreground="#444").grid(row=5, column=0, columnspan=6, sticky="w")

        def _refresh_suggestions(_=None):
            items = suggest_substances(f_entry.get(), limit=60)
            sugg.delete(0, tk.END)
            for label, _f in items: sugg.insert(tk.END, label)
            if items:
                sugg.configure(height=min(8, len(items))); sugg.grid(); scx.grid()
            else:
                sugg.grid_remove(); scx.grid_remove()
            _update_preview()

        def _accept_suggestion(_=None):
            sel = sugg.curselection()
            if not sel: return
            name = sugg.get(sel[0]).split(" — ", 1)[0]
            f_entry.delete(0, tk.END); f_entry.insert(0, name)
            sugg.grid_remove(); scx.grid_remove()
            _update_preview()

        def _update_preview():
            f = name_to_formula(f_entry.get().strip())
            p_formula.set(f"Formula: {f}")
            try:
                mm = molar_mass(f)
                p_mm.set(f"Molar mass = {mm:g} g/mol")
            except Exception:
                p_mm.set("Molar mass = (unknown)")

        f_entry.bind("<KeyRelease>", _refresh_suggestions)
        f_entry.bind("<Down>", lambda e: (sugg.focus_set(), sugg.selection_clear(0, tk.END), sugg.selection_set(0), sugg.activate(0)) if sugg.winfo_ismapped() else None)
        sugg.bind("<Return>", _accept_suggestion)
        sugg.bind("<Double-Button-1>", _accept_suggestion)
        _update_preview()

        # Dissociation map for ideal i (extend as you like)
        DISS = {
            "NaCl": {"Na+":1,"Cl-":1},
            "KNO3":{"K+":1,"NO3-":1},
            "K2SO4":{"K+":2,"SO4^2-":1},
            "CaCl2":{"Ca^2+":1,"Cl-":2},
            "Al2(SO4)3":{"Al^3+":2,"SO4^2-":3},
            "MgCl2":{"Mg^2+":1,"Cl-":2},
            "HCl": {"H+":1, "Cl-":1},
            "H2SO4": {"H+":2, "SO4^2-":1},
            # add more
        }

        hint_i = tk.StringVar(value="Non-electrolyte → i ≈ 1 (e.g., sugars)")
        ttk.Label(box, textvariable=hint_i, foreground="#555").grid(row=6, column=0, columnspan=4, sticky="w", pady=(2,0))

        def compute_M():
            # mass to moles
            try:
                mass = float(massE.get().replace(",", "."))
                vol = float(volE.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid", "Mass/volume must be numbers.", parent=win); return None

            # mass units
            if massU.get() == "mg": mass_g = mass * 1e-3
            elif massU.get() == "kg": mass_g = mass * 1e3
            else: mass_g = mass

            # volume units
            V_L = vol / 1000.0 if volU.get()=="mL" else vol

            try:
                f = name_to_formula(f_entry.get().strip())
                mm = molar_mass(f)  # g/mol
            except Exception as ex:
                messagebox.showerror("Unknown solute", f"Can't get molar mass: {ex}", parent=win); return None

            if V_L <= 0:
                messagebox.showerror("Invalid volume", "Solution volume must be > 0.", parent=win); return None

            n = mass_g / mm  # mol
            M = n / V_L      # mol/L
            ME.delete(0, "end"); ME.insert(0, f"{M:g}")

            # i helper
            ideal = DISS.get(f, None)
            if ideal:
                i_val = float(sum(ideal.values()))
                iE.delete(0, "end"); iE.insert(0, f"{i_val:g}")
                hint_i.set(f"Ideal i for {f}: {i_val:g}  ({' + '.join([f'{k}:{v}' for k,v in ideal.items()])})")
            else:
                hint_i.set("Not in electrolyte list ⇒ assume non-electrolyte (i ≈ 1)")

            return M

        ttk.Button(box, text="Compute M →", command=compute_M).grid(row=3, column=6, sticky="w", padx=(10,0))

        # ---------------- Output & explanation ----------------
        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 10, "bold")).grid(row=3, column=0, columnspan=4, sticky="w", pady=(8,2))

        steps = tk.Text(frm, height=12, wrap="word")
        steps.grid(row=4, column=0, columnspan=4, sticky="nsew")
        frm.grid_rowconfigure(4, weight=1)
        scy = ttk.Scrollbar(frm, orient="vertical", command=steps.yview)
        steps.configure(yscrollcommand=scy.set)
        scy.grid(row=4, column=4, sticky="ns")

        def compute():
            steps.configure(state="normal"); steps.delete("1.0", "end")

            # Parse core entries
            try:
                M = float(ME.get().replace(",", "."))
                i = float(iE.get().replace(",", "."))
                T_in = float(TE.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid", "M, i, and T must be numeric.", parent=win); return

            T_K = T_in if Tunit.get()=="K" else (T_in + 273.15)
            R = float(self.known_base.get("R", 8.314462618))  # J/(mol·K)
            M_m3 = M * 1000.0                                  # mol/m^3

            # π in Pa
            pi = i * M_m3 * R * T_K

            # Output summary
            out.set(f"π = {pi/101325:.5g} atm  = {pi/1e5:.5g} bar  = {pi/1e3:.5g} kPa  (= {pi:.5g} Pa)")

            # Steps text
            def g(x): return f"{x:.6g}"
            s = []
            s.append("Given / using:\n")
            s.append(f"  M = {g(M)} mol/L  → {g(M_m3)} mol·m⁻³\n")
            s.append(f"  i = {g(i)}\n")
            s.append(f"  T = {g(T_in)} {Tunit.get()}  → {g(T_K)} K\n")
            s.append(f"  R = {g(R)} J·mol⁻¹·K⁻¹\n\n")
            s.append("Compute π = i·M·R·T (SI):\n")
            s.append(f"  π = {g(i)} × {g(M_m3)} × {g(R)} × {g(T_K)} = {g(pi)} Pa\n")
            s.append(f"  = {g(pi/1e3)} kPa  = {g(pi/1e5)} bar  = {g(pi/101325)} atm\n")
            steps.insert("1.0", "".join(s))
            steps.configure(state="disabled")

            # Save to Known (base = Pa; UI kPa)
            self.known_base["π_osm"] = float(pi)
            self.known_ui["π_osm"] = {"unit": "kPa", "display_value": float(pi/1e3)}
            self._refresh_known_table(); self._refresh_equation_list()

            return pi

        # Attach the bottom compute bar to the FRAME (so it appears)
        self._add_compute_bar(win, compute)


    def _open_kb_kf_lookup(self):
        win = tk.Toplevel(self)
        win.title("Kb / Kf — common solvents")
        win.transient(self); win.grab_set(); win.geometry("560x420")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # small embedded table (values are typical textbook values)
        KBKF = {
            "water (H2O)":        {"Kb": 0.512, "Kf": 1.86, "Tb_C": 100.0, "Tf_C": 0.0},
            "benzene (C6H6)":     {"Kb": 2.53,  "Kf": 5.12, "Tb_C": 80.1,  "Tf_C": 5.5},
            "ethanol (C2H5OH)":   {"Kb": 1.22,  "Kf": 1.99, "Tb_C": 78.37, "Tf_C": -114.1},
            "chloroform (CHCl3)": {"Kb": 3.63,  "Kf": 4.68, "Tb_C": 61.2,  "Tf_C": -63.5},
            "acetic acid (CH3COOH)": {"Kb": 2.93, "Kf": 3.90, "Tb_C": 118.1, "Tf_C": 16.6},
        }

        # --- UI -----------------------------------------------------------------
        ttk.Label(frm, text="Filter:").grid(row=0, column=0, sticky="e")
        q = tk.StringVar()
        ent = ttk.Entry(frm, textvariable=q, width=28); ent.grid(row=0, column=1, sticky="w", padx=(6,8))

        lst = tk.Listbox(frm, height=8); lst.grid(row=1, column=0, columnspan=4, sticky="nsew", pady=(6,6))
        frm.grid_rowconfigure(1, weight=1); frm.grid_columnconfigure(3, weight=1)

        info = tk.StringVar(value="—")
        ttk.Label(frm, textvariable=info, font=("Segoe UI",10,"bold")).grid(row=2, column=0, columnspan=4, sticky="w", pady=(4,2))

        note = ttk.Label(frm, text="Units: K·kg/mol (molal). Values are approximate; use your course table if it differs.",
                        foreground="#555", wraplength=520, justify="left")
        note.grid(row=3, column=0, columnspan=4, sticky="w")

        # populate list
        names_all = sorted(KBKF.keys())
        def refresh():
            pat = q.get().strip().lower()
            lst.delete(0, tk.END)
            for name in names_all:
                if pat in name.lower():
                    lst.insert(tk.END, name)
            if lst.size():
                lst.selection_clear(0, tk.END); lst.selection_set(0); lst.activate(0); show()

        def current():
            sel = lst.curselection()
            if not sel: return None, None
            name = lst.get(sel[0])
            return name, KBKF[name]

        def show(_=None):
            name, data = current()
            if not data:
                info.set("—"); return
            info.set(f"{name} →  Kb = {data['Kb']}  ;  Kf = {data['Kf']}   "
                    f"(Tb ≈ {data['Tb_C']} °C, Tf ≈ {data['Tf_C']} °C)")

        ent.bind("<KeyRelease>", lambda e: refresh())
        lst.bind("<<ListboxSelect>>", show)
        refresh()

        # actions
        btns = ttk.Frame(frm); btns.grid(row=4, column=0, columnspan=4, sticky="w", pady=(10,0))
        def add_kb():
            name, d = current()
            if not d: return
            self.known_base["Kb"] = float(d["Kb"])              # base already K·kg/mol
            self.known_ui["Kb"]   = {"unit":"K·kg/mol","display_value": float(d["Kb"])}
            self._refresh_known_table(); self._refresh_equation_list()
            messagebox.showinfo("Added", f"Added Kb = {d['Kb']} K·kg/mol", parent=win)

        def add_kf():
            name, d = current()
            if not d: return
            self.known_base["Kf"] = float(d["Kf"])
            self.known_ui["Kf"]   = {"unit":"K·kg/mol","display_value": float(d["Kf"])}
            self._refresh_known_table(); self._refresh_equation_list()
            messagebox.showinfo("Added", f"Added Kf = {d['Kf']} K·kg/mol", parent=win)

        ttk.Button(btns, text="Add Kb to Known", command=add_kb).pack(side=tk.LEFT, padx=(0,8))
        ttk.Button(btns, text="Add Kf to Known", command=add_kf).pack(side=tk.LEFT)

        # a tiny no-op compute to keep UI consistent (shows current line again)
        def compute():
            show()
            return current()

        # put the bar LAST so it can’t be overlapped
        self._add_compute_bar(frm, compute, text="Show selection")


    def _open_rate_law_tool(self):
        import math
        win = tk.Toplevel(self)
        win.title("Rate laws — zero / first / second order")
        win.transient(self); win.grab_set(); win.geometry("860x640")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # --- Top: reaction & order ------------------------------------------------
        ttk.Label(frm, text="Reaction (display only):").grid(row=0, column=0, sticky="e")
        rxn = ttk.Entry(frm, width=36); rxn.grid(row=0, column=1, columnspan=3, sticky="we", padx=(6,12))
        rxn.insert(0, "A → products")

        ttk.Label(frm, text="Order n:").grid(row=1, column=0, sticky="e")
        order_var = tk.StringVar(value="1 (first-order)")
        order_box = ttk.Combobox(frm, width=20, state="readonly",
                                values=["0 (zero-order)", "1 (first-order)", "2 (second-order)"],
                                textvariable=order_var)
        order_box.grid(row=1, column=1, sticky="w", padx=(6,12))

        # k value + units (units depend on order)
        ttk.Label(frm, text="Rate constant k:").grid(row=1, column=2, sticky="e")
        k_entry = ttk.Entry(frm, width=14); k_entry.grid(row=1, column=3, sticky="w", padx=(6,6))
        k_unit = ttk.Combobox(frm, width=12, state="readonly"); k_unit.grid(row=1, column=4, sticky="w")
        k_entry.insert(0, "1.58e-3")  # handy default for first-order example

        # Prefill k if user has one in Known Variables
        for key in ("k_rate", "k", "k1"):
            if key in self.known_ui:
                try:
                    k_entry.delete(0, "end")
                    k_entry.insert(0, f"{float(self.known_ui[key]['display_value']):g}")
                except Exception:
                    pass
                break

        # --- Inputs block ---------------------------------------------------------
        ttk.Label(frm, text="[A]₀ (M):").grid(row=2, column=0, sticky="e")
        A0_entry = ttk.Entry(frm, width=12); A0_entry.grid(row=2, column=1, sticky="w", padx=(6,12))
        A0_entry.insert(0, "0.400")  # from the screenshot example

        ttk.Label(frm, text="Time t:").grid(row=2, column=2, sticky="e")
        t_entry = ttk.Entry(frm, width=12); t_entry.grid(row=2, column=3, sticky="w", padx=(6,6))
        t_unit = ttk.Combobox(frm, width=8, state="readonly",
                            values=["s", "min", "h"]); t_unit.grid(row=2, column=4, sticky="w")
        t_unit.set("s"); t_entry.insert(0, "180")  # example

        ttk.Label(frm, text="Target [A] or fraction:").grid(row=3, column=0, sticky="e")
        At_entry = ttk.Entry(frm, width=12); At_entry.grid(row=3, column=1, sticky="w", padx=(6,12))
        ttk.Label(frm, text="or % consumed:").grid(row=3, column=2, sticky="e")
        pct_entry = ttk.Entry(frm, width=12); pct_entry.grid(row=3, column=3, sticky="w", padx=(6,6))
        ttk.Label(frm, text="(e.g. 75)").grid(row=3, column=4, sticky="w")

        
                # --- Help popup ---------------------------------------------------------
        def _open_help_kinetics():
            h = tk.Toplevel(win)
            h.title("Kinetics – Help")
            h.transient(win); h.grab_set(); h.geometry("760x620")

            wrap = ttk.Frame(h, padding=12); wrap.pack(fill=tk.BOTH, expand=True)
            txt = tk.Text(wrap, wrap="word")
            txt.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
            vs = ttk.Scrollbar(wrap, orient="vertical", command=txt.yview)
            txt.configure(yscrollcommand=vs.set); vs.pack(side=tk.RIGHT, fill=tk.Y)

            help_text = (
                "Kinetics tool — what the controls do\n"
                "=====================================\n\n"
                "Order selector (0th / 1st / 2nd)\n"
                "• This chooses the integrated rate law used by all the buttons.\n"
                "• It also sets the correct units for k and the half-life formula.\n\n"
                "Cheat sheet:\n"
                "  0th order:   rate = k\n"
                "               [A] = [A]0 − k t\n"
                "               t1/2 = [A]0 / (2k)        (depends on [A]0)\n"
                "               units(k) = M·time⁻¹       linear plot: [A] vs t (slope −k)\n\n"
                "  1st order:   rate = k[A]\n"
                "               ln[A] = ln[A]0 − k t\n"
                "               t1/2 = ln 2 / k           (constant)\n"
                "               units(k) = time⁻¹         linear plot: ln[A] vs t (slope −k)\n\n"
                "  2nd order:   rate = k[A]²\n"
                "               1/[A] = 1/[A]0 + k t\n"
                "               t1/2 = 1 / (k[A]0)        (depends on [A]0)\n"
                "               units(k) = M⁻¹·time⁻¹     linear plot: 1/[A] vs t (slope +k)\n\n"
                "Buttons (what they do)\n"
                "• Predict [A](t): uses the integrated law for the chosen order to compute concentration at time t.\n"
                "• Rate at time t: returns instantaneous rate = k·[A]^n using [A](t) from the integrated law.\n"
                "• Time for fraction/percent: solves for t given a target fraction (e.g., 75% consumed).\n"
                "• k from half-life: computes k from a provided half-life (needs [A]0 for 0th/2nd order).\n"
                "• Solve k from two points: uses two (t, [A]) data points to estimate k for the chosen order.\n"
                "• Add to Known: pushes results (k, [A]t, t, etc.) to the Known Variables table for reuse elsewhere.\n\n"
                "Unit tips\n"
                "• Keep time consistent (s, min, h). The tool won’t change k units for you—choose the right unit for k.\n"
                "• Concentration is in M unless you specify otherwise; conversions are not automatic here.\n\n"
                "Worked mini-example (first-order):\n"
                "  C2H5I → C2H4 + HI,\n"
                "  k = 1.58×10⁻³ s⁻¹,  [A]0 = 0.400 M,  t = 180 s\n"
                "  ln[A] = ln[A]0 − k t → [A] = 0.400·exp(−1.58×10⁻³·180) ≈ 0.322 M\n"
                "  Rate at 180 s = k·[A] ≈ (1.58×10⁻³ s⁻¹)·(0.322 M) ≈ 5.1×10⁻⁴ M·s⁻¹\n"
                "  Time to 75% consumed (i.e., [A]/[A]0 = 0.25): t = (ln 4)/k ≈ 0.693·2 / k ≈ 877 s.\n\n"
                "Choosing the order quickly\n"
                "• Constant half-life → 1st order.  [A] vs t linear → 0th.  1/[A] vs t linear → 2nd.\n"
                "• In mixed/excess cases, you may have pseudo-first-order: k_obs = k·[excess]^m.\n"
            )
            txt.insert("1.0", help_text)
            txt.configure(state="disabled")

            btnrow = ttk.Frame(wrap); btnrow.pack(fill=tk.X, pady=(8,0))
            ttk.Button(btnrow, text="Close (Esc)", command=h.destroy).pack(side=tk.RIGHT)
            h.bind("<Escape>", lambda e: h.destroy())

        # --- Task chooser ---------------------------------------------------------
        ttk.Label(frm, text="Task:").grid(row=4, column=0, sticky="e")
        task = tk.StringVar(value="A_t")
        task_box = ttk.Combobox(
            frm, width=36, state="readonly", textvariable=task,
            values=[
                "A_t  • Find [A] at time t",
                "rate • Find rate at time t",
                "t_%  • Time for X% consumed",
                "t_A  • Time to reach target [A]",
                "k_2p • Solve k from two points",
                "t12  • Half-life (t₁/₂) from k and [A]₀",
                "k_t12 • Solve k from half-life",
            ])
        task_box.grid(row=4, column=1, columnspan=3, sticky="w", padx=(6,12))

        ttk.Label(frm, text="Half-life t₁/₂:").grid(row=5, column=0, sticky="e")
        t12_entry = ttk.Entry(frm, width=12); t12_entry.grid(row=5, column=1, sticky="w", padx=(6,12))
        t12_unit  = ttk.Combobox(frm, width=8, state="readonly", values=["s","min","h"])
        t12_unit.grid(row=5, column=2, sticky="w"); t12_unit.set("s")

        # --- Optional second-point inputs (for k from two points) -----------------
        box2p = ttk.LabelFrame(frm, text="Two-point data (only for ‘k from two points’)")
        box2p.grid(row=6, column=0, columnspan=5, sticky="we", pady=(8,6))
        for c in range(5): box2p.grid_columnconfigure(c, weight=1)

        ttk.Label(box2p, text="t₁:").grid(row=0, column=0, sticky="e")
        t1_e = ttk.Entry(box2p, width=10); t1_e.grid(row=0, column=1, sticky="w", padx=(6,6))
        t1_u = ttk.Combobox(box2p, width=6, state="readonly", values=["s","min","h"]); t1_u.grid(row=0, column=2, sticky="w"); t1_u.set("s")
        ttk.Label(box2p, text="[A]₁ (M):").grid(row=0, column=3, sticky="e")
        A1_e = ttk.Entry(box2p, width=10); A1_e.grid(row=0, column=4, sticky="w", padx=(6,6))
        # A tiny top-right help button
        frm.grid_columnconfigure(5, weight=1)  # gives some stretch so Help can sit at the right
        help_btn = ttk.Button(frm, text="❓ Help (F1)", command=_open_help_kinetics)
        help_btn.grid(row=0, column=5, sticky="e", padx=(8,0))
        win.bind("<F1>", lambda e: _open_help_kinetics())

        ttk.Label(box2p, text="t₂:").grid(row=1, column=0, sticky="e")
        t2_e = ttk.Entry(box2p, width=10); t2_e.grid(row=1, column=1, sticky="w", padx=(6,6))
        t2_u = ttk.Combobox(box2p, width=6, state="readonly", values=["s","min","h"]); t2_u.grid(row=1, column=2, sticky="w"); t2_u.set("s")
        ttk.Label(box2p, text="[A]₂ (M):").grid(row=1, column=3, sticky="e")
        A2_e = ttk.Entry(box2p, width=10); A2_e.grid(row=1, column=4, sticky="w", padx=(6,6))

        # --- Preview of equations --------------------------------------------------
        eq_var = tk.StringVar(value="")
        ttk.Label(frm, textvariable=eq_var, foreground="#444").grid(row=7, column=0, columnspan=5, sticky="w", pady=(2,6))



        def _set_k_units(*_):
            o = order_var.get()[0]  # "0", "1", or "2"
            if o == "0":
                k_unit.config(values=["M/s", "M/min", "M/h"]); k_unit.set("M/s")
                eq_var.set("0th: rate = k ;  [A]t = [A]0 − k·t ;  t1/2 = [A]0/(2k)")
            elif o == "1":
                k_unit.config(values=["s^-1", "min^-1", "h^-1"]); k_unit.set("s^-1")
                eq_var.set("1st: rate = k·[A] ;  [A]t = [A]0·e^(−k t) ;  t1/2 = ln2 / k")
            else:
                k_unit.config(values=["M^-1 s^-1", "M^-1 min^-1", "M^-1 h^-1"]); k_unit.set("M^-1 s^-1")
                eq_var.set("2nd: rate = k·[A]^2 ;  1/[A]t = 1/[A]0 + k·t ;  t1/2 = 1/(k·[A]0)")
        order_box.bind("<<ComboboxSelected>>", _set_k_units)
        _set_k_units()

        # --- Output + explanation --------------------------------------------------
        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 11, "bold")).grid(row=7, column=0, columnspan=5, sticky="w", pady=(8,4))

        explain = tk.Text(frm, height=14, wrap="word")
        explain.grid(row=9, column=0, columnspan=5, sticky="nsew")
        frm.grid_rowconfigure(8, weight=1)
        scy = ttk.Scrollbar(frm, orient="vertical", command=explain.yview); scy.grid(row=9, column=5, sticky="ns")
        explain.configure(yscrollcommand=scy.set)

        # --- Helpers ---------------------------------------------------------------
        def _num(entry):
            return float(entry.get().strip().replace(",", "."))

        def _sec_factor(unit_str: str) -> float:
            return {"s":1.0, "min":60.0, "h":3600.0}[unit_str]

        def _k_to_base(k_val: float, k_unit_txt: str) -> float:
            # Convert k to base seconds units
            if "min" in k_unit_txt:
                return k_val / 60.0
            elif "h" in k_unit_txt:
                return k_val / 3600.0
            else:
                return k_val

        # --- Compute core ----------------------------------------------------------
        def compute():
            explain.configure(state="normal"); explain.delete("1.0","end")
            try:
                A0 = _num(A0_entry)
            except Exception:
                messagebox.showerror("Need [A]₀", "Enter a numeric [A]₀ (M).", parent=win); return

            # k (if present)
            k_val = None
            if k_entry.get().strip():
                try:
                    k_val = _k_to_base(_num(k_entry), k_unit.get())
                except Exception:
                    messagebox.showerror("Bad k", "Enter a numeric k.", parent=win); return

            # time t (if present)
            t_sec = None
            if t_entry.get().strip():
                try:
                    t_sec = _num(t_entry) * _sec_factor(t_unit.get())
                except Exception:
                    messagebox.showerror("Bad t", "Enter a numeric time.", parent=win); return

            ord0 = (order_var.get().startswith("0"))
            ord1 = (order_var.get().startswith("1"))
            ord2 = (order_var.get().startswith("2"))

            # pretty format
            def g(x): return f"{x:.6g}"

            lines = []
            task_key = task.get().split()[0]  # "A_t", "rate", "t_%", "t_A", "k_2p", "t12"

            # ---- A_t (concentration after t) ----
            if task_key == "A_t":
                if k_val is None or t_sec is None:
                    messagebox.showerror("Inputs", "Need k and t for [A] at time.", parent=win); return
                if ord0:
                    At = A0 - k_val * t_sec
                    lines += [f"[0th]  [A]t = [A]0 − k t = {g(A0)} − {g(k_val)}×{g(t_sec)} = {g(At)} M"]
                elif ord1:
                    At = A0 * math.exp(-k_val * t_sec)
                    lines += [f"[1st]  [A]t = [A]0·e^(−k t) = {g(A0)}·e^(−{g(k_val)}×{g(t_sec)}) = {g(At)} M"]
                else:
                    At = A0 / (1.0 + k_val * A0 * t_sec)
                    lines += [f"[2nd]  [A]t = [A]0 /(1 + k [A]0 t) = {g(A0)}/(1+{g(k_val)}×{g(A0)}×{g(t_sec)}) = {g(At)} M"]
                rate = k_val * (At if ord1 else (At**2 if ord2 else 1.0))
                out.set(f"[A](t) = {g(At)} M    ;    rate = {g(rate)} M/s")
                lines.append(f"rate(t) = k·[A]^n = {g(k_val)} × {g(At)}^{('1' if ord1 else '2' if ord2 else '0')} = {g(rate)} M/s")

            # ---- rate at t ----
            elif task_key == "rate":
                if k_val is None or t_sec is None:
                    messagebox.showerror("Inputs", "Need k and t to find rate.", parent=win); return
                if ord0:
                    At = A0 - k_val * t_sec
                elif ord1:
                    At = A0 * math.exp(-k_val * t_sec)
                else:
                    At = A0 / (1.0 + k_val * A0 * t_sec)
                rate = k_val * (At if ord1 else (At**2 if ord2 else 1.0))
                out.set(f"rate(t) = {g(rate)} M/s   (with [A](t) = {g(At)} M)")
                lines += [f"[A](t) computed as above → rate = k·[A]^n = {g(rate)} M/s"]

            # ---- t for % consumed ----
            elif task_key == "t_%":
                if k_val is None:
                    messagebox.showerror("Inputs", "Need k for time.", parent=win); return
                try:
                    pct = _num(pct_entry) / 100.0
                except Exception:
                    messagebox.showerror("Percent", "Enter percent consumed (e.g., 75).", parent=win); return
                if pct <= 0 or pct >= 1:
                    messagebox.showerror("Percent", "Use 0–100 (exclusive).", parent=win); return
                At = A0 * (1.0 - pct)
                if ord0:
                    t_sec = (A0 - At) / k_val
                    lines += [f"[0th]  At = A0(1−f) → t = (A0−At)/k = ({g(A0)}−{g(At)})/{g(k_val)} = {g(t_sec)} s"]
                elif ord1:
                    t_sec = (1.0 / k_val) * math.log(A0 / At)
                    lines += [f"[1st]  t = (1/k) ln(A0/At) = (1/{g(k_val)}) ln({g(A0)}/{g(At)}) = {g(t_sec)} s"]
                else:
                    t_sec = (1.0 / k_val) * (1.0/At - 1.0/A0)
                    lines += [f"[2nd]  t = (1/k)(1/At − 1/A0) = (1/{g(k_val)})({g(1/At)}−{g(1/A0)}) = {g(t_sec)} s"]
                out.set(f"Time for {pct*100:.1f}% consumed:  t ≈ {g(t_sec)} s  ({g(t_sec/60)} min)")

            # ---- t to reach a target At ----
            elif task_key == "t_A":
                if k_val is None or not At_entry.get().strip():
                    messagebox.showerror("Inputs", "Need k and target [A].", parent=win); return
                try:
                    At = _num(At_entry)
                except Exception:
                    messagebox.showerror("Target", "Enter numeric target [A].", parent=win); return
                if ord0:
                    t_sec = (A0 - At) / k_val
                    lines += [f"[0th]  t = (A0−At)/k = ({g(A0)}−{g(At)})/{g(k_val)} = {g(t_sec)} s"]
                elif ord1:
                    t_sec = (1.0 / k_val) * math.log(A0 / At)
                    lines += [f"[1st]  t = (1/k) ln(A0/At) = {g(t_sec)} s"]
                else:
                    t_sec = (1.0 / k_val) * (1.0/At - 1.0/A0)
                    lines += [f"[2nd]  t = (1/k)(1/At − 1/A0) = {g(t_sec)} s"]
                out.set(f"Time to reach [A]={g(At)} M:  t ≈ {g(t_sec)} s  ({g(t_sec/60)} min)")

            # ---- solve k from two points ----
            elif task_key == "k_2p":
                try:
                    t1 = _num(t1_e) * _sec_factor(t1_u.get())
                    t2 = _num(t2_e) * _sec_factor(t2_u.get())
                    A1 = _num(A1_e); A2 = _num(A2_e)
                except Exception:
                    messagebox.showerror("Two-point data", "Enter t₁, t₂, [A]₁, [A]₂.", parent=win); return
                dt = t2 - t1
                if dt == 0:
                    messagebox.showerror("Two-point data", "t₂ must differ from t₁.", parent=win); return
                if ord0:
                    k_val = (A1 - A2) / dt
                    lines += [f"[0th]  [A]t = [A]0 − k t  ⇒  slope = −k  ⇒  k = (A1−A2)/(t2−t1) = {g(k_val)} M/s"]
                    k_unit.set("M/s")
                elif ord1:
                    k_val = (math.log(A1) - math.log(A2)) / dt
                    lines += [f"[1st]  ln[A] vs t slope = −k  ⇒  k = (ln A1 − ln A2)/Δt = {g(k_val)} s^-1"]
                    k_unit.set("s^-1")
                else:
                    k_val = (1.0/A2 - 1.0/A1) / dt
                    lines += [f"[2nd]  1/[A] vs t slope = +k  ⇒  k = (1/A2 − 1/A1)/Δt = {g(k_val)} M^-1 s^-1"]
                    k_unit.set("M^-1 s^-1")
                out.set(f"k ≈ {g(k_val)}  ({k_unit.get()})")
                k_entry.delete(0,"end"); k_entry.insert(0, g(k_val))
            
            elif task_key == "k_t12":
                if not t12_entry.get().strip():
                    messagebox.showerror("Half-life", "Enter t₁/₂.", parent=win); return
                try:
                    t12s = _num(t12_entry) * _sec_factor(t12_unit.get())   # seconds
                except Exception:
                    messagebox.showerror("Half-life", "Bad t₁/₂.", parent=win); return

                # compute k in base (per second) units
                if ord0:
                    # t1/2 = [A]0 / (2k)  →  k = [A]0 / (2 t1/2)
                    k_base = A0 / (2.0 * t12s)                 # M/s
                    target_units = {"s":"M/s", "min":"M/min", "h":"M/h"}[t12_unit.get()]
                elif ord1:
                    # t1/2 = ln 2 / k  →  k = ln 2 / t1/2
                    k_base = math.log(2.0) / t12s              # s^-1
                    target_units = {"s":"s^-1", "min":"min^-1", "h":"h^-1"}[t12_unit.get()]
                else:
                    # t1/2 = 1 / (k [A]0)  →  k = 1 / (t1/2 · [A]0)
                    k_base = 1.0 / (t12s * A0)                 # M^-1 s^-1
                    target_units = {"s":"M^-1 s^-1", "min":"M^-1 min^-1", "h":"M^-1 h^-1"}[t12_unit.get()]

                # set k_unit to match the half-life time unit
                for opt in k_unit["values"]:
                    if opt == target_units:
                        k_unit.set(opt); break

                # convert base (per second) → display unit (s/min/h)
                scale = 1.0
                if "min" in k_unit.get():   scale = 60.0
                elif "h" in k_unit.get():   scale = 3600.0
                k_disp = k_base * scale

                k_entry.delete(0, "end"); k_entry.insert(0, f"{k_disp:.6g}")
                out.set(f"k = {k_disp:.6g}  {k_unit.get()}")
                explain.insert("end",
                    f"Using half-life relation for order {order_var.get()[0]}:\n"
                    + (f"  k = [A]0/(2 t1/2) = {A0:.6g}/(2×{_num(t12_entry):.6g} {t12_unit.get()})\n" if ord0 else
                        f"  k = ln2 / t1/2 = 0.693/{_num(t12_entry):.6g} {t12_unit.get()}\n" if ord1 else
                        f"  k = 1/(t1/2 · [A]0) = 1/({_num(t12_entry):.6g} {t12_unit.get()} × {A0:.6g})\n")
                )

            # ---- half-life ----
            else:  # "t12"
                if k_val is None:
                    messagebox.showerror("Need k", "Enter k to get t₁/₂.", parent=win); return
                if ord0:
                    t12 = A0 / (2.0 * k_val)
                    lines += [f"[0th]  t1/2 = [A]0/(2k) = {g(A0)}/(2×{g(k_val)}) = {g(t12)} s"]
                elif ord1:
                    t12 = math.log(2.0) / k_val
                    lines += [f"[1st]  t1/2 = ln2 / k = {g(t12)} s"]
                else:
                    t12 = 1.0 / (k_val * A0)
                    lines += [f"[2nd]  t1/2 = 1/(k [A]0) = 1/({g(k_val)}×{g(A0)}) = {g(t12)} s"]
                out.set(f"t₁/₂ ≈ {g(t12)} s  ({g(t12/60)} min)")

            explain.insert("end", "\n".join(lines) + "\n")
            explain.configure(state="disabled")

            return  # (we don’t need to return a number; UI shows all)
        

        # Put the compute bar at the bottom (Enter triggers compute)
        self._add_compute_bar(win, compute)

        # --- Add-to-known shortcuts ----------------------------------------------
        btns = ttk.Frame(frm); btns.grid(row=8, column=0, columnspan=5, sticky="w", pady=(8,0))
        def add_k_known():
            if not k_entry.get().strip(): 
                messagebox.showinfo("k", "Compute or type k first.", parent=win); return
            try:
                k_base = _k_to_base(_num(k_entry), k_unit.get())
            except Exception:
                messagebox.showerror("k", "Bad k.", parent=win); return
            self.known_base["k_rate"] = float(k_base)   # store in base (seconds)
            self.known_ui["k_rate"]   = {"unit": k_unit.get(), "display_value": float(_num(k_entry))}
            self._refresh_known_table(); self._refresh_equation_list()
            messagebox.showinfo("Added", f"Saved k = {k_entry.get()} {k_unit.get()} as ‘k_rate’.", parent=win)

        def add_t12_known():
            # simulate pressing compute on half-life
            old = task.get(); task.set("t12  • Half-life (t₁/₂) from k and [A]₀"); compute(); task.set(old)
            txt = out.get()
            # crude parse for seconds value
            import re
            m = re.search(r"≈ ([\d.eE+-]+) s", txt)
            if not m: return
            t12 = float(m.group(1))
            self.known_base["t1/2"] = t12  # seconds
            self.known_ui["t1/2"]   = {"unit":"s","display_value":t12}
            self._refresh_known_table(); self._refresh_equation_list()
            messagebox.showinfo("Added", f"Saved t₁/₂ ≈ {t12:g} s.", parent=win)

        ttk.Button(btns, text="Add k to Known", command=add_k_known).pack(side=tk.LEFT, padx=(0,8))
        ttk.Button(btns, text="Add t₁/₂ to Known (compute first)", command=add_t12_known).pack(side=tk.LEFT)

    
    def _open_arrhenius_tool(self):
        import math
        win = tk.Toplevel(self)
        win.title("Arrhenius Equation — k, Ea, A, T")
        win.transient(self); win.grab_set(); win.geometry("820x680")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # --- Constants / R ---------------------------------------------------------
        ttk.Label(frm, text="Gas constant R:").grid(row=0, column=0, sticky="e")
        R_entry = ttk.Entry(frm, width=14); R_entry.grid(row=0, column=1, sticky="w", padx=(6,12))
        R_default = float(self.known_base.get("R", 8.314462618))
        R_entry.insert(0, f"{R_default:g}")
        ttk.Label(frm, text="J/(mol·K)").grid(row=0, column=2, sticky="w")

        # ln vs log10 checkbox (output formatting/help only; math uses ln)
        use_log10 = tk.BooleanVar(value=False)
        ttk.Checkbutton(frm, text="Show log₁₀ form too", variable=use_log10).grid(row=0, column=3, sticky="w")

        # --- Solve-for selector ----------------------------------------------------
        ttk.Label(frm, text="Solve for:").grid(row=1, column=0, sticky="e", pady=(4,0))
        solve_for = tk.StringVar(value="k")
        srow = ttk.Frame(frm); srow.grid(row=1, column=1, columnspan=3, sticky="w", pady=(2,0))
        for lab, val in [("k (rate constant)", "k"), ("Ea (activation energy)", "Ea"),
                        ("A (pre-exponential)", "A"), ("T (temperature)", "T")]:
            ttk.Radiobutton(srow, text=lab, variable=solve_for, value=val).pack(side=tk.LEFT, padx=(0,12))

        # --- Single-point inputs ---------------------------------------------------
        box = ttk.LabelFrame(frm, text="Single-point Arrhenius  k = A · exp(−Ea/(R·T))")
        box.grid(row=2, column=0, columnspan=4, sticky="we", pady=(8,8))
        for c in range(8): box.grid_columnconfigure(c, weight=1)

        # k
        ttk.Label(box, text="k:").grid(row=0, column=0, sticky="e")
        k_entry = ttk.Entry(box, width=16); k_entry.grid(row=0, column=1, sticky="w", padx=(6,8))
        ttk.Label(box, text="units of k:").grid(row=0, column=2, sticky="e")
        k_unit = ttk.Combobox(box, width=14, state="readonly",
                            values=["s⁻¹","min⁻¹","h⁻¹","M⁻¹·s⁻¹","M⁻¹·min⁻¹","custom…"])
        k_unit.grid(row=0, column=3, sticky="w"); k_unit.set("s⁻¹")

        # A
        ttk.Label(box, text="A (same units as k):").grid(row=1, column=0, sticky="e")
        A_entry = ttk.Entry(box, width=16); A_entry.grid(row=1, column=1, sticky="w", padx=(6,8))

        # Ea + unit
        ttk.Label(box, text="Ea:").grid(row=2, column=0, sticky="e")
        Ea_entry = ttk.Entry(box, width=16); Ea_entry.grid(row=2, column=1, sticky="w", padx=(6,8))
        Ea_unit = ttk.Combobox(box, width=10, state="readonly",
                            values=["J/mol","kJ/mol","cal/mol","kcal/mol"])
        Ea_unit.grid(row=2, column=2, sticky="w"); Ea_unit.set("kJ/mol")

        # Temperature + unit
        ttk.Label(box, text="T:").grid(row=3, column=0, sticky="e")
        T_entry = ttk.Entry(box, width=16); T_entry.grid(row=3, column=1, sticky="w", padx=(6,8))
        T_unit = ttk.Combobox(box, width=8, state="readonly", values=["K","°C"])
        T_unit.grid(row=3, column=2, sticky="w"); T_unit.set("K")

        # Prefill from Known if present
        if "Ea" in self.known_ui and not Ea_entry.get():
            try:
                valJ = float(self.known_base["Ea"])
                Ea_entry.insert(0, f"{valJ/1000.0:g}"); Ea_unit.set("kJ/mol")
            except Exception: pass
        if "T" in self.known_ui and not T_entry.get():
            try: T_entry.insert(0, f"{float(self.known_base['T']):g}"); T_unit.set("K")
            except Exception: pass
        if "k" in self.known_ui and not k_entry.get():
            try: k_entry.insert(0, f"{float(self.known_ui['k']['display_value']):g}")
            except Exception: pass
        if "A" in self.known_ui and not A_entry.get():
            try: A_entry.insert(0, f"{float(self.known_ui['A']['display_value']):g}")
            except Exception: pass

        # --- Two-point panel -------------------------------------------------------
        box2 = ttk.LabelFrame(frm, text="Two-point (Arrhenius)   ln(k₂/k₁) = −(Ea/R)·(1/T₂ − 1/T₁)")
        box2.grid(row=3, column=0, columnspan=4, sticky="we", pady=(4,8))
        for c in range(10): box2.grid_columnconfigure(c, weight=1)

        ttk.Label(box2, text="k₁:").grid(row=0, column=0, sticky="e"); k1E=ttk.Entry(box2,width=12); k1E.grid(row=0,column=1,sticky="w",padx=(6,8))
        ttk.Label(box2, text="T₁:").grid(row=0, column=2, sticky="e"); T1E=ttk.Entry(box2,width=12); T1E.grid(row=0,column=3,sticky="w",padx=(6,8))
        T1U=ttk.Combobox(box2,width=6, state="readonly", values=["K","°C"]); T1U.grid(row=0,column=4,sticky="w"); T1U.set("K")

        ttk.Label(box2, text="k₂:").grid(row=1, column=0, sticky="e"); k2E=ttk.Entry(box2,width=12); k2E.grid(row=1,column=1,sticky="w",padx=(6,8))
        ttk.Label(box2, text="T₂:").grid(row=1, column=2, sticky="e"); T2E=ttk.Entry(box2,width=12); T2E.grid(row=1,column=3,sticky="w",padx=(6,8))
        T2U=ttk.Combobox(box2,width=6, state="readonly", values=["K","°C"]); T2U.grid(row=1,column=4,sticky="w"); T2U.set("K")

        btns2 = ttk.Frame(box2); btns2.grid(row=2, column=0, columnspan=6, sticky="w", pady=(6,0))
        def _two_point_Ea():
            try:
                R = float(R_entry.get().replace(",", "."))
                k1 = float(k1E.get().replace(",", ".")); k2 = float(k2E.get().replace(",", "."))
                T1 = float(T1E.get().replace(",", ".")); T2 = float(T2E.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid", "Enter numbers for k₁, k₂, T₁, T₂."); return
            T1K = T1 if T1U.get()=="K" else T1 + 273.15
            T2K = T2 if T2U.get()=="K" else T2 + 273.15
            if k1<=0 or k2<=0 or T1K<=0 or T2K<=0:
                messagebox.showerror("Invalid", "k and T must be > 0."); return
            Ea = R * math.log(k1/k2) / (1.0/T2K - 1.0/T1K)  # J/mol
            Ea_entry.delete(0,"end"); Ea_entry.insert(0, f"{Ea/1000.0:g}"); Ea_unit.set("kJ/mol")
            out.set(f"Two-point: Ea = {Ea/1000.0:.6g} kJ/mol  (= {Ea:.6g} J/mol)")
            _append_steps(
                f"Two-point Ea:\n  ln(k2/k1) = -Ea/R · (1/T2 - 1/T1)\n"
                f"  Ea = R·ln(k1/k2) / (1/T2 - 1/T1)\n"
                f"  → Ea = {R:g}·ln({k1:g}/{k2:g}) / (1/{T2K:g} - 1/{T1K:g}) = {Ea/1000.0:.6g} kJ/mol\n"
            )

        def _two_point_k2():
            try:
                R = float(R_entry.get().replace(",", "."))
                k1 = float(k1E.get().replace(",", ".")); T1 = float(T1E.get().replace(",", "."))
                T2 = float(T2E.get().replace(",", "."))
                # Ea from main box (convert to J/mol)
                Ea_val = float(Ea_entry.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid", "Need k₁, T₁, T₂ and Ea (in main box)."); return
            T1K = T1 if T1U.get()=="K" else T1 + 273.15
            T2K = T2 if T2U.get()=="K" else T2 + 273.15
            factor = {"J/mol":1.0,"kJ/mol":1000.0,"cal/mol":4.184,"kcal/mol":4184.0}[Ea_unit.get()]
            EaJ = Ea_val * factor
            k2 = k1 * math.exp(-EaJ/R * (1.0/T2K - 1.0/T1K))
            k2E.delete(0,"end"); k2E.insert(0, f"{k2:.6g}")
            out.set(f"Two-point: k₂ = {k2:.6g} (units same as k₁)")
            _append_steps(
                f"Two-point k₂:\n  ln(k2/k1) = -Ea/R · (1/T2 - 1/T1)\n"
                f"  k2 = k1·exp[-Ea/R·(1/T2 - 1/T1)] = {k2:.6g}\n"
            )

        ttk.Button(btns2, text="Compute Ea from (k₁,T₁),(k₂,T₂)", command=_two_point_Ea).pack(side=tk.LEFT, padx=(0,8))
        ttk.Button(btns2, text="Predict k₂ from k₁ + Ea", command=_two_point_k2).pack(side=tk.LEFT)

        # --- Output widgets --------------------------------------------------------
        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI",10,"bold")).grid(row=4, column=0, columnspan=4, sticky="w", pady=(8,4))

        steps = tk.Text(frm, height=12, wrap="word")
        steps.grid(row=5, column=0, columnspan=4, sticky="nsew")
        frm.grid_rowconfigure(5, weight=1)
        sc = ttk.Scrollbar(frm, orient="vertical", command=steps.yview); steps.configure(yscrollcommand=sc.set)
        sc.grid(row=5, column=4, sticky="ns")

        def _append_steps(s):
            steps.configure(state="normal")
            steps.insert("end", s + ("\n" if not s.endswith("\n") else ""))
            steps.configure(state="disabled")

        # --- Single-point compute --------------------------------------------------
        def compute():
            steps.configure(state="normal"); steps.delete("1.0", "end"); steps.configure(state="disabled")
            try:
                R = float(R_entry.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid R", "Enter a numeric gas constant R."); return

            # Pull & convert inputs (some may be blank)
            k_str, A_str, Ea_str, T_str = k_entry.get().strip(), A_entry.get().strip(), Ea_entry.get().strip(), T_entry.get().strip()
            # Ea to J/mol if present
            EaJ = None
            if Ea_str:
                try:
                    eaval = float(Ea_str.replace(",", "."))
                    EaJ = eaval * {"J/mol":1.0,"kJ/mol":1000.0,"cal/mol":4.184,"kcal/mol":4184.0}[Ea_unit.get()]
                except ValueError:
                    messagebox.showerror("Invalid Ea", "Enter a numeric Ea."); return
            # T to K if present
            TK = None
            if T_str:
                try:
                    tval = float(T_str.replace(",", "."))
                    TK = tval if T_unit.get()=="K" else tval + 273.15
                except ValueError:
                    messagebox.showerror("Invalid T", "Enter a numeric temperature."); return
                if TK <= 0:
                    messagebox.showerror("Invalid T", "Absolute temperature must be > 0 K."); return

            # Parse k and A if present
            kval = None; Aval = None
            if k_str:
                try: kval = float(k_str.replace(",", "."))
                except ValueError: messagebox.showerror("Invalid k","Enter a numeric k."); return
                if kval <= 0: messagebox.showerror("Invalid k","k must be > 0."); return
            if A_str:
                try: Aval = float(A_str.replace(",", "."))
                except ValueError: messagebox.showerror("Invalid A","Enter a numeric A."); return
                if Aval <= 0: messagebox.showerror("Invalid A","A must be > 0."); return

            what = solve_for.get()

            # Helper writers
            def show_result(label, value, unit_text=""):
                if unit_text: unit_text = " " + unit_text
                out.set(f"{label} = {value:.6g}{unit_text}")

            # Compute by target
            if what == "k":
                if Aval is None or EaJ is None or TK is None:
                    messagebox.showerror("Need inputs", "Provide A, Ea and T to solve for k."); return
                k = Aval * math.exp(-EaJ/(R*TK))
                k_entry.delete(0,"end"); k_entry.insert(0, f"{k:.6g}")
                show_result("k", k, k_unit.get())
                _append_steps(
                    f"k = A·exp(-Ea/(R·T))\n"
                    f"  = {Aval:g}·exp(-{EaJ:g}/({R:g}·{TK:g})) = {k:.6g} {k_unit.get()}\n"
                    + (_log10_info(Aval, EaJ, R, TK, k) if use_log10.get() else "")
                )

            elif what == "Ea":
                if kval is None or Aval is None or TK is None:
                    messagebox.showerror("Need inputs", "Provide k, A and T to solve for Ea."); return
                EaJ = - R*TK*math.log(kval/Aval)
                Ea_entry.delete(0,"end"); Ea_entry.insert(0, f"{EaJ/1000.0:.6g}"); Ea_unit.set("kJ/mol")
                show_result("Ea", EaJ/1000.0, "kJ/mol")
                _append_steps(
                    f"ln k = ln A − Ea/(R·T)  ⇒  Ea = −R·T·ln(k/A)\n"
                    f"  = −{R:g}·{TK:g}·ln({kval:g}/{Aval:g}) = {EaJ/1000.0:.6g} kJ/mol\n"
                    + (_log10_rearranged(kval, Aval, R, TK) if use_log10.get() else "")
                )

            elif what == "A":
                if kval is None or EaJ is None or TK is None:
                    messagebox.showerror("Need inputs", "Provide k, Ea and T to solve for A."); return
                A = kval * math.exp(EaJ/(R*TK))
                A_entry.delete(0,"end"); A_entry.insert(0, f"{A:.6g}")
                show_result("A", A, f"{k_unit.get()}")
                _append_steps(
                    f"A = k·exp(Ea/(R·T))\n"
                    f"  = {kval:g}·exp({EaJ:g}/({R:g}·{TK:g})) = {A:.6g} ({k_unit.get()})\n"
                    + (_log10_A(kval, EaJ, R, TK, A) if use_log10.get() else "")
                )

            else:  # what == "T"
                if kval is None or Aval is None or EaJ is None:
                    messagebox.showerror("Need inputs", "Provide k, A and Ea to solve for T."); return
                arg = Aval / kval
                if arg <= 0 or arg == 1.0:
                    messagebox.showerror("Invalid", "A/k must be > 0 and ≠ 1."); return
                Tsol = EaJ / (R * math.log(arg))
                if T_unit.get() == "°C":
                    disp = Tsol - 273.15; unitlab = "°C"
                else:
                    disp = Tsol; unitlab = "K"
                T_entry.delete(0,"end"); T_entry.insert(0, f"{disp:.6g}")
                show_result("T", disp, unitlab)
                _append_steps(
                    f"k = A·exp(-Ea/(R·T)) ⇒ T = Ea / [R·ln(A/k)]\n"
                    f"  = {EaJ:g} / [{R:g}·ln({Aval:g}/{kval:g})] = {Tsol:.6g} K\n"
                    + (_log10_T(Aval, kval, EaJ, R, Tsol) if use_log10.get() else "")
                )

        def _log10_info(A, EaJ, R, T, k):
            return (f"\nlog₁₀ form:\n  log k = log A − Ea/(2.303·R)·(1/T)\n"
                    f"  log k = log({A:g}) − {EaJ:g}/(2.303·{R:g})·(1/{T:g}) ⇒ k ≈ {k:.6g}\n")

        def _log10_rearranged(k, A, R, T):
            return (f"\nlog₁₀ rearranged:\n  Ea = 2.303·R·T·(log A − log k)\n")

        def _log10_A(k, EaJ, R, T, A):
            return (f"\nlog₁₀ form:\n  log A = log k + Ea/(2.303·R·T) ⇒ A ≈ {A:.6g}\n")

        def _log10_T(A, k, EaJ, R, TK):
            return (f"\nlog₁₀ form:\n  1/T = [log A − log k]·(2.303·R)/Ea  ⇒  T ≈ {TK:.6g} K\n")

        # --- Compute bar at the bottom (Enter triggers) ----------------------------
        self._add_compute_bar(win, compute)

        # --- Add to Known Variables ------------------------------------------------
        btm = ttk.Frame(frm); btm.grid(row=6, column=0, columnspan=4, sticky="w", pady=(8,0))
        def _add_known():
            # Save what we have (if parseable)
            saved = []
            # T
            if T_entry.get().strip():
                try:
                    tval = float(T_entry.get().replace(",",".")); TK = tval if T_unit.get()=="K" else tval + 273.15
                    self.known_base["T"] = float(TK); self.known_ui["T"]={"unit":"K","display_value":float(TK)}; saved.append("T")
                except: pass
            # k
            if k_entry.get().strip():
                try:
                    kval = float(k_entry.get().replace(",",".")); self.known_base["k"]=float(kval)
                    self.known_ui["k"]={"unit":k_unit.get(), "display_value":float(kval)}; saved.append("k")
                except: pass
            # A
            if A_entry.get().strip():
                try:
                    Aval = float(A_entry.get().replace(",",".")); self.known_base["A"]=float(Aval)
                    self.known_ui["A"]={"unit":k_unit.get(), "display_value":float(Aval)}; saved.append("A")
                except: pass
            # Ea
            if Ea_entry.get().strip():
                try:
                    eaval = float(Ea_entry.get().replace(",",".")); fac={"J/mol":1.0,"kJ/mol":1000.0,"cal/mol":4.184,"kcal/mol":4184.0}[Ea_unit.get()]
                    EaJ = eaval*fac; self.known_base["Ea"]=float(EaJ)
                    self.known_ui["Ea"]={"unit":"kJ/mol","display_value":float(EaJ/1000.0)}; saved.append("Ea")
                except: pass

            # R
            try:
                Rv = float(R_entry.get().replace(",",".")); self.known_base["R"]=float(Rv)
                self.known_ui["R"]={"unit":"J/(mol·K)","display_value":float(Rv)}; saved.append("R")
            except: pass

            self._refresh_known_table(); self._refresh_equation_list()
            messagebox.showinfo("Saved", "Added to Known: " + (", ".join(saved) if saved else "(nothing parseable)"), parent=win)

        ttk.Button(btm, text="Add current values to Known Variables", command=_add_known).pack(side=tk.LEFT)

        # --- Help popup ------------------------------------------------------------
        def _open_help():
            h = tk.Toplevel(win); h.title("Arrhenius — Help"); h.transient(win); h.grab_set(); h.geometry("740x560")
            pad = ttk.Frame(h, padding=12); pad.pack(fill=tk.BOTH, expand=True)
            t = tk.Text(pad, wrap="word"); t.pack(fill=tk.BOTH, expand=True, side=tk.LEFT)
            vs = ttk.Scrollbar(pad, orient="vertical", command=t.yview); t.configure(yscrollcommand=vs.set); vs.pack(side=tk.RIGHT, fill=tk.Y)
            t.insert("1.0",
                    "Arrhenius Equation\n"
                    "==================\n"
                    "k = A·exp(−Ea/(R·T))  (natural log form)\n"
                    "ln k = ln A − Ea/(R·T)\n"
                    "log10 k = log10 A − Ea/(2.303·R·T)\n\n"
                    "Solve-for tips:\n"
                    "  • Solve k: needs A, Ea, T\n"
                    "  • Solve Ea: needs k, A, T\n"
                    "  • Solve A: needs k, Ea, T\n"
                    "  • Solve T: needs k, A, Ea\n\n"
                    "Two-point form (no A needed):\n"
                    "  ln(k2/k1) = −(Ea/R)·(1/T2 − 1/T1)\n"
                    "  • Use to extract Ea from two measurements, or predict k2 from k1 + Ea.\n\n"
                    "Units:\n"
                    "  • R is J/(mol·K). Ea is converted to J/mol internally (kJ/mol, cal/mol, kcal/mol supported).\n"
                    "  • Temperature may be K or °C (converted to K).\n"
                    "  • k and A share units; choose any label (s⁻¹, M⁻¹·s⁻¹, …). The tool treats them as consistent.\n\n"
                    "Workflow examples:\n"
                    "  1) Given A, Ea (kJ/mol), T(°C) → choose 'solve k' and Compute.\n"
                    "  2) From (k1,T1),(k2,T2) → press 'Compute Ea', then optionally use that Ea to predict k at another T.\n"
                    )
            t.configure(state="disabled")
            ttk.Button(pad, text="Close (Esc)", command=h.destroy).pack(anchor="e", pady=(8,0))
            h.bind("<Escape>", lambda e: h.destroy())

        help_btn = ttk.Button(frm, text="❓ Help (F1)", command=_open_help)
        help_btn.grid(row=6, column=3, sticky="e")
        win.bind("<F1>", lambda e: _open_help())

    
    # ---- tiny parser for "a A + b B <-> c C + d D" ----
    def _parse_reaction(self, s: str):
        import re
        s = s.replace("⇌", "->").replace("↔", "->").replace("⟷", "->").replace("=", "->")
        if "->" not in s:
            raise ValueError("Use an arrow like '->' or '⇌' between sides.")
        L, R = [part.strip() for part in s.split("->", 1)]

        def parse_side(side, sign):
            species = {}
            if not side:
                return species
            for term in side.split("+"):
                term = term.strip()
                if not term: 
                    continue
                m = re.match(r"^\s*(\d*\.?\d*)\s*([A-Za-z0-9()^+\-·_]+)\s*$", term)
                if not m:
                    # allow plain species (implicit 1)
                    coef, sp = "", term
                else:
                    coef, sp = m.group(1), m.group(2)
                nu = float(coef) if coef else 1.0
                species[sp] = species.get(sp, 0.0) + sign * nu
            return species

        st = {}
        for k, v in parse_side(L, -1).items(): st[k] = st.get(k,0)+v
        for k, v in parse_side(R, +1).items(): st[k] = st.get(k,0)+v
        # ensure both sides present
        if all(v>=0 for v in st.values()) or all(v<=0 for v in st.values()):
            raise ValueError("Reaction needs reactants and products on different sides.")
        return st  # dict: species -> ν (reactants negative; products positive)

    # ---- robust 1D extent solver for ICE (bisection + fallback) ----
    def _solve_extent_for_K(self, stoich, C0, K, is_pressure=False):
        import math
        # admissible interval for x s.t. all Ci(x)=C0_i + ν_i x >= 0
        lo = -1e30; hi = 1e30
        for sp, nu in stoich.items():
            if abs(nu) < 1e-16: 
                continue
            if nu > 0:
                lo = max(lo, -C0[sp]/nu)  # x >= -C0/nu
            else:
                hi = min(hi, -C0[sp]/nu)  # x <= -C0/nu
        if not (lo < hi):
            raise ValueError("No feasible extent (check initials).")

        def Q_of(x):
            num = den = 1.0
            for sp, nu in stoich.items():
                c = C0[sp] + nu * x
                if c <= 0: 
                    return None
                if nu > 0:
                    num *= c**nu
                elif nu < 0:
                    den *= (c**(-nu))
            return num/den

        # objective: f(x) = Q(x) - K
        def f(x):
            q = Q_of(x)
            return None if q is None else (q - K)

        # try bracket with coarse scan
        N = 200
        xs = [lo + (hi-lo)*i/N for i in range(N+1)]
        vals = [f(x) for x in xs]
        bracket = None
        for i in range(N):
            a,b = xs[i], xs[i+1]
            fa, fb = vals[i], vals[i+1]
            if fa is None or fb is None: 
                continue
            if fa==0: return a
            if fb==0: return b
            if fa*fb < 0:
                bracket = (a,b); break

        # bisection if bracketed
        if bracket:
            a,b = bracket
            for _ in range(100):
                m = 0.5*(a+b)
                fm = f(m)
                if fm is None: 
                    break
                if abs(fm) < 1e-12:
                    return m
                fa = f(a)
                if fa*fm <= 0:
                    b = m
                else:
                    a = m
            return 0.5*(a+b)

        # fallback: choose x with minimal |f|
        best = None; bestx = None
        for x,v in zip(xs, vals):
            if v is None: 
                continue
            if best is None or abs(v) < best:
                best = abs(v); bestx = x
        if bestx is None:
            raise ValueError("Failed to evaluate Q over feasible x (check inputs).")
        return bestx
    

    def _open_equilibrium_ice_tool(self):
        import math
        win = tk.Toplevel(self)
        win.title("ICE Solver — Kc/Kp, Q, equilibrium shift")
        win.transient(self); win.grab_set(); win.geometry("900x680")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # --- reaction + type -------------------------------------------------------
        ttk.Label(frm, text="Reaction (e.g.,  N2 + 3 H2 ⇌ 2 NH3 ):").grid(row=0, column=0, sticky="w")
        rxn = ttk.Entry(frm, width=50); rxn.grid(row=0, column=1, columnspan=3, sticky="we", padx=(6,0))
        rxn.insert(0, "N2 + 3 H2 ⇌ 2 NH3")

        kind = tk.StringVar(value="Kc")
        ttk.Radiobutton(frm, text="Use Kc (M)", variable=kind, value="Kc").grid(row=1, column=1, sticky="w", pady=(6,0))
        ttk.Radiobutton(frm, text="Use Kp (atm)", variable=kind, value="Kp").grid(row=1, column=2, sticky="w", pady=(6,0))

        # Gas constant / temperature (for Kp↔Kc hints if needed)
        ttk.Label(frm, text="T:").grid(row=1, column=3, sticky="e")
        TE = ttk.Entry(frm, width=10); TE.grid(row=1, column=4, sticky="w")
        if "T" in self.known_base:
            TE.insert(0, f"{float(self.known_base['T']):g}")
        else:
            TE.insert(0, "298.15")
        ttk.Label(frm, text="K").grid(row=1, column=5, sticky="w")

        # --- species table built after "Parse" ------------------------------------
        table = ttk.Frame(frm); table.grid(row=3, column=0, columnspan=6, sticky="nsew", pady=(8,6))
        frm.grid_rowconfigure(3, weight=1); frm.grid_columnconfigure(1, weight=1)

        species_widgets = {}  # sp -> (label, init_entry)
        stoich = {}

        def build_table():
            nonlocal stoich, species_widgets
            for w in table.winfo_children(): w.destroy()
            try:
                stoich = self._parse_reaction(rxn.get())
            except Exception as ex:
                messagebox.showerror("Parse error", str(ex), parent=win); return

            hdr = ("Species","ν (prod + / react −)","Initial (M or atm)")
            for j,h in enumerate(hdr):
                ttk.Label(table, text=h, font=("Segoe UI",10,"bold")).grid(row=0, column=j, sticky="w", padx=4, pady=2)

            species_widgets = {}
            for i,(sp,nu) in enumerate(sorted(stoich.items()), start=1):
                ttk.Label(table, text=sp).grid(row=i, column=0, sticky="w", padx=4)
                ttk.Label(table, text=f"{nu:g}").grid(row=i, column=1, sticky="w")
                e = ttk.Entry(table, width=14); e.grid(row=i, column=2, sticky="w")
                e.insert(0, "0.00")
                species_widgets[sp] = e

        ttk.Button(frm, text="Parse/Refresh species", command=build_table).grid(row=2, column=0, sticky="w", pady=(8,0))
        build_table()

        # K entry & Q button
        ttk.Label(frm, text="K value:").grid(row=4, column=0, sticky="e")
        KE = ttk.Entry(frm, width=16); KE.grid(row=4, column=1, sticky="w")
        ttk.Label(frm, textvariable=tk.StringVar(value="(enter K to solve equilibrium; leave blank to compute Q only)"),
                foreground="#666").grid(row=4, column=2, columnspan=3, sticky="w")

        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI",10,"bold")).grid(row=5, column=0, columnspan=6, sticky="w", pady=(8,4))

        steps = tk.Text(frm, height=14, wrap="word")
        steps.grid(row=6, column=0, columnspan=6, sticky="nsew")
        frm.grid_rowconfigure(6, weight=1)
        sc = ttk.Scrollbar(frm, orient="vertical", command=steps.yview); steps.configure(yscrollcommand=sc.set)
        sc.grid(row=6, column=6, sticky="ns")

        def note(s):
            steps.configure(state="normal"); steps.insert("end", s + ("\n" if not s.endswith("\n") else "")); steps.configure(state="disabled")

        def compute():
            steps.configure(state="normal"); steps.delete("1.0","end"); steps.configure(state="disabled")
            # gather initials
            C0 = {}
            try:
                for sp,e in species_widgets.items():
                    C0[sp] = float(e.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid", "Initial entries must be numbers.", parent=win); return

            # quotient from initials
            def Q_from(C):
                num = den = 1.0
                for sp,nu in stoich.items():
                    if C[sp] <= 0: 
                        return None
                    if nu>0: num *= (C[sp]**nu)
                    elif nu<0: den *= (C[sp]**(-nu))
                return num/den

            Q0 = Q_from(C0)
            if Q0 is None:
                out.set("Q not defined (zero/negative entry).")
            else:
                out.set(f"Initial reaction quotient Q = {Q0:.6g}  ({kind.get()})")
                note(f"Equilibrium expression ({kind.get()}):  K = Π(products)^ν / Π(reactants)^|ν|\n"
                    f"Initial Q with given amounts: {Q0:.6g}")
                if Q0 is not None:
                    if KE.get().strip():
                        try:
                            K_val = float(KE.get().replace(",", "."))
                            if Q0 < K_val:
                                note("Direction: Q < K → reaction will run toward products (→ right) to reach equilibrium.")
                            elif Q0 > K_val:
                                note("Direction: Q > K → reaction will run toward reactants (← left) to reach equilibrium.")
                            else:
                                note("Direction: Q = K → system is at equilibrium; no net change.")
                        except ValueError:
                            pass
                    else:
                        note("Direction: (K not entered) — cannot predict shift without K.")

            # If K given, solve equilibrium via single-extent method
            if KE.get().strip():
                try:
                    K = float(KE.get().replace(",", "."))
                except ValueError:
                    messagebox.showerror("Bad K", "K must be numeric.", parent=win); return

                try:
                    x = self._solve_extent_for_K(stoich, C0, K, is_pressure=(kind.get()=="Kp"))
                except Exception as ex:
                    messagebox.showerror("Solve failed", str(ex), parent=win); return

                # equilibrium concentrations/pressures
                Ceq = {sp: C0[sp] + stoich[sp]*x for sp in stoich}
                # post Q (≈K)
                num = den = 1.0
                for sp,nu in stoich.items():
                    c = Ceq[sp]
                    if nu>0: num *= c**nu
                    elif nu<0: den *= c**(-nu)
                qeq = num/den

                # print ICE-like summary
                note("\nICE table summary (values in M or atm per your choice):")
                note("species      I         Δ(ν·x)      E")
                for sp in sorted(stoich.keys()):
                    I = C0[sp]; d = stoich[sp]*x; E = Ceq[sp]
                    note(f"{sp:>8s}  {I:>9.5g}  {d:>10.5g}  {E:>11.5g}")

                # direction
                direction = "→ products" if Q0 is not None and Q0 < K else ("→ reactants" if Q0 is not None and Q0 > K else "already at K")
                out.set(f"Equilibrium found: x = {x:.6g} ;  K(target)={K:.6g},  Q(eq)≈{qeq:.6g}  → shift {direction}")

                # copy to a small table below
                res = tk.Toplevel(win); res.title("Equilibrium values"); res.transient(win); res.geometry("360x280")
                tfrm = ttk.Frame(res, padding=10); tfrm.pack(fill=tk.BOTH, expand=True)
                ttk.Label(tfrm, text=f"Equilibrium ({kind.get()}):", font=("Segoe UI",10,"bold")).pack(anchor="w")
                for sp in sorted(Ceq.keys()):
                    ttk.Label(tfrm, text=f"{sp} : {Ceq[sp]:.6g}").pack(anchor="w")
                def add_known():
                    # store Ceq as [sp]_eq in Known (names "[X]_eq")
                    for sp,val in Ceq.items():
                        key = f"[{sp}]_eq" if kind.get()=="Kc" else f"P_{sp}_eq"
                        self.known_base[key] = float(val)
                        self.known_ui[key]   = {"unit":"M" if kind.get()=="Kc" else "atm", "display_value": float(val)}
                    self._refresh_known_table(); self._refresh_equation_list()
                    messagebox.showinfo("Saved", "Equilibrium values added to Known Variables.", parent=res)
                ttk.Button(tfrm, text="Add equilibrium values to Known", command=add_known).pack(anchor="w", pady=(8,0))

        # compute bar
        self._add_compute_bar(win, compute)

        # quick help
        def help_popup():
            h = tk.Toplevel(win); h.title("ICE Solver — Help"); h.transient(win); h.grab_set(); h.geometry("720x520")
            txt = tk.Text(h, wrap="word"); txt.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
            txt.insert("1.0",
                "Dynamic equilibrium & ICE\n"
                "• Enter a balanced reaction. Parser understands coefficients like '2 NH3'.\n"
                "• Choose Kc (use molarity) or Kp (use partial pressures).\n"
                "• Leave K blank to compute Q only (predict direction vs. K later).\n"
                "• With K entered, the solver uses a single extent x so [i]_eq = [i]_0 + ν_i x.\n"
                "• It enforces non-negative concentrations by solving within the feasible x interval.\n"
                "\nGas equilibria (Kp vs Kc):  Kp = Kc·(R·T)^{Δn_gas}.\n"
                "Pure solids & pure liquids do not appear in K (activity ≈ 1). Omit them from inputs.\n")
            txt.configure(state="disabled")
            ttk.Button(h, text="Close", command=h.destroy).pack(pady=6)
        ttk.Button(frm, text="❓ Help", command=help_popup).grid(row=2, column=5, sticky="e")

    def _open_kp_kc_converter(self):
        import math, tkinter as tk
        from tkinter import ttk, messagebox

        win = tk.Toplevel(self); win.title("Kp ↔ Kc Converter")
        win.transient(self); win.grab_set(); win.geometry("560x280")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        ttk.Label(frm, text="Δn(gas) = Σν(products) − Σν(reactants):").grid(row=0, column=0, sticky="w")
        dnE = ttk.Entry(frm, width=8); dnE.grid(row=0, column=1, sticky="w", padx=(6,12)); dnE.insert(0, "0")

        ttk.Label(frm, text="Temperature T:").grid(row=1, column=0, sticky="e")
        TE = ttk.Entry(frm, width=10); TE.grid(row=1, column=1, sticky="w", padx=(6,6)); TE.insert(0, "298.15")
        ttk.Label(frm, text="K").grid(row=1, column=2, sticky="w")

        ttk.Label(frm, text="Kp unit:").grid(row=1, column=3, sticky="e")
        punit = tk.StringVar(value="atm")
        ttk.Combobox(frm, width=6, state="readonly", textvariable=punit,
                    values=["atm","bar","Pa(m³/mol)"]).grid(row=1, column=4, sticky="w")

        ttk.Label(frm, text="Kc:").grid(row=2, column=0, sticky="e")
        KcE = ttk.Entry(frm, width=14); KcE.grid(row=2, column=1, sticky="w", padx=(6,0))
        ttk.Label(frm, text="Kp:").grid(row=3, column=0, sticky="e")
        KpE = ttk.Entry(frm, width=14); KpE.grid(row=3, column=1, sticky="w", padx=(6,0))

        out = tk.StringVar(value=""); ttk.Label(frm, textvariable=out, font=("Segoe UI",10,"bold")).grid(row=4, column=0, columnspan=5, sticky="w", pady=(8,0))

        def R_match():
            u = punit.get()
            if u == "atm": return 0.082057338  # L·atm/(mol·K), Kc in mol/L
            if u == "bar": return 0.083144626  # L·bar/(mol·K)
            # Pa option expects Kc in mol/m^3, Kp in Pa; then R must be 8.314462618 J/(mol·K)=Pa·m^3/(mol·K)
            return 8.314462618

        def compute():
            try:
                dn = float(dnE.get().replace(",", "."))
                T  = float(TE.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid", "Δn and T must be numeric.", parent=win); return
            R = R_match()
            factor = (R * T) ** dn
            powtxt = f"(RT)^{dn:g}"

            if KcE.get().strip():
                Kc = float(KcE.get().replace(",", "."))
                Kp = Kc * factor
                KpE.delete(0,"end"); KpE.insert(0, f"{Kp:.6g}")
                out.set(f"Kp = Kc·{powtxt} = {Kp:.6g}   (R matched to {punit.get()})")
            elif KpE.get().strip():
                Kp = float(KpE.get().replace(",", "."))
                if (R*T) == 0:
                    messagebox.showerror("Oops", "RT = 0.", parent=win); return
                Kc = Kp / factor
                KcE.delete(0,"end"); KcE.insert(0, f"{Kc:.6g}")
                out.set(f"Kc = Kp/{powtxt} = {Kc:.6g}   (R matched to {punit.get()})")
            else:
                messagebox.showerror("Need a value", "Enter either Kc or Kp.", parent=win)

        self._add_compute_bar(win, compute)

    def _open_equilibrium_combine_tool(self):
        win = tk.Toplevel(self); win.title("Combine Equilibria — overall K"); win.transient(self); win.grab_set(); win.geometry("620x420")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        ttk.Label(frm, text="Enter component equilibria:  K_overall = Π K_i^{ν_i}\nUse ν_i>0 to add (multiply reaction), ν_i<0 to subtract (reverse), fractional allowed."
                ).grid(row=0,column=0,sticky="w")

        box = ttk.Frame(frm); box.grid(row=1,column=0,sticky="nsew",pady=(6,6))
        frm.grid_rowconfigure(1, weight=1); frm.grid_columnconfigure(0, weight=1)

        rows = []
        def add_row():
            r = ttk.Frame(box); r.pack(fill=tk.X, pady=2)
            Ke = ttk.Entry(r, width=16); Ke.pack(side=tk.LEFT); Ke.insert(0, "1.0")
            ttk.Label(r, text="  exponent ν:").pack(side=tk.LEFT)
            ne = ttk.Entry(r, width=10); ne.pack(side=tk.LEFT); ne.insert(0, "1")
            rows.append((Ke, ne))
        for _ in range(3): add_row()
        ttk.Button(frm, text="+ Add row", command=add_row).grid(row=2,column=0,sticky="w")

        out = tk.StringVar(value=""); ttk.Label(frm, textvariable=out, font=("Segoe UI",10,"bold")).grid(row=3,column=0,sticky="w",pady=(8,0))

        def compute():
            import math
            if not rows:
                messagebox.showerror("Empty","Add at least one K row.", parent=win); return
            try:
                logK = 0.0
                for Ke, ne in rows:
                    K = float(Ke.get().replace(",", "."))
                    n = float(ne.get().replace(",", "."))
                    if K <= 0: 
                        messagebox.showerror("Bad K","K must be > 0.", parent=win); return
                    logK += n * math.log(K)
                Koverall = math.exp(logK)
                out.set(f"K_overall = {Koverall:.6g}")
            except ValueError:
                messagebox.showerror("Invalid","Use numeric K and ν.", parent=win)

        self._add_compute_bar(win, compute)

    def _open_q_direction_tool(self):
        win = tk.Toplevel(self); win.title("Q & Le Châtelier — direction helper"); win.transient(self); win.grab_set(); win.geometry("780x560")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        ttk.Label(frm, text="Reaction:").grid(row=0,column=0,sticky="e")
        rxn = ttk.Entry(frm, width=48); rxn.grid(row=0,column=1,columnspan=3,sticky="we",padx=(6,0))
        rxn.insert(0, "N2 + 3 H2 ⇌ 2 NH3")

        ttk.Label(frm, text="K value:").grid(row=1,column=0,sticky="e")
        KE = ttk.Entry(frm, width=14); KE.grid(row=1,column=1,sticky="w",padx=(6,0))

        box = ttk.LabelFrame(frm, text="Current mixture (use M for Kc or atm for Kp; omit pure s,l)")
        box.grid(row=2,column=0,columnspan=4,sticky="we",pady=(8,6))
        tbl = ttk.Frame(box); tbl.pack(fill=tk.X, pady=6)

        entries = {}; stoich = {}
        def refresh():
            for w in tbl.winfo_children(): w.destroy()
            try:
                nonlocal stoich
                stoich = self._parse_reaction(rxn.get())
            except Exception as ex:
                messagebox.showerror("Parse error", str(ex), parent=win); return
            ttk.Label(tbl, text="Species").grid(row=0,column=0,sticky="w")
            ttk.Label(tbl, text="Amount now").grid(row=0,column=1,sticky="w")
            for i,(sp,_) in enumerate(sorted(stoich.items()),start=1):
                ttk.Label(tbl, text=sp).grid(row=i,column=0,sticky="w")
                e = ttk.Entry(tbl, width=12); e.grid(row=i,column=1,sticky="w"); e.insert(0,"0.0")
                entries[sp]=e
        refresh()

        out = tk.StringVar(value=""); ttk.Label(frm, textvariable=out, font=("Segoe UI",10,"bold")).grid(row=3,column=0,columnspan=4,sticky="w",pady=(8,4))
        txt = tk.Text(frm, height=16, wrap="word"); txt.grid(row=4,column=0,columnspan=4,sticky="nsew"); frm.grid_rowconfigure(4, weight=1)

        def compute():
            txt.configure(state="normal"); txt.delete("1.0","end"); txt.configure(state="disabled")
            try:
                K = float(KE.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Need K","Enter a numeric K.", parent=win); return
            C = {}
            try:
                for sp,e in entries.items():
                    C[sp]=float(e.get().replace(",", "."))
            except ValueError:
                messagebox.showerror("Invalid","Numbers only.", parent=win); return
            num=den=1.0
            for sp,nu in stoich.items():
                if C[sp]<=0: out.set("Q undefined (zero/negative)"); return
                if nu>0: num*=C[sp]**nu
                elif nu<0: den*=C[sp]**(-nu)
            Q=num/den
            direct = "→ products" if Q<K else ("→ reactants" if Q>K else "at equilibrium")
            out.set(f"Q = {Q:.6g}; compare to K={K:.6g}  ⇒ shift {direct}")
            txt.configure(state="normal")
            txt.insert("end",
                "Le Châtelier quick guide:\n"
                "• Add product → system consumes product (shifts left). Add reactant → shifts right.\n"
                "• Reduce volume / increase pressure (gas): shifts to side with fewer gas moles (Δn_gas<0).\n"
                "• Add inert gas at constant V: no effect on equilibrium position. At constant P: expands volume → follow Δn.\n"
                "• Temperature change acts like adding product/reactant for endothermic/exothermic:\n"
                "    – Endothermic (ΔH>0): increase T → shift right (favours products).\n"
                "    – Exothermic (ΔH<0): increase T → shift left.\n"
                "• Catalyst: speeds approach to equilibrium but does not change K, Q, or position.\n")
            txt.configure(state="disabled")

        self._add_compute_bar(win, compute)
        ttk.Button(frm, text="Parse species again", command=refresh).grid(row=2,column=3,sticky="e")

    def _open_ksp_tool(self):
        import re, math, tkinter as tk
        from tkinter import ttk, messagebox
        win = tk.Toplevel(self); win.title("Ksp & Heterogeneous Equilibria")
        win.transient(self); win.grab_set(); win.geometry("760x620")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # --- Ksp quick solver ------------------------------------------------------
        box = ttk.LabelFrame(frm, text="Ksp — molar solubility of AxBy(s) ⇌ x A^{z+} + y B^{z−}")
        box.grid(row=0, column=0, sticky="we")
        for c in range(6): box.grid_columnconfigure(c, weight=1)

        ttk.Label(box, text="Salt formula (e.g. AgCl, CaF2; complex polyatomic salts not supported here):")\
            .grid(row=0, column=0, sticky="w", columnspan=4)
        fE = ttk.Entry(box, width=16); fE.grid(row=1, column=0, sticky="w"); fE.insert(0, "CaF2")
        ttk.Label(box, text="Ksp:").grid(row=1, column=1, sticky="e")
        KspE = ttk.Entry(box, width=14); KspE.grid(row=1, column=2, sticky="w")
        ttk.Label(box, text="or  s (mol/L):").grid(row=1, column=3, sticky="e")
        sE = ttk.Entry(box, width=14); sE.grid(row=1, column=4, sticky="w")

        out1 = tk.StringVar(value="")
        ttk.Label(box, textvariable=out1, font=("Segoe UI", 10, "bold"))\
            .grid(row=2, column=0, columnspan=5, sticky="w", pady=(6,0))

        def parse_axby(formula: str):
            """
            Return (x, y) for AxBy-type salts with two ionic groups.
            Supports: Element–Element (AgCl), Element–polyatomic (CaCO3, NaNO3),
            parentheses ((NH4)2SO4, Ca(OH)2, Fe(NO3)3), and element counts (Na2SO4).
            """
            import re
            f = re.sub(r"\s+", "", formula)
            if not f:
                raise ValueError("Enter a salt formula.")

            # --- first (cation) group ---
            if f.startswith("("):
                m = re.match(r"^\(([^)]+)\)(\d*)(.+)$", f)
                if not m: raise ValueError("Could not parse the first ion group.")
                x = int(m.group(2)) if m.group(2) else 1
                rest = m.group(3)
            elif f.startswith("NH4"):                      # common polyatomic cation without ()
                x, rest = 1, f[3:]
            else:
                m = re.match(r"^([A-Z][a-z]?)(\d*)(.+)$", f)
                if not m: raise ValueError("Enter a salt with two ionic groups, e.g. CaCO3, Ca(OH)2, Fe(NO3)3, Na2SO4, AgCl.")
                x = int(m.group(2)) if m.group(2) else 1
                rest = m.group(3)

            # --- second (anion) group ---
            rest = rest.strip()
            if rest.startswith("("):
                m = re.match(r"^\(([^)]+)\)(\d*)$", rest)
                if not m:
                    m = re.match(r"^\(([^)]+)\)(\d*)", rest)  # take first group if more follows
                    if not m: raise ValueError("Could not parse the second ion group.")
                y = int(m.group(2)) if m.group(2) else 1
            else:
                # Single element like Cl2 -> y=2; multi-element like CO3/NO3 (no ()) -> y=1
                if re.fullmatch(r"[A-Z][a-z]?\d*", rest):
                    m = re.fullmatch(r"[A-Z][a-z]?(\d*)", rest)
                    y = int(m.group(1)) if m.group(1) else 1
                else:
                    y = 1
            return x, y


        def compute_ksp_block():
            try:
                x, y = parse_axby(fE.get().strip())
            except Exception as ex:
                messagebox.showerror("Formula", str(ex), parent=win); return
            if KspE.get().strip():
                try:
                    Ksp = float(KspE.get().replace(",", "."))
                except ValueError:
                    messagebox.showerror("Ksp", "Numeric.", parent=win); return
                # Ksp = (x·s)^x (y·s)^y = (x^x y^y) s^{x+y}
                s = (Ksp / ((x**x) * (y**y))) ** (1.0 / (x + y))
                sE.delete(0, "end"); sE.insert(0, f"{s:.6g}")
                out1.set(f"molar solubility s = {s:.6g} mol/L   (x={x}, y={y})")
            elif sE.get().strip():
                try:
                    s = float(sE.get().replace(",", "."))
                except ValueError:
                    messagebox.showerror("s", "Numeric.", parent=win); return
                Ksp = ((x * s) ** x) * ((y * s) ** y)
                KspE.delete(0, "end"); KspE.insert(0, f"{Ksp:.6g}")
                out1.set(f"Ksp = (x·s)^x (y·s)^y = {Ksp:.6g}")
            else:
                messagebox.showerror("Need one", "Enter either Ksp or s.", parent=win)

        # Consistent bottom bar for this first block
        self._add_compute_bar(win, compute_ksp_block)  # original behavior kept. :contentReference[oaicite:3]{index=3}

        # --- NEW: Qsp checker (mass → M or direct M) -------------------------------
        boxQ = ttk.LabelFrame(frm, text="Qsp — ‘will it precipitate?’ checker")
        boxQ.grid(row=1, column=0, sticky="we", pady=(10,0))
        for c in range(8): boxQ.grid_columnconfigure(c, weight=1)

        ttk.Label(boxQ, text="Salt formula (AxBy):").grid(row=0, column=0, sticky="e")
        qfE = ttk.Entry(boxQ, width=16); qfE.grid(row=0, column=1, sticky="w")
        qfE.insert(0, "NaCl")

        ttk.Label(boxQ, text="Mass:").grid(row=0, column=2, sticky="e")
        mE = ttk.Entry(boxQ, width=12); mE.grid(row=0, column=3, sticky="w")
        mU = ttk.Combobox(boxQ, values=("g","mg","kg"), width=6, state="readonly"); mU.set("g")
        mU.grid(row=0, column=4, sticky="w")

        ttk.Label(boxQ, text="Volume:").grid(row=0, column=5, sticky="e")
        vE = ttk.Entry(boxQ, width=12); vE.grid(row=0, column=6, sticky="w")
        vU = ttk.Combobox(boxQ, values=("mL","L"), width=6, state="readonly"); vU.set("mL")
        vU.grid(row=0, column=7, sticky="w")

        ttk.Label(boxQ, text="(or) C (M):").grid(row=1, column=0, sticky="e", pady=(6,0))
        CE = ttk.Entry(boxQ, width=14); CE.grid(row=1, column=1, sticky="w", pady=(6,0))
        ttk.Label(boxQ, text="Ksp (blank → use above):").grid(row=1, column=2, sticky="e", pady=(6,0))
        KspQE = ttk.Entry(boxQ, width=14); KspQE.grid(row=1, column=3, sticky="w", pady=(6,0))

        outQ = tk.StringVar(value="")
        ttk.Label(boxQ, textvariable=outQ, font=("Segoe UI", 10, "bold"))\
            .grid(row=2, column=0, columnspan=8, sticky="w", pady=(8,4))

        def _num(s):
            s = (s or "").strip()
            if not s: return None
            try: return float(s.replace(",", "."))
            except Exception: return None

        def _mass_g(x, u):
            if u == "mg": return x * 1e-3
            if u == "kg": return x * 1e3
            return x

        def _vol_L(x, u):
            return x / 1000.0 if u.lower() == "ml" else x

        def compute_qsp():
            # Which formula to use for stoichiometry (Qsp uses same AxBy parser)
            ftxt = (qfE.get().strip() or fE.get().strip())
            try:
                x, y = parse_axby(ftxt)
            except Exception as ex:
                messagebox.showerror("Formula", str(ex), parent=win); return

            # Concentration C (M): either direct or mass/volume via MM
            C = _num(CE.get())
            steps = []

            if C is None:
                # Need mass, volume, and molar mass
                mg = _num(mE.get()); vg = _num(vE.get())
                if mg is None or vg is None:
                    messagebox.showerror("Need amount",
                                        "Enter either C (M) directly OR mass and volume.",
                                        parent=win); return
                # Use your existing helpers to get molar mass
                try:
                    f = name_to_formula(ftxt)
                    MM = molar_mass(f)  # g/mol
                except Exception as ex:
                    messagebox.showerror("Molar mass", str(ex), parent=win); return
                n = _mass_g(mg, mU.get()) / MM
                VL = _vol_L(vg, vU.get())
                if VL <= 0:
                    messagebox.showerror("Volume", "Volume must be > 0.", parent=win); return
                C = n / VL
                steps.append(f"MM({ftxt}) = {MM:.6g} g/mol;  n = m/MM = {_mass_g(mg, mU.get()):.6g}/{MM:.6g} = {n:.6g} mol")
                steps.append(f"C = n/V = {n:.6g}/{VL:.6g} L = {C:.6g} M")
            else:
                steps.append(f"C (given) = {C:.6g} M")

            # Ion product for AxBy: Qsp = (x·C)^x (y·C)^y
            Qsp = ((x * C) ** x) * ((y * C) ** y)
            steps.append(f"Qsp = (x·C)^x (y·C)^y = ({x:g}·{C:.6g})^{x:g}({y:g}·{C:.6g})^{y:g} = {Qsp:.6g}")

            # Ksp for verdict: from above or local box
            Ktxt = (KspQE.get().strip() or KspE.get().strip())
            verdict = ""
            if Ktxt:
                try:
                    Ksp_val = float(Ktxt.replace(",", "."))
                except ValueError:
                    messagebox.showerror("Ksp", "Ksp must be numeric.", parent=win); return
                # Compare with a soft tolerance
                tol = 1e-6 * max(1.0, Ksp_val)
                if Qsp < Ksp_val - tol:
                    verdict = "→ Solution is UNSATURATED (no precipitate)."
                elif abs(Qsp - Ksp_val) <= max(tol, 0.05*Ksp_val):
                    verdict = "→ ≈ at EQUILIBRIUM (saturated)."
                else:
                    verdict = "→ SUPERSATURATED — precipitation expected."
                outQ.set(f"C ≈ {C:.6g} M;  Qsp ≈ {Qsp:.6g};  Ksp ≈ {Ksp_val:.6g}  {verdict}")
            else:
                outQ.set(f"C ≈ {C:.6g} M;  Qsp ≈ {Qsp:.6g}  (Enter Ksp for a verdict.)")

        ttk.Button(boxQ, text="Check Qsp vs Ksp", command=compute_qsp)\
            .grid(row=1, column=4, columnspan=2, sticky="w", padx=(8,0), pady=(6,0))

        noteQ = ttk.Label(
            boxQ,
            text="Notes: assumes ideal behavior (activities≈concentrations); ignores complexation/common-ion effects.",
            foreground="#555", wraplength=700, justify="left"
        )
        noteQ.grid(row=3, column=0, columnspan=8, sticky="w", pady=(2,0))

        # --- Heterogeneous expression helper --------------------------------------
        box2 = ttk.LabelFrame(frm, text="Heterogeneous equilibrium expression builder")
        box2.grid(row=2, column=0, sticky="nsew", pady=(10,0))  # moved down by one row
        for c in range(5): box2.grid_columnconfigure(c, weight=1)

        ttk.Label(
            box2,
            text="Type a reaction and mark pure s,l to omit from K. Example:\nCaCO3(s) ⇌ Ca2+(aq) + CO2(g) + O2(g)"
        ).grid(row=0, column=0, columnspan=5, sticky="w")
        rxn = ttk.Entry(box2, width=64); rxn.grid(row=1, column=0, columnspan=5, sticky="we", pady=(4,6))
        out2 = tk.StringVar("")
        ttk.Label(box2, textvariable=out2, font=("Segoe UI", 10, "bold"))\
            .grid(row=2, column=0, columnspan=5, sticky="w")

        def build_K_expr():
            try:
                st = self._parse_reaction(rxn.get())
            except Exception as ex:
                messagebox.showerror("Parse", str(ex), parent=win); return
            include = []
            for sp, nu in st.items():
                sp_l = sp.lower()
                if sp_l.endswith("(s)") or sp_l.endswith("(l)") or sp_l.endswith("(pure)"):
                    continue
                include.append((sp, nu))
                if not include:
                    out2.set("All species were pure s/l; K = 1.")
                return
            num = " · ".join([f"[{sp}]^{int(nu) if abs(nu-int(nu))<1e-12 else nu:g}" for sp,nu in include if nu>0])
            den = " · ".join([f"[{sp}]^{int(-nu) if abs(-nu-int(-nu))<1e-12 else -nu:g}" for sp,nu in include if nu<0])
            expr = f"K = ({num}) / ({den})" if den else f"K = {num}"
            out2.set(expr)

        ttk.Button(box2, text="Build K expression", command=build_K_expr).grid(row=3, column=0, sticky="w", pady=(6,8))
    

    def _open_initial_rates_tool(self):
        import math
        win = tk.Toplevel(self)
        win.title("Initial Rates — determine orders m, n")
        win.transient(self); win.grab_set(); win.geometry("820x580")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # --- Header / labels --------------------------------------------------------
        ttk.Label(frm, text="Rate law form:  rate = k · [A]^m · [B]^n",
                font=("Segoe UI", 10, "bold")).grid(row=0, column=0, columnspan=8, sticky="w", pady=(0,8))

        ttk.Label(frm, text="Name A:").grid(row=1, column=0, sticky="e")
        nameA = ttk.Entry(frm, width=12); nameA.grid(row=1, column=1, sticky="w", padx=(6,12)); nameA.insert(0, "A")

        ttk.Label(frm, text="Name B:").grid(row=1, column=2, sticky="e")
        nameB = ttk.Entry(frm, width=12); nameB.grid(row=1, column=3, sticky="w", padx=(6,12)); nameB.insert(0, "B")

        log10_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(frm, text="Use ln (natural log)  [uncheck for log10 display only]",
                        variable=log10_var).grid(row=1, column=4, columnspan=3, sticky="w")

        # --- Table of experiments ---------------------------------------------------
        # 6 rows is usually plenty; blank rows are ignored
        ttk.Label(frm, text=f"[A] ({nameA.get()})").grid(row=2, column=1, sticky="w")
        ttk.Label(frm, text=f"[B] ({nameB.get()})").grid(row=2, column=2, sticky="w")
        ttk.Label(frm, text="rate").grid(row=2, column=3, sticky="w")
        ttk.Label(frm, text="(units consistent but arbitrary)").grid(row=2, column=4, sticky="w")

        rows = []
        for i in range(6):
            r = 3 + i
            ttk.Label(frm, text=f"Exp {i+1}:").grid(row=r, column=0, sticky="e")
            aE = ttk.Entry(frm, width=12); aE.grid(row=r, column=1, sticky="w", padx=(6,6))
            bE = ttk.Entry(frm, width=12); bE.grid(row=r, column=2, sticky="w", padx=(6,6))
            vE = ttk.Entry(frm, width=12); vE.grid(row=r, column=3, sticky="w", padx=(6,6))
            rows.append((aE, bE, vE))

        # Small hint
        ttk.Label(frm, text="Leave [B] blank for single-reactant fits. Use strictly positive values.",
                foreground="#555").grid(row=9, column=0, columnspan=8, sticky="w", pady=(6,8))

        # --- Output widgets ---------------------------------------------------------
        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 10, "bold")).grid(row=10, column=0, columnspan=8, sticky="w", pady=(4,2))

        steps = tk.Text(frm, height=12, wrap="word")
        steps.grid(row=11, column=0, columnspan=8, sticky="nsew")
        frm.grid_rowconfigure(11, weight=1)
        scy = ttk.Scrollbar(frm, orient="vertical", command=steps.yview); steps.configure(yscrollcommand=scy.set)
        scy.grid(row=11, column=8, sticky="ns")

        # --- tiny linear solver (2 or 3 parameters) --------------------------------
        def solve_normal_eq(XTX, XTy):
            # Gaussian elimination for 2x2 or 3x3
            n = len(XTy)
            # build augmented matrix
            A = [XTX[i][:] + [XTy[i]] for i in range(n)]
            # elimination
            for i in range(n):
                # pivot
                pivot = i
                for r in range(i+1, n):
                    if abs(A[r][i]) > abs(A[pivot][i]): pivot = r
                if abs(A[pivot][i]) < 1e-14:
                    return None  # singular
                if pivot != i:
                    A[i], A[pivot] = A[pivot], A[i]
                # normalize
                fac = A[i][i]
                for c in range(i, n+1):
                    A[i][c] /= fac
                # eliminate
                for r in range(n):
                    if r == i: continue
                    fac = A[r][i]
                    for c in range(i, n+1):
                        A[r][c] -= fac * A[i][c]
            return [A[i][n] for i in range(n)]  # solution

        # --- core compute -----------------------------------------------------------
        def compute():
            steps.configure(state="normal"); steps.delete("1.0", "end"); steps.configure(state="disabled")
            Avals, Bvals, Rvals = [], [], []
            # gather rows
            for aE, bE, vE in rows:
                a = aE.get().strip(); b = bE.get().strip(); v = vE.get().strip()
                if not v:  # no rate → ignore whole row
                    continue
                try:
                    rate = float(v.replace(",", "."))
                except ValueError:
                    messagebox.showerror("Invalid rate", "Use numbers only for rate."); return
                if rate <= 0:
                    messagebox.showerror("Invalid rate", "Rate must be > 0 for logarithms."); return
                # A (required)
                if not a:
                    messagebox.showerror("Missing [A]", "Every used row needs [A] and rate."); return
                try:
                    Aval = float(a.replace(",", "."))
                except ValueError:
                    messagebox.showerror("Invalid [A]", "Use numbers for [A]."); return
                if Aval <= 0:
                    messagebox.showerror("Invalid [A]", "[A] must be > 0 for logarithms."); return
                # B (optional per row; treat blank as None)
                Bval = None
                if b:
                    try:
                        Bval = float(b.replace(",", "."))
                    except ValueError:
                        messagebox.showerror("Invalid [B]", "Use numbers for [B] or leave blank."); return
                    if Bval <= 0:
                        messagebox.showerror("Invalid [B]", "[B] must be > 0 for logarithms."); return
                Avals.append(Aval); Bvals.append(Bval); Rvals.append(rate)

            if len(Rvals) < 2:
                messagebox.showerror("Not enough data", "Enter at least two experiments."); return

            # Decide if B participates (at least two valid B and variance nonzero)
            Bs_present = [b for b in Bvals if b is not None]
            use_B = len(Bs_present) >= 2 and (max(Bs_present) - min(Bs_present) > 1e-12)

            # Build design matrix for ln form
            ln = math.log if log10_var.get() else math.log  # both are ln here; checkbox is just a label toggle
            X = []; y = []
            for A, B, r in zip(Avals, Bvals, Rvals):
                xi = [1.0, ln(A)]
                if use_B:
                    if B is None:
                        # if B column is used globally, missing B on a row is not allowed
                        messagebox.showerror("Missing [B]", "Some rows have [B]; all used rows must have [B] to fit m and n.")
                        return
                    xi.append(ln(B))
                X.append(xi); y.append(ln(r))

            p = len(X[0])  # 2 (intercept + lnA) or 3 (intercept + lnA + lnB)
            if len(X) < p:
                messagebox.showerror("Not enough data", f"Need at least {p} experiments for this fit."); return

            # Normal equations: (X^T X) beta = X^T y
            # compute Gram matrix and RHS
            XTX = [[0.0]*p for _ in range(p)]
            XTy = [0.0]*p
            for xi, yi in zip(X, y):
                for i in range(p):
                    XTy[i] += xi[i]*yi
                    for j in range(p):
                        XTX[i][j] += xi[i]*xi[j]

            beta = solve_normal_eq(XTX, XTy)
            if beta is None:
                messagebox.showerror("Singular system",
                                    "Data are collinear (e.g., no variation in A or B). Vary concentrations and try again.")
                return

            b0 = beta[0]
            m  = beta[1]
            n  = beta[2] if use_B else 0.0
            k  = math.exp(b0)  # because we used natural logs

            # Goodness of fit
            yhat = []
            for xi in X:
                yh = sum(b*xx for b,xx in zip(beta, xi))
                yhat.append(yh)
            ss_tot = sum((yi - sum(y)/len(y))**2 for yi in y)
            ss_res = sum((yi - yh)**2 for yi, yh in zip(y, yhat))
            R2 = 1.0 - (ss_res/ss_tot) if ss_tot > 0 else float("nan")

            # Output
            Aname = (nameA.get().strip() or "A")
            Bname = (nameB.get().strip() or "B")
            law = f"rate = {k:.6g} · [{Aname}]^{m:.4g}" + (f" · [{Bname}]^{n:.4g}" if use_B else "")
            out.set(law + f"    (R² = {R2:.4f})")

            txt = []
            txt.append("Fit details (log-linear least squares):\n")
            txt.append("  model: ln(rate) = ln k + m·ln[A] " + ("+ n·ln[B]\n" if use_B else "\n"))
            txt.append(f"  ln k = {b0:.6g}  →  k = {k:.6g}\n")
            txt.append(f"  m    = {m:.6g}\n")
            if use_B: txt.append(f"  n    = {n:.6g}\n")
            txt.append(f"  R²   = {R2:.6g}\n\n")
            if use_B: txt.append(f"  Reactionorder is givin: x+y\n")
            if use_B: txt.append(f"  {m + n}\n\n")
            txt.append("Notes:\n")
            txt.append("• Units must be consistent across experiments, but k’s units are whatever make rate’s units match.\n")
            txt.append("• Rows with blanks are ignored. Values must be > 0 to take logarithms.\n")
            if not use_B:
                txt.append("• [B] column left blank or not varied → fit uses single-reactant form (n = 0).\n")

            steps.configure(state="normal"); steps.insert("1.0", "".join(txt)); steps.configure(state="disabled")

            # Save result on the object so the Add button can use it
            win._fit_result = {"k": float(k), "m": float(m), "n": float(n), "use_B": bool(use_B), "Aname": Aname, "Bname": Bname}
            return win._fit_result

        self._add_compute_bar(win, compute)

        # --- Add-to-known -----------------------------------------------------------
        btns = ttk.Frame(frm); btns.grid(row=12, column=0, columnspan=8, sticky="w", pady=(8,0))
        def add_known():
            res = getattr(win, "_fit_result", None) or compute()
            if not res: return
            # store dimensionless m, n; store k under 'k_rate' (display unit left blank)
            self.known_base["m"] = res["m"]; self.known_ui["m"] = {"unit":"", "display_value":res["m"]}
            self.known_base["n"] = res["n"]; self.known_ui["n"] = {"unit":"", "display_value":res["n"]}
            self.known_base["k_rate"] = res["k"]; self.known_ui["k_rate"] = {"unit":"", "display_value":res["k"]}
            self._refresh_known_table(); self._refresh_equation_list()
            messagebox.showinfo("Saved", f"Saved m={res['m']:.4g}, n={res['n']:.4g}, k={res['k']:.6g} to Known Variables.", parent=win)

        ttk.Button(btns, text="Add m, n, k to Known Variables", command=add_known).pack(side=tk.LEFT)

    def _open_k_expression_builder(parent, on_done=None):
        import math, re
        """
        Open a modal window to build an equilibrium-constant expression.
        If `on_done` is provided, it is called with the final expression string.

        New features
        ------------
        • Equation parser (optional): type an equation like
            H2O(l) ⇌ H3O+(aq) + OH-(aq)
        or
            HF(aq) + H2O(l) ⇌ F-(aq) + H3O+(aq)
        then click Parse → tables auto-fill.
        • Smart helper: if the equation is water autoprotolysis, a Kw/pH panel appears:
            - enter any one of: [H+], [OH-], pH, pOH (and Kw if not 1.0e-14)
            - it computes the rest, plus fraction ionized and "1 in N" water molecules.

        Existing behavior
        -----------------
        • Choose Kc or Kp.
        • "Exclude pure liquids and solids" = activity 1 (ignored).
        • Fill Reactants (left) and Products (right).
            - Species: name (e.g., HF, H3O+, F-)
            - Coeff: stoichiometric coefficient (defaults to 1)
            - Phase: aq, g, l, s
            - Value: optional numeric value
                • For Kc: concentration in M
                • For Kp: partial pressure (same units for all gases)
        • Click Compute to generate K; Copy expression puts the symbolic K on clipboard.
        """
        win = tk.Toplevel(parent)
        win.title("Equilibrium — K-expression builder")
        win.transient(parent)
        win.grab_set()

        # --- Top controls --------------------------------------------------------
        header = ttk.Label(
            win,
            text="Build the equilibrium constant for:  Products / Reactants, ignoring pure liquids and solids by default."
        )
        header.grid(row=0, column=0, columnspan=2, sticky="w", pady=(6,8))

        mode_var = tk.StringVar(value="Kc")
        frm_mode = ttk.Frame(win)
        frm_mode.grid(row=1, column=0, columnspan=2, sticky="w", pady=(0,6))
        ttk.Radiobutton(frm_mode, text="Kc", variable=mode_var, value="Kc").grid(row=0, column=0, padx=(0,12))
        ttk.Radiobutton(frm_mode, text="Kp (use partial pressures P(X))", variable=mode_var, value="Kp").grid(row=0, column=1)

        exclude_pure_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frm_mode, text="Exclude pure liquids and solids (activity = 1)", variable=exclude_pure_var)\
            .grid(row=0, column=2, padx=(16,0))

        # --- Equation parser row -------------------------------------------------
        eq_frame = ttk.Frame(win)
        eq_frame.grid(row=2, column=0, columnspan=2, sticky="we", pady=(4,2))
        ttk.Label(eq_frame, text="Equation (optional):").grid(row=0, column=0, sticky="w", padx=(0,6))
        eq_var = tk.StringVar()
        eq_entry = ttk.Entry(eq_frame, textvariable=eq_var, width=80)
        eq_entry.grid(row=0, column=1, sticky="we")
        eq_frame.grid_columnconfigure(1, weight=1)
        parse_msg = ttk.Label(eq_frame, text="", foreground="#666")
        parse_msg.grid(row=1, column=1, sticky="w", pady=(2,0))

        # --- Tables --------------------------------------------------------------
        def make_table(parent_frame):
            headers = ["Species", "Coeff", "Phase", "Value (optional)"]
            for j, h in enumerate(headers):
                ttk.Label(parent_frame, text=h).grid(row=0, column=j, sticky="w")
            rows = []
            for i in range(6):
                r = i + 1
                e_sp = ttk.Entry(parent_frame, width=22)
                e_sp.grid(row=r, column=0, padx=(0,8), pady=2, sticky="w")

                e_cf = ttk.Entry(parent_frame, width=6)
                e_cf.insert(0, "1")
                e_cf.grid(row=r, column=1, padx=(0,8), pady=2, sticky="w")

                cb_ph = ttk.Combobox(parent_frame, values=["aq","g","l","s"], width=5, state="readonly")
                cb_ph.set("aq")
                cb_ph.grid(row=r, column=2, padx=(0,8), pady=2, sticky="w")

                e_val = ttk.Entry(parent_frame, width=14)
                e_val.grid(row=r, column=3, pady=2, sticky="w")

                rows.append((e_sp, e_cf, cb_ph, e_val))
            return rows

        left = ttk.LabelFrame(win, text="Reactants")
        left.grid(row=3, column=0, sticky="nw", padx=(0,18), pady=(4,0))
        right = ttk.LabelFrame(win, text="Products")
        right.grid(row=3, column=1, sticky="nw", pady=(4,0))

        react_rows = make_table(left)
        prod_rows = make_table(right)

        win.grid_columnconfigure(0, weight=1)
        win.grid_columnconfigure(1, weight=1)

        # --- Output areas --------------------------------------------------------
        out_lbl = ttk.Label(win, text="K expression:")
        out_lbl.grid(row=4, column=0, columnspan=2, sticky="w", pady=(10,2))

        out_expr = tk.Text(win, height=3, width=100, wrap="word")
        out_expr.grid(row=5, column=0, columnspan=2, sticky="we")
        out_expr.configure(font=("TkDefaultFont", 10))

        steps = tk.Text(win, height=10, width=100, wrap="word")
        steps.grid(row=6, column=0, columnspan=2, sticky="nsew", pady=(8,0))
        steps.configure(font=("TkDefaultFont", 9))
        win.grid_rowconfigure(6, weight=1)

        # --- Smart Helper: Kw / pH ----------------------------------------------
        kw_frame = ttk.LabelFrame(win, text="Kw / pH helper (auto-shows for water autoprotolysis)")
        # initially not gridded; call _show_kw_helper(True) to show

        def _show_kw_helper(show=True):
            if show:
                kw_frame.grid(row=7, column=0, columnspan=2, sticky="we", pady=(10,4))
            else:
                kw_frame.grid_forget()

        # inputs
        rowi = 0
        ttk.Label(kw_frame, text="Kw (@25°C default 1.0e-14):").grid(row=rowi, column=0, sticky="w")
        kw_var = tk.StringVar(value="1.0e-14")
        ttk.Entry(kw_frame, textvariable=kw_var, width=16).grid(row=rowi, column=1, sticky="w", padx=(4,12))
        ttk.Label(kw_frame, text="T (°C) [optional]:").grid(row=rowi, column=2, sticky="e")
        t_var = tk.StringVar()
        ttk.Entry(kw_frame, textvariable=t_var, width=10).grid(row=rowi, column=3, sticky="w", padx=(4,0))
        rowi += 1

        ttk.Label(kw_frame, text="[H+] (M):").grid(row=rowi, column=0, sticky="w")
        H_var = tk.StringVar()
        ttk.Entry(kw_frame, textvariable=H_var, width=16).grid(row=rowi, column=1, sticky="w", padx=(4,12))
        ttk.Label(kw_frame, text="[OH-] (M):").grid(row=rowi, column=2, sticky="e")
        OH_var = tk.StringVar()
        ttk.Entry(kw_frame, textvariable=OH_var, width=16).grid(row=rowi, column=3, sticky="w", padx=(4,0))
        rowi += 1

        ttk.Label(kw_frame, text="pH:").grid(row=rowi, column=0, sticky="w")
        pH_var = tk.StringVar()
        ttk.Entry(kw_frame, textvariable=pH_var, width=16).grid(row=rowi, column=1, sticky="w", padx=(4,12))
        ttk.Label(kw_frame, text="pOH:").grid(row=rowi, column=2, sticky="e")
        pOH_var = tk.StringVar()
        ttk.Entry(kw_frame, textvariable=pOH_var, width=16).grid(row=rowi, column=3, sticky="w", padx=(4,0))
        rowi += 1

        kw_info = ttk.Label(kw_frame, text="", foreground="#555")
        kw_info.grid(row=rowi, column=0, columnspan=4, sticky="w", pady=(6,4))
        rowi += 1

        def _parse_float(s):
            s = (s or "").strip()
            if not s: return None
            try:
                return float(s)
            except Exception:
                try:
                    return float(s.replace("×", "e").replace("·", "").replace("^", "**"))
                except Exception:
                    return None

        def solve_kw():
            Kw = _parse_float(kw_var.get()) or 1.0e-14
            H = _parse_float(H_var.get())
            OH = _parse_float(OH_var.get())
            pH = _parse_float(pH_var.get())
            pOH = _parse_float(pOH_var.get())

            # choose a starting value
            if H is None and OH is None and pH is None and pOH is None:
                H = math.sqrt(Kw)
            if H is None and pH is not None:
                H = 10.0**(-pH)
            if OH is None and pOH is not None:
                OH = 10.0**(-pOH)
            if H is None and OH is not None:
                H = Kw / OH
            if OH is None and H is not None:
                OH = Kw / H

            # now compute pH/pOH
            if H is not None and H > 0:
                pH = -math.log10(H)
            if OH is not None and OH > 0:
                pOH = -math.log10(OH)

            # display back
            H_var.set("" if H is None else f"{H:.6g}")
            OH_var.set("" if OH is None else f"{OH:.6g}")
            pH_var.set("" if pH is None else f"{pH:.4f}")
            pOH_var.set("" if pOH is None else f"{pOH:.4f}")

            # fraction ionized of water
            info_txt = []
            info_txt.append(f"Using Kw = {Kw:.3g}.")
            if H is not None:
                try:
                    alpha = H / 55.5  # 55.5 M water in pure water
                    if alpha > 0:
                        N = 1.0 / alpha
                        info_txt.append(f"Fraction of water ionized ≈ {alpha:.3e}  (≈ 1 in {N:.2e} molecules).")
                except Exception:
                    pass
            if (t_var.get() or "").strip():
                info_txt.append("Note: T entered does not auto-adjust Kw; supply temperature-correct Kw if needed.")
            kw_info.config(text="  ".join(info_txt))

        ttk.Button(kw_frame, text="Solve Kw / pH", command=solve_kw).grid(row=rowi, column=0, sticky="w", pady=(2,8))
        rowi += 1

        # --- Helpers / internal logic -------------------------------------------
        def parse_coeff(s):
            try:
                v = float(s)
                return v
            except Exception:
                return 1.0

        def parse_value(s):
            if not s.strip():
                return None
            try:
                return float(s)
            except Exception:
                return None

        def collect(rows):
            """Return list of dicts for non-empty species rows."""
            items = []
            for e_sp, e_cf, cb_ph, e_val in rows:
                name = e_sp.get().strip()
                if not name:
                    continue
                items.append({
                    "name": name,
                    "coeff": parse_coeff(e_cf.get()),
                    "phase": cb_ph.get(),
                    "val": parse_value(e_val.get())
                })
            return items

        def term_symbol(name, coeff, phase, mode):
            if mode == "Kp":
                base = f"P({name})"
            else:  # Kc
                base = f"[{name}]"
            expo = "" if abs(coeff) == 1 else f"^{int(abs(coeff)) if abs(coeff).is_integer() else f'^{abs(coeff)}'}"
            return base + expo

        def include_in_K(item, mode, exclude_pure):
            ph = item["phase"]
            if mode == "Kp":
                return ph == "g"  # only gases in Kp
            if ph in ("aq", "g"):
                return True
            if ph in ("l", "s"):
                return not exclude_pure
            return True

        # --- Equation string parser ---------------------------------------------
        ARROW_RX = re.compile(r"(⇌|↔|⟺|<=>|<=|=>|=|→|<-*->)")

        def parse_species_token(tok):
            """Return dict with coeff, core_formula (no phase), phase, charge, atoms."""
            s = tok.strip()

            # leading stoichiometric coefficient
            m = re.match(r"^(\d+)\s*(.*)$", s)
            coeff = 1
            if m:
                coeff = int(m.group(1))
                s = m.group(2).strip()

            # strip phase from the END first so the charge is truly at the end
            phase = None
            m = re.search(r"\((s|l|g|aq)\)\s*$", s, flags=re.I)
            if m:
                phase = m.group(1).lower()
                s = s[:m.start()].strip()  # remove '(aq)' etc. before reading charge

            # read charge (now at end if present)
            core, charge = _read_charge(s)   # supports Pt^2+, Pt2+, Ag+, etc.

            # tokenize core to count atoms
            atoms = _tokenize(_strip_phase_and_ws(core))
            return {"coeff": coeff, "core": core, "phase": phase, "charge": charge, "atoms": atoms}



        def _fill_table(rows, items):
            for i, (e_sp, e_cf, cb_ph, e_val) in enumerate(rows):
                if i < len(items):
                    it = items[i]
                    e_sp.delete(0, "end"); e_sp.insert(0, it["name"])
                    e_cf.delete(0, "end"); e_cf.insert(0, str(it["coeff"]))
                    cb_ph.set(it["phase"])
                    e_val.delete(0, "end")
                else:
                    e_sp.delete(0, "end")
                    e_cf.delete(0, "end"); e_cf.insert(0, "1")
                    cb_ph.set("aq")
                    e_val.delete(0, "end")

        def _names(side):
            return [re.sub(r"\s+", "", d["name"]).upper() for d in side]

        def _is_water_autoprotolysis(reacts, prods):
            r = _names(reacts); p = _names(prods)
            has_H2O = any(x == "H2O" for x in r)
            has_H3O = any(x == "H3O+" for x in p) or any(x == "H+" for x in p)
            has_OH  = any(x in ("OH-", "HO-") for x in p)
            return has_H2O and has_H3O and has_OH

        def parse_equation_into_tables():
            s = eq_var.get().strip()
            if not s:
                parse_msg.config(text="Type an equation first.")
                return
            parts = ARROW_RX.split(s, maxsplit=1)
            if len(parts) < 3:
                parse_msg.config(text="Could not find an equilibrium arrow (⇌, ↔, <=>, =>, =, →).")
                return
            lhs, _, rhs = parts[0], parts[1], parts[2]
            react_tokens = [t for t in lhs.split('+') if t.strip()]
            prod_tokens  = [t for t in rhs.split('+') if t.strip()]
            reacts = [_parse_species_token(t) for t in react_tokens]
            prods  = [_parse_species_token(t) for t in prod_tokens]

            if len(reacts) > len(react_rows) or len(prods) > len(prod_rows):
                messagebox.showwarning("Too many species",
                                    "This builder shows up to 6 species per side. Extra items will be ignored.")
            _fill_table(react_rows, reacts[:len(react_rows)])
            _fill_table(prod_rows,  prods[:len(prod_rows)])

            # Smart helper detection
            if _is_water_autoprotolysis(reacts, prods):
                _show_kw_helper(True)
                parse_msg.config(text="Recognized water autoprotolysis: Kw/pH helper opened.")
            else:
                _show_kw_helper(False)
                parse_msg.config(text="Equation parsed into tables.")

        ttk.Button(eq_frame, text="Parse", command=parse_equation_into_tables)\
            .grid(row=0, column=2, padx=(8,0), sticky="w")

        # --- Compute / Copy actions ---------------------------------------------
        def build_expression():
            mode = mode_var.get()
            exclude_pure = exclude_pure_var.get()

            prods = collect(prod_rows)
            reacts = collect(react_rows)

            if not prods or not reacts:
                messagebox.showwarning("Missing data", "Enter at least one Reactant and one Product.")
                return

            num_syms, den_syms = [], []
            num_val, den_val = 1.0, 1.0
            have_full_numeric = True  # becomes False if any included term lacks a numeric value

            step_lines = []
            step_lines.append(f"Mode: {mode}  |  Exclude pure l/s: {exclude_pure}")
            step_lines.append("Included species and exponents:")

            # Products (positive exponents)
            for it in prods:
                if include_in_K(it, mode, exclude_pure):
                    num_syms.append(term_symbol(it["name"], it["coeff"], it["phase"], mode))
                    if it["val"] is not None:
                        num_val *= (it["val"] ** it["coeff"])
                    else:
                        have_full_numeric = False
                    step_lines.append(f"  + {it['name']} ({it['phase']}): exponent {it['coeff']}"
                                    + (f", value {it['val']}" if it['val'] is not None else ", value —"))
                else:
                    step_lines.append(f"  · {it['name']} ({it['phase']}) excluded")

            # Reactants (negative exponents)
            for it in reacts:
                if include_in_K(it, mode, exclude_pure):
                    den_syms.append(term_symbol(it["name"], it["coeff"], it["phase"], mode))
                    if it["val"] is not None:
                        den_val *= (it["val"] ** it["coeff"])
                    else:
                        have_full_numeric = False
                    step_lines.append(f"  − {it['name']} ({it['phase']}): exponent {it['coeff']}"
                                    + (f", value {it['val']}" if it['val'] is not None else ", value —"))
                else:
                    step_lines.append(f"  · {it['name']} ({it['phase']}) excluded")

            num_str = " · ".join(num_syms) if num_syms else "1"
            den_str = " · ".join(den_syms) if den_syms else "1"

            if den_str == "1":
                expr = f"{mode} = {num_str}"
            else:
                expr = f"{mode} = ({num_str}) / ({den_str})"

            out_expr.config(state="normal")
            out_expr.delete("1.0", "end")
            out_expr.insert("1.0", expr + "\n")
            if have_full_numeric:
                try:
                    k_val = num_val / den_val if den_val != 0 else float("inf")
                    out_expr.insert("end", f"{mode} ≈ {k_val:g}\n")
                    step_lines.append(f"\nNumeric evaluation: numerator={num_val:g}, denominator={den_val:g}, {mode}={k_val:g}")
                except Exception as e:
                    step_lines.append(f"\nNumeric evaluation failed: {e}")
            else:
                step_lines.append("\nNumeric K not computed (one or more included species missing a value).")
            out_expr.config(state="disabled")

            steps.config(state="normal")
            steps.delete("1.0", "end")
            steps.insert("1.0", "\n".join(step_lines))
            steps.config(state="disabled")

            # Callback for host app if provided
            if on_done is not None:
                on_done(expr)

        def copy_expression():
            text = out_expr.get("1.0", "end").strip().splitlines()
            if not text:
                return
            expr_line = text[0]
            win.clipboard_clear()
            win.clipboard_append(expr_line)
            messagebox.showinfo("Copied", "Symbolic K expression copied to clipboard.")

        # --- Buttons -------------------------------------------------------------
        btns = ttk.Frame(win)
        btns.grid(row=8, column=0, columnspan=2, sticky="w", pady=(10,10))
        ttk.Button(btns, text="Compute", command=build_expression).grid(row=0, column=0, padx=(0,8))
        ttk.Button(btns, text="Copy expression", command=copy_expression).grid(row=0, column=1)
        ttk.Button(btns, text="Show/Hide Kw helper", command=lambda: _show_kw_helper(not kw_frame.winfo_ismapped()))\
            .grid(row=0, column=2, padx=(8,0))

        # Close behavior
        def _close(*_):
            win.grab_release()
            win.destroy()
        win.protocol("WM_DELETE_WINDOW", _close)
        win.bind("<Escape>", _close)

        # Initial focus
        eq_entry.focus_set()


    def _open_ph_modal(parent):
        """
        pH / pOH calculator (modal).
        - Enter ANY one (or more) of: [H+], [OH-], pH, pOH and Kw (defaults to 1.0e-14 @25 °C).
        - Click Solve → computes all related quantities.
        - Strong-acid/base shortcut: give concentration and stoichiometric factor (n),
        it pre-fills [H+] = n*C (acid) or [OH-] = n*C (base) before Solve.
        - Quick parser: paste a question sentence like "[H+] = 0.0030 M" or "pH = 3.52".
        """

        win = tk.Toplevel(parent)
        win.title("pH / pOH calculator")
        win.transient(parent)
        win.grab_set()

        # ------------------------------------------------------------------ header
        ttk.Label(win, text="Compute pH, pOH, [H+], [OH-] from any one input. Kw default is 1.0e-14 (25 °C).")\
            .grid(row=0, column=0, columnspan=4, sticky="w", pady=(8,6))
        # Kw + temperature note
        ttk.Label(win, text="Kw:").grid(row=1, column=0, sticky="e")
        kw_var = tk.StringVar(value="1.0e-14")
        ttk.Entry(win, textvariable=kw_var, width=16).grid(row=1, column=1, sticky="w")
            # ------------------------------------------------------------------ Kw helper
        ttk.Label(win, text="Temperature (°C):").grid(row=1, column=3, sticky="e")
        T_var = tk.StringVar()
        ttk.Entry(win, textvariable=T_var, width=10).grid(row=1, column=4, sticky="w")

        def temp_to_kw():
            try:
                T = float(T_var.get())
            except:
                messagebox.showwarning("Invalid input", "Enter a valid temperature in °C.")
                return

            # rough Kw lookup table (°C: Kw)
            table = {
                0: 1.1e-15,
                10: 2.9e-15,
                25: 1.0e-14,
                50: 5.5e-14,
                100: 5.1e-13
            }

            # exact match
            if T in table:
                kw_val = table[T]
            else:
                # linear interpolate between closest two points
                keys = sorted(table.keys())
                for i in range(len(keys)-1):
                    if keys[i] <= T <= keys[i+1]:
                        T1, T2 = keys[i], keys[i+1]
                        Kw1, Kw2 = table[T1], table[T2]
                        # log-linear interpolation (Kw ~ exponential with T)
                        logKw = ( (T2-T)*math.log10(Kw1) + (T-T1)*math.log10(Kw2) ) / (T2-T1)
                        kw_val = 10**logKw
                        break
                else:
                    messagebox.showwarning("Out of range", "Temperature must be between 0 and 100 °C.")
                    return

            kw_var.set(f"{kw_val:.2e}")

        ttk.Button(win, text="Convert", command=temp_to_kw).grid(row=1, column=5, padx=(6,0))

        ttk.Label(win, text="(enter a different Kw for other temperatures)").grid(row=1, column=2, columnspan=2, sticky="w")

        # ------------------------------------------------------------------ inputs
        ttk.Separator(win).grid(row=2, column=0, columnspan=4, sticky="we", pady=(8,8))

        ttk.Label(win, text="[H+] (M):").grid(row=3, column=0, sticky="e")
        H_var = tk.StringVar()
        ttk.Entry(win, textvariable=H_var, width=18).grid(row=3, column=1, sticky="w")

        ttk.Label(win, text="[OH-] (M):").grid(row=3, column=2, sticky="e")
        OH_var = tk.StringVar()
        ttk.Entry(win, textvariable=OH_var, width=18).grid(row=3, column=3, sticky="w")

        ttk.Label(win, text="pH:").grid(row=4, column=0, sticky="e")
        pH_var = tk.StringVar()
        ttk.Entry(win, textvariable=pH_var, width=18).grid(row=4, column=1, sticky="w")

        ttk.Label(win, text="pOH:").grid(row=4, column=2, sticky="e")
        pOH_var = tk.StringVar()
        ttk.Entry(win, textvariable=pOH_var, width=18).grid(row=4, column=3, sticky="w")

        # --------------------------------------------------- strong acid/base aid
        aid = ttk.LabelFrame(win, text="Strong acid/base shortcut (optional)")
        aid.grid(row=5, column=0, columnspan=4, sticky="we", pady=(10,0))
        kind_var = tk.StringVar(value="acid")
        ttk.Radiobutton(aid, text="Strong acid → [H+] = n·C", variable=kind_var, value="acid").grid(row=0, column=0, sticky="w")
        ttk.Radiobutton(aid, text="Strong base → [OH-] = n·C", variable=kind_var, value="base").grid(row=0, column=1, sticky="w")
        ttk.Label(aid, text="C (M):").grid(row=1, column=0, sticky="e", padx=(0,4))
        C_var = tk.StringVar()
        ttk.Entry(aid, textvariable=C_var, width=14).grid(row=1, column=1, sticky="w")
        ttk.Label(aid, text="n (stoichiometric factor):").grid(row=1, column=2, sticky="e", padx=(8,4))
        n_var = tk.StringVar(value="1")
        ttk.Entry(aid, textvariable=n_var, width=10).grid(row=1, column=3, sticky="w")

        def apply_strong():
            def _pf(s):
                s = (s or "").strip()
                if not s: return None
                try:
                    return float(s)
                except Exception:
                    try:
                        return float(s.replace("×","e").replace("^","**"))
                    except Exception:
                        return None
            C = _pf(C_var.get())
            n = _pf(n_var.get()) or 1.0
            if C is None or C < 0:
                messagebox.showwarning("Input needed", "Enter a valid concentration C.")
                return
            if kind_var.get()=="acid":
                H_var.set(f"{(n*C):.6g}")
            else:
                OH_var.set(f"{(n*C):.6g}")
        ttk.Button(aid, text="Apply", command=apply_strong).grid(row=1, column=4, padx=(8,0))

        # ---------------------------------------------------------- text parser
        par = ttk.LabelFrame(win, text="Quick text parser (paste a sentence and Detect)")
        par.grid(row=6, column=0, columnspan=4, sticky="we", pady=(10,0))
        text_var = tk.StringVar()
        ttk.Entry(par, textvariable=text_var, width=92).grid(row=0, column=0, columnspan=3, sticky="we")
        detect_msg = ttk.Label(par, text="", foreground="#666")
        detect_msg.grid(row=1, column=0, columnspan=3, sticky="w")

        def detect_values():
            s = (text_var.get() or "").strip()
            if not s:
                detect_msg.config(text="Paste something like: [H+] = 0.0030 M  or  pH = 3.52")
                return
            # patterns for [H+], [OH-], pH, pOH
            def _grab(rx):
                m = re.search(rx, s, flags=re.I)
                if m: 
                    raw = m.group(1)
                    raw = raw.replace("×10^","e").replace("×10","e").replace("×","e").replace("^","**")
                    try:
                        return float(raw)
                    except Exception:
                        try:
                            return float(eval(raw))
                        except Exception:
                            return None
                return None
            H = _grab(r"\[?\s*H\+\s*\]?\s*=\s*([0-9eE\.\-+×^]+)")
            OH = _grab(r"\[?\s*OH-\s*\]?\s*=\s*([0-9eE\.\-+×^]+)")
            pH = _grab(r"\bpH\s*=\s*([0-9eE\.\-+×^]+)")
            pOH = _grab(r"\bpOH\s*=\s*([0-9eE\.\-+×^]+)")
            taken = []
            if H is not None: H_var.set(f"{H:.6g}"); taken.append("[H+]")
            if OH is not None: OH_var.set(f"{OH:.6g}"); taken.append("[OH-]")
            if pH is not None: pH_var.set(f"{pH:.6g}"); taken.append("pH")
            if pOH is not None: pOH_var.set(f"{pOH:.6g}"); taken.append("pOH")
            detect_msg.config(text=("Detected: " + ", ".join(taken)) if taken else "No numbers detected.")
        ttk.Button(par, text="Detect", command=detect_values).grid(row=0, column=3, padx=(8,0), sticky="w")

        # --------------------------------------------------------------- results
        ttk.Separator(win).grid(row=7, column=0, columnspan=4, sticky="we", pady=(10,6))
        out = tk.Text(win, height=9, width=100, wrap="word")
        out.grid(row=8, column=0, columnspan=4, sticky="nsew")
        out.configure(font=("TkDefaultFont", 10))
        win.grid_rowconfigure(8, weight=1)

        # ------------------------------------------------------------ utilities
        def _pf(s):
            s = (s or "").strip()
            if not s: return None
            try:
                return float(s)
            except Exception:
                try:
                    # tolerate "3.0×10^-3" etc.
                    s2 = s.replace("×10^","e").replace("×10","e").replace("×","e").replace("^","**")
                    return float(eval(s2))
                except Exception:
                    return None

        # --------------------------------------------------------------- solver
        def solve():
            Kw = _pf(kw_var.get()) or 1.0e-14
            H  = _pf(H_var.get())
            OH = _pf(OH_var.get())
            pH = _pf(pH_var.get())
            pOH= _pf(pOH_var.get())

            steps = []
            steps.append(f"Using Kw = {Kw:.6g}")

            # derive missing items
            if H is None and pH is not None:
                H = 10.0**(-pH); steps.append(f"H from pH: [H+] = 10^(-pH) = {H:.6g}")
            if OH is None and pOH is not None:
                OH = 10.0**(-pOH); steps.append(f"OH from pOH: [OH-] = 10^(-pOH) = {OH:.6g}")
            if H is None and OH is not None:
                H = Kw / OH; steps.append(f"H from Kw/[OH-]: [H+] = Kw/[OH-] = {H:.6g}")
            if OH is None and H is not None:
                OH = Kw / H; steps.append(f"OH from Kw/[H+]: [OH-] = Kw/[H+] = {OH:.6g}")
            if pH is None and H is not None and H>0:
                pH = -math.log10(H); steps.append(f"pH from [H+]: pH = -log10([H+]) = {pH:.4f}")
            if pOH is None and OH is not None and OH>0:
                pOH = -math.log10(OH); steps.append(f"pOH from [OH-]: pOH = -log10([OH-]) = {pOH:.4f}")

            # sanity & classification
            flavor = ""
            if pH is not None:
                # Neutral pH when pH = 0.5 * pKw
                pKw = -math.log10(Kw)
                neutral_pH = 0.5 * pKw
                if abs(pH - neutral_pH) < 1e-6:
                    flavor = "neutral"
                elif pH < neutral_pH:
                    flavor = "acidic"
                else:
                    flavor = "basic"
                steps.append(f"pKw = {pKw:.4f} → neutral pH = {neutral_pH:.4f} (solution is {flavor}).")

            # fraction of water ionized (only a fun fact)
            if H is not None and H>0:
                try:
                    alpha = H / 55.5
                    N = (1.0/alpha) if alpha>0 else float("inf")
                    steps.append(f"Water ionization fraction ≈ {alpha:.3e} (≈ 1 in {N:.2e} molecules).")
                except Exception:
                    pass

            # show results
            out.config(state="normal"); out.delete("1.0","end")
            def _fmt(v, d=6): return "" if v is None else f"{v:.6g}"
            def _fmtp(v): return "" if v is None else f"{v:.4f}"
            out.insert("end", f"[H+] = {_fmt(H)} M\n[OH-] = {_fmt(OH)} M\npH = {_fmtp(pH)}\npOH = {_fmtp(pOH)}\n")
            out.insert("end", "\nSteps / notes:\n" + "\n".join(steps))
            out.config(state="disabled")

            # push the canonical numbers back to inputs
            if H is not None:   H_var.set(f"{H:.6g}")
            if OH is not None: OH_var.set(f"{OH:.6g}")
            if pH is not None: pH_var.set(f"{pH:.4f}")
            if pOH is not None:pOH_var.set(f"{pOH:.4f}")

        def clear_all():
            for v in (H_var, OH_var, pH_var, pOH_var, C_var, n_var, text_var):
                v.set("")
            kw_var.set("1.0e-14")
            out.config(state="normal"); out.delete("1.0","end"); out.config(state="disabled")

        # --------------------------------------------------------------- buttons
        btns = ttk.Frame(win); btns.grid(row=9, column=0, columnspan=4, sticky="w", pady=(8,10))
        ttk.Button(btns, text="Solve", command=solve).grid(row=0, column=0, padx=(0,8))
        ttk.Button(btns, text="Clear", command=clear_all).grid(row=0, column=1)

        # close behavior
        def _close(*_):
            win.grab_release(); win.destroy()
        win.protocol("WM_DELETE_WINDOW", _close)
        win.bind("<Escape>", _close)

        # focus
        H_var.set("0.0030")  # example-friendly; remove if you prefer blank start

    def _open_weak_acid_base_modal(parent):
        """
        Weak acid/base modal — solve pH from Ka/Kb *or* solve Ka/Kb from pH (ICE).
        + Buffer (Henderson–Hasselbalch) block with Ka/pKa/Kb/pKb selector.
        + Salt helper (mass → C0) for salts of weak acids/bases.
        Accepts forgiving inputs like: 0.30 M, 5.6×10^-4, 25 °C, 3.50, 3,0e-3, etc.
        """
        import math, tkinter as tk
        from tkinter import ttk, messagebox

        win = tk.Toplevel(parent)
        win.title("Weak acid/base — pH or Ka/Kb (ICE)")
        win.transient(parent); win.grab_set()

        # ----------------------------- tolerant parser -----------------------------
        def _pf(s: str):
            s = (s or "").strip()
            if not s:
                return None
            s2 = (s.replace(",", ".")
                    .replace("×10^","e").replace("x10^","e")
                    .replace("×10","e").replace("x10","e")
                    .replace("×","e").replace("^","**")
                    .replace("°C","").replace("°c","")
                    .replace(" mol/L","").replace("mol/L","")
                    .replace(" M","").replace("m",""))
            try:
                return float(s2)
            except Exception:
                try:
                    return float(eval(s2))
                except Exception:
                    return None

        def _kw_from_T_table(T):
            """Kw via small log-linear table (0–100 °C)."""
            table = {0:1.1e-15, 10:2.9e-15, 25:1.0e-14, 50:5.5e-14, 100:5.1e-13}
            keys = sorted(table)
            if T in table: return table[T]
            if T < keys[0] or T > keys[-1]: return None
            for i in range(len(keys)-1):
                T1, T2 = keys[i], keys[i+1]
                if T1 <= T <= T2:
                    logKw = ((T2-T)*math.log10(table[T1]) + (T-T1)*math.log10(table[T2]))/(T2-T1)
                    return 10**logKw

        # -------------------------------- Kw + T row --------------------------------
        ttk.Label(win, text="Kw:").grid(row=0, column=0, sticky="e", padx=(6,2), pady=(8,4))
        kw_var = tk.StringVar(value="1.0e-14")
        ttk.Entry(win, textvariable=kw_var, width=16).grid(row=0, column=1, sticky="w", pady=(8,4))

        ttk.Label(win, text="Temperature (°C):").grid(row=0, column=2, sticky="e", padx=(12,2))
        T_var = tk.StringVar(value="25")
        ttk.Entry(win, textvariable=T_var, width=10).grid(row=0, column=3, sticky="w")

        def convert_kw():
            T = _pf(T_var.get())
            if T is None:
                messagebox.showwarning("Invalid", "Enter a valid temperature in °C.")
                return
            Kw = _kw_from_T_table(T)
            if Kw is None:
                messagebox.showwarning("Out of range", "T must be between 0 and 100 °C.")
                return
            kw_var.set(f"{Kw:.2e}")
        ttk.Button(win, text="Convert Kw(T)", command=convert_kw)\
            .grid(row=0, column=4, sticky="w", padx=(8,6))

        ttk.Separator(win).grid(row=1, column=0, columnspan=5, sticky="we", pady=(6,8))

        # --------------------------------- Inputs (ICE) -----------------------------
        frm = ttk.LabelFrame(win, text="ICE: pH from Ka/Kb  ⇄  Ka/Kb from pH")
        frm.grid(row=2, column=0, columnspan=5, sticky="we", padx=6)

        mode_var = tk.StringVar(value="acid")
        ttk.Radiobutton(frm, text="Weak acid (Ka)", variable=mode_var, value="acid").grid(row=0, column=0, padx=(6,8), pady=6, sticky="w")
        ttk.Radiobutton(frm, text="Weak base (Kb)", variable=mode_var, value="base").grid(row=0, column=1, padx=(0,8), pady=6, sticky="w")

        solve_var = tk.StringVar(value="ph")  # 'ph' or 'K'
        ttk.Radiobutton(frm, text="Solve for pH",   variable=solve_var, value="ph").grid(row=0, column=2, padx=(12,8), pady=6, sticky="w")
        ttk.Radiobutton(frm, text="Solve for Ka/Kb", variable=solve_var, value="K").grid(row=0, column=3, padx=(0,8),  pady=6, sticky="w")

        # C0
        ttk.Label(frm, text="Formal concentration C₀ (M):").grid(row=1, column=0, sticky="e", padx=(6,4))
        C0_var = tk.StringVar(value="")
        ttk.Entry(frm, textvariable=C0_var, width=18).grid(row=1, column=1, sticky="w", padx=(0,8), pady=4)

        # Ka/Kb with dropdown kind
        ttk.Label(frm, text="Equilibrium constant:").grid(row=1, column=2, sticky="e")
        K_var = tk.StringVar(value="")
        ttk.Entry(frm, textvariable=K_var, width=18).grid(row=1, column=3, sticky="w", padx=(4,6))
        const_kind = tk.StringVar(value="Ka")  # Ka, pKa, Kb, pKb
        ttk.Combobox(frm, textvariable=const_kind, values=("Ka","pKa","Kb","pKb"),
                    state="readonly", width=6).grid(row=1, column=4, sticky="w", padx=(0,6))

        approx_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(frm, text="Show √(K·C₀) approximation when solving pH", variable=approx_var)\
            .grid(row=2, column=0, columnspan=5, sticky="w", padx=6, pady=(0,6))

        # ----------------------- Salt helper (mass → C0) ---------------------------
        salt = ttk.LabelFrame(win, text="Salt helper (mass → C₀) — salts of weak acids/bases")
        salt.grid(row=3, column=0, columnspan=5, sticky="we", padx=6)

        salt_type = tk.StringVar(value="A-")  # A- (salt of HA) or BH+ (salt of B)
        ttk.Radiobutton(salt, text="A⁻ (salt of weak acid HA → behaves as weak base)", variable=salt_type, value="A-")\
            .grid(row=0, column=0, columnspan=3, sticky="w", padx=6, pady=(6,2))
        ttk.Radiobutton(salt, text="BH⁺ (salt of weak base B → behaves as weak acid)", variable=salt_type, value="BH+")\
            .grid(row=1, column=0, columnspan=3, sticky="w", padx=6)

        ttk.Label(salt, text="Mass of salt (g):").grid(row=0, column=3, sticky="e", padx=(8,4))
        mass_var = tk.StringVar(); ttk.Entry(salt, textvariable=mass_var, width=12).grid(row=0, column=4, sticky="w")
        ttk.Label(salt, text="Molar mass (g/mol):").grid(row=0, column=5, sticky="e", padx=(8,4))
        mw_var = tk.StringVar(); ttk.Entry(salt, textvariable=mw_var, width=12).grid(row=0, column=6, sticky="w")
        ttk.Label(salt, text="Solution volume:").grid(row=1, column=3, sticky="e", padx=(8,4))
        vol_var = tk.StringVar(); ttk.Entry(salt, textvariable=vol_var, width=12).grid(row=1, column=4, sticky="w")
        vol_unit = tk.StringVar(value="mL")
        ttk.Combobox(salt, textvariable=vol_unit, values=("mL","L"), state="readonly", width=5)\
            .grid(row=1, column=5, sticky="w", padx=(4,0))

        ttk.Label(salt, text="Parent constant (for HA or B):").grid(row=2, column=0, sticky="e", padx=(6,4), pady=(6,6))
        parentK_var = tk.StringVar(); ttk.Entry(salt, textvariable=parentK_var, width=14).grid(row=2, column=1, sticky="w", pady=(6,6))
        parent_kind = tk.StringVar(value="pKa")
        parent_box = ttk.Combobox(salt, textvariable=parent_kind, values=("pKa","Ka"), state="readonly", width=6)
        parent_box.grid(row=2, column=2, sticky="w", pady=(6,6))

        def _sync_parent_kind(*_):
            if salt_type.get() == "A-":
                parent_box.configure(values=("pKa","Ka"))
                if parent_kind.get() not in ("pKa","Ka"): parent_kind.set("pKa")
            else:
                parent_box.configure(values=("pKb","Kb"))
                if parent_kind.get() not in ("pKb","Kb"): parent_kind.set("pKb")
        _sync_parent_kind()
        salt_type.trace_add("write", _sync_parent_kind)

        def apply_salt():
            m = _pf(mass_var.get()); MW = _pf(mw_var.get()); V = _pf(vol_var.get())
            if any(v is None or v <= 0 for v in (m, MW, V)):
                messagebox.showwarning("Salt helper", "Enter mass, molar mass, and volume (>0).")
                return
            VL = V/1000.0 if vol_unit.get().lower()=="ml" else V
            C0 = (m/MW)/VL
            C0_var.set(f"{C0:.6g}")

            Kw = _pf(kw_var.get()) or 1.0e-14
            pKpar = _pf(parentK_var.get())
            if pKpar is None:
                messagebox.showwarning("Salt helper", "Enter the parent constant (pKa/Ka or pKb/Kb).")
                return

            if salt_type.get()=="A-":
                Ka = 10.0**(-pKpar) if parent_kind.get()=="pKa" else pKpar
                Kb = Kw/Ka
                mode_var.set("base"); const_kind.set("Kb"); K_var.set(f"{Kb:.6g}"); solve_var.set("ph")
            else:
                Kb = 10.0**(-pKpar) if parent_kind.get()=="pKb" else pKpar
                Ka = Kw/Kb
                mode_var.set("acid"); const_kind.set("Ka"); K_var.set(f"{Ka:.6g}"); solve_var.set("ph")

            messagebox.showinfo("Salt helper",
                                "Applied to Inputs:\n"
                                f"C₀ = {C0:.6g} M\n"
                                f"Mode = {'Weak base (Kb)' if salt_type.get()=='A-' else 'Weak acid (Ka)'}\n"
                                f"K = {K_var.get()} ({const_kind.get()})")
        ttk.Button(salt, text="Apply to Inputs", command=apply_salt)\
            .grid(row=2, column=3, columnspan=2, sticky="w", padx=(8,0), pady=(4,6))

        # --------------------- Buffer (Henderson–Hasselbalch) ----------------------
        buf = ttk.LabelFrame(win, text="Buffer (Henderson–Hasselbalch)")
        buf.grid(row=4, column=0, columnspan=5, sticky="we", padx=6, pady=(4,0))

        buf_type = tk.StringVar(value="acid")  # 'acid' → HA/A-, 'base' → B/BH+
        ttk.Radiobutton(buf, text="HA/A⁻ buffer (use Ka or pKa)", variable=buf_type, value="acid")\
            .grid(row=0, column=0, columnspan=2, sticky="w", padx=6, pady=(6,2))
        ttk.Radiobutton(buf, text="B/BH⁺ buffer (use Kb or pKb)", variable=buf_type, value="base")\
            .grid(row=0, column=2, columnspan=2, sticky="w", padx=6, pady=(6,2))

        ttk.Label(buf, text="Constant:").grid(row=1, column=0, sticky="e", padx=(6,4))
        bufK_var = tk.StringVar()
        ttk.Entry(buf, textvariable=bufK_var, width=14).grid(row=1, column=1, sticky="w")
        bufK_kind = tk.StringVar(value="pKa")
        bufK_box = ttk.Combobox(buf, textvariable=bufK_kind, values=("pKa","Ka"), state="readonly", width=6)
        bufK_box.grid(row=1, column=2, sticky="w", padx=(4,6))

        def _sync_buf_kind(*_):
            if buf_type.get()=="acid":
                bufK_box.configure(values=("pKa","Ka"))
                if bufK_kind.get() not in ("pKa","Ka"): bufK_kind.set("pKa")
            else:
                bufK_box.configure(values=("pKb","Kb"))
                if bufK_kind.get() not in ("pKb","Kb"): bufK_kind.set("pKb")
        _sync_buf_kind()
        buf_type.trace_add("write", _sync_buf_kind)

        # concentrations
        ttk.Label(buf, text="[HA] (M):").grid(row=2, column=0, sticky="e", padx=(6,4))
        ha_var = tk.StringVar(); ttk.Entry(buf, textvariable=ha_var, width=12).grid(row=2, column=1, sticky="w")
        ttk.Label(buf, text="[A⁻] (M):").grid(row=2, column=2, sticky="e", padx=(6,4))
        a_var = tk.StringVar();  ttk.Entry(buf, textvariable=a_var,  width=12).grid(row=2, column=3, sticky="w")

        ttk.Label(buf, text="[B] (M):").grid(row=3, column=0, sticky="e", padx=(6,4))
        b_var = tk.StringVar();  ttk.Entry(buf, textvariable=b_var,  width=12).grid(row=3, column=1, sticky="w")
        ttk.Label(buf, text="[BH⁺] (M):").grid(row=3, column=2, sticky="e", padx=(6,4))
        bh_var = tk.StringVar(); ttk.Entry(buf, textvariable=bh_var, width=12).grid(row=3, column=3, sticky="w")

        def compute_buffer():
            Kw = _pf(kw_var.get()) or 1.0e-14
            pKw = -math.log10(Kw)

            # constant
            K_in = _pf(bufK_var.get())
            if K_in is None or K_in <= 0:
                messagebox.showwarning("Buffer", "Enter Ka/Kb or pKa/pKb.")
                return

            if buf_type.get()=="acid":
                Ka = 10.0**(-K_in) if bufK_kind.get()=="pKa" else K_in
                pKa = -math.log10(Ka)
                HA = _pf(ha_var.get()); A = _pf(a_var.get())
                if HA is None or A is None or HA<=0 or A<=0:
                    messagebox.showwarning("Buffer", "Enter [HA] and [A⁻] (>0).")
                    return
                pH = pKa + math.log10(A/HA)
                msg = (f"Buffer (acid): pH = pKa + log([A⁻]/[HA])\n"
                    f"pKa = {pKa:.4f}, [A⁻]/[HA] = {A/HA:.4g} → pH = {pH:.4f}")
            else:
                Kb = 10.0**(-K_in) if bufK_kind.get()=="pKb" else K_in
                pKb = -math.log10(Kb)
                B = _pf(b_var.get()); BH = _pf(bh_var.get())
                if B is None or BH is None or B<=0 or BH<=0:
                    messagebox.showwarning("Buffer", "Enter [B] and [BH⁺] (>0).")
                    return
                pOH = pKb + math.log10(B/BH)
                pH  = pKw - pOH
                msg = (f"Buffer (base): pOH = pKb + log([B]/[BH⁺])\n"
                    f"pKb = {pKb:.4f}, [B]/[BH⁺] = {B/BH:.4g} → pOH = {pOH:.4f},  pH = {pH:.4f}")

            out.config(state="normal"); out.insert("end", "\n" + msg + "\n"); out.config(state="disabled")
            # Push pH back into ICE section as a convenience when user wants to continue
            # (does not change C0/K there; this is only a readout).
        ttk.Button(buf, text="Compute buffer pH", command=compute_buffer)\
            .grid(row=4, column=0, columnspan=2, sticky="w", padx=6, pady=(4,6))

        ttk.Separator(win).grid(row=5, column=0, columnspan=5, sticky="we", pady=(6,6))

        # --- Buffer recipe helper (put this right under your Buffer block) -------------
        recipe = ttk.LabelFrame(win, text="Buffer recipe helper (fills [HA]/[A⁻])")
        recipe.grid(row=5, column=0, columnspan=5, sticky="we", padx=6, pady=(4,0))

        rec_mode = tk.StringVar(value="HA + A-")
        ttk.Label(recipe, text="Recipe:").grid(row=0, column=0, sticky="e", padx=(6,4))
        ttk.Combobox(recipe, textvariable=rec_mode,
                    values=("HA + A− (acid + acetate salt)",
                            "A− + HCl (neutralize conjugate base)"),
                    state="readonly", width=28).grid(row=0, column=1, columnspan=3, sticky="w")

        # Inputs used by the two recipes (unused fields can be left blank)
        ttk.Label(recipe, text="Acid C (M):").grid(row=1, column=0, sticky="e", padx=(6,4))
        rec_Cacid = tk.StringVar(); ttk.Entry(recipe, textvariable=rec_Cacid, width=10).grid(row=1, column=1, sticky="w")
        ttk.Label(recipe, text="Acid V (mL):").grid(row=1, column=2, sticky="e", padx=(6,4))
        rec_Vacid = tk.StringVar(); ttk.Entry(recipe, textvariable=rec_Vacid, width=10).grid(row=1, column=3, sticky="w")

        ttk.Label(recipe, text="Salt mass (g):").grid(row=2, column=0, sticky="e", padx=(6,4))
        rec_msg = tk.StringVar(); ttk.Entry(recipe, textvariable=rec_msg, width=10).grid(row=2, column=1, sticky="w")
        ttk.Label(recipe, text="Salt MW (g/mol):").grid(row=2, column=2, sticky="e", padx=(6,4))
        rec_mw  = tk.StringVar(); ttk.Entry(recipe, textvariable=rec_mw, width=10).grid(row=2, column=3, sticky="w")

        ttk.Label(recipe, text="Strong acid C (M):").grid(row=3, column=0, sticky="e", padx=(6,4))
        rec_Csa = tk.StringVar(); ttk.Entry(recipe, textvariable=rec_Csa, width=10).grid(row=3, column=1, sticky="w")
        ttk.Label(recipe, text="Strong acid V (mL):").grid(row=3, column=2, sticky="e", padx=(6,4))
        rec_Vsa = tk.StringVar(); ttk.Entry(recipe, textvariable=rec_Vsa, width=10).grid(row=3, column=3, sticky="w")

        ttk.Label(recipe, text="Final volume (mL):").grid(row=4, column=0, sticky="e", padx=(6,4))
        rec_Vfin = tk.StringVar(); ttk.Entry(recipe, textvariable=rec_Vfin, width=10).grid(row=4, column=1, sticky="w")

        def _pf(s):
            s = (s or "").strip()
            if not s: return None
            s2 = (s.replace(",", ".")
                    .replace("×10^","e").replace("×10","e").replace("×","e").replace("^","**"))
            try: return float(s2)
            except Exception:
                try: return float(eval(s2))
                except Exception: return None

        def fill_buffer_from_recipe():
            mode = rec_mode.get()
            VL = _pf(rec_Vfin.get())
            if VL is None or VL <= 0:
                messagebox.showwarning("Recipe", "Enter final volume (mL) > 0.")
                return
            VL /= 1000.0  # liters

            # output placeholders
            HA = A = None

            if mode.startswith("HA + A−"):
                Cacid = _pf(rec_Cacid.get()); Vac = _pf(rec_Vacid.get())
                msg = _pf(rec_msg.get()); MW = _pf(rec_mw.get())
                if None in (Cacid, Vac, msg, MW):
                    messagebox.showwarning("Recipe", "Fill Acid C/V and Salt mass & MW.")
                    return
                n_HA = Cacid * (Vac/1000.0)
                n_A  = (msg / MW)
                HA, A = n_HA/VL, n_A/VL

            else:  # "A− + HCl"
                msg = _pf(rec_msg.get()); MW = _pf(rec_mw.get())
                Csa = _pf(rec_Csa.get()); Vsa = _pf(rec_Vsa.get())
                if None in (msg, MW, Csa, Vsa):
                    messagebox.showwarning("Recipe", "Fill Salt mass & MW and strong acid C/V.")
                    return
                n_A0 = msg / MW
                n_H  = Csa * (Vsa/1000.0)
                # A− + H+ → HA
                n_A  = max(n_A0 - n_H, 0.0)
                n_HA = min(n_H, n_A0)
                HA, A = n_HA/VL, n_A/VL

            # push into the Buffer block (acid-type buffer)
            buf_type.set("acid")
            ha_var.set(f"{HA:.6g}")
            a_var.set(f"{A:.6g}")
            messagebox.showinfo("Recipe", f"Filled Buffer fields:\n[HA] = {HA:.6g} M\n[A⁻] = {A:.6g} M\nSelect pKa/pKb above and click 'Compute buffer pH'.")

        ttk.Button(recipe, text="Fill Buffer fields", command=fill_buffer_from_recipe)\
            .grid(row=4, column=2, columnspan=2, sticky="w", padx=(8,0))

        # --------------------------------- Output ----------------------------------
        out = tk.Text(win, height=12, width=96, wrap="word")
        out.grid(row=6, column=0, columnspan=5, sticky="nsew", padx=6, pady=(0,6))
        win.grid_rowconfigure(6, weight=1)
        out.configure(font=("TkDefaultFont", 10))

        # ------------------------------ constant reader ----------------------------
        def _read_K(mode, rawK, kind, Kw):
            """Return numeric Ka (mode='acid') or Kb (mode='base') from the entry + dropdown."""
            K = _pf(rawK)
            if K is None: return None
            if mode == "acid":
                if kind == "Ka":   return K
                if kind == "pKa":  return 10.0**(-K)
                if kind == "Kb":   return Kw / K
                if kind == "pKb":  return Kw / (10.0**(-K))
            else:
                if kind == "Kb":   return K
                if kind == "pKb":  return 10.0**(-K)
                if kind == "Ka":   return Kw / K
                if kind == "pKa":  return Kw / (10.0**(-K))
            return None

        # --------------------------------- compute (ICE) ---------------------------
        def compute_ice():
            C0 = _pf(C0_var.get())
            Kw = _pf(kw_var.get()) or 1.0e-14
            mode = mode_var.get()
            goal = solve_var.get()

            if C0 is None or C0 <= 0:
                messagebox.showwarning("Input", "Enter a valid C₀ > 0.")
                return

            lines = [f"Kw = {Kw:.6g}",
                    f"Mode: {'Weak acid (Ka)' if mode=='acid' else 'Weak base (Kb)'}",
                    f"Goal: {'pH' if goal=='ph' else ('Ka' if mode=='acid' else 'Kb')}",
                    f"C₀ = {C0:.6g} M"]

            if goal == "ph":
                Knum = _read_K(mode, K_var.get(), const_kind.get(), Kw)
                if Knum is None or Knum <= 0:
                    messagebox.showwarning("Input", "Enter a valid Ka/Kb (choose Ka/pKa/Kb/pKb from the dropdown).")
                    return
                lines.append(f"K (numeric for this mode) = {Knum:.6g}\n")

                # exact quadratic: x = (-K + sqrt(K^2 + 4KC0))/2
                disc = Knum*Knum + 4*Knum*C0
                x = (-Knum + math.sqrt(disc)) / 2.0

                if mode == "acid":
                    H  = x; OH = Kw / H
                    pH  = -math.log10(H); pOH = -math.log10(OH)
                    lines += ["Reaction: HA ⇌ H⁺ + A⁻",
                            f"x = [-Ka + √(Ka² + 4KaC₀)]/2 = {x:.6g} M",
                            f"[H⁺] = {H:.6g} M,  [OH⁻] = {OH:.6g} M",
                            f"pH = {pH:.4f},  pOH = {pOH:.4f}",
                            f"% ionization = {100*x/C0:.3g}%"]
                    if approx_var.get():
                        x_est = math.sqrt(Knum*C0)
                        err = abs((x - x_est)/x)*100 if x>0 else float("inf")
                        lines += [f"Approx: √(Ka·C₀) = {x_est:.6g} M  (error ≈ {err:.3g}%)"]
                    lines += [f"Conjugate: Kb = Kw/Ka = {Kw/Knum:.6g}"]
                else:
                    OH = x; H = Kw / OH
                    pOH = -math.log10(OH); pH = -math.log10(H)
                    lines += ["Reaction: B + H₂O ⇌ BH⁺ + OH⁻",
                            f"x = [-Kb + √(Kb² + 4KbC₀)]/2 = {x:.6g} M",
                            f"[OH⁻] = {OH:.6g} M,  [H⁺] = {H:.6g} M",
                            f"pOH = {pOH:.4f},  pH = {pH:.4f}",
                            f"% ionization = {100*x/C0:.3g}%"]
                    if approx_var.get():
                        x_est = math.sqrt(Knum*C0)
                        err = abs((x - x_est)/x)*100 if x>0 else float("inf")
                        lines += [f"Approx: √(Kb·C₀) = {x_est:.6g} M  (error ≈ {err:.3g}%)"]
                    lines += [f"Conjugate: Ka = Kw/Kb = {Kw/Knum:.6g}"]

            else:  # goal == "K"
                messagebox.showinfo("Ka/Kb from pH", "For Ka/Kb from pH, give C₀ and measured pH/[H⁺]/[OH⁻] in your pH modal, then compute Ka or Kb with the ICE relations.\n(Kept minimal here to keep the UI compact.)")
                return

            out.config(state="normal"); out.delete("1.0","end"); out.insert("1.0", "\n".join(lines)); out.config(state="disabled")

        # -------------------------------- buttons ----------------------------------
        btns = ttk.Frame(win); btns.grid(row=7, column=0, columnspan=5, sticky="w", padx=6, pady=(0,8))
        ttk.Button(btns, text="Solve (ICE)", command=compute_ice).grid(row=0, column=0, padx=(0,8))
        ttk.Button(btns, text="Apply to Inputs", command=apply_salt).grid(row=0, column=1, padx=(0,8))

        def clear_all():
            for v in (C0_var, K_var, kw_var, T_var,
                    mass_var, mw_var, vol_var, parentK_var,
                    bufK_var, ha_var, a_var, b_var, bh_var):
                v.set("")
            kw_var.set("1.0e-14"); const_kind.set("Ka"); parent_kind.set("pKa"); bufK_kind.set("pKa")
            mode_var.set("acid"); solve_var.set("ph"); buf_type.set("acid")
            out.config(state="normal"); out.delete("1.0","end"); out.config(state="disabled")
        ttk.Button(btns, text="Clear", command=clear_all).grid(row=0, column=2)

        def _close(*_):
            win.grab_release(); win.destroy()
        win.protocol("WM_DELETE_WINDOW", _close)
        win.bind("<Escape>", _close)


    def _open_command_palette(self):
        """
        Global 'BIG search' that suggests which modal/tool to open based on
        English/Danish buzzwords. Opens on Ctrl+K / ⌘K and via menu.
        """
        import difflib, unicodedata, tkinter as tk
        from tkinter import ttk

        # --- helpers -----------------------------------------------------------
        def norm(s: str) -> str:
            # fold to ascii, lowercase, trim
            s = unicodedata.normalize("NFKD", s)
            s = s.encode("ascii", "ignore").decode("ascii")
            return " ".join(s.lower().split())

        def score(query: str, hay: str) -> float:
            # token hits + fuzzy fallback
            q = norm(query); h = norm(hay)
            if not q:
                return 0.0
            tokens = q.split()
            hits = sum(1 for t in tokens if t in h)
            fuzz = difflib.SequenceMatcher(None, q, h).ratio()
            return hits * 2.0 + fuzz  # tune weighting as desired

        # --- catalog of tools (label, callable, keywords, category) ------------
        # NOTE: all handlers already exist in your GUI. :contentReference[oaicite:2]{index=2}
        tools = [
            # Equilibrium / solubility
            ("Ksp & Heterogeneous Helper", self._open_ksp_tool,
             ["ksp", "solubility product", "molar solubility", "molær opløselighed",
              "molære opløselighed", "opløselighedsprodukt", "qsp", "precipitation",
              "utfældning", "saturated", "unsaturated", "umættet", "overmættet"],
             "Equilibrium"),
            ("ICE Solver (Kc/Kp, Q, shift)", self._open_equilibrium_ice_tool,
             ["ice", "equilibrium", "ligevægt", "extent", "solve K", "solve Q",
              "predict direction", "Kc", "Kp"],
             "Equilibrium"),
            ("Q & Le Châtelier Helper", self._open_q_direction_tool,
             ["Q", "direction", "le chatelier", "tryk", "pressure", "volume",
              "add product", "add reactant", "skift", "shift"],
             "Equilibrium"),
            ("Kp ↔ Kc Converter", self._open_kp_kc_converter,
             ["kp", "kc", "delta n", "Δn_gas", "gas", "tryk", "koncentration"],
             "Equilibrium"),
            ("K-expression (builder)", self._open_k_expression_builder,
             ["mass action", "equilibrium expression", "heterogeneous",
              "omit solids", "fast stof", "ren vaeske"],
             "Equilibrium"),

            # Solutions
            ("Percent & Mixing", self._open_percent_mixer_tool,
             ["percent", "pct", "mixer", "blanding", "dilution", "fortynding",
              "ppm", "w/w", "w/v"],
             "Solutions"),
            ("Henry’s Law", self._open_henry_tool,
             ["henry", "gas solubility", "opløselighed af gas", "kh", "partial pressure"],
             "Solutions"),
            ("Raoult’s Law (2 comp)", self._open_raoult_tool,
             ["raoult", "vapor pressure", "damptryk", "ideal solution"],
             "Solutions"),
            ("Ion Concentrations / van ’t Hoff", self._open_ion_vant_hoff_tool,
             ["van 't hoff", "vant hoff", "i factor", "ion factor", "osmolarity",
              "electrolyte", "dissociation"],
             "Solutions"),
            ("Colligative ΔT (Kb/Kf)", self._open_colligative_tool,
             ["boiling elevation", "freezing depression", "kogepunkt", "frysepunkt",
              "Kb", "Kf", "molality", "molalitet"],
             "Solutions"),
            ("Osmotic Pressure", self._open_osmotic_tool,
             ["osmotic", "osmotisk tryk", "pi", "Π=iMRT", "semipermeable"],
             "Solutions"),
            ("Kb / Kf (common solvents) — lookup", self._open_kb_kf_lookup,
             ["kb", "kf", "ebulioskopi", "kryoskopi", "solvent table"],
             "Lookup"),

            # Acid/base
            ("pH (strong acids/bases)", self._open_ph_modal,
             ["ph", "strong acid", "stærk syre", "stærk base", "[h+]", "[oh-]"],
             "Acid & Base"),
            ("pH (weak / buffer, Ka/Kb, Henderson–Hasselbalch)", self._open_weak_acid_base_modal,
             ["weak acid", "svag syre", "buffer", "henderson", "hasselbalch", "pka",
              "ka", "kb", "pkb", "salt of weak", "konjugeret"],
             "Acid & Base"),

            # Kinetics
            ("Rate laws (0/1/2 order)", self._open_rate_law_tool,
             ["rate law", "reaktionsorden", "half-life", "halveringstid", "integrated"],
             "Kinetics"),
            ("Arrhenius (k, Ea, A, T)", self._open_arrhenius_tool,
             ["arrhenius", "activation energy", "aktiveringsenergi", "ln k vs 1/t"],
             "Kinetics"),
            ("Initial Rates (find m, n)", self._open_initial_rates_tool,
             ["initial rates", "begin", "m", "n", "order exponents", "metode for begyndelseshastigheder"],
             "Kinetics"),

            # Conversions / thermo / misc
            ("Mass ↔ Moles (element/formula)", self._open_mol_converter,
             ["molar mass", "molarmasse", "stofmængde", "g→mol", "mol→g", "n", "M"],
             "Convert"),
            ("Pressure units", self._open_pressure_converter,
             ["pressure", "tryk", "atm", "kpa", "mmhg", "torr", "pa"],
             "Convert"),
            ("Thermo lookup (ΔH°f, ΔG°f, S°)", self._open_thermo_lookup,
             ["thermo table", "enthalpy of formation", "standard entropy", "gibbs", "termokemi"],
             "Thermo"),
            ("Reaction builder (ΔH°, ΔS°, ΔG°)", self._open_reaction_builder,
             ["reaction builder", "balance", "enthalpy change", "gibbs", "entropy"],
             "Thermo"),
            ("Photon color from E/λ", self._open_photon_color_tool,
             ["photon", "wavelength", "energy", "lambda", "nm", "color", "farve"],
             "Other"),
            ("Heating/Cooling (phase changes)", self._open_phase_change_tool,
             ["specific heat", "latent heat", "q=mcΔt", "smelte", "fordampe", "fase"],
             "Other"),
        ]

        # make a flattened searchable string for each
        indexed = []
        for label, fn, kws, cat in tools:
            hay = " ".join([label] + kws + [cat])
            indexed.append({"label": label, "open": fn, "cat": cat, "hay": hay})

        # --- UI ----------------------------------------------------------------
        win = tk.Toplevel(self)
        win.title("Search tools")
        win.transient(self); win.grab_set(); win.geometry("720x520")

        wrap = ttk.Frame(win, padding=12); wrap.pack(fill=tk.BOTH, expand=True)
        q = tk.StringVar()
        ent = ttk.Entry(wrap, textvariable=q, font=("Segoe UI", 14))
        ent.pack(fill=tk.X, pady=(0,8)); ent.focus_set()

        info = tk.StringVar(value="Type keywords (English / Dansk). ↑/↓ to pick, Enter to open.")
        ttk.Label(wrap, textvariable=info, foreground="#555").pack(anchor="w", pady=(0,6))

        lst = tk.Listbox(wrap, height=12)
        lst.pack(fill=tk.BOTH, expand=True)

        detail = tk.StringVar(value="—")
        ttk.Label(wrap, textvariable=detail, font=("Segoe UI", 10, "bold")).pack(anchor="w", pady=(8,0))

        def refresh(*_):
            lst.delete(0, tk.END)
            query = q.get()
            scored = [(score(query, item["hay"]), item) for item in indexed]
            scored.sort(key=lambda t: t[0], reverse=True)
            top = [it for s,it in scored[:16] if s > 0] or [it for s,it in scored[:16]]
            for item in top:
                lst.insert(tk.END, f"{item['label']}   —   {item['cat']}")
            if lst.size():
                lst.selection_clear(0, tk.END)
                lst.selection_set(0); lst.activate(0)
                detail.set(top[0]["label"])
            else:
                detail.set("—")

        def choose(*_):
            sel = lst.curselection()
            if not sel: 
                win.destroy(); return
            text = lst.get(sel[0]).split("   —   ", 1)[0]
            for item in indexed:
                if item["label"] == text:
                    # open the modal/tool
                    item["open"]()
                    break
            try:
                win.grab_release()
            except Exception:
                pass
            win.destroy()

        def on_move(_):
            sel = lst.curselection()
            if sel:
                detail.set(lst.get(sel[0]))

        ent.bind("<KeyRelease>", refresh)
        lst.bind("<<ListboxSelect>>", on_move)
        lst.bind("<Double-Button-1>", choose)
        lst.bind("<Return>", choose)
        win.bind("<Return>", choose)
        win.bind("<Escape>", lambda e: (win.grab_release(), win.destroy()))

        refresh()

    def _open_vant_hoff_tool(self):
        import tkinter as tk
        from tkinter import ttk, messagebox
        import math

        win = tk.Toplevel(self)
        win.title("van ’t Hoff — ΔH° (and ΔS°) from K vs T")
        win.transient(self); win.grab_set(); win.geometry("600x360")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # Inputs
        ttk.Label(frm, text="T₁ (K):").grid(row=0, column=0, sticky="w")
        t1 = ttk.Entry(frm, width=10); t1.grid(row=0, column=1, sticky="w", padx=(6,12)); t1.insert(0,"298")
        ttk.Label(frm, text="K₁ (dimensionless):").grid(row=0, column=2, sticky="w")
        k1 = ttk.Entry(frm, width=12); k1.grid(row=0, column=3, sticky="w", padx=(6,0)); k1.insert(0,"3.10")

        ttk.Label(frm, text="T₂ (K):").grid(row=1, column=0, sticky="w")
        t2 = ttk.Entry(frm, width=10); t2.grid(row=1, column=1, sticky="w", padx=(6,12)); t2.insert(0,"323")
        ttk.Label(frm, text="K₂ (dimensionless):").grid(row=1, column=2, sticky="w")
        k2 = ttk.Entry(frm, width=12); k2.grid(row=1, column=3, sticky="w", padx=(6,0)); k2.insert(0,"0.550")

        # Optional prediction at T3
        ttk.Label(frm, text="Predict K at T₃ (K) [optional]:").grid(row=2, column=0, columnspan=2, sticky="w", pady=(8,0))
        t3 = ttk.Entry(frm, width=10); t3.grid(row=2, column=2, sticky="w", padx=(6,0)); t3.insert(0,"")

        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI",10,"bold"), justify="left").grid(
            row=3, column=0, columnspan=4, sticky="w", pady=(12,4)
        )

        expl = tk.Text(frm, height=7, wrap="word"); expl.grid(row=4, column=0, columnspan=4, sticky="nsew")
        frm.grid_rowconfigure(4, weight=1)
        sv = ttk.Scrollbar(frm, orient="vertical", command=expl.yview)
        expl.configure(yscrollcommand=sv.set); sv.grid(row=4, column=4, sticky="ns")
        expl.configure(state="disabled")

        def write(lines):
            expl.configure(state="normal"); expl.delete("1.0","end")
            expl.insert("1.0","\n".join(lines)); expl.configure(state="disabled")

        def compute():
            try:
                T1 = float(t1.get().replace(",",".")); T2 = float(t2.get().replace(",",".")); 
                K1 = float(k1.get().replace(",",".")); K2v = float(k2.get().replace(",","."))
                if T1<=0 or T2<=0 or K1<=0 or K2v<=0: raise ValueError
            except ValueError:
                messagebox.showerror("Bad input","T must be > 0 K and K must be > 0.", parent=win); return

            R = 8.314462618  # J mol^-1 K^-1
            dH = -R * math.log(K2v/K1) / (1.0/T2 - 1.0/T1)     # J/mol
            # Solve ΔS° from either point: ΔS° = R ln K + ΔH°/T
            dS = R*math.log(K1) + dH/T1                        # J/(mol·K)
            heat = "exothermic (ΔH° < 0)" if dH < 0 else "endothermic (ΔH° > 0)" if dH > 0 else "thermally neutral"

            # Optional prediction at T3
            pred = ""
            t3txt = (t3.get() or "").strip()
            if t3txt:
                try:
                    T3 = float(t3txt.replace(",",".")); 
                    if T3 <= 0: raise ValueError
                    lnK3 = (-dH/(R*T3)) + (dS/R)
                    K3 = math.exp(lnK3)
                    pred = f"\nPredicted K at T₃={T3:g} K → K ≈ {K3:.4g}"
                except ValueError:
                    messagebox.showerror("Bad T₃","T₃ must be > 0 K.", parent=win); return

            out.set(f"ΔH° = {dH/1000.0:.3g} kJ/mol   ({heat})\nΔS° = {dS:.3g} J/(mol·K){pred}")

            lines = [
                "van ’t Hoff (two-point) with ΔH°, ΔS° approx. constant:",
                "  ln(K₂/K₁) = −(ΔH°/R)·(1/T₂ − 1/T₁)",
                f"  T₁={T1:g} K,  K₁={K1:g};  T₂={T2:g} K,  K₂={K2v:g}",
                f"  ΔH° = −R·ln(K₂/K₁) / (1/T₂ − 1/T₁) = {dH/1000.0:.6g} kJ/mol",
                "  From  ln K = −ΔH°/(R T) + ΔS°/R:",
                f"  ΔS° = R ln K₁ + ΔH°/T₁ = {dS:.6g} J/(mol·K)",
                "Notes: assumes ΔH°, ΔS° are T-independent over this range; K is dimensionless.",
            ]
            write(lines)
            return dH, dS

        def add_known():
            res = compute()
            if not res: return
            dH, dS = res  # J/mol and J/(mol·K)
            self.known_base["ΔH°rxn"] = float(dH)
            self.known_ui["ΔH°rxn"] = {"unit":"kJ/mol","display_value": dH/1000.0}
            self.known_base["ΔS°rxn"] = float(dS)
            self.known_ui["ΔS°rxn"] = {"unit":"J/(mol·K)","display_value": dS}
            self._refresh_known_table(); self._refresh_equation_list()
            messagebox.showinfo("Added","ΔH°, ΔS° added to Known Variables.", parent=win)

        btns = ttk.Frame(frm); btns.grid(row=5, column=0, columnspan=4, sticky="w", pady=(8,0))
        ttk.Button(btns, text="Compute", command=compute).pack(side=tk.LEFT, padx=(0,8))
        ttk.Button(btns, text="Add to Known Variables", command=add_known).pack(side=tk.LEFT)


    def _open_oxidation_number_tool(self):
        import re, math
        import tkinter as tk
        from tkinter import ttk, messagebox

        win = tk.Toplevel(self)
        win.title("Oxidation number — single species")
        win.transient(self); win.grab_set(); win.geometry("680x460")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # ---------- Inputs ----------
        ttk.Label(frm, text="Species (e.g., S^2-, SOCl2, HSO4^-, Fe(CN)6^3-):").grid(row=0, column=0, sticky="w")
        ent_spec = ttk.Entry(frm, width=36); ent_spec.grid(row=0, column=1, sticky="w", padx=(6,8))
        ent_spec.insert(0, "HSO4^-")

        ttk.Label(frm, text="Element to solve (symbol):").grid(row=1, column=0, sticky="w", pady=(8,0))
        ent_elem = ttk.Entry(frm, width=10); ent_elem.grid(row=1, column=1, sticky="w", padx=(6,0), pady=(8,0))
        ent_elem.insert(0, "S")

        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 10, "bold"), justify="left").grid(
            row=2, column=0, columnspan=2, sticky="w", pady=(12,4)
        )

        expl = tk.Text(frm, height=12, wrap="word")
        expl.grid(row=3, column=0, columnspan=2, sticky="nsew")
        frm.grid_rowconfigure(3, weight=1)
        sv = ttk.Scrollbar(frm, orient="vertical", command=expl.yview)
        expl.configure(yscrollcommand=sv.set); sv.grid(row=3, column=2, sticky="ns")
        expl.configure(state="disabled")

        # Quick examples
        ex = ttk.Frame(frm); ex.grid(row=4, column=0, columnspan=2, sticky="w", pady=(8,0))
        def _fill(x, e):
            ent_spec.delete(0, tk.END); ent_spec.insert(0, x)
            ent_elem.delete(0, tk.END); ent_elem.insert(0, e)
        ttk.Label(ex, text="Examples: ").pack(side=tk.LEFT)
        ttk.Button(ex, text="S^2−",        command=lambda:_fill("S^2-","S")).pack(side=tk.LEFT, padx=4)
        ttk.Button(ex, text="SOCl2",      command=lambda:_fill("SOCl2","S")).pack(side=tk.LEFT, padx=4)
        ttk.Button(ex, text="HSO4^−",      command=lambda:_fill("HSO4^-","S")).pack(side=tk.LEFT, padx=4)
        ttk.Button(ex, text="Fe(CN)6^3−", command=lambda:_fill("Fe(CN)6^3-","Fe")).pack(side=tk.LEFT, padx=4)

        # ---------- Chemistry helpers ----------
        GROUP1 = {"Li","Na","K","Rb","Cs","Fr"}
        GROUP2 = {"Be","Mg","Ca","Sr","Ba","Ra"}
        METALS  = GROUP1 | GROUP2 | {"Al","Zn","Fe","Cu","Ag","Au","Co","Ni","Sn","Pb","Hg","Cr","Mn","Ti","V","Cd","Pt","Pd"}
        HALOGENS = {"F","Cl","Br","I"}

        def _strip_phase_and_ws(s):
            s = s.strip()
            s = re.sub(r"\([a-z]{1,3}\)\s*", "", s, flags=re.I)  # remove (s),(l),(g),(aq)
            return s.replace("·","").replace(" ", "")

        def _read_charge(s):
            # Accept trailing forms: ^3-, ^2+, 3-, -, +, 2+, etc.
            # Return (clean_formula, charge_int)
            m = re.search(r"(?:\^?\s*([+-]?\d+)?\s*([+-]))\s*$", s)
            if m:
                num = m.group(1)
                sign = m.group(2)
                mag = int(num) if num not in (None,"","+","-") else 1
                chg = mag if sign=="+" else -mag
                core = s[:m.start()].strip()
                return core, chg
            # also allow something like "S2-" (where the trailing '2-' is charge and not stoich)
            m2 = re.search(r"(\d+)([+-])\s*$", s)
            if m2 and re.search(r"[A-Za-z)]$", s[:m2.start()] or ""):
                mag = int(m2.group(1)); chg = mag if m2.group(2)=="+" else -mag
                core = s[:m2.start()].strip()
                return core, chg
            # If no explicit charge, neutral
            return s, 0

        def _tokenize(formula):
            # returns list of (token, count) without charge; supports nested ()n
            # e.g. Fe(CN)6 -> {"Fe":1,"C":6,"N":6}
            i, n = 0, len(formula)
            stack = [dict()]
            def push_count(elem, cnt):
                if not elem: return
                d = stack[-1]
                d[elem] = d.get(elem, 0) + cnt

            while i < n:
                if formula[i] == '(':
                    stack.append(dict()); i += 1
                elif formula[i] == ')':
                    i += 1
                    # read multiplier
                    j = i
                    while j < n and formula[j].isdigit():
                        j += 1
                    mult = int(formula[i:j] or "1")
                    grp = stack.pop()
                    for k,v in grp.items():
                        push_count(k, v*mult)
                    i = j
                else:
                    # element symbol
                    if not formula[i].isalpha():
                        raise ValueError(f"Bad token at {formula[i:]}")
                    sym = formula[i]
                    i += 1
                    if i<n and formula[i].islower():
                        sym += formula[i]; i += 1
                    # optional count
                    j = i
                    while j < n and formula[j].isdigit():
                        j += 1
                    cnt = int(formula[i:j] or "1")
                    push_count(sym, cnt)
                    i = j
            if len(stack) != 1:
                raise ValueError("Unmatched parentheses")
            return stack[0]

        def parse_species(spec):
            s = _strip_phase_and_ws(spec)
            core, charge = _read_charge(s)
            atoms = _tokenize(core)
            return atoms, charge

        def oxidation_number(species, target=None):
            """
            Return:
            ox_map: dict element -> oxidation number (float) for all elements we can determine
            target_sym: the symbol we solved for (if any)
            lines: list of explanation lines including the equation solved
            If there is exactly one unknown element after applying rules, it will be solved as x.
            """
            atoms, charge = parse_species(species)
            if not atoms:
                raise ValueError("Could not parse any atoms from the species.")

            # normalize target symbol if provided
            target_sym = None
            if target:
                t = target.strip().capitalize()
                if t not in atoms:
                    raise ValueError(f"Element {t} not found in {species}")
                target_sym = t

            # --- Rules (same logic as before) ---
            GROUP1 = {"Li","Na","K","Rb","Cs","Fr"}
            GROUP2 = {"Be","Mg","Ca","Sr","Ba","Ra"}
            METALS  = GROUP1 | GROUP2 | {"Al","Zn","Fe","Cu","Ag","Au","Co","Ni","Sn","Pb","Hg","Cr","Mn","Ti","V","Cd","Pt","Pd"}
            HALOGENS = {"F","Cl","Br","I"}

            notes = []
            known = {}

            # F always −1
            if "F" in atoms:
                known["F"] = -1; notes.append("Rule: F is −1 in compounds.")
            # Group 1/2, Al
            for e in atoms.keys() & GROUP1:
                known[e] = +1; notes.append(f"Rule: {e} (Group 1) is +1.")
            for e in atoms.keys() & GROUP2:
                known[e] = +2; notes.append(f"Rule: {e} (Group 2) is +2.")
            if "Al" in atoms:
                known["Al"] = +3; notes.append("Rule: Al is +3 in compounds.")

            # Oxygen (common cases + quick exceptions)
            if "O" in atoms and "O" not in known:
                formula = _strip_phase_and_ws(species)
                # OF2 (O+2) and O2F2 (O+1)
                if set(atoms.keys()) == {"O","F"}:
                    if atoms.get("O",0) == 1 and atoms.get("F",0) == 2:
                        known["O"] = +2; notes.append("Exception: O is +2 in OF₂.")
                    elif atoms.get("O",0) == 2 and atoms.get("F",0) == 2:
                        known["O"] = +1; notes.append("Exception: O is +1 in O₂F₂.")
                    else:
                        known["O"] = -2; notes.append("Rule: O is usually −2.")
                elif ("O2^2-" in formula) or re.search(r"\(O2\)\s*2-\s*$", formula):
                    known["O"] = -1; notes.append("Exception: peroxide detected → O is −1.")
                elif ("O2^-" in formula) or re.search(r"\(O2\)\s*-\s*$", formula):
                    known["O"] = -0.5; notes.append("Exception: superoxide detected → O is −½.")
                else:
                    known["O"] = -2; notes.append("Rule: O is usually −2.")

            # Hydrogen (+1 with nonmetals, −1 with metals)
            if "H" in atoms and "H" not in known:
                if (atoms.keys() & METALS):
                    known["H"] = -1; notes.append("Rule: H is −1 in metal hydrides.")
                else:
                    known["H"] = +1; notes.append("Rule: H is +1 with nonmetals.")

            # Halogens (Cl/Br/I usually −1 unless with O or F)
            for X in (atoms.keys() & (HALOGENS - {"F"})):
                if "O" not in atoms and "F" not in atoms and X not in known:
                    known[X] = -1; notes.append(f"Rule: {X} is −1 (no O/F present).")

            # If user didn't specify a target, and exactly one element remains unknown, use that
            unknown = [e for e in atoms if e not in known]
            if target_sym is None and len(unknown) == 1:
                target_sym = unknown[0]

            # If there is exactly one unknown element, solve it
            ox_map = dict(known)  # seed with known values
            lines = []
            fmt = lambda v: f"{v:+g}"

            if target_sym and target_sym not in ox_map:
                n_target = atoms[target_sym]
                sum_known = sum(atoms[e]*ox_map[e] for e in ox_map)
                # Equation: n_target·x + sum_known = charge
                # Build pretty equation string like: +1 + x + 4(−2) = −1
                term_strs = []
                # known terms
                for e, n in atoms.items():
                    if e == target_sym: 
                        continue
                    if e in ox_map:
                        v = ox_map[e]
                        if n == 1:
                            term_strs.append(fmt(v))
                        else:
                            term_strs.append(f"{n}({fmt(v)})")
                # unknown term
                if n_target == 1:
                    term_strs.insert(0, "x")
                else:
                    term_strs.insert(0, f"{n_target}(x)")

                eqn_pretty = " + ".join(term_strs).replace("+ -", " - ")
                # Solve x
                x = (charge - sum_known) / n_target
                ox_map[target_sym] = x

                lines.append(f"Let x = oxidation number of {target_sym}.")
                lines.append(f"Sum of oxidation numbers = overall charge {fmt(charge)}.")
                lines.append(f"Equation: {eqn_pretty} = {fmt(charge)}")
                lines.append(f"⇒ x = (charge − Σknown)/n({target_sym}) = ({fmt(charge)} − {fmt(sum_known)})/{n_target:g} = {fmt(x)}")
            else:
                # nothing to solve (everything assigned by rules, or multiple unknowns)
                pass

            # Per-element summary (like labels above each element)
            lines.append("\nPer-element oxidation numbers:")
            for e in sorted(atoms.keys()):
                v = ox_map.get(e, None)
                lines.append(f"  {e}: {fmt(v) if v is not None else '(undetermined)'} (count n={atoms[e]})")

            if notes:
                lines.append("\nRules/exceptions used:")
                lines += [f"  • {t}" for t in notes]

            return ox_map, target_sym, lines
        
        def write(lines):
            expl.configure(state="normal"); expl.delete("1.0","end")
            expl.insert("1.0","\n".join(lines)); expl.configure(state="disabled")


        def compute():
            spec = ent_spec.get().strip()
            elem = ent_elem.get().strip()
            # Allow blank target; we’ll pick the sole unknown automatically
            target = elem if elem else None
            try:
                ox_map, tgt, steps = oxidation_number(spec, target)
            except Exception as e:
                messagebox.showerror("Error", str(e), parent=win); return

            # Headline + compact list (H:+1, S:+6, O:−2)
            def fmt(v): 
                return f"{v:+g}"
            listing = ", ".join(f"{k}:{fmt(ox_map[k])}" for k in sorted(ox_map))
            title = f"Oxidation numbers in {spec}:  {listing}"
            if tgt:
                title += f"   (solved: {tgt})"

            out.set(title)
            write(steps)


        # Buttons
        btns = ttk.Frame(frm); btns.grid(row=5, column=0, columnspan=2, sticky="w", pady=(10,0))
        ttk.Button(btns, text="Compute", command=compute).pack(side=tk.LEFT, padx=(0,8))

    
    def _open_redox_analysis_tool(self):
        import re, math
        import tkinter as tk
        from tkinter import ttk, messagebox

        def _split_species_list(s: str):
            """Split a side like 'Ag + Pt^2+ + NO3-' into tokens without breaking charges."""
            import re
            # only split on plus signs that are surrounded by whitespace (the separators)
            return [t.strip() for t in re.split(r"\s+\+\s+", s.strip()) if t.strip()]


        # Reuse the species parser from the single-species modal via inner copy
        # --- parsing helpers (inside _open_redox_analysis_tool) ---
        def _strip_phase_and_ws(s):
            import re
            s = s.strip()
            s = re.sub(r"\([a-z]{1,3}\)\s*", "", s, flags=re.I)  # remove (s),(l),(g),(aq)
            return s.replace("·", "").replace(" ", "")

        def _read_charge(s):
            # accept: ^3-, ^2+, 3-, -, +, etc.
            import re
            m = re.search(r"(?:\^?\s*([+-]?\d+)?\s*([+-]))\s*$", s)
            if m:
                num = m.group(1); sign = m.group(2)
                mag = int(num) if num not in (None,"","+","-") else 1
                chg = mag if sign == "+" else -mag
                core = s[:m.start()].strip()
                return core, chg
            return s, 0

        def _tokenize(formula):
            i, n = 0, len(formula)
            stack = [dict()]
            def push(E, c):
                d = stack[-1]; d[E] = d.get(E, 0) + c
            while i < n:
                ch = formula[i]
                if ch == '(':
                    stack.append({}); i += 1
                elif ch == ')':
                    i += 1
                    j = i
                    while j < n and formula[j].isdigit(): j += 1
                    mult = int(formula[i:j] or "1")
                    grp = stack.pop()
                    for k, v in grp.items(): push(k, v*mult)
                    i = j
                else:
                    if not ch.isalpha():
                        raise ValueError(f"Unexpected token near '{formula[i:]}'")
                    sym = ch; i += 1
                    if i < n and formula[i].islower():
                        sym += formula[i]; i += 1
                    j = i
                    while j < n and formula[j].isdigit(): j += 1
                    cnt = int(formula[i:j] or "1"); push(sym, cnt); i = j
            if len(stack) != 1:
                raise ValueError("Unmatched parentheses")
            return stack[0]

        def parse_species(spec):
            import re
            s = _strip_phase_and_ws(spec)
            # leading stoichiometric coefficient (ALWAYS accept digits if present)
            m = re.match(r"^(\d+)(.*)$", s)
            coeff = 1
            if m:
                coeff = int(m.group(1))
                s = m.group(2)
            # trailing charge, then tokenize the core
            core, chg = _read_charge(s)
            atoms = _tokenize(core)
            return coeff, atoms, chg


        # very light rules to get oxidation numbers (same as in single-species)
        GROUP1={"Li","Na","K","Rb","Cs","Fr"}; GROUP2={"Be","Mg","Ca","Sr","Ba","Ra"}
        METALS  = GROUP1|GROUP2|{"Al","Zn","Fe","Cu","Ag","Au","Co","Ni","Sn","Pb","Hg","Cr","Mn","Ti","V","Cd","Pt","Pd"}
        HALO={"F","Cl","Br","I"}

        def guess_ox(atoms, charge):
            # 0) Elemental species (O2, H2, N2, Cl2, P4, S8, Na, etc.)
            if len(atoms) == 1 and charge == 0:
                el = next(iter(atoms))
                return {el: 0.0}

            # 1) Rules
            known = {}
            if "F" in atoms: known["F"] = -1
            for g1 in atoms.keys() & GROUP1: known[g1] = +1
            for g2 in atoms.keys() & GROUP2: known[g2] = +2
            if "Al" in atoms: known["Al"] = +3

            if "O" in atoms and "O" not in known:
                known["O"] = -2  # default (handled elemental O2 above)

            if "H" in atoms and "H" not in known:
                if (atoms.keys() & METALS):
                    known["H"] = -1   # metal hydrides
                else:
                    known["H"] = +1

            for X in (atoms.keys() & (HALO - {"F"})):
                if "O" not in atoms and "F" not in atoms and X not in known:
                    known[X] = -1

            # 2) If exactly one element still unknown, solve from Σ(n·ox#)=charge
            unknown = [e for e in atoms if e not in known]
            if len(unknown) == 1:
                e = unknown[0]
                s = sum(atoms[k] * known[k] for k in known)
                known[e] = (charge - s) / atoms[e]
            return known


        def analyze(side_str):
            # split by + and collect average ox# per element weighted by stoich
            parts = _split_species_list(side_str)
            el2avg={}
            for p in parts:
                coeff, atoms, chg = parse_species(p)
                ox = guess_ox(atoms, chg)
                for el,n in atoms.items():
                    if el in ox:
                        el2avg.setdefault(el,0.0); el2avg[el]+= coeff*n*ox[el]
            # normalize by total atoms of each element
            counts={}
            for p in parts:
                coeff, atoms, _ = parse_species(p)
                for el,n in atoms.items():
                    counts[el]=counts.get(el,0)+coeff*n
            for el in list(el2avg):
                el2avg[el]/=counts[el]
            return el2avg

        win = tk.Toplevel(self)
        win.title("Redox analysis — reaction")
        win.transient(self); win.grab_set(); win.geometry("720x420")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        ttk.Label(frm, text="Reactants (e.g., 2Na + 2H2O):").grid(row=0, column=0, sticky="w")
        eL = ttk.Entry(frm, width=50); eL.grid(row=0, column=1, sticky="w", padx=(6,8)); eL.insert(0, "2Na + 2H2O")
        ttk.Label(frm, text="Products (e.g., 2Na+ + 2OH- + H2):").grid(row=1, column=0, sticky="w")
        eR = ttk.Entry(frm, width=50); eR.grid(row=1, column=1, sticky="w", padx=(6,8)); eR.insert(0, "2Na+ + 2OH- + H2")

        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 10, "bold")).grid(row=2, column=0, columnspan=2, sticky="w", pady=(10,4))

        box = tk.Text(frm, height=12, wrap="word"); box.grid(row=3, column=0, columnspan=2, sticky="nsew")
        frm.grid_rowconfigure(3, weight=1)
        sv = ttk.Scrollbar(frm, orient="vertical", command=box.yview)
        box.configure(yscrollcommand=sv.set); sv.grid(row=3, column=2, sticky="ns")
        box.configure(state="disabled")

        def write(lines):
            box.configure(state="normal"); box.delete("1.0","end"); box.insert("1.0","\n".join(lines)); box.configure(state="disabled")

        def compute():
            try:
                L = analyze(eL.get()); R = analyze(eR.get())
            except Exception as err:
                messagebox.showerror("Parse error", str(err), parent=win); return
            elems = sorted(set(L)|set(R))
            lines = ["Average oxidation numbers per element:"]
            for el in elems:
                l = L.get(el, None); r = R.get(el, None)
                if l is None or r is None: 
                    lines.append(f"  {el}: (missing on one side)"); continue
                change = r - l
                tag = "oxidized (↑)" if change>0 else ("reduced (↓)" if change<0 else "no change")
                lines.append(f"  {el}: left {l:.3g} → right {r:.3g}   Δ = {change:+.3g}  → {tag}")
            out.set("Oxidation change summary:")
            def per_species_breakdown(side_str, tag):
                parts = _split_species_list(side_str)
                lines = [f"{tag}:"]
                for p in parts:
                    coeff, atoms, chg = parse_species(p)
                    ox = guess_ox(atoms, chg)
                    listing = ", ".join(f"{el}:{ox[el]:+g}" for el in sorted(atoms) if el in ox)
                    lines.append(f"  {p}: {listing}")
                return lines

            # after you compute L and R (your current averages), append:
            lines += per_species_breakdown(eL.get(), "Per-species (reactants)")
            lines += per_species_breakdown(eR.get(), "Per-species (products)")
            write(lines)

        ttk.Button(frm, text="Analyze", command=compute).grid(row=4, column=0, sticky="w", pady=(8,0))


    def _open_net_ionic_tool(self):
        import re, math
        from fractions import Fraction
        import tkinter as tk
        from tkinter import ttk, messagebox

        win = tk.Toplevel(self)
        win.title("Net-ionic equation")
        win.transient(self); win.grab_set(); win.geometry("840x560")
        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # ---------- Inputs ----------
        ttk.Label(frm, text="Molecular Reactants (e.g., 2Ag(s) + Pt(NO3)2(aq))").grid(row=0, column=0, sticky="w")
        eL = ttk.Entry(frm, width=80); eL.grid(row=1, column=0, sticky="we", pady=(2,8))
        ttk.Label(frm, text="Molecular Products (e.g., 2AgNO3(aq) + Pt(s))").grid(row=2, column=0, sticky="w")
        eR = ttk.Entry(frm, width=80); eR.grid(row=3, column=0, sticky="we", pady=(2,6))
        frm.grid_columnconfigure(0, weight=1)

        do_balance = tk.BooleanVar(value=True)
        ttk.Checkbutton(frm, text="Attempt to balance the molecular equation first",
                        variable=do_balance).grid(row=4, column=0, sticky="w")

        already_ionic = tk.BooleanVar(value=False)
        ttk.Checkbutton(frm, text="Inputs are already ionic (skip dissociation)",
                        variable=already_ionic).grid(row=5, column=0, sticky="w", pady=(2,0))

        # Output widgets (pushed down one row)
        out = tk.Text(frm, height=20, wrap="word"); out.grid(row=6, column=0, sticky="nsew", pady=(10,0))
        frm.grid_rowconfigure(6, weight=1)
        sv = ttk.Scrollbar(frm, orient="vertical", command=out.yview)
        sv.grid(row=6, column=1, sticky="ns")

        def write(lines):
            out.configure(state="normal"); out.delete("1.0","end")
            out.insert("1.0","\n".join(lines)); out.configure(state="disabled")



        # ---------- Helpers: parsing ----------
        def _strip_phase_and_ws(s):
            s = s.strip()
            s = re.sub(r"\([a-z]{1,3}\)\s*", "", s, flags=re.I)   # remove (s),(l),(g),(aq) for counting
            return s.replace("·","").replace(" ", "")

        def _read_phase(s):
            m = re.search(r"\((s|l|g|aq)\)\s*$", s, flags=re.I)
            return (m.group(1).lower() if m else None)

        def _read_charge(s):
            m = re.search(r"(?:\^?\s*([+-]?\d+)?\s*([+-]))\s*$", s)
            if m:
                num = m.group(1); sign = m.group(2)
                mag = int(num) if num not in (None,"","+","-") else 1
                chg = mag if sign=="+" else -mag
                core = s[:m.start()].strip()
                return core, chg
            return s, 0

        def _tokenize(formula):
            i, n = 0, len(formula)
            stack = [dict()]
            def push(E,c):
                d = stack[-1]; d[E] = d.get(E,0) + c
            while i < n:
                ch = formula[i]
                if ch == '(':
                    stack.append({}); i += 1
                elif ch == ')':
                    i += 1
                    j = i
                    while j < n and formula[j].isdigit(): j += 1
                    mult = int(formula[i:j] or "1")
                    grp = stack.pop()
                    for k,v in grp.items(): push(k, v*mult)
                    i = j
                else:
                    if not ch.isalpha(): raise ValueError(f"Unexpected token near '{formula[i:]}'")
                    sym = ch; i += 1
                    if i<n and formula[i].islower(): sym += formula[i]; i += 1
                    j = i
                    while j<n and formula[j].isdigit(): j += 1
                    cnt = int(formula[i:j] or "1"); push(sym, cnt); i = j
            if len(stack)!=1: raise ValueError("Unmatched parentheses")
            return stack[0]

        def _split_species_list(s: str):
            # split only on + that act as separators (have spaces around)
            return [t.strip() for t in re.split(r"\s+\+\s+", s.strip()) if t.strip()]

        def parse_species_token(tok):
            s = tok.strip()
            m = re.match(r"^(\d+)\s*(.*)$", s)
            coeff = 1
            if m:
                coeff = int(m.group(1)); s = m.group(2).strip()

            phase = None
            m = re.search(r"\((s|l|g|aq)\)\s*$", s, flags=re.I)
            if m:
                phase = m.group(1).lower()
                s = s[:m.start()].strip()  # remove phase FIRST

            core, charge = _read_charge(s)   # now sees ...^2+ at end
            atoms = _tokenize(_strip_phase_and_ws(core))
            return {"coeff": coeff, "core": core, "phase": phase, "charge": charge, "atoms": atoms}



        # ---------- Dissociation rules ----------
        POLY = {
            "NO3": -1, "NO2": -1, "ClO4": -1, "ClO3": -1, "ClO": -1, "CN": -1,
            "OH": -1, "SCN": -1, "C2H3O2": -1, "CH3COO": -1,
            "SO4": -2, "SO3": -2, "CO3": -2, "S": -2, "Cr2O7": -2, "HCO3": -1, "HSO4": -1,
            "PO4": -3, "HPO4": -2, "H2PO4": -1, "MnO4": -1
        }
        # strong acids (fully dissociate)
        STRONG_ACIDS = {
            "HCl": ("H^+", "Cl^-"),
            "HBr": ("H^+", "Br^-"),
            "HI":  ("H^+", "I^-"),
            "HNO3":("H^+", "NO3^-"),
            "HClO3":("H^+", "ClO3^-"),
            "HClO4":("H^+", "ClO4^-"),
            "H2SO4":("H^+", "H^+", "SO4^2-"),   # treat as 2 H+ + SO4^2-
        }
        STRONG_BASES = {
            "NaOH":("Na^+", "OH^-"), "KOH":("K^+","OH^-"),
            "LiOH":("Li^+","OH^-"), "RbOH":("Rb^+","OH^-"), "CsOH":("Cs^+","OH^-"),
            "Ba(OH)2":("Ba^2+","OH^-","OH^-"), "Sr(OH)2":("Sr^2+","OH^-","OH^-"), "Ca(OH)2":("Ca^2+","OH^-","OH^-")
        }
        GROUP1 = {"Li","Na","K","Rb","Cs","Fr"}
        GROUP2 = {"Be","Mg","Ca","Sr","Ba","Ra"}
        NH4 = "NH4"

        def as_ion_str(elem, z):
            sign = "+" if z>0 else "-"
            mag = abs(int(z))
            return f"{elem}{'^'+str(mag) if mag>1 else ''}{sign}"

        def try_dissociate(spec):
            """Return list of (ion_str, count) or None if left as molecule."""
            # keep non-aqueous intact
            if spec["phase"] != "aq":
                return None

            core = spec["core"]
            bare = _strip_phase_and_ws(core)

            # Already an ion like "Pt^2+"? Keep as is (stoich applied later).
            if spec["charge"] != 0:
                ion = as_ion_str(bare, spec["charge"])
                return [(ion, 1)]

            # strong acids / bases (full dissociation)
            if bare in STRONG_ACIDS:
                ions = STRONG_ACIDS[bare]
                return [(i, ions.count(i)) for i in sorted(set(ions), key=ions.index)]
            if bare in STRONG_BASES:
                ions = STRONG_BASES[bare]
                return [(i, ions.count(i)) for i in sorted(set(ions), key=ions.index)]

            # Generic ionic salts that end with a known polyatomic anion, e.g. "Pt(NO3)2", "AgNO3"
            for an, z_an in POLY.items():
                # does the formula end with (an)n or ann ?
                if re.search(rf"(?:\({an}\)|{an})\d*$", bare):
                    # how many anions?
                    m = re.search(rf"(?:\({an}\)|{an})(\d*)$", bare)
                    count_an = int(m.group(1) or "1")

                    # remove the trailing anion block; then strip any dangling parentheses
                    cat_block = re.sub(rf"(?:\({an}\)|{an})\d*$", "", bare)
                    cat_block = re.sub(r"[()]+$", "", cat_block)

                    if not cat_block:
                        continue

                    # Identify cation unit and charge
                    if cat_block == NH4:
                        cat_unit, z_cat_unit = NH4, +1
                    else:
                        toks = list(_tokenize(cat_block).keys())
                        if len(toks) != 1:
                            # complex/unknown cation → give up on splitting
                            continue
                        cat_unit = toks[0]
                        if cat_unit in GROUP1:
                            z_cat_unit = +1
                        elif cat_unit in GROUP2:
                            z_cat_unit = +2
                        else:
                            z_cat_unit = None  # will deduce from neutrality below

                    total_an = count_an * z_an
                    if z_cat_unit is None:
                        # assume one cation unit per formula unit, deduce charge by neutrality
                        z_cat_unit = -total_an
                        if z_cat_unit == 0:
                            continue

                    cat_ion = as_ion_str(cat_unit, z_cat_unit)
                    an_ion  = as_ion_str(an, z_an)

                    # how many cations per formula unit?  n_c*z_cat + count_an*z_an = 0
                    n_c = abs(total_an) // abs(z_cat_unit)
                    if n_c == 0:
                        n_c = 1

                    return [(cat_ion, n_c), (an_ion, count_an)]

            # Simple binary chlorides like NaCl, CaCl2, etc. (quick convenience rule)
            m = re.match(r"^([A-Z][a-z]?)(\d*)Cl(\d*)$", bare)
            if m:
                el = m.group(1)
                n1 = int(m.group(2) or "1")
                n2 = int(m.group(3) or "1")
                z_an = -1 * n2
                if el in GROUP1:
                    z_cat = +1
                elif el in GROUP2:
                    z_cat = +2
                else:
                    z_cat = (n2 // max(n1,1)) or 1
                return [(as_ion_str(el, z_cat), n1), ("Cl^-", n2)]

            # otherwise keep as molecule (weak electrolyte or unknown)
            return None


        # ---------- Equation handling ----------
        def format_side(specs):
            def show(sp):
                # rebuild label including charge if present
                if sp["charge"]:
                    # e.g. "Pt^2+" or "SO4^2-"
                    mag = abs(int(sp["charge"]))
                    sign = "+" if sp["charge"] > 0 else "-"
                    core = f"{sp['core']}{'^'+str(mag) if mag>1 else ''}{sign}"
                else:
                    core = sp["core"]
                if sp["phase"]:
                    core += f"({sp['phase']})"
                return (f"{sp['coeff']} " if sp['coeff'] != 1 else "") + core

            return " + ".join(show(s) for s in specs) or "—"


        def parse_side(s):
            return [parse_species_token(tok) for tok in _split_species_list(s)]

        def gather_elements(specs):
            els = set()
            for sp in specs:
                els |= set(sp["atoms"].keys())
            return sorted(els)

        def balance_molecular(left, right):
            # Build stoichiometric matrix A * coeffs = 0 over elements
            # unknowns: all species, coefficients positive; we compute nullspace rational vector
            specs = left + right
            els = sorted(set().union(*[s["atoms"].keys() for s in specs]))
            m = []
            for el in els:
                row = []
                for s in left:
                    row.append(s["atoms"].get(el,0))
                for s in right:
                    row.append(-s["atoms"].get(el,0))
                m.append(row)

            # Gauss-Jordan rational nullspace (small sizes typical)
            R = len(m); C = len(m[0]) if m else 0
            A = [[Fraction(x) for x in row] for row in m]
            # augment with zeros
            piv_col = [-1]*R
            r = c = 0
            while r<R and c<C:
                # find non-zero pivot
                piv = None
                for i in range(r, R):
                    if A[i][c] != 0:
                        piv = i; break
                if piv is None:
                    c += 1; continue
                A[r], A[piv] = A[piv], A[r]
                piv_col[r] = c
                # normalize
                pivval = A[r][c]
                A[r] = [v/pivval for v in A[r]]
                # eliminate others
                for i in range(R):
                    if i==r: continue
                    factor = A[i][c]
                    if factor != 0:
                        A[i] = [A[i][j] - factor*A[r][j] for j in range(C)]
                r += 1; c += 1

            # free variables → set 1; back-substitute
            sol = [Fraction(0) for _ in range(C)]
            free_cols = [j for j in range(C) if j not in piv_col]
            if not free_cols:
                free_cols = [C-1]
            for j in free_cols:
                sol[j] = Fraction(1)
            for i in reversed(range(R)):
                pc = piv_col[i]
                if pc == -1: continue
                s = -sum(A[i][j]*sol[j] for j in range(pc+1, C))
                sol[pc] = s

            # convert to smallest integers
            den = 1
            for v in sol: den = math.lcm(den, v.denominator)
            ints = [int(v*den) for v in sol]
            # make them positive
            sign = -1 if any(x<0 for x in ints) else 1
            ints = [sign*abs(x) for x in ints]
            # apply to copies
            L = [dict(sp) for sp in left]
            Rg = [dict(sp) for sp in right]
            for i,coef in enumerate(ints[:len(L)]):
                L[i]["coeff"] = coef
            for j,coef in enumerate(ints[len(L):]):
                Rg[j]["coeff"] = coef
            return L, Rg

        def multiply_counts(dic, k):
            return {k2: v2*k for k2,v2 in dic.items()}

        def count_ionic_side(side):
            """"side" is list of (ion_str,count) pairs already multiplied by stoich; combine like terms."""
            count = {}
            for ion, n in side:
                count[ion] = count.get(ion, 0) + n
            return count

        def dissociate_side(specs):
            ionic = []
            for sp in specs:
                if already_ionic.get():
                    # keep ions as ions; keep s/l/g neutral species intact
                    if sp["charge"] != 0:
                        ion = as_ion_str(_strip_phase_and_ws(sp["core"]), sp["charge"])
                        ionic.append((ion, sp["coeff"]))
                    else:
                        ionic.append((sp["core"] + ("("+sp["phase"]+")" if sp["phase"] else ""), sp["coeff"]))
                else:
                    ds = try_dissociate(sp)
                    if ds is None:
                        ionic.append((sp["core"] + ("("+sp["phase"]+")" if sp["phase"] else ""), sp["coeff"]))
                    else:
                        for ion, n in ds:
                            ionic.append((ion, sp["coeff"]*n))
            return ionic


        def format_ionic(counts):
            if not counts: return "—"
            parts = []
            for ion, n in counts.items():
                if n == 0: continue
                parts.append((n, ion))
            if not parts: return "—"
            return " + ".join([ (f"{n} " if n!=1 else "") + ion for n,ion in parts ])


        def _atoms_and_charge(label: str):
            # label like "Ag^+", "Pt^2+", "SO4^2-", "Ag(s)"
            s = label.strip()
            s, z = _read_charge(s)          # remove ^2+, get charge
            atoms = _tokenize(_strip_phase_and_ws(s))
            return atoms, int(z)

        def _balance_net_ionic(netL: dict, netR: dict):
            """Return scaled (netL, netR) with integer coefficients that balance atoms and charge.
            If something goes wrong, return inputs unchanged."""
            # collect species with side flag
            specs = []
            for sp, n in netL.items():
                atoms, z = _atoms_and_charge(sp)
                specs.append({"side":"L","label":sp,"atoms":atoms,"z":z})
            for sp, n in netR.items():
                atoms, z = _atoms_and_charge(sp)
                specs.append({"side":"R","label":sp,"atoms":atoms,"z":z})

            if not specs:
                return netL, netR

            elements = sorted(set().union(*[s["atoms"].keys() for s in specs]))
            rows = []
            # element conservation rows
            for el in elements:
                row = []
                for s in specs:
                    coef = s["atoms"].get(el,0)
                    row.append(coef if s["side"]=="L" else -coef)
                rows.append(row)
            # charge conservation row
            row = []
            for s in specs:
                coef = s["z"]
                row.append(coef if s["side"]=="L" else -coef)
            rows.append(row)

            # Gauss-Jordan nullspace (like your molecular balancer)
            from fractions import Fraction
            R = len(rows)
            C = len(rows[0]) if rows else 0
            A = [[Fraction(v) for v in r] for r in rows]
            piv = [-1]*R
            r = c = 0
            while r<R and c<C:
                k = next((i for i in range(r,R) if A[i][c]!=0), None)
                if k is None: c += 1; continue
                A[r], A[k] = A[k], A[r]
                piv[r] = c
                pv = A[r][c]
                A[r] = [x/pv for x in A[r]]
                for i in range(R):
                    if i==r: continue
                    f = A[i][c]
                    if f!=0:
                        A[i] = [A[i][j]-f*A[r][j] for j in range(C)]
                r += 1; c += 1

            free = [j for j in range(C) if j not in piv] or [C-1]
            sol = [Fraction(0) for _ in range(C)]
            for j in free: sol[j] = Fraction(1)
            for i in reversed(range(R)):
                pc = piv[i]
                if pc==-1: continue
                s = -sum(A[i][j]*sol[j] for j in range(pc+1,C))
                sol[pc] = s

            # scale to smallest positive integers
            den = 1
            for v in sol: den = math.lcm(den, v.denominator)
            ints = [int(v*den) for v in sol]
            # make positive
            if any(x<0 for x in ints):
                ints = [abs(x) for x in ints]

            # apply back to netL/netR
            i = 0
            outL, outR = {}, {}
            for sp, _ in netL.items():
                outL[sp] = ints[i]; i += 1
            for sp, _ in netR.items():
                outR[sp] = ints[i]; i += 1
            return outL, outR


        # ---------- Compute button ----------
        def go():
            arrow_pat = r"\s*(?:->|→|=>|=)\s*"
            L_raw = eL.get().strip()
            R_raw = eR.get().strip()

            if not R_raw and re.search(arrow_pat, L_raw):
                parts = re.split(arrow_pat, L_raw, maxsplit=1)
                if len(parts) == 2:
                    L_raw, R_raw = parts[0], parts[1]

            try:
                L = parse_side(L_raw)
                R = parse_side(R_raw)
            except Exception as err:
                messagebox.showerror("Parse error", str(err), parent=win); return

            lines = []

            # --- decide whether to balance ---
            has_ionic = any(sp["charge"] != 0 for sp in (L + R))
            do_balance_now = bool(do_balance.get()) and not (already_ionic.get() or has_ionic)

            if do_balance_now:
                try:
                    L, R = balance_molecular(L, R)
                    lines.append("Balanced molecular equation:")
                except Exception:
                    lines.append("(Balancing failed; proceeding as entered)")
            else:
                if already_ionic.get() or has_ionic:
                    lines.append("(Detected ionic species — skipping balancing.)")
                else:
                    lines.append("Molecular equation (as entered):")
            lines.append("  " + format_side(L) + "  →  " + format_side(R))

            # Total ionic (dissociate aq strong electrolytes; keep s, l, g and weak/unknown intact)
            L_ionic = dissociate_side(L)
            R_ionic = dissociate_side(R)
            Lc = count_ionic_side(L_ionic)
            Rc = count_ionic_side(R_ionic)

            lines.append("\nTotal ionic equation (after dissociation):")
            lines.append("  " + format_ionic(Lc) + "  →  " + format_ionic(Rc))

            # Cancel spectators
            spectators = []
            for ion in list(Lc.keys() | Rc.keys()):
                nL = Lc.get(ion, 0); nR = Rc.get(ion, 0)
                m = min(nL, nR)
                if m > 0:
                    spectators.append((ion, m))
                    Lc[ion] = nL - m
                    Rc[ion] = nR - m
            if spectators:
                      lines.append("\nSpectators cancelled: " + ", ".join([ (f"{n} {ion}" if n!=1 else ion) for ion,n in spectators]))
            else:
                lines.append("\nSpectators cancelled: (none)")

            netL = {k:v for k,v in Lc.items() if v!=0}
            netR = {k:v for k,v in Rc.items() if v!=0}

            if netL or netR:
                try:
                    netL, netR = _balance_net_ionic(netL, netR)
                except Exception:
                    pass  # if anything fails, we just print the unscaled counts

            if not netL and not netR:
                lines.append("\nNet-ionic: NR (no reaction)")
            else:
                lines.append("\nNet-ionic equation:")
                lines.append("  " + format_ionic(netL) + "  →  " + format_ionic(netR))

            write(lines)

        ttk.Button(frm, text="Compute net-ionic", command=go).grid(row=7, column=0, sticky="w", pady=(10,0))

        # Prefill with your Ag/Pt(II) example
        eL.insert(0, "Ag(s) + Pt(NO3)2(aq)")
        eR.insert(0, "AgNO3(aq) + Pt(s)")
        out.configure(state="disabled")


    def _open_functional_group_finder(self):
        import re, tkinter as tk
        from tkinter import ttk

        win = tk.Toplevel(self)
        win.title("Functional-group finder (SMILES or condensed)")
        win.transient(self); win.grab_set(); win.geometry("920x560")

        frm = ttk.Frame(win, padding=12); frm.pack(fill="both", expand=True)

        ttk.Label(frm, text="Enter a structure (SMILES or condensed formula):").grid(row=0, column=0, sticky="w")
        ent = ttk.Entry(frm, width=88)
        ent.grid(row=1, column=0, sticky="we", pady=(2,8))
        frm.grid_columnconfigure(0, weight=1)

        # Example picker
        examples = [
            ("p-aminobenzoic acid (Ph–NH2 & –COOH)", "Nc1ccc(cc1)C(=O)O"),
            ("p-fluoronitrobenzene (Ph–F & –NO2)",   "Fc1ccc(cc1)[N+](=O)[O-]"),
            ("aspirin (ester + acid + aryl)",         "CC(=O)OC1=CC=CC=C1C(=O)O"),
            ("ethyl acetate (ester)",                 "CCOC(=O)C"),
            ("acetonitrile (nitrile)",                "CC#N"),
            ("acetone (ketone)",                      "CC(=O)C"),
            ("ethanol (alcohol)",                     "CCO"),
            ("aniline (amine + aryl)",                "Nc1ccccc1"),
            ("nitroethane (nitro)",                   "CC[N+](=O)[O-]"),
            ("phenylacetic acid (Ph + –COOH)",        "c1ccccc1CH2C(=O)O"),
        ]
        pick = ttk.Combobox(frm, values=[name for name,_ in examples], state="readonly", width=56)
        pick.grid(row=2, column=0, sticky="w")
        def load_example(_=None):
            i = pick.current()
            if i >= 0:
                ent.delete(0, tk.END); ent.insert(0, examples[i][1])
                run()
        pick.bind("<<ComboboxSelected>>", load_example)

        # Output pane
        out = tk.Text(frm, height=24, wrap="word")
        out.grid(row=3, column=0, sticky="nsew", pady=(10,0))
        frm.grid_rowconfigure(3, weight=1)
        sv = ttk.Scrollbar(frm, orient="vertical", command=out.yview)
        out.configure(yscrollcommand=sv.set); sv.grid(row=3, column=1, sticky="ns")
        out.configure(state="disabled")

        def write(lines):
            out.configure(state="normal"); out.delete("1.0", "end")
            out.insert("1.0", "\n".join(lines)); out.configure(state="disabled")

        # --- Patterns (SMILES and "condensed" tokens). Very lightweight/heuristic. ---
        # Each entry: (name, [regex_alternatives], hint)
        FG = [
            ("carboxylic acid (–COOH)",
            [r"C\(=O\)O[Hh]\b", r"COOH\b", r"CO2H\b", r"COOH(?=[^A-Za-z]|$)"],
            "protonated carboxyl; acids often written COOH / CO2H"),
            ("carboxylate (–COO⁻)",
            [r"C\(=O\)\[O-]", r"COO-\b", r"CO2-"],
            "deprotonated carboxylic acid"),
            ("ester (–COOR)",
            [r"C\(=O\)O[A-Za-z0-9]", r"COO[A-Za-z]"],
            "acyl–oxygen single bond"),
            ("amide (–CONH– / –CONR–)",
            [r"C\(=O\)N"],
            "peptide/amide linkage"),
            ("aldehyde (–CHO)",
            [r"CHO\b", r"(?<![A-Za-z0-9])O=CH(?![A-Za-z0-9])", r"C=O(?=$)"],
            "terminal carbonyl with –H"),
            ("ketone (–C(=O)–)",
            [r"[A-Za-z0-9]C\(=O\)[A-Za-z0-9]"],
            "internal carbonyl"),
            ("alcohol (–OH / ROH)",
            [r"(?<![A-Za-z])OH(?![A-Za-z])", r"O[Hh]\b", r"\[OH]"],
            "hydroxy group"),
            ("phenyl / aryl ring",
            [r"c1ccccc1", r"(?<![A-Za-z])Ph(?![A-Za-z])", r"(?i)\bphenyl\b", r"C6H5"],
            "benzene ring or Ph"),
            ("ether (–O–)",
            [r"[A-Za-z0-9]O[A-Za-z0-9]"],
            "alkoxy linkage (excludes acid/ester by more specific matches)"),
            ("amine (–NH2 / –NHR / –NR2)",
            [r"(?<!N)\bNH2\b", r"(?<!\[)N(?![+=\[])"],   # crude: avoid NO2, [N+], etc.
            "basic nitrogen (primary/secondary/tertiary)"),
            ("nitrile (–C≡N)",
            [r"C#N", r"C≡N"],
            "cyano group"),
            ("nitro (–NO2 / –N(=O)O–)",
            [r"N\(=O\)O", r"\[N\+\]\(=O\)\[O-]", r"NO2"],
            "nitro substituent"),
            ("sulfonic acid (–SO3H)",
            [r"S\(=O\)\(=O\)O", r"SO3H"],
            "strong acid group"),
            ("thiol (–SH)",
            [r"(?<![A-Za-z])SH(?![A-Za-z])", r"S[Hh]\b"],
            "mercapto group"),
            ("thioether (–S–)",
            [r"[A-Za-z0-9]S[A-Za-z0-9]"],
            "sulfide linkage"),
            ("alkene (C=C)",
            [r"C=C"],
            "double bond"),
            ("alkyne (C≡C)",
            [r"C#C", r"C≡C"],
            "triple bond"),
        ]

        # Halogens counted separately to show X type
        HALO = [
            ("fluoro", r"(?<![a-z])F(?![a-z])"),
            ("chloro", r"Cl"),
            ("bromo",  r"Br"),
            ("iodo",   r"(?<![a-z])I(?![a-z])"),
        ]

        def count_any(patterns, s):
            return sum(len(re.findall(p, s)) for p in patterns)

        def detect(struct):
            s = struct.strip()
            s_nos = s.replace(" ", "")
            s_low = s.lower()

            lines = []
            found = []

            # 1) Specific functional groups
            for name, pats, hint in FG:
                n = count_any(pats, s_nos)
                # Heuristic: avoid double-counting ether when an ester/acid was matched in same region.
                if name.startswith("ether"):
                    # If an obvious ester or acid matched, ether count can be misleading; keep only if
                    # there is an O between carbons and no C(=O)O nearby in the input at all.
                    if re.search(r"C\(=O\)O|COO", s_nos):
                        continue
                if n > 0:
                    found.append((name, n, hint))

            # 2) Halogens
            halo_hits = []
            for label, pat in HALO:
                m = len(re.findall(pat, s_nos))
                if m:
                    halo_hits.append((label, m))
            if halo_hits:
                # present as one line
                parts = [f"{lab} ×{m}" if m>1 else lab for lab,m in halo_hits]
                found.append(("halogen(s)", sum(m for _,m in halo_hits),
                            "substituted halides: " + ", ".join(parts)))

            # Print nicely
            if not found:
                lines.append("No common patterns recognized.\n")
                lines.append("Tips:")
                lines.append("• SMILES works best (e.g., Nc1ccc(cc1)C(=O)O).")
                lines.append("• For condensed input, use tokens like Ph, NO2, COOH, NH2, OH, Cl/Br/F/I, C=C, C#C.")
            else:
                lines.append("Functional groups detected:")
                for (name, n, hint) in found:
                    lines.append(f"  • {name}" + (f"  ×{n}" if n>1 else "") + f"   — {hint}")

            # Small legend
            lines.append("\nLegend (what I look for):")
            lines.append("  –COOH / CO2H → carboxylic acid;  C(=O)O–R → ester;  C(=O)N → amide")
            lines.append("  –CHO (end) → aldehyde;  –C(=O)– → ketone;  –OH → alcohol;  Ph / c1ccccc1 → aryl")
            lines.append("  –NO2 / N(=O)O → nitro;  –C#N → nitrile;  –O– / –S– bridges → ether / thioether")
            lines.append("  C=C → alkene;  C#C → alkyne;  F/Cl/Br/I → halo")

            return lines

        def run():
            txt = ent.get()
            write(detect(txt))

        ttk.Button(frm, text="Detect functional groups", command=run).grid(row=4, column=0, sticky="w", pady=(6,0))
        out.configure(state="disabled")

        # Prefill a useful demo (Ph–NH2, –F, –NO2, –COOH)
        ent.insert(0, "O2N-Ph-COOH with F and NH2  →  Fc1ccc(cc1)N and O2N-Ph-CO2H")

    def _open_structure_drawer(self):
        import math, tkinter as tk
        from tkinter import ttk

        # ---------- window ----------
        win = tk.Toplevel(self)
        win.title("Structure sketch + functional-group detect")
        win.transient(self); win.grab_set(); win.geometry("1100x720")

        # ---------- state ----------
        atom_id_seq = [0]
        atoms = {}   # id -> {"elem": "C","x":float,"y":float,"charge":0,"oval":canvas_id,"label":canvas_id}
        bonds = []   # {"a":id,"b":id,"order":1|2|3,"aromatic":False,"lines":[canvas_ids]}
        selected = {"atom": None}
        R = 16    # atom radius for hit-test/drawing

        # ---------- UI ----------
        root = ttk.Frame(win, padding=10); root.pack(fill="both", expand=True)
        left = ttk.Frame(root); left.pack(side="left", fill="y")
        right = ttk.Frame(root); right.pack(side="right", fill="both", expand=True)

        # tools
        tool = tk.StringVar(value="atom")
        ttk.Label(left, text="Tool").pack(anchor="w")
        ttk.Radiobutton(left, text="Atom",  variable=tool, value="atom").pack(anchor="w")
        ttk.Radiobutton(left, text="Bond",  variable=tool, value="bond").pack(anchor="w")

        ttk.Separator(left).pack(fill="x", pady=6)

        elem = tk.StringVar(value="C")
        ttk.Label(left, text="Element").pack(anchor="w")
        ttk.Combobox(left, textvariable=elem, state="readonly",
                    values=["C","O","N","S","F","Cl","Br","I","H"], width=6).pack(anchor="w")

        ttk.Separator(left).pack(fill="x", pady=6)

        bond_order = tk.IntVar(value=1)
        ttk.Label(left, text="Bond order").pack(anchor="w")
        for n in (1,2,3):
            ttk.Radiobutton(left, text=str(n), variable=bond_order, value=n).pack(anchor="w")
        aromatic = tk.BooleanVar(value=False)
        ttk.Checkbutton(left, text="Aromatic (ring)", variable=aromatic).pack(anchor="w")
        exam_mode = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            left,
            text="Exam-style names (phenol→alcohol, aromatic ring→benzene)",
            variable=exam_mode
        ).pack(anchor="w")

        ttk.Separator(left).pack(fill="x", pady=6)

        ttk.Label(left, text="Templates / fragments").pack(anchor="w")

        # ---------- classifier ----------
        def classify_functional_groups(atoms_in, bonds_in):
            from collections import defaultdict

            # ---------- graph helpers ----------
            adj = defaultdict(list)
            bond_order_map = {}
            aro_edges = set()  # edges explicitly marked aromatic in the drawer

            for i, b in enumerate(bonds_in):
                a, c = b["a"], b["b"]
                o = b.get("order", 1)
                if isinstance(o, str) and o.lower().startswith("arom"):
                    aro_edges.add(frozenset((a, c)))
                    o = 1  # treat as single for generic rules; aromatic handled separately
                adj[a].append((c, i)); adj[c].append((a, i))
                bond_order_map[frozenset((a, c))] = o

            def neighbors(v, el=None):
                for w, idx in adj[v]:
                    if el is None or atoms_in[w]["el"] == el:
                        yield w, idx

            def order(a, b):
                return bond_order_map.get(frozenset((a, b)), 1)

            # ---------- 1) aromatic ring (benzenoid) ----------
            in_aromatic = set()

            def six_cycles():
                # tiny DFS enumerator (good enough for small sketches)
                for start in atoms_in:
                    stack = [start]
                    seen = {start}
                    def dfs(v):
                        if len(stack) == 6:
                            # can we return to start?
                            if any(w == start for (w, _) in adj[v]):
                                yield stack + [start]
                            return
                        for w, _ in adj[v]:
                            if w not in seen and len(stack) < 6:
                                stack.append(w); seen.add(w)
                                yield from dfs(w)
                                seen.remove(w); stack.pop()
                    yield from dfs(start)

            rings_found = []  # store each benzenoid ring as a frozenset of nodes
            for cyc in six_cycles():
                if len(cyc) != 7:
                    continue
                edges = list(zip(cyc[:-1], cyc[1:]))
                dbl = 0; aro = 0; ok = True
                for a, b in edges:
                    if not (atoms_in[a]["el"] == "C" and atoms_in[b]["el"] == "C"):
                        ok = False; break
                    if order(a, b) == 2: dbl += 1
                    if frozenset((a, b)) in aro_edges: aro += 1
                # accept classic 3 doubles OR all 6 edges flagged aromatic
                if ok and (dbl == 3 or aro == 6):
                    ring_nodes = frozenset(cyc[:-1])
                    rings_found.append(ring_nodes)
                    in_aromatic.update(ring_nodes)

            # ---------- collectors ----------
            found = defaultdict(list)
            used_alkene_edges = set()

            def is_carbonyl_carbon(c):
                return (atoms_in[c]["el"] == "C" and
                        any(atoms_in[w]["el"] == "O" and order(c, w) == 2 for w, _ in neighbors(c)))

            # ---------- 2) carbonyl family: acid / ester / amide / aldehyde / ketone ----------
            for c_id, a in atoms_in.items():
                if a["el"] != "C":
                    continue
                O_dbl = [w for w, _ in neighbors(c_id, "O") if order(c_id, w) == 2]
                if len(O_dbl) != 1:
                    continue
                o = O_dbl[0]
                nbrs = [w for w, _ in neighbors(c_id) if w != o]
                single_Os = [w for w in nbrs if atoms_in[w]["el"] == "O" and order(c_id, w) == 1]
                single_Ns = [w for w in nbrs if atoms_in[w]["el"] == "N" and order(c_id, w) == 1]
                has_H     = any(atoms_in[w]["el"] == "H" for w in nbrs)
                carbon_ns = [w for w in nbrs if atoms_in[w]["el"] == "C"]

                if single_Os:
                    o_single = single_Os[0]
                    # O–H?  → carboxylic acid
                    if any(atoms_in[w]["el"] == "H" and order(o_single, w) == 1 for w, _ in neighbors(o_single)):
                        found["carboxylic acid"].append(((c_id, o, o_single),
                                                        f"C{c_id}(=O{o})–O{o_single}–H"))
                        continue
                    # O–C (not the carbonyl carbon)? → ester
                    if any(atoms_in[w]["el"] == "C" and w != c_id for w, _ in neighbors(o_single)):
                        found["ester"].append(((c_id, o, o_single),
                                            f"C{c_id}(=O{o})–O{o_single}–C"))
                        continue

                if single_Ns:
                    found["amide"].append(((c_id, o, single_Ns[0]),
                                        f"C{c_id}(=O{o})–N{single_Ns[0]}"))
                    continue

                # aldehyde vs ketone
                if has_H:
                    found["aldehyde"].append(((c_id, o),
                                            f"C{c_id}=O{o} with H on/implicit at C{c_id}"))
                elif len(carbon_ns) >= 2:
                    found["ketone"].append(((c_id, o),
                                            f"internal C{c_id}=O{o} (two carbon neighbors)"))

            # ---------- 3) ether: R–O–R' (exclude ester context) ----------
            def carbonyl_on(c):
                return is_carbonyl_carbon(c)

            for o_id, a in atoms_in.items():
                if a["el"] != "O": continue
                neigh = list(neighbors(o_id))
                if len(neigh) != 2: continue
                (n1, _), (n2, _) = neigh
                if atoms_in[n1]["el"] == "C" and atoms_in[n2]["el"] == "C":
                    if order(o_id, n1) == 1 and order(o_id, n2) == 1:
                        if not (carbonyl_on(n1) or carbonyl_on(n2)):
                            found["ether"].append(((o_id, n1, n2),
                                                f"O{o_id} single-bonded to C{n1} and C{n2}"))

            # --- 4) ALCOHOL / PHENOL: –O–H  (skip carboxylic –OH) ---
            for o_id, a in atoms_in.items():
                if a["el"] != "O":
                    continue
                c_nb = None; h_nb = None; ok = True
                for w, _ in neighbors(o_id):
                    if order(o_id, w) != 1: ok = False; break
                    if atoms_in[w]["el"] == "H": h_nb = w
                    elif atoms_in[w]["el"] == "C": c_nb = w
                if not ok or not c_nb or not h_nb:
                    continue
                # don't miscount –C(=O)–OH as alcohol (it's a carboxylic acid)
                if (atoms_in[c_nb]["el"] == "C" and
                    any(atoms_in[w]["el"] == "O" and order(c_nb, w) == 2 and w != o_id
                        for w, _ in neighbors(c_nb))):
                    continue
                if c_nb in in_aromatic:
                    found["phenol"].append(((o_id, c_nb, h_nb),
                                            f"Ar–O{o_id}–H{h_nb} (C{c_nb} aromatic)"))
                else:
                    found["alcohol"].append(((o_id, c_nb, h_nb),
                                            f"C{c_nb}–O{o_id}–H{h_nb}"))

            # ---------- 5) amine (exclude amide) ----------
            for n_id, a in atoms_in.items():
                if a["el"] != "N": 
                    continue
                # attached to any carbonyl carbon? then it's an amide (already counted)
                if any(is_carbonyl_carbon(c) and order(n_id, c) == 1 for c, _ in neighbors(n_id, "C")):
                    continue
                # simple amine if N–C single bonds exist
                if any(atoms_in[w]["el"] == "C" and order(n_id, w) == 1 for w, _ in neighbors(n_id)):
                    found["amine"].append(((n_id,),
                                        f"N{n_id} bound to carbon(s) and not to a carbonyl"))

            # ---------- 6) halogens (F/Cl/Br/I) ----------
            halo_names = {"F":"fluoro","Cl":"chloro","Br":"bromo","I":"iodo"}
            halo_counts = defaultdict(int)
            for x_id, a in atoms_in.items():
                el = a["el"]
                if el in halo_names:
                    # count when single-bonded to a carbon (alkyl/vinyl/aryl halide)
                    if any(atoms_in[w]["el"]=="C" and order(x_id, w)==1 for w, _ in neighbors(x_id)):
                        halo_counts[halo_names[el]] += 1
            for label, cnt in halo_counts.items():
                found[label].append(((label,), f"{label} substituent(s)"))

            # ---------- 7) alkene outside aromatic ----------
            for key, o in bond_order_map.items():
                if o != 2: 
                    continue
                a, b = tuple(key)
                if a in in_aromatic and b in in_aromatic:
                    continue
                if atoms_in[a]["el"] == "C" and atoms_in[b]["el"] == "C":
                    used_alkene_edges.add((a, b))
            for a, b in sorted(used_alkene_edges):
                found["alkene"].append(((a, b), f"C{a}=C{b} outside aromatic ring"))

            # ---------- 8) phenyl substituent(s) (aryl attachment) ----------
            if in_aromatic:
                # count ring carbons that connect by a single bond to an external carbon
                attachments = set()
                for r in in_aromatic:
                    for w, _ in neighbors(r, "C"):
                        if w not in in_aromatic and order(r, w) == 1:
                            attachments.add((min(r, w), max(r, w)))
                if attachments:
                    found["phenyl substituent"].append((tuple(sorted(list(a for pair in attachments for a in pair))),
                                                        f"{len(attachments)} aryl (phenyl) attachment point(s)"))
                # also report the ring itself once
                found["aromatic ring"].append((tuple(sorted(in_aromatic)),
                                            "6-membered benzenoid ring"))

            # ---------- pretty report ----------
            lines = ["Functional groups detected:"]
            display_order = [
                "aromatic ring", "phenyl substituent",
                "carboxylic acid", "ester", "amide",
                "aldehyde", "ketone",
                "phenol", "alcohol", "ether", "amine",
                "alkene",
                "fluoro","chloro","bromo","iodo",
            ]
            any_found = False
            for k in display_order:
                if k in found:
                    lines.append(f"  • {k} ×{len(found[k])}")
                    any_found = True
            if not any_found:
                lines.append("  • none")

            lines.append("\nWhy we classified them:")
            for k in display_order:
                if k in found:
                    for atoms_tuple, msg in found[k]:
                        lines.append(f"  – {k}: {msg}")

            return lines, found
        
        
        def attach_Ph():
            if selected["atom"] is None: return
            a = selected["atom"]
            x0, y0 = atoms[a]["x"], atoms[a]["y"]
            # small ring to the right of the selected atom
            r = 40
            import math
            ang0 = math.radians(0)
            pts = [(x0 + 70 + r*math.cos(ang0 + k*math.pi/3),
                    y0 + r*math.sin(ang0 + k*math.pi/3)) for k in range(6)]
            ring = [add_atom("C", x, y) for (x,y) in pts]
            # alternating bonds; mark as aromatic so it’s recognized robustly
            for i in range(6):
                add_bond(ring[i], ring[(i+1)%6], 1 if i%2 else 2, aromatic=True)
            # attach ring carbon 0 back to the selected atom
            add_bond(a, ring[0], 1)
        ttk.Button(left, text="Attach –Ph (phenyl)", command=attach_Ph).pack(anchor="w", pady=(4,8))



                # ---------- naming: acyclic hydrocarbons (alkane/alkene/alkyne) ----------
        def _name_hydrocarbon_from_canvas():
            # Build carbon-only graph, keep bond orders (1/2/3). Ignore aromatic edges.
            carbons = [i for i,a in atoms.items() if a["elem"]=="C"]
            if not carbons:
                say(["Draw some carbons first."]); return

            adj = {i: [] for i in carbons}     # neighbors (carbon only)
            ordmap = {}                        # frozenset({a,b}) -> 1/2/3
            hetero = False
            aromatic_present = any(bd.get("aromatic") for bd in bonds)

            for bd in bonds:
                a,b = bd["a"], bd["b"]
                ea, eb = atoms[a]["elem"], atoms[b]["elem"]
                # ignore bonds to H; reject hetero atoms
                if ea=="C" and eb=="C":
                    if not bd.get("aromatic"):
                        o = int(bd["order"]) if isinstance(bd["order"], (int,float)) else 1
                        adj[a].append(b); adj[b].append(a); ordmap[frozenset((a,b))]=int(o)
                else:
                    if ea not in ("C","H") or eb not in ("C","H"):
                        hetero = True

            if hetero or aromatic_present:
                say(["Naming supports **acyclic hydrocarbons only** (no hetero atoms, no aromatics)."])
                return

            # ---- choose component with most carbons
            seen=set(); comps=[]
            for s in carbons:
                if s in seen: continue
                stack=[s]; seen.add(s); comp=[]
                while stack:
                    v=stack.pop(); comp.append(v)
                    for w in adj[v]:
                        if w not in seen:
                            seen.add(w); stack.append(w)
                comps.append(comp)
            nodes = max(comps, key=len)

            # reject cycles (we only do trees right now)
            edge_count = sum(len(adj[i]) for i in nodes)//2
            if edge_count != len(nodes)-1:
                say(["Cyclo structure detected — cyclo naming not implemented (yet)."])
                return

            # ---- longest weighted path in a tree (prefer multiple bonds)
            def edge_weight(a,b):
                o = ordmap.get(frozenset((a,b)),1)
                return (1000 if o>1 else 0) + 1  # prioritize unsaturations, then length

            def farthest_from(src):
                best_len = 0
                best_node = src
                parent = {src: None}
                def dfs(u, p, acc):
                    nonlocal best_len, best_node
                    if acc > best_len:
                        best_len, best_node = acc, u
                    for v in adj[u]:
                        if v == p: continue
                        parent[v] = u
                        dfs(v, u, acc + edge_weight(u,v))
                dfs(src, None, 0)
                return best_node, parent

            u,_parent0 = farthest_from(nodes[0])
            v,parent = farthest_from(u)

            # rebuild the chain u..v via parent map
            chain=[v]
            while chain[-1] != u:
                nxt = parent.get(chain[-1])
                if nxt is None: break  # defensive; shouldn't happen in a tree
                chain.append(nxt)
            chain=list(reversed(chain))
            n=len(chain)
            chain_set=set(chain)

            # ---- collect unsaturations (double/triple) along the chain
            dbl_pos=[]; trp_pos=[]
            for i in range(n-1):
                a,b = chain[i], chain[i+1]
                o = ordmap.get(frozenset((a,b)),1)
                if o==2: dbl_pos.append(i+1)   # locant = lower index (1-based)
                elif o==3: trp_pos.append(i+1)

            # ---- collect simple alkyl substituents (sizes off the chain)
            from collections import deque, defaultdict
            def branch_size(start, parent_node):
                q=deque([start]); visited={start,parent_node}; size=0; is_simple=True
                while q:
                    x=q.popleft(); size+=1
                    off=[w for w in adj[x] if w not in visited and w not in chain_set]
                    if len(off)>1: is_simple=False
                    for w in off:
                        visited.add(w); q.append(w)
                return size, is_simple

            subs=[]; simple_ok=True
            for idx,c in enumerate(chain):
                for w in adj[c]:
                    if w not in chain_set:
                        sz, ok = branch_size(w, c)
                        simple_ok &= ok
                        subs.append((idx+1, sz))  # (locant, size)

            # ---- choose numbering: lowest locants for multiple bonds; if tie -> substituents
            def rev_locs(locs): return sorted(n - p for p in locs)
            A = (sorted(dbl_pos), sorted(trp_pos))
            B = (rev_locs(dbl_pos), rev_locs(trp_pos))
            use_reverse = False
            if B < A:
                use_reverse=True
            elif B == A:
                # tiebreak by substituent locants
                subsA = sorted(p for p,_ in subs)
                subsB = sorted(n + 1 - p for p,_ in subs)
                if subsB < subsA:
                    use_reverse = True

            if use_reverse:
                chain=list(reversed(chain))
                chain_set=set(chain)
                dbl_pos = rev_locs(dbl_pos)
                trp_pos = rev_locs(trp_pos)
                subs = [(n+1-p, sz) for (p,sz) in subs]

            # ---- build the name
            root = {1:"meth",2:"eth",3:"prop",4:"but",5:"pent",6:"hex",7:"hept",8:"oct",
                    9:"non",10:"dec",11:"undec",12:"dodec"}.get(n, f"{n}-carbon")
            mult = {2:"di",3:"tri",4:"tetra",5:"penta",6:"hexa",7:"hepta",8:"octa"}
            alkyl = {1:"methyl",2:"ethyl",3:"propyl",4:"butyl",5:"pentyl",6:"hexyl",
                    7:"heptyl",8:"octyl",9:"nonyl",10:"decyl"}

            # substituents
            groups=defaultdict(list)
            for p,sz in subs:
                nm = alkyl.get(sz, f"{sz}-alkyl")
                groups[nm].append(p)
            sub_parts=[]
            for nm in sorted(groups):  # alphabetical
                locs=",".join(str(x) for x in sorted(groups[nm]))
                cnt=len(groups[nm]); prefix = mult.get(cnt,"")
                sub_parts.append(f"{locs}-" + (prefix+nm if prefix else nm))
            sub_prefix = ("-".join(sub_parts)+"-") if sub_parts else ""

            # unsaturation suffix/infix
            end_en=end_yn=""
            if dbl_pos:
                en_suffix = "ene" if len(dbl_pos)==1 else mult[len(dbl_pos)] + "ene"
                end_en = f"{'-'.join(map(str,sorted(dbl_pos)))}-{en_suffix}"
            if trp_pos:
                yn_suffix = "yne" if len(trp_pos)==1 else mult[len(trp_pos)] + "yne"
                end_yn = f"{'-'.join(map(str,sorted(trp_pos)))}-{yn_suffix}"

            if dbl_pos and trp_pos:
                base = f"{root}-{end_en}-{end_yn}"
                base_da = base.replace("ane","an").replace("ene","en").replace("yne","yn")
            elif dbl_pos:
                base = f"{root}-{end_en}"; base_da = base.replace("ene","en")
            elif trp_pos:
                base = f"{root}-{end_yn}"; base_da = base.replace("yne","yn")
            else:
                base = root + "ane"; base_da = root + "an"

            final_name = sub_prefix + base
            say([f"Hydrocarbon (experimental): {final_name}",
                f"Dansk: {sub_prefix + base_da}",
                "",
                f"Main chain length = {n}",
                ("Double-bond locants: " + ", ".join(map(str,sorted(dbl_pos)))) if dbl_pos else "No C=C",
                ("Triple-bond locants: " + ", ".join(map(str,sorted(trp_pos)))) if trp_pos else "No C≡C",
                ("Substituents: " + ", ".join(sub_parts)) if sub_parts else "No substituents",
                *(["(Substituent shapes seem branched — named as simple alkyls only.)"] if not simple_ok else [])
            ])


        ttk.Button(right, text="Name hydrocarbon (experimental)",
                   command=_name_hydrocarbon_from_canvas).pack(pady=(0,6))


        # ---- template buttons ----
        def add_benzene():
            cx = last_xy[0] if last_xy[0] is not None else (canvas.winfo_width()/2)
            cy = last_xy[1] if last_xy[1] is not None else (canvas.winfo_height()/2)
            r = 70
            ang0 = math.radians(90)
            pts = [(cx + r*math.cos(ang0 + k*math.pi/3), cy - r*math.sin(ang0 + k*math.pi/3)) for k in range(6)]
            atom_ids = [add_atom("C", x, y) for (x,y) in pts]
            # alternating double; also mark aromatic
            for i in range(6):
                a, b = atom_ids[i], atom_ids[(i+1)%6]
                add_bond(a, b, 1 if i%2 else 2, aromatic=True)

        ttk.Button(left, text="Insert benzene ring", command=add_benzene).pack(anchor="w", pady=(0,4))

        def attach_OH():
            if selected["atom"] is None: return
            a = selected["atom"]; x,y = atoms[a]["x"], atoms[a]["y"]
            o = add_atom("O", x+36, y-10)
            h = add_atom("H", x+66, y-10)
            add_bond(a, o, 1); add_bond(o, h, 1)

        def attach_NH2():
            if selected["atom"] is None: return
            a = selected["atom"]; x,y = atoms[a]["x"], atoms[a]["y"]
            n = add_atom("N", x+36, y+6)
            h1 = add_atom("H", x+64, y-6); h2 = add_atom("H", x+64, y+18)
            add_bond(a, n, 1); add_bond(n, h1, 1); add_bond(n, h2, 1)

        def attach_COOH():
            if selected["atom"] is None: return
            a = selected["atom"]; x,y = atoms[a]["x"], atoms[a]["y"]
            c = add_atom("C", x+38, y)
            o1 = add_atom("O", x+70, y-16)
            o2 = add_atom("O", x+70, y+16)
            h  = add_atom("H", x+94, y+16)
            add_bond(a, c, 1); add_bond(c, o1, 2); add_bond(c, o2, 1); add_bond(o2, h, 1)

        ttk.Button(left, text="Attach –OH",   command=attach_OH).pack(anchor="w", pady=(4,0))
        ttk.Button(left, text="Attach –NH₂",  command=attach_NH2).pack(anchor="w", pady=(4,0))
        ttk.Button(left, text="Attach –COOH", command=attach_COOH).pack(anchor="w", pady=(4,8))

        ttk.Separator(left).pack(fill="x", pady=6)

        def clear_all():
            canvas.delete("all"); atoms.clear(); bonds.clear()
            selected["atom"] = None
        ttk.Button(left, text="Clear", command=clear_all).pack(anchor="w")

        ttk.Separator(left).pack(fill="x", pady=8)

        ttk.Label(left, text="How to use").pack(anchor="w", pady=(0,2))
        tips = tk.Text(left, width=34, height=14, wrap="word")
        tips.insert("1.0",
            "• Atom tool: left-click to place the selected element.\n"
            "• Bond tool: click first atom, then second atom (set order 1/2/3).\n"
            "• ‘Aromatic’ toggles dashed bonds used to flag benzene rings.\n"
            "• Right-click an atom to delete it with its bonds.\n"
            "• Use fragments or draw freely, then click ‘Detect functional groups’.\n"
            "• For aldehydes, draw the H on the carbonyl carbon (–CHO).\n"
        )
        tips.configure(state="disabled"); tips.pack(anchor="w")

        # canvas + output
        right_top = ttk.Frame(right); right_top.pack(side="top", fill="both", expand=True)
        canvas = tk.Canvas(right_top, bg="#ffffff", height=470); canvas.pack(fill="both", expand=True)

        right_bottom = ttk.Frame(right); right_bottom.pack(side="bottom", fill="x")
        out = tk.Text(right_bottom, height=10, wrap="word"); out.pack(side="left", fill="both", expand=True, pady=(6,0))
        sc = ttk.Scrollbar(right_bottom, orient="vertical", command=out.yview); sc.pack(side="right", fill="y", pady=(6,0))
        out.configure(yscrollcommand=sc.set, state="disabled")

        def say(lines):
            out.configure(state="normal"); out.delete("1.0","end")
            out.insert("1.0","\n".join(lines)); out.configure(state="disabled")

        # ---------- drawing helpers ----------
        def add_atom(sym, x, y, charge=0):
            aid = atom_id_seq[0]; atom_id_seq[0] += 1
            o = canvas.create_oval(x-R, y-R, x+R, y+R, fill="#f7f7f7", outline="#444", width=1.5)
            t = canvas.create_text(x, y, text=sym, font=("Segoe UI", 11, "bold"))
            atoms[aid] = {"elem":sym, "x":x, "y":y, "charge":charge, "oval":o, "label":t}
            canvas.tag_bind(o, "<Button-1>", lambda e, i=aid: select_atom(i))
            canvas.tag_bind(t, "<Button-1>", lambda e, i=aid: select_atom(i))
            canvas.tag_bind(o, "<Button-3>", lambda e, i=aid: delete_atom(i))
            canvas.tag_bind(t, "<Button-3>", lambda e, i=aid: delete_atom(i))
            return aid

        def add_bond(a, b, order=1, aromatic=False):
            if a==b: return
            # prevent duplicate
            for bd in bonds:
                if {bd["a"], bd["b"]} == {a,b}:
                    bd["order"] = order; bd["aromatic"]=aromatic
                    redraw_bond(bd); return
            bd = {"a":a,"b":b,"order":order,"aromatic":aromatic,"lines":[]}
            bonds.append(bd); redraw_bond(bd)

        def redraw_bond(bd):
            for i in bd["lines"]: canvas.delete(i)
            bd["lines"].clear()
            ax, ay = atoms[bd["a"]]["x"], atoms[bd["a"]]["y"]
            bx, by = atoms[bd["b"]]["x"], atoms[bd["b"]]["y"]
            vx, vy = bx-ax, by-ay
            L = math.hypot(vx,vy) or 1.0
            nx, ny = -vy/L, vx/L      # perpendicular
            offset = 3.5
            orders = bd["order"]
            if bd["aromatic"]:
                line = canvas.create_line(ax, ay, bx, by, width=2, dash=(4,3), fill="#333")
                bd["lines"].append(line)
                orders = max(1, orders)
            if orders == 1:
                line = canvas.create_line(ax, ay, bx, by, width=2, fill="#333")
                bd["lines"].append(line)
            elif orders == 2:
                for s in (-offset, +offset):
                    line = canvas.create_line(ax+nx*s, ay+ny*s, bx+nx*s, by+ny*s, width=2, fill="#333")
                    bd["lines"].append(line)
            else:  # 3
                line = canvas.create_line(ax, ay, bx, by, width=2, fill="#333")
                bd["lines"].append(line)
                for s in (-2*offset, +2*offset):
                    line = canvas.create_line(ax+nx*s, ay+ny*s, bx+nx*s, by+ny*s, width=2, fill="#333")
                    bd["lines"].append(line)

        def select_atom(i):
            if selected["atom"] is not None:
                canvas.itemconfigure(atoms[selected["atom"]]["oval"], outline="#444", width=1.5)
            selected["atom"] = i
            canvas.itemconfigure(atoms[i]["oval"], outline="#2b7", width=3)

        def delete_atom(i):
            for bd in list(bonds):
                if bd["a"]==i or bd["b"]==i:
                    for k in bd["lines"]: canvas.delete(k)
                    bonds.remove(bd)
            canvas.delete(atoms[i]["oval"]); canvas.delete(atoms[i]["label"])
            atoms.pop(i, None)
            if selected["atom"] == i: selected["atom"] = None

        def nearest_atom(x, y, tol=R+6):
            best, dmin = None, tol
            for i,a in atoms.items():
                d = math.hypot(a["x"]-x, a["y"]-y)
                if d < dmin:
                    best, dmin = i, d
            return best

        # canvas interactions
        last_xy = [None, None]
        bond_first = {"atom": None}
        def on_click(e):
            last_xy[0], last_xy[1] = e.x, e.y
            if tool.get()=="atom":
                aid = add_atom(elem.get(), e.x, e.y); select_atom(aid)
            else:
                i = nearest_atom(e.x, e.y)
                if i is None: return
                if bond_first["atom"] is None:
                    bond_first["atom"] = i; select_atom(i)
                else:
                    a, b = bond_first["atom"], i
                    add_bond(a, b, bond_order.get(), aromatic.get())
                    bond_first["atom"] = None
        canvas.bind("<Button-1>", on_click)

        # ---------- detection ----------
        def detect_groups():
            # adapt canvas data -> classifier input
            atoms2 = {i: {"el": a["elem"]} for i, a in atoms.items()}
            bonds2 = [
                {"a": bd["a"], "b": bd["b"],
                "order": "aromatic" if bd.get("aromatic") else bd["order"]}
                for bd in bonds
            ]

            # run classifier
            lines_detailed, found = classify_functional_groups(atoms2, bonds2)

            if exam_mode.get():
                # ---- exam-style summary (merge + rename) ----
                # copy so we can mutate
                show = {k: list(v) for k, v in found.items()}

                # phenol counted as alcohol
                if "phenol" in show:
                    show.setdefault("alcohol", []).extend(show.pop("phenol"))

                # aromatic ring labeled as benzene
                if "aromatic ring" in show:
                    show["benzene"] = show.pop("aromatic ring")

                # optional: keep only the buckets expected by the answer key
                order = ("ether", "alcohol", "benzene", "aldehyde")

                lines = ["Functional groups detected:"]
                any_hit = False
                for k in order:
                    if k in show and show[k]:
                        lines.append(f"  • {k} ×{len(show[k])}")
                        any_hit = True
                if not any_hit:
                    lines.append("  • none")

                # still show the detailed reasoning so you can verify the hits
                lines.append("\nWhy we classified them:")
                for k, L in found.items():
                    for _ids, msg in L:
                        lines.append(f"  – {k}: {msg}")

                say(lines)

            else:
                # original detailed output
                say(lines_detailed)


        ttk.Button(right, text="Detect functional groups", command=detect_groups).pack(pady=6)

        # initial canvas size update
        win.update_idletasks()

    
    def _open_nanoparticle_atoms_tool(self):
        import tkinter as tk
        from tkinter import ttk, messagebox
        import math

        win = tk.Toplevel(self)
        win.title("Atoms in nanoparticle — size & density")
        win.transient(self); win.grab_set()
        win.geometry("760x540")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        # ---------- Substance with suggestions (reuses your pattern) ----------
        ttk.Label(frm, text="Substance (name/symbol/formula):").grid(row=0, column=0, sticky="w")
        ent = ttk.Entry(frm, width=36); ent.grid(row=0, column=1, columnspan=3, sticky="we", padx=(6,8))
        ent.insert(0, "platinum")

        sugg = tk.Listbox(frm, height=0); sugg.grid(row=1, column=0, columnspan=4, sticky="we")
        scx = ttk.Scrollbar(frm, orient="horizontal", command=sugg.xview)
        sugg.configure(xscrollcommand=scx.set); scx.grid(row=2, column=0, columnspan=4, sticky="we", pady=(0,6))
        sugg.grid_remove(); scx.grid_remove()

        # Live preview (formula & molar mass)
        preview = tk.StringVar(value="")
        mm_info = tk.StringVar(value="")
        ttk.Label(frm, textvariable=preview, font=("Segoe UI", 10, "bold")).grid(row=3, column=0, columnspan=4, sticky="w", pady=(2,2))
        ttk.Label(frm, textvariable=mm_info, foreground="#444").grid(row=4, column=0, columnspan=4, sticky="w")

        def _refresh_suggestions(_=None):
            from molar_masses import suggest_substances, name_to_formula, molar_mass
            items = suggest_substances(ent.get(), limit=60)
            sugg.delete(0, tk.END)
            for label, f in items:
                sugg.insert(tk.END, label)
            if items:
                sugg.configure(height=min(8, len(items))); sugg.grid(); scx.grid()
            else:
                sugg.grid_remove(); scx.grid_remove()
            # live preview
            f = name_to_formula(ent.get().strip())
            preview.set(f"Formula: {f}")
            try:
                mm = molar_mass(f)
                mm_info.set(f"Molar mass = {mm:g} g/mol")
            except Exception:
                mm_info.set("Molar mass = (unknown — check name/formula)")

        def _accept_suggestion(_=None):
            sel = sugg.curselection()
            if not sel: return
            name = sugg.get(sel[0]).split(" — ", 1)[0]
            ent.delete(0, tk.END); ent.insert(0, name)
            sugg.grid_remove(); scx.grid_remove()
            _refresh_suggestions()

        ent.bind("<KeyRelease>", _refresh_suggestions)
        ent.bind("<Down>", lambda e: (sugg.focus_set(), sugg.selection_clear(0, tk.END), sugg.selection_set(0), sugg.activate(0)) if sugg.winfo_ismapped() else None)
        sugg.bind("<Return>", _accept_suggestion)
        sugg.bind("<Double-Button-1>", _accept_suggestion)
        sugg.bind("<Escape>", lambda e: (sugg.grid_remove(), scx.grid_remove()))

        # ---------- Geometry / size ----------
        size_box = ttk.LabelFrame(frm, text="Particle geometry & size")
        size_box.grid(row=5, column=0, columnspan=4, sticky="we", pady=(10,6))
        for c in range(4): size_box.grid_columnconfigure(c, weight=1)

        shape = tk.StringVar(value="sphere")
        ttk.Radiobutton(size_box, text="Sphere", variable=shape, value="sphere").grid(row=0, column=0, sticky="w")

        mode = tk.StringVar(value="diam")  # diam or rad
        ttk.Radiobutton(size_box, text="Diameter", variable=mode, value="diam").grid(row=1, column=0, sticky="w")
        ttk.Radiobutton(size_box, text="Radius",   variable=mode, value="rad").grid(row=1, column=1, sticky="w")

        ttk.Label(size_box, text="Size value:").grid(row=2, column=0, sticky="w", pady=(4,0))
        val_size = ttk.Entry(size_box, width=12); val_size.grid(row=2, column=1, sticky="w", padx=(6,12), pady=(4,0))
        val_size.insert(0, "5")  # nice default

        unit_size = ttk.Combobox(size_box, width=8, state="readonly",
                                values=["nm", "Å", "pm", "µm", "m"])
        unit_size.grid(row=2, column=2, sticky="w", pady=(4,0))
        unit_size.set("nm")

        # ---------- Density ----------
        dens_box = ttk.LabelFrame(frm, text="Density")
        dens_box.grid(row=6, column=0, columnspan=4, sticky="we", pady=(6,6))
        ttk.Label(dens_box, text="Value:").grid(row=0, column=0, sticky="w")
        dens_val = ttk.Entry(dens_box, width=12); dens_val.grid(row=0, column=1, sticky="w", padx=(6,12))
        dens_val.insert(0, "21.43")  # Pt example

        dens_unit = ttk.Combobox(dens_box, width=10, state="readonly",
                                values=["g/cm³", "kg/m³"])
        dens_unit.grid(row=0, column=2, sticky="w")
        dens_unit.set("g/cm³")

        ttk.Label(dens_box, text="(Tip: for Pt use 21.43 g/cm³)").grid(row=0, column=3, sticky="w", padx=(10,0))

        # ---------- Output ----------
        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 11, "bold")).grid(row=7, column=0, columnspan=4, sticky="w", pady=(10,4))

        steps = tk.Text(frm, height=8, wrap="word")
        steps.grid(row=8, column=0, columnspan=4, sticky="nsew", pady=(4,0))
        frm.grid_rowconfigure(8, weight=1)
        sv = ttk.Scrollbar(frm, orient="vertical", command=steps.yview)
        steps.configure(yscrollcommand=sv.set); sv.grid(row=8, column=4, sticky="ns")
        steps.configure(state="disabled")

        def _write_steps(lines):
            steps.configure(state="normal")
            steps.delete("1.0", "end")
            steps.insert("1.0", "\n".join(lines))
            steps.configure(state="disabled")

        # ---------- Compute ----------
        def compute():
            from molar_masses import name_to_formula, molar_mass
            from constants import CONSTANTS

            # parse molar mass
            f = name_to_formula(ent.get().strip())
            try:
                MM = molar_mass(f)  # g/mol
            except Exception as ex:
                messagebox.showerror("Invalid substance", f"{ex}", parent=win); return

            # parse size
            try:
                s = float(val_size.get().strip().replace(",", "."))
            except ValueError:
                messagebox.showerror("Bad size", "Enter a numeric size.", parent=win); return
            if s <= 0:
                messagebox.showerror("Bad size", "Size must be positive.", parent=win); return

            # to centimeters (because density default is g/cm³)
            unit = unit_size.get()
            to_cm = {"nm": 1e-7, "Å": 1e-8, "pm": 1e-10, "µm": 1e-4, "m": 100.0}
            if unit not in to_cm:
                messagebox.showerror("Unit", "Choose a valid size unit.", parent=win); return
            length_cm = s * to_cm[unit]
            r_cm = (length_cm/2.0) if mode.get()=="diam" else length_cm

            # volume (sphere)
            V_cm3 = (4.0/3.0) * math.pi * (r_cm**3)

            # density to g/cm³
            try:
                rho = float(dens_val.get().strip().replace(",", "."))
            except ValueError:
                messagebox.showerror("Bad density", "Enter a numeric density.", parent=win); return
            rho_g_cm3 = rho if dens_unit.get()=="g/cm³" else rho * 0.001  # 1 kg/m³ = 0.001 g/cm³

            # mass (g), moles, atoms
            m_g = rho_g_cm3 * V_cm3
            n_mol = m_g / MM if MM > 0 else 0.0
            N_A = CONSTANTS.get("N_A", 6.02214076e23)
            atoms = n_mol * N_A

            out.set(f"≈ {atoms:.3g} atoms")
            _write_steps([
                f"Given: {f}   MM = {MM:g} g/mol",
                f"Size: {'diameter' if mode.get()=='diam' else 'radius'} = {s:g} {unit}  →  r = {r_cm:.6g} cm",
                f"Volume (sphere): V = 4/3·π·r³ = {V_cm3:.6g} cm³",
                f"Density: ρ = {rho_g_cm3:g} g/cm³   → mass m = ρV = {m_g:.6g} g",
                f"Moles: n = m/MM = {n_mol:.6g} mol",
                f"Atoms: N = n·N_A = {n_mol:.6g}·{N_A:.6g} = {atoms:.6g}"
            ])
            return atoms

        # Buttons
        btns = ttk.Frame(frm); btns.grid(row=9, column=0, columnspan=4, sticky="w", pady=(10,0))
        ttk.Button(btns, text="Compute", command=compute).pack(side=tk.LEFT, padx=(0,8))
        ttk.Button(btns, text="Reset",
                command=lambda: (ent.delete(0, tk.END), ent.insert(0, "platinum"),
                                    val_size.delete(0, tk.END), val_size.insert(0, "5"),
                                    unit_size.set("nm"), mode.set("diam"),
                                    dens_val.delete(0, tk.END), dens_val.insert(0, "21.43"),
                                    dens_unit.set("g/cm³"), out.set(""), _write_steps([]))
                ).pack(side=tk.LEFT)

        # Prepopulate preview and suggestions
        _refresh_suggestions()

    def _open_unpaired_d_electrons_tool(self):
        import re, tkinter as tk
        from tkinter import ttk, messagebox

        # ---------- små tabeller ----------
        # --- Transition metal "group numbers" (bruges til d^n = GROUP_NUM[M] - oxidation) ---
        GROUP_NUM = {
            # 3d
            "Sc":3, "Ti":4, "V":5, "Cr":6, "Mn":7, "Fe":8, "Co":9, "Ni":10, "Cu":11, "Zn":12,
            # 4d
            "Y":3, "Zr":4, "Nb":5, "Mo":6, "Tc":7, "Ru":8, "Rh":9, "Pd":10, "Ag":11, "Cd":12,
            # 5d
            "La":3, "Hf":4, "Ta":5, "W":6, "Re":7, "Os":8, "Ir":9, "Pt":10, "Au":11, "Hg":12,
            # valgfrit (ofte sat i gruppe 3; god at have med)
            "Lu":3
        }

        # --- Ligand-aliaser (normaliser brugerinput) ---
        LIGAND_ALIAS = {
            # navne ↔ formler
            "aqua":"H2O", "water":"H2O",
            "ammine":"NH3", "amine":"NH3",
            "pyridine":"py", "pyr":"py",
            "bipy":"bpy", "bpy":"bpy", "phen":"phen", "terpy":"terpy",
            "acac-":"acac", "acetylacetonate":"acac",
            "acetate":"OAc", "ch3coo":"OAc", "ac-":"OAc", "oac":"OAc",
            "oxalate":"ox", "c2o4":"ox",
            "nitrate":"NO3", "carbonate":"CO3", "sulfate":"SO4", "perchlorate":"ClO4",
            "thiocyanate":"SCN", "ncs":"SCN", "scn":"SCN",
            "cyanate":"OCN", "nco":"NCO", "ocn":"OCN",
            "azide":"N3", "hydride":"H",
            "acetonitrile":"MeCN", "ch3cn":"MeCN",
            "dimethylsulfoxide":"DMSO", "dmso":"DMSO",
        }

        def normalize_ligand(L: str) -> str:
            key = L.strip()
            key = key.replace(" ", "")
            key_u = key.upper()
            # bevar case for “py”, “en”, “bpy”, “phen” osv.
            return LIGAND_ALIAS.get(key.lower(), LIGAND_ALIAS.get(key_u, key))

        # --- Typiske ligandladninger (per koordinationsenhed) ---
        # Bemærk: NO (nitrosyl) kan være NO+ (lineær, stærk-felt) eller NO− (bøjet, svagere).
        # Her tager vi kun sikre standarder; tvivlsomme systemer håndteres som “borderline”.
        LIGAND_CHARGE = {
            # halider & O-donorer
            "F":-1, "Cl":-1, "Br":-1, "I":-1,
            "OH":-1, "O":-2,      # oxo
            "H2O":0, "MeOH":0, "EtOH":0,
            "NO3":-1, "CO3":-2, "SO4":-2, "ClO4":-1,
            "OAc":-1, "acac":-1, "ox":-2, "gly":-1,      # acetat, acetylacetonat, oxalat, glycinat
            # N/C-donorer
            "NH3":0, "py":0, "Im":0,                     # imidazol
            "CN":-1, "NCO":-1, "OCN":-1, "N3":-1, "SCN":-1,
            "MeCN":0,
            # P-donorer
            "PPh3":0, "PR3":0, "PF3":0,
            # S-donorer
            "DMSO":0, "S2":-2, "SH":-1,                  # simplificeret: thiolat SH−
            # stærke π-acceptorer
            "CO":0,
            # polypyridyl & chelater
            "en":0, "dien":0, "terpy":0, "bpy":0, "phen":0,
            "EDTA":-4
            # (tilføj nemt flere efter behov)
        }

        # --- Spektrokemisk feltstyrke (meget grov auto-klassifikation) ---
        # Bruges kun til at vælge HIGH/LOW spin automatisk. Blandinger → "ambiguous".
        WEAK = {
            "I","Br","Cl","F","OH","H2O","S2","SH","SCN","NO3","CO3","SO4","ClO4","OAc"
        }
        BORDERLINE = {
            # typisk midt i serien eller følsom for metal/ox-stadie
            "NH3","py","Im","MeCN","DMSO","acac","ox","NCO","OCN","N3"
        }
        STRONG = {
            # stærk-felt/π-acceptorer & chelater der ofte giver lav-spin (især for 3d)
            "NO2","ONO","CN","CO","PF3","PR3","PPh3","bpy","phen","terpy","en","dien"
        }
                
        # ---------- helpers ----------
        def parse_complex(txt):
            """
            Returnerer: metal, {ligand: antal}, total_ladning (int)
            Accepterer fx [FeCl6]^4-, [Co(CN)6]4-, [Fe(NH3)4(H2O)2]^2+, [Co(en)3]3+, ...
            """
            import re

            s = txt.strip().replace(" ", "")
            m = re.search(r"\[([A-Za-z][a-z]?)(.+?)\](?:\^?\s*([0-9]+)?\s*([+-]))?$", s)
            if not m:
                raise ValueError("Kan ikke parse. Brug fx [FeCl6]^4- eller [Co(CN)6]4-.")
            metal, inner, mag, sign = m.group(1), m.group(2), m.group(3), m.group(4)

            # total ladning (tillad også ...]4- uden '^')
            if mag and sign:
                q = int(mag) * (+1 if sign == "+" else -1)
            else:
                m2 = re.search(r"\](\d+)([+-])$", s)
                q = int(m2.group(1)) * (+1 if m2.group(2) == "+" else -1) if m2 else 0

            # Tokeniser ligander:
            #  1) (Lig)k         → plig, pn
            #  2) Ligk (uden ()) → lig,  n
            #  Lookahead stopper før næste Uppercase, '(' eller strengslut,
            #  så 'Cl6' → lig='Cl', n='6' og 'NO26' → lig='NO2', n='6'.
            pattern = re.compile(r"""
                \(
                    (?P<plig>[A-Za-z0-9]+)
                \)
                    (?P<pn>\d*)
                |
                (?P<lig>[A-Za-z][A-Za-z0-9]*?)
                    (?P<n>\d*)
                    (?=(?:[A-Z]|\(|$))
            """, re.VERBOSE)

            ligs = {}
            i = 0
            for m in pattern.finditer(inner):
                L = m.group("plig") or m.group("lig")
                n = m.group("pn") or m.group("n") or "1"

                # normaliser aliaser/typiske navne
                L = L.replace("ammine", "NH3").replace("aqua", "H2O")
                L = normalize_ligand(L)  # brug din helper

                # pæn case for halider mv.
                if L.upper() in {"F","CL","BR","I","OH","CN","NO2","ONO"}:
                    L = L.upper()
                    if L in {"F","CL","BR","I"}:
                        L = L.title()  # F, Cl, Br, I

                # gem
                ligs[L] = ligs.get(L, 0) + int(n)

                # (optionelt) detekter huller i match-kæden til debugging:
                i = m.end()

            # sikkerhed: alt skal være matchet
            if i != len(inner):
                # kan fx ske ved uventede tegn
                rest = inner[i:]
                raise ValueError(f"Kunne ikke parse rest: '{rest}' i '{inner}'")

            return metal, ligs, q


        def oxidation_state(metal, ligs, q):
            s = 0
            for L, n in ligs.items():
                charge = LIGAND_CHARGE.get(L, LIGAND_CHARGE.get(normalize_ligand(L)))
                if charge is None:
                    raise ValueError(f"Ukendt ligand/ladning: {L} (tilføj i tabellen).")
                s += n * charge
            return q - s


        def d_count(metal, ox):
            if metal not in GROUP_NUM:
                raise ValueError(f"Metal {metal} ikke i hurtig-tabellen.")
            return GROUP_NUM[metal] - ox

        def auto_spin_choice(ligand_dict: dict[str,int]) -> str:
            """Returner 'high', 'low' eller 'ambiguous' ud fra ligandtyperne."""
            s = w = b = 0
            for L, n in ligand_dict.items():
                L = normalize_ligand(L)
                if L in STRONG:     s += n
                elif L in WEAK:     w += n
                elif L in BORDERLINE: b += n
                else:                b += n  # ukendt → behandles som borderline
            # ren stærk eller ren svag → tydelig afgørelse
            if s > 0 and w == 0 and b == 0: return "low"
            if w > 0 and s == 0 and b == 0: return "high"
            # blandinger eller borderline → tvetydigt
            if s > 0 and w == 0 and b > 0:  return "low"        # stærk + borderline → typisk lav-spin
            if w > 0 and s == 0 and b > 0:  return "high"       # svag + borderline → typisk høj-spin
            if s > 0 and w > 0:             return "ambiguous"  # stærk + svag i samme kompleks
            return "ambiguous"

        def fill_oct_unpaired(d, mode):
            """Returnerer (n_unpaired, t2g_e, eg_e) for octahedral."""
            # orbitaler: 3 t2g + 2 eg
            occ = [0,0,0, 0,0]
            def place(seq, eleft):
                for i in seq:
                    if eleft<=0: break
                    if occ[i] < 2:
                        occ[i] += 1
                        eleft  -= 1
                return eleft
            e = d
            if mode=="low":
                # t2g: singly → pair, derefter eg: singly → pair
                e = place([0,1,2], e)
                e = place([0,1,2], e)
                e = place([3,4],   e)
                e = place([3,4],   e)
            else:  # "high"
                # t2g singly → eg singly → t2g pair → eg pair
                e = place([0,1,2], e)
                e = place([3,4],   e)
                e = place([0,1,2], e)
                e = place([3,4],   e)
            unp = sum(1 for x in occ if x==1)
            return unp, sum(occ[:3]), sum(occ[3:])

        # ---------- UI ----------
        win = tk.Toplevel(self)
        win.title("Unpaired d-electrons — octahedral complexes")
        win.transient(self); win.grab_set(); win.geometry("820x600")

        frm = ttk.Frame(win, padding=12); frm.pack(fill=tk.BOTH, expand=True)

        ttk.Label(frm, text="Complex (e.g. [FeCl6]^4-, [Co(CN)6]^4-)").grid(row=0,column=0,sticky="w")
        eComplex = ttk.Entry(frm, width=38); eComplex.grid(row=1,column=0,columnspan=3,sticky="we",padx=(0,10))
        eComplex.insert(0, "[FeCl6]^4-")

        ttk.Label(frm, text="Spin mode:").grid(row=0,column=3,sticky="e")
        spin = tk.StringVar(value="auto")
        spin_box = ttk.Combobox(frm, width=18, state="readonly",
                                values=["auto (by ligands)","force high-spin","force low-spin"])
        spin_box.grid(row=1,column=3,sticky="w")
        spin_box.current(0)

        # Resultater
        out = tk.StringVar(value="")
        ttk.Label(frm, textvariable=out, font=("Segoe UI", 11, "bold")).grid(row=2,column=0,columnspan=4,sticky="w",pady=(10,4))

        steps = tk.Text(frm, height=20, wrap="word"); steps.grid(row=3,column=0,columnspan=4,sticky="nsew")
        frm.grid_rowconfigure(3, weight=1)
        sv = ttk.Scrollbar(frm, orient="vertical", command=steps.yview); steps.configure(yscrollcommand=sv.set)
        sv.grid(row=3, column=4, sticky="ns")

        def say(lines):
            steps.configure(state="normal"); steps.delete("1.0","end")
            steps.insert("1.0", "\n".join(lines)); steps.configure(state="disabled")

        def compute():
            txt = eComplex.get().strip()
            try:
                metal, ligs, q = parse_complex(txt)
                ox = oxidation_state(metal, ligs, q)
                d  = d_count(metal, ox)
            except Exception as ex:
                messagebox.showerror("Parse/chem error", str(ex), parent=win); return

            if not (0 <= d <= 10):
                messagebox.showerror("Out of range", f"d-count = {d} (udenfor 0–10).", parent=win); return

            # vælg spin
            chosen = spin_box.get()
            if chosen.startswith("force high"): use = "high"
            elif chosen.startswith("force low"): use = "low"
            else:
                auto = auto_spin_choice(ligs)
                use = {"ambiguous":"high"}.get(auto, auto)  # default til high ved tvivl
            unp, t2, eg = fill_oct_unpaired(d, use)

            # hvis auto var tvivlsomt, vis begge
            lines = []
            lines.append(f"Parsed: metal = {metal}, ligands = {', '.join(f'{L}×{n}' for L,n in ligs.items())}, overall charge = {q:+d}")
            lines.append(f"Oxidation state: {metal} = {ox:+d}")
            lines.append(f"d-electron count: d^{d} (group({metal}) − ox)")
            auto = auto_spin_choice(ligs)
            lines.append(f"Spin assessment: auto = {auto}; using = {use}")
            lines.append(f"Octahedral filling → t₂g^{t2} e_g^{eg}")
            lines.append(f"Unpaired electrons = {unp}")

            # hvis ambivalent: giv også alternativt udfald
            if chosen.startswith("auto") and auto=="ambiguous":
                unp_hi, t2_hi, eg_hi = fill_oct_unpaired(d, "high")
                unp_lo, t2_lo, eg_lo = fill_oct_unpaired(d, "low")
                lines.append("")
                lines.append("Ambiguous ligands (fx NH3/H2O blandet med halider).")
                lines.append(f"  • If HIGH-spin:  t₂g^{t2_hi} e_g^{eg_hi}  →  {unp_hi} unpaired")
                lines.append(f"  • If LOW-spin:   t₂g^{t2_lo} e_g^{eg_lo}  →  {unp_lo} unpaired")

            out.set(f"Result: {unp} unpaired electron(s)")
            say(lines)
            return unp

        # Compute-bar i bunden (samme som dine andre modals)
        self._add_compute_bar(win, compute)



    



if __name__ == "__main__":
    app = ChemGUI()
    app.mainloop()
