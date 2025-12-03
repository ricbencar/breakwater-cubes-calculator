import math
import sys
import os

# ======================================================================================
# PROGRAM DESCRIPTION & METHODOLOGY
# ======================================================================================
#
# 1. PURPOSE:
#    This software calculates the required size and weight of artificial concrete 
#    armor units (specifically Cubes and Antifer blocks) for rubble mound breakwaters.
#    It determines requirements for two distinct sections of the structure:
#    - The Trunk (the main longitudinal section).
#    - The Head (the exposed tip/end of the breakwater).
#
# 2. METHODOLOGY:
#    The calculator utilizes empirical stability formulas derived from hydraulic 
#    model testing:
#    - For Antifer Cubes: Uses Chegini & Aghtouman (2006) power-law formulas.
#    - For Simple Cubes: Uses Van der Meer (1988) stability formulas.
#    - For Underlayers: maps calculated requirements to EN 13383 Standard Rock Grades.
#
# 3. DESIGN LOGIC & STRATEGY:
#    a. Trunk Design: 
#       Calculates the Stability Number (Ns) based on wave height, duration, and damage 
#       allowance. This determines the nominal diameter (Dn) and Weight (W).
#
#    b. Head Design (High-Density Strategy): 
#       The breakwater head is subjected to more severe turbulence than the trunk. 
#       Standard practice usually requires larger blocks. 
#       However, this specific script uses a "Constant Geometry, Variable Density" logic:
#       - It enforces a Safety Factor (Kd Ratio = 1.5).
#       - It keeps the Block Size (Dn) identical to the Trunk (simplifying construction).
#       - It increases the Concrete Density (Wc) to achieve the required stability 
#         against the higher forces at the head.
#
# 4. TECHNICAL BIBLIOGRAPHY & REFERENCES:
#
# 1. Van der Meer, J.W. (1988). "Stability of Cubes, Tetrapods and Accropode."
#    Proceedings of the Conference Breakwaters '88, Eastbourne, Thomas Telford.
#
# 2. Van der Meer, J.W. (1998). "Application and stability criteria for rock and 
#    artificial units." Chapter 11 in: Dikes and Revetments. Balkema.
#
# 3. Chegini, V. & Aghtouman, P. (2006). "An Investigation on the Stability of Rubble 
#    Mound Breakwaters with Armour Layers of Antifer Cubes." 
#    Journal of Marine Engineering.
#
# 4. USACE (2006). "Coastal Engineering Manual" (CEM), Chapter VI-5.
#
# 5. EN 13383-1. "Armourstone - Part 1: Specification." (Standard Rock Grades).
#
# ======================================================================================

class BreakwaterCalculator:
    """
    Coastal Engineering Calculator for armor block sizing.
    
    Structure:
    1. Trunk Armor (Van der Meer / Chegini-Aghtouman)
    2. Trunk Underlayer (Selected from EN 13383 Standard Grades)
    3. Head Armor (High Density, Fixed Kd Ratio 1.5)
    4. Head Underlayer (Selected from EN 13383 Standard Grades)
    """

    def __init__(self):
        # ----------------------------------------------------------------------
        # PHYSICAL CONSTANTS
        # ----------------------------------------------------------------------
        self.g = 9.80665  # Acceleration due to gravity (m/s^2)
        self.W_rock_spec = 26.5 # Standard Specific weight for underlayer rock (kN/m3) (Basalt/Granite)
        
        # ----------------------------------------------------------------------
        # LAYER CHARACTERISTICS
        # ----------------------------------------------------------------------
        # Porosity (P) is the percentage of void space in the armor layer.
        self.P_cubes = 0.40  # Porosity of armor layers (Standard for double-layer Cubes)
        self.P_rock = 0.25   # Porosity of rock underlayers (Standard approximation)
        
        # ----------------------------------------------------------------------
        # HEAD TO TRUNK TRANSFER RATIO
        # ----------------------------------------------------------------------
        # The head of a breakwater experiences higher hydraulic loads (3D effects).
        # We define a fixed ratio between the stability coefficient (Kd) of the Trunk vs Head.
        # Ratio = Kd_trunk / Kd_head = 1.5
        # This implies the Head needs to be ~1.5x more stable than the Trunk.
        self.KD_RATIO_FIXED = 1.5

        # ----------------------------------------------------------------------
        # FORMULA DATABASE
        # ----------------------------------------------------------------------
        # These coefficients (k1-k5) correspond to the empirical power-law formulas:
        # Ns = (k1 * (Nod^k2 / Nz^k3) + k4) * s0m^-k5
        
        self.formulas = {
            1: {
                "name": "Van Der Meer (1988a) - Cubes (Slope 2.0:1)",
                "type": "Cubes", "slope_ratio": 2.0, 
                "k1": 6.7, "k2": 0.4, "k3": 0.3, "k4": 1.0, "k5": 0.1
            },
            2: {
                "name": "Van Der Meer (1988a) - Cubes (Slope 1.5:1)",
                "type": "Cubes", "slope_ratio": 1.5,
                "k1": 6.7, "k2": 0.4, "k3": 0.3, "k4": 1.0, "k5": 0.1
            },
            3: {
                "name": "Chegini-Aghtouman (2006) - Antifer (Slope 2:1)",
                "type": "Antifer", "slope_ratio": 2.0,
                "k1": 6.138, "k2": 0.443, "k3": 0.276, "k4": 1.164, "k5": 0.07
            },
            4: {
                "name": "Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)",
                "type": "Antifer", "slope_ratio": 1.5,
                "k1": 6.951, "k2": 0.443, "k3": 0.291, "k4": 1.082, "k5": 0.082
            }
        }

        # ----------------------------------------------------------------------
        # ROCK GRADING DATABASE (EN 13383)
        # ----------------------------------------------------------------------
        # Mass ranges converted to weight (kN) dynamically during calculation.
        # These represent standard commercial availability of rock.
        self.standard_gradings = [
            # Light Mass (LMA)
            {"name": "EN 13383 LMA 5-40 kg",     "min_kg": 5,    "max_kg": 40},
            {"name": "EN 13383 LMA 10-60 kg",    "min_kg": 10,   "max_kg": 60},
            {"name": "EN 13383 LMA 40-200 kg",   "min_kg": 40,   "max_kg": 200},
            {"name": "EN 13383 LMA 60-300 kg",   "min_kg": 60,   "max_kg": 300},
            # Heavy Mass (HMA)
            {"name": "EN 13383 HMA 300-1000 kg", "min_kg": 300,  "max_kg": 1000},
            {"name": "EN 13383 HMA 1-3 ton",     "min_kg": 1000, "max_kg": 3000},
            {"name": "EN 13383 HMA 3-6 ton",     "min_kg": 3000, "max_kg": 6000},
            {"name": "EN 13383 HMA 6-10 ton",    "min_kg": 6000, "max_kg": 10000},
        ]

        # ----------------------------------------------------------------------
        # DEFAULT INPUTS
        # ----------------------------------------------------------------------
        self.defaults = {
            "Hs": 10.0,                # Significant Wave Height (m)
            "Tm": 13.0,                # Mean Wave Period (s)
            "Storm_Duration_hr": 12.0, # Duration of the design storm
            "Nod": 1.0,                # Damage Number (Nod=0 to 1 implies no damage/start of damage)
            "Wc": 24.0,                # Specific weight of concrete (kN/m3)
            "Ww": 10.05,               # Specific weight of seawater (kN/m3)
            "Formula_ID": 1            # Default is now VdM 2.0 (ID 1)
        }

    def calculate_L0(self, Tm):
        """Calculates Deepwater Wavelength L0 = g * Tm^2 / (2 * pi)"""
        return (self.g * Tm**2) / (2 * math.pi)

    def calculate_underlayer_params(self, W_armor):
        """
        Helper to calculate underlayer parameters.
        Logic: The underlayer rock is typically 1/10th the weight of the armor unit.
        This function calculates that target, then finds the closest EN 13383 
        standard rock grade.
        """
        target_weight = W_armor / 10.0
        
        selected_grading = None
        min_diff = float('inf')
        
        # Iterate through standard gradings to find the closest match
        for grading in self.standard_gradings:
            # Convert Mass (kg) to Weight (kN)
            # W = m * g / 1000
            w_min = grading["min_kg"] * self.g / 1000.0
            w_max = grading["max_kg"] * self.g / 1000.0
            
            # Calculate Mean Weight of the range
            w_mean_range = (w_min + w_max) / 2.0
            
            # Find closest to target (W_armor/10)
            diff = abs(w_mean_range - target_weight)
            
            if diff < min_diff:
                min_diff = diff
                selected_grading = {
                    "grading_name": grading["name"],
                    "W_mean": w_mean_range, # Use mean of the range for thickness calculations
                    "W1": w_min,
                    "W2": w_max,
                    "target_W": target_weight
                }
        
        # Use the selected grading's mean weight for hydraulic calculations
        W_mean = selected_grading["W_mean"]
        
        # Calculate nominal diameter of the rock: Dn = (Weight / SpecificWeight)^(1/3)
        Dn_rock = (W_mean / self.W_rock_spec) ** (1/3.0)
        
        # Thickness r2 (Double layer rock)
        r2 = 2.0 * Dn_rock
        
        # Packing Density f2 (Rocks per 100m2)
        # Formula: 100 * n_layers * k_layer * (1 - Porosity) / Dn^2
        # Here n=2, k=1.0 approx.
        f2 = 100.0 * 2.0 * 1.0 * (1.0 - self.P_rock) / (Dn_rock**2)
        
        # Return merged dictionary containing all hydraulic and physical properties
        return {
            "grading_name": selected_grading["grading_name"],
            "target_W": selected_grading["target_W"],
            "W_mean": W_mean, 
            "W1": selected_grading["W1"], 
            "W2": selected_grading["W2"],
            "Dn_rock": Dn_rock, 
            "r2": r2, 
            "f2": f2,
            "W_rock_spec": self.W_rock_spec
        }

    # --------------------------------------------------------------------------
    # CORE CALCULATION METHOD
    # --------------------------------------------------------------------------
    def solve(self, formula_id, user_inputs=None):
        # 1. Load Parameters (Merge defaults with any user overrides)
        params = self.defaults.copy()
        if user_inputs:
            params.update(user_inputs)
        
        # 2. Load Coefficients from the Formula Database
        if formula_id not in self.formulas:
            raise ValueError(f"Invalid Formula ID: {formula_id}")
        
        coeffs = self.formulas[formula_id]
        k1, k2, k3, k4, k5 = coeffs['k1'], coeffs['k2'], coeffs['k3'], coeffs['k4'], coeffs['k5']

        # 3. Preliminary Hydraulic Calculations
        # Calculate deep water wavelength and wave steepness (s0m)
        L0 = self.calculate_L0(params['Tm'])
        k0 = (2 * math.pi) / L0
        s0m = params['Hs'] / L0
        sm = s0m 

        # Calculate Number of Waves (Nz) during the storm
        Nz = (params['Storm_Duration_hr'] * 3600) / params['Tm']
        
        # Calculate Relative Density Delta = (Rho_concrete / Rho_water) - 1
        delta_trunk = (params['Wc'] / params['Ww']) - 1

        # 4. Algorithmic Core (Chegini-Aghtouman / Van der Meer) - TRUNK
        # This calculates the Stability Number (Ns)
        
        # Term 1: Damage / Waves relationship
        term_damage = params['Nod'] ** k2
        term_waves = Nz ** k3
        damage_wave_ratio = term_damage / term_waves
        
        # Combine terms based on power-law formula
        scaled_term = k1 * damage_wave_ratio
        inv_f = scaled_term + k4
        f_stab = 1.0 / inv_f 
        
        # Wave steepness influence
        steepness_factor = s0m ** (-k5)
        
        # The logic below adjusts Ns for VdM Slope 2.0 specifically if needed.
        # This accounts for slight geometric variations in the original VdM 1988 dataset.
        if formula_id == 1:
            Ns_trunk = inv_f * steepness_factor * (2.0/1.5)**(1/3)
        else:
            Ns_trunk = inv_f * steepness_factor

        # 5. Block Sizing (Armor) - TRUNK
        # Main Stability Equation: Dn = Hs / (Delta * Ns)
        Dn = params['Hs'] / (delta_trunk * Ns_trunk)
        
        # Weight = Specific_Weight * Volume (Dn^3)
        W_trunk = params['Wc'] * (Dn ** 3)
        
        # --- PACKING DENSITY CALCULATION (TRUNK) ---
        # Formula: 100 * Layers(2) * Coeff(1) * (1 - Porosity) / Dn^2
        packing_density_trunk = 100.0 * 2.0 * 1.0 * (1.0 - self.P_cubes) / (Dn**2)

        # 6. Hudson Comparative Calculation (Kd Trunk)
        # Calculates the equivalent Hudson Stability Coefficient (Kd) for reference.
        slope = coeffs['slope_ratio']
        kd_trunk_equiv = (params['Wc'] * (params['Hs']**3)) / (W_trunk * (delta_trunk**3) * slope)

        # 7. UNDERLAYER - TRUNK
        # Independently calculated based on Trunk Armor Weight
        ul_trunk = self.calculate_underlayer_params(W_trunk)

        # 8. HEAD CALCULATION (FIXED RATIO 1.5)
        # -------------------------------------------------------------
        # STRATEGY: ISO-GEOMETRIC HIGH DENSITY HEAD
        # Instead of increasing the block size (Dn) for the head, we keep
        # the size constant to match the Trunk cube size.
        # To achieve stability, we increase the density of the concrete.
        # -------------------------------------------------------------
        
        # Logic: Kd_trunk / Kd_head = 1.5
        kd_ratio = self.KD_RATIO_FIXED
        kd_head_derived = kd_trunk_equiv / kd_ratio
        
        # Calculate Delta Required for Head (to keep Dn same as trunk)
        # From Hudson: Kd ~ Delta^3. Therefore Delta_head = Delta_trunk * (Ratio)^(1/3)
        delta_head = delta_trunk * (kd_ratio**(1/3.0))
        
        # Convert required Delta back to Concrete Specific Weight (Wc)
        Wc_head = params['Ww'] * (delta_head + 1)
        
        # Head Weight: W_head = W_trunk * (Wc_head / Wc_trunk)
        # Note: W_head is usually > W_trunk due to higher density required
        W_head = W_trunk * (Wc_head / params['Wc'])

        # Calculate Stability Number for Head
        Ns_head = params['Hs'] / (delta_head * Dn)
        
        # --- PACKING DENSITY CALCULATION (HEAD) ---
        # Since Dn is maintained constant between Head and Trunk, 
        # the packing density (units/100m2) will actually be identical.
        # We perform the calculation again for verification.
        packing_density_head = 100.0 * 2.0 * 1.0 * (1.0 - self.P_cubes) / (Dn**2)

        # 9. UNDERLAYER - HEAD
        # Independently calculated based on Head Armor Weight (W_head)
        # Since W_head differs from W_trunk (heavier), this may result in a different grading.
        ul_head = self.calculate_underlayer_params(W_head)

        # 10. Armor Layer Details (Common)
        # r1 is the theoretical thickness of the double armor layer
        r1 = 2.0 * Dn

        # --- CALCULATE CUBE DIMENSIONS (H, A, B) ---
        # These shape factors are specific to the unit type (Cubes vs Antifer).
        # Note: 1.0247 is a shape factor approximation used here for volume.
        
        # Trunk Dimensions
        vol_trunk = W_trunk / params['Wc']
        h_trunk = (vol_trunk / 1.0247)**(1/3.0)
        a_trunk = 1.086 * h_trunk
        b_trunk = 1.005 * h_trunk

        # Head Dimensions
        vol_head = W_head / Wc_head
        h_head = (vol_head / 1.0247)**(1/3.0)
        a_head = 1.086 * h_head
        b_head = 1.005 * h_head

        # 11. Compile Results into a Dictionary
        results = {
            "inputs": params,
            "coefficients": coeffs,
            "intermediate": {
                "L0": L0, 
                "k0": k0,    # Added to results
                "s0m": s0m,  # Added to results
                "sm": sm, 
                "Nz": Nz, 
                "delta": delta_trunk, 
                "Ns_trunk": Ns_trunk 
            },
            "final_trunk": {
                "Dn": Dn,
                "W": W_trunk,
                "Mass_tonnes": W_trunk / self.g,
                "Kd_Equivalent": kd_trunk_equiv,
                "r1": r1,
                "packing_density": packing_density_trunk,
                "dims": {
                    "H": h_trunk,
                    "A": a_trunk,
                    "B": b_trunk
                }
            },
            "underlayer_trunk": ul_trunk,
            "final_head": {
                "Ns_head": Ns_head, 
                "Kd_Derived": kd_head_derived,
                "Kd_Ratio": kd_ratio,
                "Delta_Required": delta_head,
                "Wc_Required": Wc_head,
                "W": W_head,
                "Mass_tonnes": W_head / self.g,
                "packing_density": packing_density_head,
                "dims": {
                    "H": h_head,
                    "A": a_head,
                    "B": b_head
                }
            },
            "underlayer_head": ul_head,
            "constants": {
                "P_rock": self.P_rock,
                "P_cubes": self.P_cubes
            }
        }
        return results

    # --------------------------------------------------------------------------
    # REPORT GENERATION
    # --------------------------------------------------------------------------
    def generate_report_file(self, results, filepath="output.txt"):
        # Unpack dictionary for easier string formatting
        p = results["inputs"]
        c = results["coefficients"]
        i = results["intermediate"]
        ft = results["final_trunk"]
        ut = results["underlayer_trunk"]
        fh = results["final_head"]
        uh = results["underlayer_head"]
        const = results["constants"]

        lines = []
        lines.append("================================================================================")
        lines.append("    TECHNICAL REPORT: BREAKWATER ARMOR & UNDERLAYER DESIGN                        ")
        lines.append("================================================================================")
        lines.append(f"Methodology: {c['name']}")
        lines.append("-" * 80)
        
        # --- INPUT PARAMETERS ---
        lines.append("1. INPUT PARAMETERS")
        lines.append(f"   Hs (Sigificant Wave Height)         : {p['Hs']:.2f} m")
        lines.append(f"   Tm (Mean Wave Period)               : {p['Tm']:.2f} s")
        lines.append(f"   Storm Duration (h)                  : {p['Storm_Duration_hr']:.2f} h")
        lines.append(f"   Nod (Damage)                        : {p['Nod']:.2f}")
        lines.append(f"   Wc Trunk (Concrete Spec. Weight)    : {p['Wc']:.2f} kN/m3")
        lines.append(f"   Ww (Water Specific Weight)          : {p['Ww']:.2f} kN/m3")
        lines.append(f"   Relative Density D=(Wc/Ww)-1        : {i['delta']:.5f}")
        lines.append(f"   Structure Slope (TRUNK & HEAD)      : {c['slope_ratio']}:1")
        lines.append(f"   Porosity (Cubes)                    : {const['P_cubes']*100:.0f}%")
        lines.append(f"   Porosity (Rock Layer)               : {const['P_rock']*100:.0f}%")
        lines.append("-" * 80)
        
        # --- INTERMEDIATE PARAMETERS ---
        lines.append("2. INTERMEDIATE PARAMETERS")
        lines.append(f"   Wave Length (L0)                    : {i['L0']:.2f} m")
        lines.append(f"   wave number (k0 = 2*pi/L0)          : {i['k0']:.5f}")
        lines.append(f"   wave steepness (s0m = Hs/L0)        : {i['s0m']:.5f}")
        lines.append(f"   Number of waves (Nz)                : {i['Nz']:.0f}")
        lines.append(f"   Stability Number TRUNK (Ns)         : {i['Ns_trunk']:.5f}")
        lines.append(f"   Stability Number HEAD (Ns)          : {fh['Ns_head']:.5f}")
        lines.append("-" * 80)
        
        # --- TRUNK SECTION ---
        lines.append("3. ARMOR LAYER RESULTS - TRUNK")
        lines.append(f"   BLOCK WEIGHT (W)                    : {ft['W']:.2f} kN")
        lines.append(f"   Mass (ton)                          : {ft['Mass_tonnes']:.2f} t")
        lines.append(f"   Nominal Diameter (Dn)               : {ft['Dn']:.3f} m")
        lines.append(f"   Cube Height (H)                     : {ft['dims']['H']:.3f} m")
        lines.append(f"   Cube Top Width (B)                  : {ft['dims']['B']:.3f} m")
        lines.append(f"   Cube Base Width (A)                 : {ft['dims']['A']:.3f} m")
        lines.append(f"   KD_TRUNK (Equivalent)               : {ft['Kd_Equivalent']:.2f}")
        lines.append(f"   Double Layer Thickness (r1)         : {ft['r1']:.2f} m")
        lines.append(f"   Packing Density, d [units/100m2]    : {ft['packing_density']:.2f}")
        lines.append("")
        
        # --- UNDERLAYER TRUNK ---
        lines.append("4. UNDERLAYER RESULTS - TRUNK")
        lines.append(f"   Theoretical Target (W/10)           : {ut['target_W']:.2f} kN")
        lines.append(f"   Adopted rock grading                : {ut['grading_name']}")
        lines.append(f"   Grading Min (Lower Limit)           : {ut['W1']:.2f} kN")
        lines.append(f"   Grading Max (Upper Limit)           : {ut['W2']:.2f} kN")
        lines.append(f"   Mean Weight (Used for thickness)    : {ut['W_mean']:.2f} kN")
        lines.append(f"   Nominal Diameter (Dn_rock)          : {ut['Dn_rock']:.3f} m")
        lines.append(f"   Double Layer Thickness (r2)         : {ut['r2']:.2f} m")
        lines.append(f"   Packing Density, f2 [rocks/100m2]   : {ut['f2']:.2f}")
        lines.append("-" * 80)

        # --- HEAD SECTION ---
        lines.append("5. ARMOR LAYER RESULTS - HEAD (High Density)")
        lines.append("   *Maintains same Dn and Slope as Trunk*")
        lines.append(f"   Stability Ratio (Kd_T/Kd_H)         : {fh['Kd_Ratio']:.2f}")
        lines.append(f"   Nominal Diameter (Dn)               : {ft['Dn']:.3f} m")
        lines.append(f"   Cube Height (H)                     : {fh['dims']['H']:.3f} m")
        lines.append(f"   Cube Top width (B)                  : {fh['dims']['B']:.3f} m")
        lines.append(f"   Cube Base Width (A)                 : {fh['dims']['A']:.3f} m")
        lines.append(f"   KD_HEAD (Equivalent)                : {fh['Kd_Derived']:.2f}")
        lines.append(f"   Required Concrete Density (Wc)      : {fh['Wc_Required']:.2f} kN/m3")
        lines.append(f"   BLOCK WEIGHT (W)                    : {fh['W']:.2f} kN")
        lines.append(f"   Mass (ton)                          : {fh['Mass_tonnes']:.2f} t")
        lines.append(f"   Packing Density, d [units/100m2]    : {fh['packing_density']:.2f}")
        lines.append("")
        
        # --- UNDERLAYER HEAD ---
        lines.append("6. UNDERLAYER RESULTS - HEAD")
        lines.append(f"   Theoretical Target (W/10)           : {uh['target_W']:.2f} kN")
        lines.append(f"   Adopted rock grading                : {uh['grading_name']}")
        lines.append(f"   Grading Min (Lower Limit)           : {uh['W1']:.2f} kN")
        lines.append(f"   Grading Max (Upper Limit)           : {uh['W2']:.2f} kN")
        lines.append(f"   Mean Weight (Used for thickness)    : {uh['W_mean']:.2f} kN")
        lines.append(f"   Nominal Diameter (Dn_rock)          : {uh['Dn_rock']:.3f} m")
        lines.append(f"   Double Layer Thickness (r2)         : {uh['r2']:.2f} m")
        lines.append(f"   Packing Density, f2 [rocks/100m2]   : {uh['f2']:.2f}")
        lines.append("================================================================================")
        
        try:
            with open(filepath, "w", encoding="utf-8") as file:
                file.write("\n".join(lines))
            print(f"\n Report generated successfully: {os.path.abspath(filepath)}")
            print("Report content:\n")
            print("\n".join(lines))
        except IOError as e:
            print(f"\n Error saving file: {e}")

# ==============================================================================
# MAIN EXECUTION BLOCK
# ==============================================================================
def main():
    calc = BreakwaterCalculator()
    
    print("\n--- COASTAL PROTECTION BLOCK CALCULATOR (TRUNK & HEAD) ---")
    # REORDERED MENU
    print("1. Van Der Meer (1988) - Cubes (Slope 2.0:1)")
    print("2. Van Der Meer (1988) - Cubes (Slope 1.5:1)")
    print("3. Chegini-Aghtouman (2006) - Antifer (Slope 2:1)")
    print("4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)")
    
    try:
        selection = input("\nOption [1-4]: ").strip()
    except KeyboardInterrupt:
        sys.exit(0)

    formula_id = 1  # Default set to 1 (Now VdM 2.0)
    user_inputs = None

    if selection in ['1', '2', '3', '4']:
        formula_id = int(selection)
        print("\n--- Enter Parameters (Press ENTER for Default) ---")
        
        def get_param(prompt, default_val):
            val = input(f"{prompt} [{default_val}]: ").strip()
            return float(val) if val else default_val

        try:
            defaults = calc.defaults
            user_inputs = {
                "Hs": get_param("Hs (m)", defaults["Hs"]),
                "Tm": get_param("Tm (s)", defaults["Tm"]),
                "Nod": get_param("Nod (Damage)", defaults["Nod"]),
                "Storm_Duration_hr": get_param("Duration (h)", defaults["Storm_Duration_hr"]),
                "Wc": get_param("Concrete Weight Trunk (kN/m3)", defaults["Wc"])
            }
        except ValueError:
            print("\n Error: Non-numeric value.")
            sys.exit(1)

    try:
        results = calc.solve(formula_id, user_inputs)
        calc.generate_report_file(results, "output.txt")
    except Exception as e:
        print(f"\n Calculation Error: {e}")

if __name__ == "__main__":
    main()