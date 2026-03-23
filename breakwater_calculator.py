# ======================================================================================
# PROGRAM DESCRIPTION & METHODOLOGY
# ======================================================================================
#
# 1. PURPOSE:
#    This software performs preliminary hydraulic sizing of rubble-mound breakwater
#    armor made with artificial concrete units and the associated rock underlayers.
#    The implemented armor-unit families are:
#    - Simple Cubes
#    - Antifer Blocks
#
#    The calculator dimensions two structural zones:
#    - The Trunk (main longitudinal section)
#    - The Head (roundhead / exposed end section)
#
#    The implemented design philosophy is iso-geometric but not iso-weight:
#    the Head keeps the same nominal unit size and geometry as the Trunk, while
#    the required concrete specific weight is increased to satisfy the adopted
#    head-stability transfer rule.
#
# 2. METHODOLOGY:
#    The calculator uses empirical hydraulic-stability formulae and an internal
#    EN 13383 rock-grading database:
#
#    - Simple Cubes:
#      Van der Meer (1988a), with embedded coefficient sets for slopes 2.0:1
#      and 1.5:1.
#
#    - Antifer Blocks:
#      Chegini-Aghtouman (2006), with embedded coefficient sets for slopes
#      2.0:1 and 1.5:1.
#
#    - Underlayers:
#      Automatic selection from an internal EN 13383 grading table using a
#      target underlayer weight equal to one tenth of armor-unit weight.
#
# 3. DESIGN LOGIC & COMPUTATIONAL STRATEGY:
#
#    a. Input philosophy:
#       The exposed hydraulic storm input is the number of waves, Nz.
#       Storm duration is not entered directly; it is calculated internally as:
#
#           Storm_Duration_hr = Nz * Tm / 3600
#
#    b. Trunk stability calculation:
#       For the selected formula ID, the program evaluates the stability number:
#
#           Ns = (k1 * (Nod^k2 / Nz^k3) + k4) * s0m^(-k5)
#
#       where:
#       - s0m = Hs / L0
#       - L0  = g * Tm^2 / (2 * pi)
#
#       The trunk armor nominal diameter is then obtained from:
#
#           Dn = Hs / (Delta_trunk * Ns)
#
#       with:
#
#           Delta_trunk = (Wc_trunk / Ww) - 1
#
#       The trunk armor-unit weight is:
#
#           W_trunk = Wc_trunk * Dn^3
#
#       The equivalent Hudson stability coefficient for the trunk is also
#       calculated for reference.
#
#    c. Head design strategy:
#       The Head is not resized by increasing nominal diameter Dn.
#       Instead, the program keeps the same Dn and the same slope as the Trunk,
#       and transfers the design to the Head through a fixed Hudson-coefficient
#       ratio:
#
#           Kd_trunk / Kd_head = 1.5
#
#       Assuming fixed geometry, the required relative-density transfer becomes:
#
#           Delta_head = Delta_trunk * (1.5)^(1/3)
#
#       and the required head concrete specific weight is:
#
#           Wc_head = Ww * (Delta_head + 1)
#
#       Because unit volume is kept constant, head unit weight increases only
#       due to the higher concrete specific weight:
#
#           W_head = W_trunk * (Wc_head / Wc_trunk)
#
#       Therefore:
#       - same molds may be used at Trunk and Head;
#       - same nominal diameter is preserved;
#       - same geometric dimensions are preserved;
#       - Head units are heavier because density is increased.
#
#    d. Underlayer selection:
#       The target underlayer weight is taken as:
#
#           W_target = W_armor / 10
#
#       This target is converted to mass and matched against the internal
#       EN 13383 nominal grading limits. The adopted grading is the one that
#       strictly contains the target mass and, if more than one grading matches,
#       the tightest grading range is selected. If no grading strictly contains
#       the target mass, the code falls back to the first grading in the
#       database.
#
# 4. DEFAULTS, CONSTANTS, AND CURRENT CODE SETTINGS:
#
#    a. Default design inputs currently used in the code:
#       - Hs         = 11.0 m
#       - Tm         = 11.9 s
#       - Nz         = 3000
#       - Nod        = 0.5
#       - Wc_trunk   = 27.48 kN/m3
#       - Ww         = 10.05 kN/m3
#       - Formula_ID = 1
#
#       Default Formula_ID = 1 corresponds to:
#       Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)
#
#    b. Physical and layer constants currently used in the code:
#       - g            = 9.80665 m/s2
#       - W_rock_spec  = 26.5 kN/m3
#       - P_cubes      = 0.40
#       - P_rock       = 0.25
#       - Kd ratio     = 1.5
#
#    c. Formula sets currently embedded in the code:
#       1. Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)
#       2. Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)
#       3. Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)
#       4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)
#
# 5. OUTPUTS PRODUCED:
#    The calculator reports, for both Trunk and Head:
#    - Stability number
#    - Nominal diameter
#    - Armor-unit weight and mass
#    - Equivalent Kd value
#    - Armor-layer thickness
#    - Packing density
#    - Underlayer grading and derived rock parameters
#
#    For Cubes, displayed volume is reported as:
#
#        V = Dn^3
#
#    For Antifer units, displayed volume is reported using the geometric
#    shape factor:
#
#        V = 1.0247 * H^3
#
# 6. TECHNICAL BIBLIOGRAPHY & REFERENCES:
#
#    1. Van der Meer, J.W. (1988). "Rock Slopes and Gravel Beaches Under Wave Attack."
#       Doctoral Thesis, Delft University of Technology.
#
#    2. Van der Meer, J.W. (1988). "Stability of Cubes, Tetrapods and Accropode."
#       Proceedings of the Conference Breakwaters '88, Eastbourne, Thomas Telford.
#
#    3. Chegini, V., & Aghtouman, P. (2006). "An Investigation on the Stability of
#       Rubble Mound Breakwaters with Armour Layers of Antifer Cubes."
#
#    4. USACE (2006). "Coastal Engineering Manual (CEM)", Chapter VI-5.
#
#    5. CEN (2002). "EN 13383-1: Armourstone - Part 1: Specification."
#
# ======================================================================================

import math
import sys
import os

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
                "name": "Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)",
                "type": "Cubes", "slope_ratio": 2.0,
                "k1": 7.374304189198, "k2": 0.4, "k3": 0.3, "k4": 1.100642416298, "k5": 0.1
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
            # Coarse / Light Gradings
            {"name": "CP32/90",        "NLL_kg": 0.868,  "NUL_kg": 19.319},
            {"name": "CP45/125",       "NLL_kg": 2.415,  "NUL_kg": 51.758},
            {"name": "CP63/180",       "NLL_kg": 6.626,  "NUL_kg": 154.548},
            {"name": "CP90/250",       "NLL_kg": 19.319, "NUL_kg": 414.063},
            {"name": "CP45/180",       "NLL_kg": 2.415,  "NUL_kg": 154.548},
            {"name": "CP90/180",       "NLL_kg": 19.319, "NUL_kg": 154.548},
            # Light Mass Armourstone (LMA)
            {"name": "LMA5-40",        "NLL_kg": 5,      "NUL_kg": 40},
            {"name": "LMA10-60",       "NLL_kg": 10,     "NUL_kg": 60},
            {"name": "LMA15-120",      "NLL_kg": 15,     "NUL_kg": 120},
            {"name": "LMA40-200",      "NLL_kg": 40,     "NUL_kg": 200},
            {"name": "LMA60-300",      "NLL_kg": 60,     "NUL_kg": 300},
            {"name": "LMA15-300",      "NLL_kg": 15,     "NUL_kg": 300},
            # Heavy Mass Armourstone (HMA)
            {"name": "HMA300-1000",    "NLL_kg": 300,    "NUL_kg": 1000},
            {"name": "HMA1000-3000",   "NLL_kg": 1000,   "NUL_kg": 3000},
            {"name": "HMA3000-6000",   "NLL_kg": 3000,   "NUL_kg": 6000},
            {"name": "HMA6000-10000",  "NLL_kg": 6000,   "NUL_kg": 10000},
            {"name": "HMA10000-15000", "NLL_kg": 10000,  "NUL_kg": 15000},
        ]
		
        # ----------------------------------------------------------------------
        # DEFAULT INPUTS
        # ----------------------------------------------------------------------
        self.defaults = {
            "Hs": 11.0,                # Significant Wave Height (m)
            "Tm": 11.9,                # Mean Wave Period (s)
            "Nz": 3000.0,              # Number of waves in the design storm
            "Nod": 0.5,                # Damage Number (Nod=0 to 1 implies no damage/start of damage)
            "Wc": 27.48,               # Specific weight of concrete (kN/m3)
            "Ww": 10.05,               # Specific weight of seawater (kN/m3)
            "Formula_ID": 1            # Default is now VdM 2.0 (ID 1)
        }

    def calculate_L0(self, Tm):
        """Calculates Deepwater Wavelength L0 = g * Tm^2 / (2 * pi)"""
        return (self.g * Tm**2) / (2 * math.pi)

    def calculate_underlayer_params(self, W_armor):
        """
        Calculates underlayer parameters.
        Logic: Target rock is 1/10th the weight of armor. Finds the EN 13383 grading
        that strictly contains the target mass, breaking ties with the tightest range.
        """
        target_weight = W_armor / 10.0
        target_mass_kg = (target_weight * 1000.0) / self.g
        
        selected_grading = None
        min_range_width = float('inf')
        
        # --- CONTAINMENT & TIGHTEST RANGE LOGIC ---
        for grading in self.standard_gradings:
            # Check containment: Target must be strictly inside nominal limits
            if grading["NLL_kg"] < target_mass_kg < grading["NUL_kg"]:
                current_range = grading["NUL_kg"] - grading["NLL_kg"]
                
                # Update if this is the first match OR if this range is tighter (smaller)
                if current_range < min_range_width:
                    min_range_width = current_range
                    selected_grading = grading
                    
        # Fallback if no grading strictly contains the target mass
        if not selected_grading and self.standard_gradings:
            selected_grading = self.standard_gradings[0]

        final_NLL = selected_grading["NLL_kg"]
        final_NUL = selected_grading["NUL_kg"]
        final_M50 = 0.5 * (final_NLL + final_NUL)
        
        # --- EN 13383 LIMIT CALCULATIONS ---
        ELL = 0.7 * final_NLL
        EUL = 1.5 * final_NUL
        W_mean_kn = (final_M50 * self.g) / 1000.0
        
        Dn_rock = (W_mean_kn / self.W_rock_spec) ** (1/3.0)
        r2 = 2.0 * Dn_rock
        f2 = 100.0 * 2.0 * 1.0 * (1.0 - self.P_rock) / (Dn_rock**2)
        
        return {
            "grading_name": selected_grading["name"],
            "target_W": target_weight,
            "target_M50_kg": target_mass_kg,
            "M50_kg": final_M50,
            "NLL_kg": final_NLL,
            "NUL_kg": final_NUL,
            "ELL_kg": ELL,
            "EUL_kg": EUL,
            "W_mean_kn": W_mean_kn, 
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

        # Calculate storm duration from the number of waves
        Nz = params['Nz']
        storm_duration_hr = (Nz * params['Tm']) / 3600.0
        
        # Calculate Relative Density Delta = (Rho_concrete / Rho_water) - 1
        delta_trunk = (params['Wc'] / params['Ww']) - 1

        # 4. Algorithmic Core (Chegini-Aghtouman / Van der Meer) - TRUNK
        # This calculates the Stability Number (Ns)
        
        # Term 1: Damage / number-of-waves relationship
        term_damage = params['Nod'] ** k2
        term_waves = Nz ** k3
        damage_wave_ratio = term_damage / term_waves
        
        # Combine terms based on power-law formula
        scaled_term = k1 * damage_wave_ratio
        inv_f = scaled_term + k4

        # Wave steepness influence
        steepness_factor = s0m ** (-k5)
        Ns_trunk = inv_f * steepness_factor

        # 5. Block Sizing (Armor) - TRUNK
        # Main Stability Equation: Dn = Hs / (Delta * Ns)
        Dn = params['Hs'] / (delta_trunk * Ns_trunk)
        
        # Weight = Specific_Weight * Volume (Dn^3)
        W_trunk = params['Wc'] * (Dn ** 3)
        
        # --- PACKING DENSITY CALCULATION (TRUNK) ---
        # Formula: 100 * Layers(2) * Coeff(1) * (1 - Porosity) / Dn^2
        packing_density_trunk = 100.0 * 2.0 * 1.1 * (1.0 - self.P_cubes) / (Dn**2)

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
        packing_density_head = 100.0 * 2.0 * 1.1 * (1.0 - self.P_cubes) / (Dn**2)

        # 9. UNDERLAYER - HEAD
        # Independently calculated based on Head Armor Weight (W_head)
        # Since W_head differs from W_trunk (heavier), this may result in a different grading.
        ul_head = self.calculate_underlayer_params(W_head)

        # 10. Armor Layer Details (Common)
        # r1 is the theoretical thickness of the double armor layer
        r1 = 2.0 * 1.1 * Dn

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
                "Storm_Duration_hr": storm_duration_hr,
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
                "Dn": Dn,
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

        vol_trunk_display = ft["Dn"] ** 3 if c["type"] == "Cubes" else 1.0247 * (ft["dims"]["H"] ** 3)
        vol_head_display = fh["Dn"] ** 3 if c["type"] == "Cubes" else 1.0247 * (fh["dims"]["H"] ** 3)

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
        lines.append(f"   Number of waves (Nz)                : {p['Nz']:.0f}")
        lines.append(f"   Nod (Damage)                        : {p['Nod']:.2f}")
        lines.append(f"   Wc Trunk (Concrete Spec. Weight)    : {p['Wc']:.2f} kN/m3")
        lines.append(f"   Ww (Water Specific Weight)          : {p['Ww']:.2f} kN/m3")
        lines.append(f"   Relative Density D=(Wc/Ww)-1        : {i['delta']:.4f}")
        lines.append(f"   Structure Slope (TRUNK & HEAD)      : {c['slope_ratio']}:1")
        lines.append(f"   Porosity (Cubes)                    : {const['P_cubes']*100:.0f}%")
        lines.append(f"   Porosity (Rock Layer)               : {const['P_rock']*100:.0f}%")
        lines.append("-" * 80)
        
        # --- INTERMEDIATE PARAMETERS ---
        lines.append("2. INTERMEDIATE PARAMETERS")
        lines.append(f"   Wave Length (L0)                    : {i['L0']:.2f} m")
        lines.append(f"   wave number (k0 = 2*pi/L0)          : {i['k0']:.4f}")
        lines.append(f"   wave steepness (s0m = Hs/L0)        : {i['s0m']:.4f}")
        lines.append(f"   Calculated Storm Duration (h)       : {i['Storm_Duration_hr']:.3f} h")
        lines.append(f"   Stability Number TRUNK (Ns)         : {i['Ns_trunk']:.4f}")
        lines.append(f"   Stability Number HEAD (Ns)          : {fh['Ns_head']:.4f}")
        lines.append("-" * 80)
        
        # --- TRUNK SECTION ---
        lines.append("3. ARMOR LAYER RESULTS - TRUNK")
        lines.append(f"   BLOCK WEIGHT (W)                    : {ft['W']:.2f} kN")
        lines.append(f"   Mass (ton)                          : {ft['Mass_tonnes']:.2f} t")
        lines.append(f"   Nominal Diameter (Dn)               : {ft['Dn']:.3f} m")
        if c["type"] == "Cubes":
            lines.append(f"   Volume = Dn^3 (V)                   : {vol_trunk_display:.3f} m3")
        else:
            lines.append(f"   Volume = 1.0247 * H^3 (V)           : {vol_trunk_display:.3f} m3")
            lines.append(f"   Cube Height (H)                     : {ft['dims']['H']:.3f} m")
            lines.append(f"   Cube Top Width (B)                  : {ft['dims']['B']:.3f} m")
            lines.append(f"   Cube Base Width (A)                 : {ft['dims']['A']:.3f} m")
        lines.append(f"   KD_TRUNK (Equivalent)               : {ft['Kd_Equivalent']:.2f}")
        lines.append(f"   Double Layer Thickness (r1)         : {ft['r1']:.2f} m")
        lines.append(f"   Packing Density, d [units/100m2]    : {ft['packing_density']:.2f}")
        lines.append("")
        
        # --- UNDERLAYER TRUNK ---
        lines.append("4. UNDERLAYER RESULTS - TRUNK")
        lines.append(f"   Theoretical Target (W/10)           : {ut['target_W']:.2f} kN ({ut['target_M50_kg']:.1f} kg)")
        lines.append(f"   Adopted rock grading                : {ut['grading_name']}")
        lines.append(f"   Representative M50                  : {ut['M50_kg']:.1f} kg")
        lines.append(f"   Nominal lower limit (NLL)           : {ut['NLL_kg']:.1f} kg")
        lines.append(f"   Nominal upper limit (NUL)           : {ut['NUL_kg']:.1f} kg")
        lines.append(f"   Extreme lower limit (ELL)           : {ut['ELL_kg']:.1f} kg")
        lines.append(f"   Extreme upper limit (EUL)           : {ut['EUL_kg']:.1f} kg")
        lines.append(f"   Nominal Diameter (Dn_rock)          : {ut['Dn_rock']:.3f} m")
        lines.append(f"   Double Layer Thickness (r2)         : {ut['r2']:.2f} m")
        lines.append(f"   Packing Density, f2 [rocks/100m2]   : {ut['f2']:.2f}")
        lines.append("-" * 80)

        # --- HEAD SECTION ---
        lines.append("5. ARMOR LAYER RESULTS - HEAD (High Density)")
        lines.append("   *Maintains same Dn and Slope as Trunk*")
        lines.append(f"   Stability Ratio (Kd_T/Kd_H)         : {fh['Kd_Ratio']:.2f}")
        lines.append(f"   Nominal Diameter (Dn)               : {fh['Dn']:.3f} m")
        if c["type"] == "Cubes":
            lines.append(f"   Volume = Dn^3 (V)                   : {vol_head_display:.3f} m3")
        else:
            lines.append(f"   Volume = 1.0247 * H^3 (V)           : {vol_head_display:.3f} m3")
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
        lines.append(f"   Theoretical Target (W/10)           : {uh['target_W']:.2f} kN ({uh['target_M50_kg']:.1f} kg)")
        lines.append(f"   Adopted rock grading                : {uh['grading_name']}")
        lines.append(f"   Representative M50                  : {uh['M50_kg']:.1f} kg")
        lines.append(f"   Nominal lower limit (NLL)           : {uh['NLL_kg']:.1f} kg")
        lines.append(f"   Nominal upper limit (NUL)           : {uh['NUL_kg']:.1f} kg")
        lines.append(f"   Extreme lower limit (ELL)           : {uh['ELL_kg']:.1f} kg")
        lines.append(f"   Extreme upper limit (EUL)           : {uh['EUL_kg']:.1f} kg")
        lines.append(f"   Nominal Diameter (Dn_rock)          : {uh['Dn_rock']:.3f} m")
        lines.append(f"   Double Layer Thickness (r2)         : {uh['r2']:.2f} m")
        lines.append(f"   Packing Density, f2 [rocks/100m2]   : {uh['f2']:.2f}")
        lines.append("=" * 80)
        lines.append("")
        
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
    print("1. Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)")
    print("2. Van Der Meer (1988a) - Simple Cubes (Slope 1.5:1)")
    print("3. Chegini-Aghtouman (2006) - Antifer (Slope 2:1)")
    print("4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)")
    
    try:
        selection = input("\nOption [1-4]: ").strip()
    except KeyboardInterrupt:
        sys.exit(0)

    # Set formula ID, default to 1 if user just pressed ENTER or typed nonsense
    if selection in ['1', '2', '3', '4']:
        formula_id = int(selection)
    else:
        formula_id = 1

    # Unconditionally ask for parameters (matching Fortran/C++ behavior)
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
            "Nz": get_param("Number of waves (Nz)", defaults["Nz"]),
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