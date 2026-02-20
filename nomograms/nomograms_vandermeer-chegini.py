# --- Requirements ---
# numpy==2.3.1
# PyNomo==0.3.3
# pyx==0.15
# scipy==1.16.0

# pip install PyNomo==0.3.3
# pip install numpy pyx scipy --upgrade

# --------------------
# The 'arange' function was removed from scipy 1.8.0+.
# This patch adds it back by borrowing it from numpy.
import numpy
import scipy
try:
    scipy.arange
except AttributeError:
    scipy.arange = numpy.arange
# --------------------

import sys
import importlib.metadata

# --- PyNomo Version Check ---
try:
    pynomo_version = importlib.metadata.version("PyNomo")
    
    if pynomo_version != "0.3.3":
        print(f"Error: Incorrect PyNomo version detected ({pynomo_version}).")
        print("This script requires exactly version 0.3.3.")
        sys.exit(1)

except importlib.metadata.PackageNotFoundError:
    print("Error: PyNomo is not installed.")
    sys.exit(1)
# ----------------------------

import math
import pyx

# --- Configure PyX to use LaTeX engine Correctly ---
pyx.text.set(pyx.text.LatexEngine)
pyx.text.preamble(r"\usepackage{amsmath}")
pyx.text.preamble(r"\usepackage{amsfonts}")
pyx.text.preamble(r"\usepackage{amssymb}")
# ----------------------------------------------

sys.path.insert(0, "..")
from pynomo.nomographer import Nomographer

# --- Optional parameter to control isopleth plotting ---
PLOT_ISOPLETHS = True

# --- Constants from breakwater_calculator.py ---
Ww = 10.05    # Specific weight of seawater (kN/m3)
g = 9.80665   # Gravity

# --- Helper for Van Der Meer Slope Correction ---
# In breakwater_calculator.py, VdM Slope 2.0 uses a correction factor C:
# Ns = (k1*X + k4) * C * ...
# Therefore, to integrate C into the coefficients p1 and p4 for the nomogram:
# Ns = (k1*C)*X + (k4*C)
vdm_slope_correction = (2.0 / 1.5)**(1.0 / 3.0)

# --- Configuration for all nomograms to be generated ---
NOMOGRAM_CONFIGS = [
    # 1. Van Der Meer (1988) - Cubes (Slope 2.0:1)
    {
        "filename": "nomogram_vandermeer_cubes_slope_2.0.pdf",
        # target_W removed; will be calculated dynamically
        "params": {
            "p1": 6.7 * vdm_slope_correction,  # k1
            "p2": 0.4,                         # k2
            "p3": 0.3,                         # k3
            "p4": 1.0 * vdm_slope_correction,  # k4
            "p5": 0.1                          # k5
        },
        "title": (
            r"1. Van Der Meer (1988a) - Cubes (Slope 2.0:1)" + "\n\n" +
            r"$N_s = \frac{H_s}{\Delta D_n} = (k_1 N_{od}^{k_2}/N_z^{k_3} + k_4) s_{om}^{-k_5} (2.0/1.5)^{1/3}$" + "\n" +
            r"$k_1=7.374; k_2=0.400; k_3=0.300; k_4=1.101; k_5=0.100$"
        )
    },
    
    # 2. Van Der Meer (1988) - Cubes (Slope 1.5:1)
    {
        "filename": "nomogram_vandermeer_cubes_slope_1.5.pdf",
        "params": {
            "p1": 6.700, 
            "p2": 0.400, 
            "p3": 0.300, 
            "p4": 1.000, 
            "p5": 0.100
        },
        "title": (
            r"2. Van Der Meer (1988a) - Cubes (Slope 1.5:1)" + "\n\n" +
            r"$N_s = \frac{H_s}{\Delta D_n} = (k_1 N_{od}^{k_2}/N_z^{k_3} + k_4) s_{om}^{-k_5}$" + "\n" +
            r"$k_1=6.700; k_2=0.400; k_3=0.300; k_4=1.000; k_5=0.100$"
        )
    },

    # 3. Chegini-Aghtouman (2006) - Antifer (Slope 2:1)
    {
        "filename": "nomogram_chegini_antifer_slope_2.0.pdf",
        "params": {
            "p1": 6.138, 
            "p2": 0.443, 
            "p3": 0.276, 
            "p4": 1.164, 
            "p5": 0.070
        },
        "title": (
            r"3. Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)" + "\n\n" +
            r"$N_s = \frac{H_s}{\Delta D_n} = (k_1 N_{od}^{k_2}/N_z^{k_3} + k_4) s_{om}^{-k_5}$" + "\n" +
            r"$k_1=6.138; k_2=0.443; k_3=0.276; k_4=1.164; k_5=0.07$"
        )
    },

    # 4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)
    {
        "filename": "nomogram_chegini_antifer_slope_1.5.pdf",
        "params": {
            "p1": 6.951, 
            "p2": 0.443, 
            "p3": 0.291, 
            "p4": 1.082, 
            "p5": 0.082
        },
        "title": (
            r"4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)" + "\n\n" +
            r"$N_s = \frac{H_s}{\Delta D_n} = (k_1 N_{od}^{k_2}/N_z^{k_3} + k_4) s_{om}^{-k_5}$" + "\n" +
            r"$k_1=6.951; k_2=0.443; k_3=0.291; k_4=1.082; k_5=0.082$"
        )
    }
]

# --- Main loop to generate all nomograms ---
for config in NOMOGRAM_CONFIGS:
    print(f"Generating {config['filename']}...")

    # Unpack parameters for the current nomogram
    param_1 = config["params"]["p1"] # k1 (adjusted)
    param_2 = config["params"]["p2"] # k2
    param_3 = config["params"]["p3"] # k3
    param_4 = config["params"]["p4"] # k4 (adjusted)
    param_5 = config["params"]["p5"] # k5
    title = config["title"]
    filename = config["filename"]

    # --- Axis and Scale Definitions ---
    
    # Axis: f (Pivot 1)
    # UPDATED: Matches breakwater_calculator.py 'direct' logic: Ns_base = k1*X + k4
    # Therefore X = (u - k4) / k1.
    param_f1 = {
        'u_min': 1.2, 'u_max': 4.0,
        'function': lambda u: math.log10((u - param_4) / param_1),
        'scale_type': 'linear smart', 'title': '$f(N_s)$',
        'tick_levels': 3, 'tick_text_levels': 2, 'tick_side': 'left', 'tag': 'f',
    }
    
    # Axis: f (Pivot 2 - Mirror)
    # UPDATED: Weight W is proportional to Ns^-3 (W ~ 1/Ns^3).
    param_f2 = {
        'u_min': 1.2, 'u_max': 4.0,
        'function': lambda u: 3 * math.log10(u),
        'scale_type': 'linear smart', 'title': '',
        'tick_levels': 0, 'tick_text_levels': 0, 'tick_side': 'left', 'tag': 'f',
    }
    
    # Axis: Tm (Mean Wave Period) - Left Block
    # log(X) = k2 log(Nod) - k3 log(Nz)
    param_Tm1 = {
        'u_min': 5.0, 'u_max': 20.0,
        'function': lambda u: -param_3 * math.log10(u),
        'scale_type': 'linear smart', 'title': '$T_m$[s]',
        'tick_levels': 1, 'tick_text_levels': 1, 'tick_side': 'left',
    }
    
    # Axis: Duration (hours)
    param_dur = {
        'u_min': 3, 'u_max': 7*24,
        'function': lambda u: param_3 * math.log10(3600 * u),
        'scale_type': 'linear smart', 'title': 'dur[h]',
        'tick_levels': 2, 'tick_text_levels': 1, 'tick_side': 'left',
    }
    
    # Axis: Nod (Damage Number)
    param_Nod = {
        'u_min': 0.1, 'u_max': 2.0,
        'function': lambda u: -param_2 * math.log10(u),
        'scale_type': 'linear smart', 'title': '$N_{od}$[ ]',
        'tick_levels': 2, 'tick_text_levels': 1, 'tick_side': 'left',
    }
    
    # Axis: Hs (Significant Wave Height)
    # Matches W ~ Hs^(3 + 3k5)
    param_Hs = {
        'u_min': 3.0, 'u_max': 12.0,
        'function': lambda u: -(3 + 3 * param_5) * math.log10(u),
        'scale_type': 'linear', 'title': '$H_s$[m]',
        'tick_levels': 4, 'tick_text_levels': 2, 'tick_side': 'left',
    }
    
    # Axis: Wc (Concrete Density)
    # Matches W ~ Delta^-3
    param_Wc = {
        'u_min': 23.0, 'u_max': 32.0,
        'function': lambda u: -math.log10(u * (2 * math.pi / g)**(3 * param_5) / (u / Ww - 1)**3),
        'title': '$W_c$[kN/m3]', 'tick_levels': 2, 'tick_text_levels': 1,
        'tick_side': 'left', 'scale_type': 'linear',
    }
    
    # Axis: Tm (Mean Wave Period) - Right Block
    # Matches W ~ Tm^(-6k5)
    param_Tm2 = {
        'u_min': 5.0, 'u_max': 20.0,
        'function': lambda u: (6 * param_5) * math.log10(u),
        'scale_type': 'linear', 'title': '$T_m$[s]',
        'tick_levels': 1, 'tick_text_levels': 1, 'tick_side': 'left',
    }
    
    # Axis: W (Block Weight)
    param_W = {
        'u_min': 50.0, 'u_max': 1000.0,
        'function': lambda u: math.log10(u),
        'title': '$W$[kN]', 'tick_levels': 4, 'tick_text_levels': 2, 'scale_type': 'linear',
        'extra_params': [
            {'u_min': 50.0, 'u_max': 100.0, 'tick_levels': 4, 'tick_text_levels': 4},
            {'u_min': 100.0, 'u_max': 200.0, 'tick_levels': 4, 'tick_text_levels': 3}
        ]
    }

    # --- Block Definitions ---
    block1_params = {
        'block_type': 'type_3', 'width': 20.0, 'height': 30.0,
        'f_params': [param_Nod, param_dur, param_Tm1, param_f1],
    }
    
    block2_params = {
        'block_type': 'type_3', 'width': 25.0, 'height': 30.0, 'mirror_x': True,
        'f_params': [param_f2, param_Hs, param_Tm2, param_Wc, param_W], 
    }

    # --- Main Nomographer Parameters ---
    main_params = {
        'filename': filename,
        'paper_width': 20.0, 'paper_height': 30.0,
        'block_params': [block1_params, block2_params],
        'transformations': [('rotate', 0), ('scale paper',)],
        'title_str': title,
        'title_y': 28.0,
    }

    # --- Conditionally add isopleths ---
    if PLOT_ISOPLETHS:
        # Define Isopleth Input Values
        iso_Nod = 1.0
        iso_dur = 12.0 # hours
        iso_Tm = 13.0 # seconds
        iso_Hs = 10.0 # meters
        iso_Wc = 24.0 # kN/m3
        
        # --- MATH: Calculate Resulting W for the Isopleth ---
        # 1. Number of Waves (Nz)
        iso_Nz = (iso_dur * 3600) / iso_Tm
        
        # 2. Stability Base Term: k1 * Nod^k2 * Nz^-k3
        iso_base = param_1 * (iso_Nod**param_2) * (iso_Nz**-param_3)
        
        # 3. Wave Steepness (som)
        iso_som = (2 * math.pi * iso_Hs) / (g * iso_Tm**2)
        
        # 4. Stability Number (Ns) = (Base + k4) * som^-k5
        iso_Ns = (iso_base + param_4) * (iso_som**-param_5)
        
        # 5. Nominal Diameter (Dn)
        iso_Delta = (iso_Wc / Ww) - 1
        iso_Dn = iso_Hs / (iso_Ns * iso_Delta)
        
        # 6. Weight (W) = Wc * Dn^3
        iso_W_calculated = iso_Wc * (iso_Dn**3)

       # Set Isopleth Values
        # Block 1 calculates 'x' (f) based on these inputs
        block1_params['isopleth_values'] = [[iso_Nod, iso_dur, iso_Tm, 'x']]
        
        # Block 2 links 'x' to inputs and the CALCULATED W
        block2_params['isopleth_values'] = [['x', iso_Hs, iso_Tm, iso_Wc, iso_W_calculated]]
        
        main_params['isopleth_params'] = [
            {'color': 'Red', 'linewidth': 'THick', 'linestyle': 'solid', 'circle_size': 0.2},
        ]

    # --- Generate the Nomogram ---
    Nomographer(main_params)

print("\nAll nomograms generated successfully.")