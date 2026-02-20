import sys
import math
import numpy as np
import scipy

# --- FIX FOR SCIPY ATTRIBUTE ERROR ---
# Older libraries like pynomo expect scipy.arange, which was removed in new versions.
# This line restores it by pointing it to numpy.arange.
scipy.arange = np.arange
# ------------------------------------

sys.path.insert(0, "..")
from pynomo.nomographer import Nomographer

# --- Optional parameter to control isopleth plotting ---
PLOT_ISOPLETHS = True

# --- Axis Definitions ---
param_f1 = {
    'u_min': 0.2,
    'u_max': 4,
    'function': lambda u: math.log10(u),
    'title': 'f',
    'tick_levels': 3,
    'tick_text_levels': 2,
    'tick_side': 'left',
    'tag': 'f',
    'text_format': r'\Large{%g}',
}
param_f2 = {
    'u_min': 0.2,
    'u_max': 4,
    'function': lambda u: -math.log10(u),
    'title': '',
    'tick_levels': 0,
    'tick_text_levels': 0,
    'tick_side': 'left',
    'tag': 'f',
    'text_format': r'\Large{%g}',
}
param_Gtrunk = {
    'u_min': 23.0,
    'u_max': 32.0,
    'function': lambda u: math.log10(u/(u-10.05)**3),
    'title': 'Gtrunk[kN/m3]',
    'tick_levels': 2,
    'tick_text_levels': 1,
    'tick_side': 'left',
    'text_format': r'\Large{%g}',
}
param_Ghead = {
    'u_min': 23.0,
    'u_max': 32.0,
    'function': lambda u: -math.log10(u/(u-10.05)**3),
    'title': 'Ghead[kN/m3]',
    'tick_levels': 2,
    'tick_text_levels': 1,
    'tick_side': 'left',
    'text_format': r'\Large{%g}',
}
param_rKd = {
    'u_min': 0.1,
    'u_max': 5.0,
    'function': lambda u: -math.log10(u),
    'title': 'Kdtrunk/Kdhead',
    'tick_levels': 4,
    'tick_text_levels': 2,
    'tick_side': 'left',
    'text_format': r'\Large{%g}',
}
param_ralfa = {
    'u_min': 0.5,
    'u_max': 2.0,
    'function': lambda u: -math.log10(u),
    'title': 'Atrunk/Ahead',
    'scale_type':'linear smart',
    'tick_levels': 3,
    'tick_text_levels': 1,
    'tick_side': 'left',
    'text_format': r'\Large{%g}',
}
param_rW = {
    'u_min': 0.5,
    'u_max': 5.0,
    'function': lambda u: math.log10(u),
    'title': 'Whead/Wtrunk',
    'tick_levels': 3,
    'tick_text_levels': 2,
    'tick_side': 'left',
    'text_format': r'\Large{%g}',
}

# --- Block Definitions (Initial) ---
block1_params = {
    'block_type': 'type_1',
    'width': 20.0,
    'height': 30.0,
    'f1_params': param_Gtrunk,
    'f2_params': param_Ghead,
    'f3_params': param_f1,
    # isopleth_values added conditionally below
}

block2_params = {
    'block_type': 'type_3',
    'width': 28.0,
    'height': 29.0,
    'mirror_x': True,
    'f_params': [param_f2, param_rKd, param_ralfa, param_rW],
    # isopleth_values added conditionally below
}

# --- Main Configuration ---
main_params = {
    'filename': sys.argv[0].split('.')[0]+'.pdf',
    'paper_width': 21.0,
    'paper_height': 29.7,    
    'block_params': [block1_params, block2_params],
    'transformations': [('rotate', 0), ('scale paper',)],
    'title_str': 'HUDSON TRUNK TO HEAD NOMOGRAM',
    'title_y': 28.0,
    # isopleth_params added conditionally below
}

# --- Conditionally Add Isopleths ---
if PLOT_ISOPLETHS:
    # 1. Set Isopleth Values (Coordinates)
    block1_params['isopleth_values'] = [[24.0, 24.0, 'x']]
    block2_params['isopleth_values'] = [['x', 7.5/5.0, 1, 'x']]
    
    # 2. Set Isopleth Styling (Red, Thick, Solid)
    main_params['isopleth_params'] = [
        {'color': 'Red', 'linewidth': 'THick', 'linestyle': 'solid', 'circle_size': 0.2},
    ]

# --- Generate Nomogram ---
Nomographer(main_params)