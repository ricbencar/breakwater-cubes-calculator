# Hydraulic Stability Calculator for Breakwater Armor Units
## Concrete Cubes & Antifer Blocks

## 1. Abstract

This repository contains a comprehensive computational suite designed for the hydraulic design of coastal breakwater armor layers. The software implements advanced empirical formulae to dimension artificial concrete units (**Simple Cubes** and **Antifer Blocks**).

A distinguishing feature of this tool is its **Iso-Geometric Design Philosophy** for the breakwater head. Rather than increasing the nominal diameter ($D_n$) of armor units at the roundhead—which necessitates different casting moulds and complex logistics—this calculator solves for the required **increase in concrete density** ($\rho_c$). This allows the trunk and the head to be constructed using geometrically identical units, optimizing construction efficiency while meeting safety factors via material density adjustments.

## 2. Theoretical Framework & Methodology

The stability of rubble mound breakwaters is governed by the interaction between hydraulic loads (waves) and the structural resistance (gravity, friction, interlocking). The core parameter defining this balance is the **Stability Number ($N_s$)**:

```math
N_s = \frac{H_s}{\Delta D_n}
```

Where:
* $H_s$: Significant wave height at the toe of the structure (m).
* $\Delta$: Relative buoyant density, defined as $\frac{\rho_c}{\rho_w} - 1$.
* $D_n$: Nominal diameter (side length of the cube), defined as $V^{1/3}$.

### 2.1 Trunk Stability: Empirical Power Laws

For the trunk section (longitudinal main body), the software abandons the simplified Hudson formula in favor of more advanced power-law formulae that account for wave period, storm duration, and damage levels.

The general form of the stability equation implemented is:

```math
N_s = \left( k_1 \frac{N_{od}^{k_2}}{N_z^{k_3}} + k_4 \right) s_{om}^{-k_5}
```

Where:
* $N_{od}$: Damage Number (Number of displaced units relative to a strip of width $D_n$).
* $N_z$: Number of waves during the storm duration ($N_z = \text{Duration} / T_m$).
* $s_{om}$: Deep water wave steepness ($s_{om} = H_s / L_0$).
* $k_{1..5}$: Empirical coefficients derived from hydraulic model tests.

**Implemented Coefficients:**

* **Simple Cubes (Van der Meer, 1988):** Relies mostly on friction.
    * Coefficients: $k_1=6.7, k_2=0.4, k_3=0.3, k_4=1.0, k_5=0.1$.
* **Antifer Blocks (Chegini & Aghtouman, 2006):** Grooved cubes with higher interlocking.
    * Slope 2:1 Coefficients: $k_1=6.138, k_2=0.443, k_3=0.276, k_4=1.164, k_5=0.07$.

### 2.2 The Head Design: The Hudson Transfer Function

The breakwater roundhead is subjected to 3D wave attack, refraction, and diffraction, resulting in significantly higher instability than the trunk. Standard design (SPM/CEM) usually dictates calculating a larger $D_n$ for the head.

**The Methodology of this Tool:**

To maintain a single mould size for the entire project, we use the Hudson Formula as a Transfer Function.

**Step A: Calculate Equivalent Trunk Stability**
We first calculate the $N_s$ required for the Trunk using the advanced formulas in section 2.1. We then reverse-solve the Hudson formula to find the "Effective $K_D$" that the trunk is experiencing:

```math
K_{D,trunk} = \frac{(N_{s,trunk})^3}{\cot\alpha}
```

**Step B: Apply Safety Factor for Head**
The stability coefficient $K_D$ represents the degree of interlocking and friction. The head requires a lower $K_D$ (higher safety margin). We apply a safety ratio ($R \approx 1.5$):

```math
K_{D,head} = \frac{K_{D,trunk}}{R}
```

**Step C: Iso-Geometric Density Solving**
We fix the geometry ($D_n$) as constant ($D_{n,head} = D_{n,trunk}$). Instead, we solve for the required relative density at the head ($\Delta_{head}$):

```math
\Delta_{head} = \left( \frac{H_s}{D_n \cdot (K_{D,head} \cot\alpha)^{1/3}} \right)
```

**Step D: Required Concrete Density**
Finally, the required density for the head units is:

```math
\rho_{c,head} = \rho_w (\Delta_{head} + 1)
```

**Result:** The trunk might use normal concrete (2.4 t/m³), while the head might require high-density concrete (e.g., 2.65 t/m³) using heavy aggregates (barytes, magnetite, or iron ore), allowing the use of the same casting forms.

### 2.3 Underlayer Design (EN 13383)

The software automatically selects the appropriate underlayer rock based on the European Standard **EN 13383**.
* **Target Weight:** $W_{underlayer} \approx W_{armor} / 10$.
* **Algorithm:** The software iterates through a database of standard grading classes (e.g., HMA 300-1000kg, HMA 1-3t). It calculates the logarithmic mean of the grading limits and selects the class where the target weight is closest to this mean.

## 3. File Descriptions & Usage

The repository includes three variations of the calculator.

### 3.1 `breakwater-cubes-calculator.py`
A Python script for rapid, interactive analysis.
* **Dependencies:** Standard Python 3 libraries (`math`, `sys`, `os`).
* **Usage:** `python3 breakwater-cubes-calculator.py`
* **Features:** Interactive prompts, instant TTY output, easy to modify for scripting/batch processing.

### 3.2 `breakwater_calculator_cli.cpp`
A high-performance Command Line Interface written in C++. Ideal for integration into pipelines or running on headless servers.
* **Compilation (GCC/Clang):** `g++ -o bw_calc breakwater_calculator_cli.cpp` then `./bw_calc`
* **Input Method:** Interactive text prompts. Inputs are validated for type safety.

### 3.3 `breakwater_calculator_gui.cpp`
A native Windows GUI application built using the pure Win32 API. This provides a user-friendly form for data entry and result visualization.
* **Requirements:** Windows OS.
* **Compilation (MinGW / G++ on Windows):** `g++ -o BreakwaterTool.exe breakwater_calculator_gui.cpp -lgdi32 -luser32`
* **Compilation (Visual Studio Developer Command Prompt):** `cl.exe breakwater_calculator_gui.cpp user32.lib gdi32.lib`
* **Features:**
    * Dropdown selection for formulas.
    * Form-based input fields.
    * "Calculate" button triggers the logic engine.
    * Results displayed in a read-only text window.
    * Visual separation of Trunk and Head results.

## 4. Input Parameter Definitions

| Parameter | Symbol | Unit | Description |
| :--- | :--- | :--- | :--- |
| **Sig. Wave Height** | $H_s$ | meters | The average height of the highest 1/3 of waves at the toe of the structure. |
| **Mean Period** | $T_m$ | seconds | The average wave period. Used to calculate $L_0$ and $N_z$. |
| **Damage Number** | $N_{od}$ | - | Allowable damage. Start of damage $\approx 0.5$. Moderate damage $\approx 2.0$. |
| **Storm Duration** | $t$ | hours | Duration of the peak storm event. Typically 6 to 24 hours. |
| **Concrete Density** | $\rho_c$ | t/m³ | Density of the concrete used in the Trunk (usually 2.4). |
| **Slope** | $\cot\alpha$ | - | The cotangent of the breakwater slope angle (e.g., 1.5 or 2.0). |

## 5. Bibliographic References

This software implementation is grounded in the following academic and technical standards:

* **Van der Meer, J.W. (1988).** *Rock Slopes and Gravel Beaches Under Wave Attack.* Doctoral Thesis, Delft University of Technology.
    * Primary source for the permeability and spectral stability formulas for cubes.
* **Chegini, V., & Aghtouman, P. (2006).** "An Investigation on the Stability of Rubble Mound Breakwaters with Armour Layers of Antifer Cubes." *Journal of Marine Engineering*.
    * Source of the regression coefficients ($k_1$ to $k_5$) for grooved Antifer units.
* **CEN (European Committee for Standardization). (2002).** *EN 13383-1: Armourstone - Part 1: Specification.* Brussels.
    * Source for the standard Heavy Mass Armourstone (HMA) grading classes.
* **USACE (U.S. Army Corps of Engineers). (1984).** *Shore Protection Manual (SPM).* Waterways Experiment Station.
    * Source for the Hudson Formula logic used here specifically for the Head/Trunk transfer function.
* **USACE (2002).** *Coastal Engineering Manual (CEM).* EM 1110-2-1100.
    * General guidance on safety factors and wave steepness calculations.
* **Gómez-Martín, M.E., & Medina, J.R. (2014).** "Heterogeneous packing of concrete armor units." *Coastal Engineering*, 88, 27-38.
    * Reference for packing density coefficients ($k_\Delta$) used in layer thickness calculations.

## 6. License

This project is open-source. Engineering judgment must be exercised when using these results for construction. Physical model testing is recommended for final design verification.