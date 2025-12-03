# Hydraulic Stability Calculator for Breakwater Armor Units
### Concrete Cubes & Antifer Blocks

## 1. Abstract

[cite_start]This repository contains a comprehensive computational suite designed for the hydraulic design of coastal breakwater armor layers[cite: 4]. [cite_start]The software implements advanced empirical formulae to dimension artificial concrete units, specifically **Simple Cubes** (Van der Meer) and **Antifer Blocks** (Chegini & Aghtouman)[cite: 5].

[cite_start]A distinguishing feature of this tool is its **Iso-Geometric Design Philosophy** for the breakwater head[cite: 6]. [cite_start]Rather than increasing the nominal diameter ($D_n$) of armor units at the roundhead—which necessitates different casting moulds, storage logistics, and crane requirements—this calculator solves for the required **increase in concrete density** ($\rho_c$)[cite: 7]. [cite_start]This allows the Trunk and the Head to be constructed using geometrically identical units (same mould), optimizing construction efficiency while meeting strict safety factors via material density adjustments[cite: 8].

---

## 2. Theoretical Framework & Methodology

[cite_start]The software calculates design parameters based on a synthesis of hydraulic model test data and established coastal engineering standards[cite: 10].

### 2.1 Wave Mechanics & Geometric Parameters

[cite_start]Before determining armor stability, the software computes the fundamental wave and structure properties[cite: 11].

**Deep Water Wavelength ($L_0$):**
[cite_start]Calculated based on the mean wave period ($T_m$) using linear wave theory approximation for deep water[cite: 12, 13]:
$$\Large L_0 = \frac{g T_m^2}{2\pi}$$

**Wave Steepness ($s_{0m}$):**
[cite_start]A critical parameter influencing breaker type and stability[cite: 15]:
$$\Large s_{0m} = \frac{H_s}{L_0}$$

**Number of Waves ($N_z$):**
[cite_start]The total number of waves attacking the structure during the design storm duration ($t$ in hours)[cite: 17, 18]:
$$\Large N_z = \frac{t \times 3600}{T_m}$$

**Relative Buoyant Density ($\Delta$):**
[cite_start]The dimensionless density of the concrete relative to seawater[cite: 20]:
$$\Large \Delta = \frac{\rho_c}{\rho_w} - 1$$

### 2.2 Trunk Stability: Empirical Power Laws

The core of the calculator utilizes power-law formulas derived from extensive hydraulic model testing. [cite_start]The Stability Number ($N_s$) is the governing dimensionless parameter[cite: 22, 23].

**General Stability Formula:**
$$\Large N_s = \left( k_1 \frac{N_{od}^{k_2}}{N_z^{k_3}} + k_4 \right) s_{0m}^{-k_5}$$

*Where:*
* $N_{od}$ is the Damage Number (normalized damage level).
* $k_{1..5}$ are empirical coefficients specific to the block type and slope.



[cite_start]**Note on Calibration:** For Van der Meer Cubes (Slope 2.0:1), the script applies a specific calibration factor adjustment to the result found in the source code[cite: 26]:
$$\Large N_{s,adjusted} = N_s \times \left(\frac{2.0}{1.5}\right)^{1/3}$$

**Nominal Diameter ($D_n$):**
[cite_start]Once $N_s$ is determined, the required characteristic size of the block is[cite: 28]:
$$\Large D_n = \frac{H_s}{\Delta \cdot N_s}$$

**Armor Unit Weight ($W$):**
$$\Large W = \rho_c \cdot D_n^3$$

#### Empirical Coefficients ($k$)

[cite_start]The software utilizes the following database of coefficients [cite: 32-37]:

| Method | Block Type | Slope | $k_1$ | $k_2$ | $k_3$ | $k_4$ | $k_5$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Van der Meer (1988)** | Simple Cubes | 2.0:1 | 6.7 | 0.4 | 0.3 | 1.0 | 0.1 |
| **Van der Meer (1988)** | Simple Cubes | 1.5:1 | 6.7 | 0.4 | 0.3 | 1.0 | 0.1 |
| **Chegini (2006)** | Antifer | 2.0:1 | 6.138 | 0.443 | 0.276 | 1.164 | 0.07 |
| **Chegini (2006)** | Antifer | 1.5:1 | 6.951 | 0.443 | 0.291 | 1.082 | 0.082 |

### 2.3 Layer Properties & Dimensions

**Packing Density ($d$) [Units/100m²]:**
[cite_start]For a double layer ($n=2$) with packing coefficient $k=1.0$ and porosity $P=0.40$[cite: 39, 40]:
$$\Large d = \frac{100 \cdot n \cdot k \cdot (1 - P)}{D_n^2}$$

**Layer Thickness ($r_1$):**
$$\Large r_1 = n \cdot D_n$$

**Unit Dimensions (Approximate Shape Factors):**
[cite_start]For volume $V = W/\rho_c$[cite: 44, 45]:
$$\Large Height (H) = \left( \frac{V}{1.0247} \right)^{1/3}$$
$$\Large A = 1.086 \cdot H \quad ; \quad B = 1.005 \cdot H$$

### 2.4 The Head Design Strategy: Iso-Geometric Transfer

The breakwater head is subjected to severe 3D turbulence and wave breaking, significantly higher than the trunk. Standard practice suggests increasing the block weight. [cite_start]However, this tool uses a "Transfer Function" approach to maintain constant geometry ($D_{n,head} = D_{n,trunk}$)[cite: 48, 49].

**The Logic:**

1.  [cite_start]**Calculate Equivalent Hudson $K_D$:** We first determine what the Hudson Stability Coefficient ($K_D$) would be for the calculated Trunk armor[cite: 50].
    $$\Large K_{D,trunk} = \frac{\rho_c H_s^3}{W \Delta^3 \cot\alpha}$$

2.  [cite_start]**Apply Safety Ratio ($R$):** We enforce that the Head must be 1.5 times more stable than the Trunk (standard engineering judgment)[cite: 52, 53].
    $$\Large K_{D,head} = \frac{K_{D,trunk}}{1.5}$$

3.  **Solve for Density Increase:** In the Hudson formula, stability is proportional to $\Delta^3$. [cite_start]To maintain the same weight $W$ and size $D_n$ while increasing stability, we must increase $\Delta$[cite: 55].
    $$\Large \Delta_{head} = \Delta_{trunk} \cdot (1.5)^{1/3}$$

4.  [cite_start]**Result:** The software outputs a required **High Density Concrete** ($\rho_{c,head}$) for the head units[cite: 57].
    $$\Large \rho_{c,head} = \rho_w (\Delta_{head} + 1)$$

### 2.5 Underlayer Geotechnical Design (EN 13383)

[cite_start]The software automatically sizes the underlayer rock based on the **EN 13383** standard European rock gradings[cite: 59].

**Target Mean Weight:**
$$\Large W_{ul} \approx \frac{W_{armor}}{10}$$

**Rock Nominal Diameter:**
$$\Large D_{n,rock} = \left( \frac{W_{mean}}{\rho_{rock}} \right)^{1/3}$$

**Underlayer Packing Density ($f_2$) [Rocks/100m²]:**
[cite_start]For a double layer ($n=2$) and porosity $P_{rock}=0.25$[cite: 64, 65]:
$$\Large f_2 = \frac{100 \cdot n \cdot 1 \cdot (1 - P_{rock})}{D_{n,rock}^2}$$

---

## 3. Repository Contents & Compilation

[cite_start]This repository includes three implementation variations of the calculator to suit different engineering workflows[cite: 67].

### 3.1 Python Script (`breakwater-cubes-calculator.py`)
[cite_start]A rapid prototyping script ideal for quick checks and education[cite: 69, 70].
* [cite_start]**Dependencies:** Python 3.x (Standard libraries: `math`, `sys`, `os`)[cite: 71].
* **Execution:**
    ```bash
    python3 breakwater-cubes-calculator.py
    ```

### 3.2 C++ CLI (`breakwater_calculator_cli.cpp`)
A high-performance command-line tool designed for batch processing or integration into pipelines. [cite_start]It utilizes static linking for portability[cite: 73, 74].

* [cite_start]**Compilation (GCC/MinGW)[cite: 76]:**
    ```bash
    g++ -O3 -march=native -std=c++17 -Wall -Wextra \
    -static -static-libgcc -static-libstdc++ \
    -o breakwater_calculator_cli breakwater_calculator_cli.cpp
    ```
* [cite_start]**Usage (Interactive):** `./breakwater_calculator_cli` [cite: 78]
* **Usage (Arguments):**
    ```bash
    ./breakwater_calculator_cli [Hs] [Tm] [Duration] [Nod] [Wc] [FormulaID]
    # Example:
    ./breakwater_calculator_cli 10.5 13.0 12.0 1.0 24.0 1
    ```

### 3.3 C++ GUI (`breakwater_calculator_gui.cpp`)
A standalone native Windows application using the Win32 API. [cite_start]It provides a visual interface for inputting parameters and viewing generated reports[cite: 80, 81].

* [cite_start]**Compilation (MinGW on Windows)[cite: 83]:**
    ```bash
    g++ -O3 -march=native -std=c++17 -municode \
    breakwater_calculator_gui.cpp -o breakwater_calculator_gui \
    -mwindows -static -static-libgcc -static-libstdc++
    ```
* [cite_start]**Features[cite: 84]:**
    * Dropdown selection for formulas.
    * Form-based input fields.
    * Visual separation of Trunk and Head results.

---

## 4. Input Parameter Definitions

| Parameter | Symbol | Unit | Description |
| :--- | :---: | :---: | :--- |
| **Sig. Wave Height** | $H_s$ | meters | [cite_start]The average height of the highest 1/3 of waves at the toe of the structure. [cite: 89] |
| **Mean Period** | $T_m$ | seconds | [cite_start]The average wave period. [cite: 89] |
| **Damage Number** | $N_{od}$ | - | Allowable damage level. [cite_start]<br> $N_{od} \approx 0.5$: Start of damage. [cite: 89] |
| **Storm Duration** | $t$ | hours | [cite_start]Duration of the peak design storm event (determines $N_z$). [cite: 89] |
| **Concrete Density** | $\rho_c$ | kN/m³ | [cite_start]Specific weight of the concrete used in the Trunk. [cite: 89] |
| **Slope** | $\cot\alpha$ | - | [cite_start]The cotangent of the slope angle (e.g., 1.5 or 2.0). [cite: 89] |

---

## 5. Bibliographic References

1.  **Van der Meer, J.W. (1988)[cite_start].** *Rock Slopes and Gravel Beaches Under Wave Attack.* Doctoral Thesis. [cite: 91]
2.  **Van der Meer, J.W. (1988).** "Stability of Cubes, Tetrapods and Accropode." [cite_start]*Proceedings of the Conference Breakwaters '88*, Eastbourne, Thomas Telford. [cite: 5]
3.  **Chegini, V., & Aghtouman, P. (2006).** "An Investigation on the Stability of Rubble Mound Breakwaters with Armour Layers of Antifer Cubes." [cite_start]*Journal of Marine Engineering*. [cite: 92, 93]
4.  [cite_start]**CEN (2002).** *EN 13383-1: Armourstone - Part 1: Specification*. [cite: 94]
5.  [cite_start]**USACE (1984/2002).** *Shore Protection Manual* & *Coastal Engineering Manual*. [cite: 95]
6.  **Gómez-Martín, M.E., & Medina, J.R. (2014).** "Heterogeneous packing of concrete armor units." [cite_start]*Coastal Engineering*, 88, 27-38. [cite: 96]

---

## 6. License

This project is open-source. Engineering judgment must be exercised when using these results for construction. [cite_start]Physical model testing is recommended for final design verification [cite: 97-99].