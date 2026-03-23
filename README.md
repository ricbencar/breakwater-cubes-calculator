# Hydraulic Stability Calculator for Breakwater Armor Units
## Simple Cubes and Antifer Blocks

## 1. Purpose

This repository contains a family of calculators for preliminary hydraulic
sizing of rubble-mound breakwater armor made with **simple concrete cubes**
and **Antifer blocks**, together with the associated **rock underlayers**.

The current scripts implement one common engineering workflow:

- size the **trunk armor** with empirical stability formulae;
- derive an equivalent **Hudson stability coefficient** for the trunk;
- transfer that trunk design to the **head** through a fixed stability ratio;
- keep the **same nominal armor size** at trunk and head;
- increase the **required concrete specific weight** at the head instead of
  increasing block size;
- select the **underlayer grading** from an internal **EN 13383** database.

In other words, the implemented design philosophy is **iso-geometric** but
**not iso-weight**:

- same nominal diameter at trunk and head;
- same armor-unit geometry at trunk and head;
- same packing-density expression at trunk and head;
- higher required concrete specific weight at the head;
- higher head unit weight because density increases while volume is kept.

---

## 2. Repository contents

The scripts covered by this README are:

- `breakwater_calculator.py` — Python command-line implementation;
- `breakwater_calculator_cli.cpp` — C++ command-line implementation;
- `breakwater_calculator_gui.cpp` — Win32 C++ graphical interface;
- `breakwater_calculator_gui_pt.cpp` — Portuguese Win32 C++ graphical
  interface;
- `breakwater_calculator_cli.f90` — Fortran command-line implementation.

### 2.1 Reference logic

The intended reference implementation is **`breakwater_calculator_gui.cpp`**.
The currently aligned version uses:

- **number of waves** `Nz` as the exposed storm input;
- **storm duration** computed internally from `Nz` and `Tm`;
- embedded precise coefficients for the **2.0:1 Van der Meer simple-cube**
  case;
- a fixed **head-to-trunk Hudson ratio** of `1.5`;
- the same **EN 13383 grading database** in all aligned scripts.

The Portuguese GUI is intended to be a language-only counterpart of the
English GUI, with the same numerical results and only translated visible text.

---

## 3. Inputs, defaults, and constants

### 3.1 Exposed inputs

| Symbol | Meaning | Unit |
| --- | --- | --- |
| `Hs` | Significant wave height | m |
| `Tm` | Mean wave period | s |
| `Nz` | Number of waves in design storm | – |
| `Nod` | Damage number | – |
| `Wc` | Trunk concrete specific weight | kN/m³ |
| `Formula_ID` | Formula selector | – |

### 3.2 Internal water specific weight

| Symbol | Meaning | Unit | Default |
| --- | --- | --- | ---: |
| `Ww` | Water specific weight | kN/m³ | 10.05 |

### 3.3 Defaults currently used by the aligned scripts

| Parameter | Default |
| --- | ---: |
| `Hs` | 11.0 m |
| `Tm` | 11.9 s |
| `Nz` | 3000 |
| `Nod` | 0.5 |
| `Wc` | 27.48 kN/m³ |
| `Ww` | 10.05 kN/m³ |
| `Formula_ID` | 1 |

These values match the current GUI defaults.

### 3.4 Physical and layer constants

The scripts use:

$$
\displaystyle g = 9.80665\ \text{m/s}^2
$$

$$
\displaystyle W_{\text{rock}} = 26.5\ \text{kN/m}^3
$$

$$
\displaystyle P_{\text{cubes}} = 0.40
$$

$$
\displaystyle P_{\text{rock}} = 0.25
$$

$$
\displaystyle \frac{K_{D,\text{trunk}}}{K_{D,\text{head}}} = 1.5
$$

---

## 4. Formula sets and coefficient database

The trunk stability number is evaluated with:

$$
\displaystyle N_s = \left(k_1\frac{N_{od}^{k_2}}{N_z^{k_3}} + k_4\right)
 s_{0m}^{-k_5}
$$

The current coefficient sets are:

| ID | Method | Unit | Slope | `k1` | `k2` | `k3` | `k4` | `k5` |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | Van der Meer (1988a) | Simple Cubes | 2.0:1 | 7.374304 | 0.400 | 0.300 | 1.100642 | 0.100 |
| 2 | Van der Meer (1988a) | Cubes | 1.5:1 | 6.700 | 0.400 | 0.300 | 1.000 | 0.100 |
| 3 | Chegini-Aghtouman (2006) | Antifer | 2.0:1 | 6.138 | 0.443 | 0.276 | 1.164 | 0.070 |
| 4 | Chegini-Aghtouman (2006) | Antifer | 1.5:1 | 6.951 | 0.443 | 0.291 | 1.082 | 0.082 |

These values are the ones used in the current source files.

### 4.1 Embedded precise coefficients for simple cubes at 2.0:1

For formula `1`, the slope-correction multiplier is absorbed into
`k1` and `k4`:

$$
\displaystyle k_{1,\,2.0:1} = 6.7\left(\frac{2.0}{1.5}\right)^{1/3}
 = 7.374304189198
$$

$$
\displaystyle k_{4,\,2.0:1} = 1.0\left(\frac{2.0}{1.5}\right)^{1/3}
 = 1.100642416298
$$

This removes the need for a special correction branch in the code.

---

## 5. Calculation sequence implemented by the scripts

### 5.1 Deep-water wave descriptors

The deep-water wavelength is:

$$
\displaystyle L_0 = \frac{gT_m^2}{2\pi}
$$

The deep-water wave number reported by the scripts is:

$$
\displaystyle k_0 = \frac{2\pi}{L_0}
$$

The deep-water wave steepness is:

$$
\displaystyle s_{0m} = \frac{H_s}{L_0}
$$

The scripts also store:

$$
\displaystyle s_m = s_{0m}
$$

### 5.2 Storm duration from number of waves

The aligned scripts use **number of waves** as input and compute storm duration
internally as:

$$
\displaystyle t_{\text{storm}} = \frac{N_z T_m}{3600}
$$

where `Tm` is in seconds and the result is in hours.

### 5.3 Relative density of the trunk armor unit

The relative buoyant density is:

$$
\displaystyle \Delta_{\text{trunk}} = \frac{W_{c,\text{trunk}}}{W_w} - 1
$$

Therefore:

$$
\displaystyle W_{c,\text{trunk}} = W_w\left(\Delta_{\text{trunk}} + 1\right)
$$

---

## 6. Trunk calculations

### 6.1 Expanded stability expression

The scripts evaluate the trunk stability expression through the following
intermediate terms:

$$
\displaystyle \text{term}_{\text{damage}} = N_{od}^{k_2}
$$

$$
\displaystyle \text{term}_{\text{waves}} = N_z^{k_3}
$$

$$
\displaystyle R_{dw} = \frac{N_{od}^{k_2}}{N_z^{k_3}}
$$

$$
\displaystyle \text{inv\_f} = k_1 R_{dw} + k_4
$$

$$
\displaystyle F_s = s_{0m}^{-k_5}
$$

So the trunk stability number is:

$$
\displaystyle N_{s,\text{trunk}} = \left(k_1\frac{N_{od}^{k_2}}
{N_z^{k_3}} + k_4\right)s_{0m}^{-k_5}
$$

### 6.2 Nominal diameter and unit weight

The trunk nominal diameter is:

$$
\displaystyle D_n = \frac{H_s}{\Delta_{\text{trunk}} N_{s,\text{trunk}}}
$$

The trunk unit weight is:

$$
\displaystyle W_{\text{trunk}} = W_{c,\text{trunk}} D_n^3
$$

The corresponding block mass is reported as:

$$
\displaystyle M_{\text{trunk}} = \frac{W_{\text{trunk}}}{g}
$$

### 6.3 Layer thickness and packing density

The double-layer armor thickness is:

$$
\displaystyle r_1 = 2 \times 1.1 \times D_n
$$

The armor placement density is:

$$
\displaystyle d_{\text{trunk}} = \frac{100 \times 2 \times 1.1
 \times (1 - P_{\text{cubes}})}{D_n^2}
$$

### 6.4 Equivalent Hudson coefficient for the trunk

The scripts derive an equivalent trunk Hudson coefficient as:

$$
\displaystyle K_{D,\text{trunk}} =
\frac{W_{c,\text{trunk}} H_s^3}
{W_{\text{trunk}}\,\Delta_{\text{trunk}}^3\,S}
$$

where `S` is the formula slope ratio, equal to `2.0` or `1.5`.

---

## 7. Head transfer strategy

### 7.1 Core design rule

The scripts do **not** size the head by increasing `Dn`. They keep the same
nominal diameter and transfer the trunk design through the fixed ratio:

$$
\displaystyle \frac{K_{D,\text{trunk}}}{K_{D,\text{head}}} = 1.5
$$

Therefore:

$$
\displaystyle K_{D,\text{head}} = \frac{K_{D,\text{trunk}}}{1.5}
$$

### 7.2 Relation between trunk and head relative density

With geometry kept fixed, the scripts use:

$$
\displaystyle \Delta_{\text{head}} = \Delta_{\text{trunk}} (1.5)^{1/3}
$$

This is the key transfer equation used in the code.

### 7.3 Relation between `Wc_trunk`, `Wc_head`, and the deltas

The code does **not** choose `Wc_head` independently. It is obtained directly
from the transferred head relative density.

For the trunk:

$$
\displaystyle \Delta_{\text{trunk}} = \frac{W_{c,\text{trunk}}}{W_w} - 1
$$

so:

$$
\displaystyle W_{c,\text{trunk}} = W_w\left(\Delta_{\text{trunk}} + 1\right)
$$

For the head:

$$
\displaystyle \Delta_{\text{head}} = \frac{W_{c,\text{head}}}{W_w} - 1
$$

so:

$$
\displaystyle W_{c,\text{head}} = W_w\left(\Delta_{\text{head}} + 1\right)
$$

Substituting the transfer rule gives:

$$
\displaystyle W_{c,\text{head}} =
W_w\left[\Delta_{\text{trunk}}(1.5)^{1/3} + 1\right]
$$

This makes the head-specific weight relation explicit:

- `Wc_trunk` defines `Δ_trunk`;
- `Δ_trunk` is scaled to `Δ_head`;
- `Δ_head` immediately defines `Wc_head`.

### 7.4 Relation between `W_trunk` and `W_head`

Because geometry is preserved, the head and trunk keep the same block volume:

$$
\displaystyle V_{\text{head}} = V_{\text{trunk}}
$$

Therefore the block weight scales only with concrete specific weight:

$$
\displaystyle W_{\text{head}} =
W_{\text{trunk}}\left(\frac{W_{c,\text{head}}}{W_{c,\text{trunk}}}\right)
$$

Using the delta definitions, the same relation is:

$$
\displaystyle W_{\text{head}} =
W_{\text{trunk}}\left(\frac{\Delta_{\text{head}} + 1}
{\Delta_{\text{trunk}} + 1}\right)
$$

This is the exact physical interpretation implemented in the scripts:

- the block **volume stays constant**;
- the concrete **specific weight increases**;
- the block **weight increases in the same proportion**.

### 7.5 Head stability number and packing density

With `Dn` unchanged, the reported head stability number is:

$$
\displaystyle N_{s,\text{head}} = \frac{H_s}{\Delta_{\text{head}} D_n}
$$

The head placement density is computed with the same expression as the trunk:

$$
\displaystyle d_{\text{head}} = \frac{100 \times 2 \times 1.1
 \times (1 - P_{\text{cubes}})}{D_n^2}
$$

Since the nominal diameter is unchanged:

$$
\displaystyle d_{\text{head}} = d_{\text{trunk}}
$$

---

## 8. Volume and geometric dimensions reported by the scripts

### 8.1 Volume for simple cubes

$$
\displaystyle V = D_n^3
$$

### 8.2 Volume for Antifer units (Pita, 1986)

$$
\displaystyle V = 1.0247\,H^3
$$

### 8.3 Geometric dimensions reported for all units

The scripts then report:

$$
\displaystyle H = \left(\frac{V}{1.0247}\right)^{1/3}
$$

$$
\displaystyle A = 1.086\,H
$$

$$
\displaystyle B = 1.005\,H
$$

The head transfer keeps geometry fixed, so trunk and head dimensions are
expected to be identical apart from rounding.

---

## 9. Underlayer selection from EN 13383 gradings

The repository uses **EN 13383** grading terminology and an internal mass-based
lookup table implemented in the source code. The code values below are embedded in the calculators.

### 9.1 Theoretical target weight and target mass

For any armor layer weight `W_armor`, the theoretical underlayer target is:

$$
\displaystyle W_{\text{target}} = \frac{W_{\text{armor}}}{10}
$$

The corresponding target mass is:

$$
\displaystyle M_{\text{target}} = \frac{1000\,W_{\text{target}}}{g}
$$

### 9.2 Selection rule used by the scripts

The scripts search the internal grading list for the grading that satisfies:

$$
\displaystyle \text{NLL} < M_{\text{target}} < \text{NUL}
$$

If more than one grading contains the target mass, the scripts choose the one
with the smallest nominal width:

$$
\displaystyle \Delta M_{\text{range}} = \text{NUL} - \text{NLL}
$$

If no grading strictly contains the target mass, the scripts fall back to the
first grading in the database.

### 9.3 Secondary underlayer quantities

Once a grading is selected, the scripts compute:

$$
\displaystyle M_{50} = \frac{\text{NLL} + \text{NUL}}{2}
$$

$$
\displaystyle \text{ELL} = 0.7\,\text{NLL}
$$

$$
\displaystyle \text{EUL} = 1.5\,\text{NUL}
$$

$$
\displaystyle W_{\text{mean, rock}} = \frac{M_{50} g}{1000}
$$

$$
\displaystyle D_{n,\text{rock}} =
\left(\frac{W_{\text{mean, rock}}}{W_{\text{rock}}}\right)^{1/3}
$$

$$
\displaystyle r_2 = 2 D_{n,\text{rock}}
$$

$$
\displaystyle f_2 = \frac{100 \times 2 \times 1 \times (1 - P_{\text{rock}})}
{D_{n,\text{rock}}^2}
$$

With `P_rock = 0.25`, this becomes:

$$
\displaystyle f_2 = \frac{150}{D_{n,\text{rock}}^2}
$$

### 9.4 Implemented EN 13383 grading tables

The following tables reproduce the grading values currently implemented in the
scripts.

#### 9.4.3 HMA gradings

| Class | NLL<br>(kg) | NUL<br>(kg) | M50<br>(kg) | ELL<br>(kg) | EUL<br>(kg) |
| --- | ---: | ---: | ---: | ---: | ---: |
| HMA300-1000   | 300.0   | 1000.0  | 650.0   | 210.0  | 1500.0  |
| HMA1000-3000  | 1000.0  | 3000.0  | 2000.0  | 700.0  | 4500.0  |
| HMA3000-6000  | 3000.0  | 6000.0  | 4500.0  | 2100.0 | 9000.0  |
| HMA6000-10000 | 6000.0  | 10000.0 | 8000.0  | 4200.0 | 15000.0 |
| HMA10000-15000 | 10000.0 | 15000.0 | 12500.0 | 7000.0 | 22500.0 |

#### 9.4.2 LMA gradings

| Class | NLL<br>(kg) | NUL<br>(kg) | M50<br>(kg) | ELL<br>(kg) | EUL<br>(kg) |
| --- | ---: | ---: | ---: | ---: | ---: |
| LMA5-40   | 5.0  | 40.0  | 22.5  | 3.5  | 60.0  |
| LMA10-60  | 10.0 | 60.0  | 35.0  | 7.0  | 90.0  |
| LMA15-120 | 15.0 | 120.0 | 67.5  | 10.5 | 180.0 |
| LMA40-200 | 40.0 | 200.0 | 120.0 | 28.0 | 300.0 |
| LMA60-300 | 60.0 | 300.0 | 180.0 | 42.0 | 450.0 |
| LMA15-300 | 15.0 | 300.0 | 157.5 | 10.5 | 450.0 |

#### 9.4.3 CP gradings

| Class | NLL<br>(kg) | NUL<br>(kg) | M50<br>(kg) | ELL<br>(kg) | EUL<br>(kg) |
| --- | ---: | ---: | ---: | ---: | ---: |
| CP32/90  | 0.868  | 19.319  | 10.094  | 0.608  | 28.979 |
| CP45/125 | 2.415  | 51.758  | 27.087  | 1.691  | 77.637 |
| CP63/180 | 6.626  | 154.548 | 80.587  | 4.638  | 231.822 |
| CP90/250 | 19.319 | 414.063 | 216.691 | 13.523 | 621.095 |
| CP45/180 | 2.415  | 154.548 | 78.482  | 1.691  | 231.822 |
| CP90/180 | 19.319 | 154.548 | 86.934  | 13.523 | 231.822 |

The nominal lower and upper limits shown above are the values used directly by
`breakwater_calculator_gui.cpp` and the aligned code family.

---

## 10. Output report structure

The generated text report is organised as:

1. **Input parameters**
2. **Intermediate parameters**
3. **Armor layer results — trunk**
4. **Underlayer results — trunk**
5. **Armor layer results — head**
6. **Underlayer results — head**

The aligned outputs report:

- `Number of waves (Nz)` as an **input**;
- `Storm Duration (h)` as an **intermediate calculated value**;
- trunk and head stability numbers;
- trunk and head volumes;
- trunk and head equivalent `Kd` values;
- trunk and head adopted underlayer classes and derived rock quantities.

---

## 11. Compilation and execution instructions

### 11.1 Recommended toolchains

For the C++ files:

- **GCC / g++** with C++17 support;
- on Windows, **MinGW-w64** is the simplest option for direct command-line
  compilation.

For the Fortran file:

- **GNU Fortran (`gfortran`)** with Fortran 2008 support.

For the Python file:

- Python 3.x;
- no third-party packages are required.

### 11.2 General preparation steps

1. Put all source files in one working folder.
2. Open a terminal in that folder.
3. Check the compiler versions before building.
4. Compile one target at a time.
5. Run the executable from the same folder if you want `output.txt` to be
   written there.

Useful version checks:

```bash
g++ --version
gfortran --version
python --version
```

### 11.3 Python command-line version

#### 11.3.1 Run directly

```bash
python breakwater_calculator.py
```

#### 11.3.2 Typical interactive workflow

1. Run the script.
2. Enter `Hs`, `Tm`, `Nz`, `Nod`, `Wc`, and `Formula_ID` when prompted.
3. Review the report printed on screen.
4. Check `output.txt` in the same folder.

### 11.4 C++ command-line version

#### 11.4.1 Recommended GCC build

```bash
g++ \
  -O3 \
  -march=native \
  -std=c++17 \
  -Wall \
  -Wextra \
  -static \
  -static-libgcc \
  -static-libstdc++ \
  -o breakwater_calculator_cli \
  breakwater_calculator_cli.cpp
```

#### 11.4.2 What the flags mean

| Flag | Meaning |
| --- | --- |
| `-O3` | high optimization level |
| `-march=native` | optimize for the current CPU |
| `-std=c++17` | compile as C++17 |
| `-Wall -Wextra` | enable common warnings |
| `-static` | prefer static linking |
| `-static-libgcc` | statically link GCC runtime |
| `-static-libstdc++` | statically link C++ runtime |
| `-o ...` | define output executable name |

#### 11.4.3 Run interactively

```bash
./breakwater_calculator_cli
```

#### 11.4.4 Run with explicit arguments

```bash
./breakwater_calculator_cli \
  11.0 \
  11.9 \
  3000 \
  0.5 \
  27.48 \
  1
```

Argument order:

```text
[Hs] [Tm] [NumberOfWaves] [Nod] [Wc] [FormulaID]
```

### 11.5 C++ GUI, English version

#### 11.5.1 Recommended build with MinGW-w64

```bash
g++ \
  -O3 \
  -march=native \
  -std=c++17 \
  -municode \
  -mwindows \
  -static \
  -static-libgcc \
  -static-libstdc++ \
  -o breakwater_calculator_gui.exe \
  breakwater_calculator_gui.cpp
```

#### 11.5.2 Notes on the GUI build

- `-municode` is needed because the GUI uses wide-character Win32 calls.
- `-mwindows` builds a windowed application instead of a console program.
- the output name ends with `.exe` because this target is intended for Windows.

#### 11.5.3 Running the GUI

1. Start `breakwater_calculator_gui.exe`.
2. Enter the design inputs.
3. Select the formula.
4. Click **Calculate Armor**.
5. Review the report in the GUI and in `output.txt`.

### 11.6 C++ GUI, Portuguese version

#### 11.6.1 Recommended build with MinGW-w64

```bash
g++ \
  -O3 \
  -march=native \
  -std=c++17 \
  -municode \
  -mwindows \
  -static \
  -static-libgcc \
  -static-libstdc++ \
  -o breakwater_calculator_gui_pt.exe \
  breakwater_calculator_gui_pt.cpp
```

#### 11.6.2 Running the Portuguese GUI

1. Start `breakwater_calculator_gui_pt.exe`.
2. Enter the parameters in Portuguese.
3. Choose the formula.
4. Click **Calcular Manto**.
5. Review the report in the window and in `output.txt`.

### 11.7 Fortran command-line version

#### 11.7.1 Recommended GNU Fortran build

```bash
gfortran \
  -O3 \
  -march=native \
  -std=f2008 \
  -Wall \
  -Wextra \
  -pedantic \
  -Wconversion \
  -static \
  -static-libgfortran \
  -static-libgcc \
  -o breakwater_calculator_cli \
  breakwater_calculator_cli.f90
```

#### 11.7.2 What the main Fortran flags mean

| Flag | Meaning |
| --- | --- |
| `-O3` | high optimization level |
| `-march=native` | optimize for the current CPU |
| `-std=f2008` | compile as Fortran 2008 |
| `-Wall -Wextra` | enable common warnings |
| `-pedantic` | request stricter standards checking |
| `-Wconversion` | warn about implicit numeric conversions |
| `-static` | prefer static linking |
| `-static-libgfortran` | statically link Fortran runtime |
| `-static-libgcc` | statically link GCC runtime |

#### 11.7.3 Run interactively

```bash
./breakwater_calculator_cli
```

#### 11.7.4 Run with explicit arguments

```bash
./breakwater_calculator_cli \
  11.0 \
  11.9 \
  3000 \
  0.5 \
  27.48 \
  1
```

Argument order:

```text
[Hs] [Tm] [NumberOfWaves] [Nod] [Wc] [FormulaID]
```

### 11.8 If static linking fails

Some systems do not provide all static runtime libraries by default. If a fully
static build fails, try the same command **without** the static-linking flags.

For C++:

```bash
g++ \
  -O3 \
  -march=native \
  -std=c++17 \
  -Wall \
  -Wextra \
  -o breakwater_calculator_cli \
  breakwater_calculator_cli.cpp
```

For Fortran:

```bash
gfortran \
  -O3 \
  -march=native \
  -std=f2008 \
  -Wall \
  -Wextra \
  -pedantic \
  -Wconversion \
  -o breakwater_calculator_cli \
  breakwater_calculator_cli.f90
```

### 11.9 Notes on Windows paths and execution

If you compile on Windows with MinGW-w64 inside `cmd.exe`, the executable is
usually run as:

```bat
breakwater_calculator_cli.exe
breakwater_calculator_gui.exe
breakwater_calculator_gui_pt.exe
```

If you compile inside PowerShell, you may need:

```powershell
.\breakwater_calculator_cli.exe
.\breakwater_calculator_gui.exe
.\breakwater_calculator_gui_pt.exe
```

---

## 12. Output-file behavior by implementation

The implementations do not all treat `output.txt` the same way.

### 12.1 Overwrite behavior

These implementations overwrite `output.txt` each time they run:

- `breakwater_calculator.py`
- `breakwater_calculator_cli.cpp`

### 12.2 Append behavior

These implementations append each new report to `output.txt`:

- `breakwater_calculator_gui.cpp`
- `breakwater_calculator_gui_pt.cpp`
- `breakwater_calculator_cli.f90`

---

## 13. Notes on interpretation

1. This repository is a **preliminary hydraulic sizing tool**, not a complete
   breakwater design system.
2. The scripts do **not** perform overtopping analysis, toe stability design,
   geotechnical verification, diffraction analysis, or physical-model
   validation.
3. The slope term in the Hudson-equivalent expression is used exactly as stored
   in the script formula database.
4. Some report labels still say **Cube Height**, **Cube Top Width**, and
   **Cube Base Width** even for Antifer calculations. Those labels are part of
   the current script outputs.
5. The head-transfer model is intentionally direct:
   **same geometry, higher density, higher block weight**.

---

## 14. Compact summary of the implemented model

The engineering chain implemented by the aligned scripts is:

$$
\displaystyle
(H_s, T_m, N_z, N_{od}, W_{c,\text{trunk}}, \text{Formula ID})
\rightarrow N_{s,\text{trunk}}
\rightarrow D_n
\rightarrow W_{\text{trunk}}
\rightarrow K_{D,\text{trunk}}
$$

$$
\displaystyle
\rightarrow \Delta_{\text{head}}
\rightarrow W_{c,\text{head}}
\rightarrow W_{\text{head}}
\rightarrow \text{underlayer grading}
$$

This is the exact conceptual sequence implemented by the source code.
