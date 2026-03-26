# Hydraulic Stability Calculator for Breakwater Armor Units
## Simple Cubes and Antifer Blocks

## 1. Purpose and scope

This repository contains a family of calculators for the **preliminary hydraulic sizing** of rubble-mound breakwater armor made with **simple concrete cubes** and **Antifer blocks**, together with the associated **rock underlayers**. The repository is intentionally focused on the hydraulic-stability stage of conceptual and preliminary structural definition. It is not a substitute for full breakwater design, which must also address wave transformation, toe stability, overtopping, crest and crown-wall performance, geotechnical resistance, construction method, settlement, durability, and physical-model validation for important structures (Van der Meer, 1995; Van der Meer, 1998; USACE, 1995; USACE, 2006).

The common engineering workflow implemented by the aligned calculators is:

- size the **trunk armor** with empirical hydraulic-stability formulae for concrete armor units;
- derive an **equivalent Hudson stability coefficient** for the trunk for comparison and transfer purposes;
- transfer the trunk design to the **head** through a fixed stability-ratio rule expressed in Hudson terms;
- keep the **same nominal armor size** at trunk and head;
- increase the **required concrete specific weight** at the head instead of increasing the nominal size;
- size the **rock underlayers** from an embedded EN 13383 grading database;
- optionally derive **custom interpolated underlayer gradings** by family (`AUTO`, `HMA`, `LMA`, or `CP`) when strict EN 13383 selection is not used or is not possible.

The design philosophy implemented in the code is therefore **iso-geometric** but **not iso-weight**:

- the nominal armor-unit size is kept constant between trunk and head;
- the unit geometry is kept constant between trunk and head;
- the slope is kept constant between trunk and head;
- the head is stabilized by increasing the required concrete specific weight;
- the head unit becomes heavier because density increases while volume is kept constant.

This is a practical design model for concept definition, cost screening, and early option studies. As emphasized repeatedly in the attached technical sources, formula-based design of rubble-mound armor remains fundamentally empirical and test-dependent, and major projects should normally be confirmed by physical modeling before final design is frozen (Van der Meer, 1995; Van der Meer, 1998; USACE, 1995; USACE, 2006).

---

## 2. What is covered by this README

This README documents the verified implementations reviewed in the current codebase:

- `breakwater_calculator_cli.cpp` — C++ command-line reference;
- `breakwater_calculator.py` — Python command-line version aligned to the C++ CLI;
- `breakwater_calculator_cli.f90` — Fortran command-line version aligned to the C++ CLI;
- `breakwater_calculator_gui.cpp` — Win32 GUI version using the same hydraulic kernel;
- `breakwater_calculator.xlsx` — spreadsheet implementation of the same engineering chain;
- `breakwater_calculator.ipynb` — notebook implementation;
- `nomograms_vandermeer-chegini.py` — nomograms for the trunk stability formulae;
- `nomogram_hudson_trunk_to_head.py` — nomogram for the trunk-to-head Hudson transfer.

The most explicit computational reference is the current **`breakwater_calculator_cli.cpp`**, because it exposes the formula database, the underlayer-selection logic, the custom-grading interpolation rules, and the full report structure in a single source file. The Python and Fortran command-line versions were aligned to that behavior, and the GUI uses the same hydraulic logic behind a graphical interface.

---

## 3. Engineering basis of the calculator family

### 3.1 Conceptual-design character of the implemented formulae

The theoretical basis of the repository is the classical empirical design framework for rubble-mound armor. That framework begins with Hudson-type stability-number thinking, progresses to the more refined Van der Meer formulations for rock, and then extends to special concrete-armour-unit relationships expressed in terms of wave steepness, damage, and storm length. The attached sources explicitly stress that these formulae are derived from laboratory tests on schematized structures and must therefore be used as **conceptual and preliminary design tools**, not as stand-alone final design substitutes (Van der Meer, 1995; Van der Meer, 1998; USACE, 1995; USACE, 2006).

### 3.2 Why the repository uses concrete-unit formulae rather than generic rock formulae

The repository sizes the **armor layer** with specialized empirical equations for **double-layer simple cubes** and **Antifer blocks**. It does **not** directly implement the general Van der Meer rock equations for armorstone sizing, even though those equations are essential to understand the theoretical evolution of stability design. This distinction is important:

- the **general rock formulae** use wave type, damage level, storm duration, and a notional permeability factor explicitly (Van der Meer, 1998);
- the **concrete-unit formulae** implemented here are reduced relationships calibrated for specific armor families and specific slope regimes, so they are expressed in terms of `Nod`, `Nz`, and `som`, with slope effects embedded through the selected coefficient set rather than through an explicit `ξm` transition function (Van der Meer, 1999; Chegini and Aghtouman, 2006).

### 3.3 Role of Hudson in the repository

Hudson is used in this repository in a **derived** and **comparative** way, not as the primary trunk design formula for cubes or Antifer units. More specifically:

1. the trunk armor is first sized by the selected Van der Meer- or Chegini-type formula;
2. the code then derives an **equivalent** `KD_trunk` by back-substituting the resulting unit size into the Hudson form;
3. the head is transferred from the trunk through the adopted fixed ratio
   `KD_trunk / KD_head = 1.5`.

This means Hudson serves as a **transfer metric** and as a familiar interpretive quantity, which is consistent with the long-standing use of `KD` values in SPM, CEM, and the Rock Manual for comparative armor performance and for engineering judgment in preliminary design (CERC, 1984; USACE, 2006; CIRIA et al., 2007).

---

## 4. Governing hydraulic and structural parameters

### 4.1 Incident wave height at the structure

In the attached Van der Meer material, the governing wave height for armor design is the **incident** wave height at the structure toe, usually expressed as `Hs` or `Hm0`, depending on the analysis basis. In deep water these are often close, but in shallow or depth-limited conditions they may diverge, and the selected characteristic value matters (Van der Meer, 1995; Van der Meer, 1998; Van der Meer, 2011). The current code uses the user-entered `Hs` directly as the design wave height attacking the structure.

### 4.2 Mean period and wave steepness

The period used by the repository is the **mean wave period** `Tm`. Following Van der Meer, the deep-water wave steepness based on `Tm` is:

$$
\displaystyle s_{0m} = \frac{2\pi H_s}{gT_m^2} = \frac{H_s}{L_0}
$$

with the corresponding deep-water wavelength:

$$
\displaystyle L_0 = \frac{gT_m^2}{2\pi}
$$

The quantity `som` is a dimensionless period descriptor. It is not merely a geometric steepness; it carries information about the incident wave climate and strongly affects hydraulic stability in both the generalized Van der Meer framework and the specialized concrete-unit relationships used in this repository (Van der Meer, 1995; Van der Meer, 1998; Van der Meer, 1999).

### 4.3 Breaker parameter and its role in the theory

Van der Meer defines the breaker or surf-similarity parameter as:

$$
\displaystyle \xi = \frac{\tan\alpha}{\sqrt{s}}
$$

where `α` is the slope angle and `s` is a suitable wave steepness. In the rock formulae this parameter controls the transition between plunging and surging behavior and therefore changes the structure of the stability equation itself (Van der Meer, 1998). The repository does not use `ξm` directly in the implemented cube and Antifer sizing equations, because those formula families were reduced to `Nod`, `Nz`, and `som` forms for the specific tested cross-sections (Van der Meer, 1999).

### 4.4 Relative buoyant density

Van der Meer expresses the relative buoyant density as:

$$
\displaystyle \Delta = \frac{\rho_r-\rho_w}{\rho_w}
$$

The code works with **specific weights** rather than mass densities, but the ratio is equivalent in form, so the implemented expression is:

$$
\displaystyle \Delta = \frac{W_c}{W_w} - 1
$$

This is a key variable because the stability number is defined as `Hs/(ΔDn)` for concrete units or `Hs/(ΔDn50)` for graded rock (Van der Meer, 1998).

### 4.5 Nominal diameter

For rock, the nominal diameter is classically defined from the representative mass:

$$
\displaystyle D_{n50} = \left(\frac{M_{50}}{\rho_r}\right)^{1/3}
$$

For concrete armor units of uniform size, the notation usually becomes simply `Dn`, because the layer is not described by a grading curve in the same way as rock (Van der Meer, 1998). The current software uses `Dn` as the characteristic unit size for cubes and Antifer blocks.

### 4.6 Damage descriptors

The repository uses `Nod`, following the concrete-unit damage convention described by Van der Meer. In that convention, `Nod` is the number of displaced units normalized over a strip of width one nominal diameter `Dn` along the structure axis. It is therefore a relative hydraulic-damage measure rather than a percentage directly. If the number of units in a one-`Dn`-wide strip is known, `Nod` can be related to a damage percentage by simple ratio (Van der Meer, 1998; Van der Meer, 1999).

This is different from the rock damage parameter `S`, which is defined from erosion area around still-water level:

$$
\displaystyle S = \frac{A_e}{D_{n50}^2}
$$

The distinction matters. The code uses `Nod` for concrete units and does not mix it with the rock damage definition `S`.

---

## 5. Hudson framework and its role in interpretation

### 5.1 Classical Hudson equation

The classical Hudson equation for a nearly uniform armor layer can be written in weight form as:

$$
\displaystyle W = \frac{W_r H^3}{K_D \Delta^3 \cot\alpha}
$$

where `W` is the required armor-unit weight, `Wr` is the unit specific weight, `H` is the design wave height at the structure, `KD` is the stability coefficient, `Δ` is the relative buoyant density, and `cot α` is the cotangent of the slope angle measured from the horizontal (CERC, 1984; USACE, 1995; USACE, 2006).

The same relationship can be rewritten as a stability-number expression:

$$
\displaystyle N_s = \frac{H}{\Delta D_n} = (K_D\cot\alpha)^{1/3}
$$

For graded riprap, the same idea is applied to the 50% mass, with the caution that grading effects become important and the rock class must be described by a mass range rather than a single block weight (CERC, 1984).

### 5.2 What Hudson captures well

Hudson remains attractive because it is compact, easy to interpret, and linked to a wide body of KD values for different armor types and placements. It is still extensively used in manuals for preliminary sizing and for comparative interpretation of different unit families (CERC, 1984; USACE, 2006; CIRIA et al., 2007).

### 5.3 What Hudson does not capture well

The attached sources also emphasize the limitations of Hudson:

- it was originally tied to **regular-wave** logic;
- it does not explicitly include **wave period**;
- it does not explicitly include **storm duration**;
- it does not describe **damage development** in a graded way;
- it is most naturally associated with **non-overtopped, permeable-core** situations;
- the simple `KD cot α` combination does not always capture slope effects optimally (Van der Meer, 1998; CERC, 1984; USACE, 2006).

This is exactly why the repository uses Van der Meer- and Chegini-type relations for trunk sizing and keeps Hudson as an equivalent interpretation variable rather than the primary sizing law.

### 5.4 KD values for cubes and Antifer in the manuals

The Rock Manual excerpts attached to this repository show the classical double-layer Hudson `KD` guidance for concrete units. For **double-layer cubes**, values of approximately `KD = 6.5` for breaking-wave trunk conditions and `KD = 7.5` for non-breaking-wave trunk conditions are reported; for **Antifer cubes**, the same table indicates `KD = 7` for breaking waves and `KD = 8` for non-breaking waves on the trunk, for typical slopes in the 1:1.5 to 1:3 range for cubes and around 1:2 for Antifer. The same excerpts stress that “breaking waves” in this context means breaking on the **foreshore approaching the structure**, not breaking induced by the armor slope itself (CIRIA et al., 2007; CERC, 1984; USACE, 2006).

Those KD values are not used directly by the software to size the trunk armor. Instead, they provide context for the equivalent KD values derived after the empirical trunk formula has already produced a unit size.

---

## 6. Van der Meer theoretical framework

### 6.1 General rock formulae

The general Van der Meer deep-water rock framework distinguishes **plunging** and **surging** wave attack. In the attached chapter, the plunging-wave equation is given in the familiar form:

$$
\displaystyle \frac{H_s}{\Delta D_{n50}} = 6.2\,P^{0.18}\left(\frac{S}{\sqrt{N}}\right)^{0.2}\xi_m^{-0.5}
$$

and the surging-wave formulation is given as a separate equation involving the same governing variables together with slope influence through `cot α` and the breaker parameter `ξm` (Van der Meer, 1998).

The key theoretical point is not merely the exact exponents. It is that Van der Meer replaces Hudson’s single-coefficient logic with a richer description in which stability depends on:

- **wave height** through `Hs`;
- **unit density and size** through `ΔDn50`;
- **storm length** through `N`;
- **damage level** through `S`;
- **wave period and breaker type** through `som` and `ξm`;
- **structural permeability** through the notional factor `P`.

This is a major conceptual advance relative to Hudson and explains why Van der Meer-type formulations are widely regarded as more precise for rubble-mound stability assessment (Van der Meer, 1998).

### 6.2 Damage level `S`

In the rock framework, the damage level is:

$$
\displaystyle S = \frac{A_e}{D_{n50}^2}
$$

where `Ae` is the erosion area around still-water level. Van der Meer explains that this definition is scale-independent in a useful engineering sense and that it incorporates both settlement and displacement. For two-diameter-thick rock armor layers, the attached tables indicate that **start of damage** is typically around `S = 2–3`, while larger `S` values correspond to intermediate damage and eventually underlayer exposure, with the exact interpretation depending strongly on slope (Van der Meer, 1998).

### 6.3 Notional permeability factor `P`

The notional permeability factor `P` is not a porosity measure. It is an engineering descriptor of how permeable the supporting cross-section is from a hydraulic-stability point of view. Van der Meer indicates that very impermeable arrangements may correspond to `P ≈ 0.1`, while very permeable, almost homogeneous rock structures may approach `P ≈ 0.6` (Van der Meer, 1998). The attached material is explicit that higher permeability generally increases stability because more water can penetrate the structure during run-up and run-down, reducing destabilizing forces on the armor layer.

### 6.4 Why this matters even though the repository does not directly size rock armor with Van der Meer

The underlayers in this repository are **not** hydraulically sized by the general Van der Meer rock equations. They are selected by practical grading criteria anchored to the armor-unit weight. Even so, the Van der Meer theory remains highly relevant because it explains why permeability, slope regime, storm duration, and damage definition matter in rubble-mound behavior. Those concepts are essential background for interpreting when a simple `W/10` underlayer rule is acceptable and when more detailed stability or filter studies would be warranted (Van der Meer, 1995; Van der Meer, 1998; USACE, 1995).

### 6.5 Reliability and uncertainty

The attached Van der Meer chapter also discusses reliability. It notes that the Van der Meer coefficients may be treated statistically and that the resulting formulae exhibit lower scatter than a simple Hudson treatment based on a mean `KD`, while still requiring sensitivity analysis and engineering judgment. This is directly relevant to software use: the repository produces deterministic outputs, but those outputs should be interpreted as one point in an uncertainty band, not as an exact truth (Van der Meer, 1998).

---

## 7. Van der Meer concrete-unit formula for simple cubes

### 7.1 The published two-layer cube relationship

For double-layer cubes, Van der Meer gives the stability-number relation:

$$
\displaystyle \frac{H_s}{\Delta D_n} = \left(6.7\,\frac{N_{od}^{0.4}}{N^{0.3}} + 1.0\right)s_{om}^{-0.1}
$$

for the tested steep-slope cube configuration. This is the core published relationship behind the cube implementation in the repository (Van der Meer, 1999; CIRIA et al., 2007).

The equation says that cube stability increases when:

- `Nod` is allowed to be larger;
- the number of waves `N` is lower;
- the wave steepness `som` is lower.

In other words, the same cube size can withstand a more severe sea state if the designer accepts more displacement damage, a shorter storm, or a milder wave steepness.

### 7.2 Why the formula uses `som` and not `ξm`

Van der Meer explains that for the concrete-unit research in question only one slope angle was investigated for each armor family, so the influence of wave period should not be expressed through `ξm`, because that parameter mixes slope and wave steepness. The period effect was therefore carried explicitly by `som` in the published equations (Van der Meer, 1998; Van der Meer, 1999). This is exactly the theoretical logic preserved in the code.

### 7.3 Formula form implemented by the code

The software implements the cube relation in the generalized five-coefficient form:

$$
\displaystyle N_s = \left(k_1\frac{N_{od}^{k_2}}{N_z^{k_3}} + k_4\right)s_{0m}^{-k_5}
$$

with the following coefficient sets:

| Formula ID | Description | Slope | `k1` | `k2` | `k3` | `k4` | `k5` |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | Van der Meer (1988a) | 2.0:1 | 7.374304 | 0.400 | 0.300 | 1.100642 | 0.100 |
| 2 | Van der Meer (1988a) | 1.5:1 | 6.700 | 0.400 | 0.300 | 1.000 | 0.100 |

### 7.4 Exact interpretation of the 2.0:1 cube coefficients

For the `2.0:1` cube case, the repository does not introduce a separate correction branch after evaluating the base 1.5:1 formula. Instead, it embeds the slope transfer directly into the coefficients:

$$
\displaystyle k_1 = 6.7\left(\frac{2.0}{1.5}\right)^{1/3} = 7.374304189198
$$

$$
\displaystyle k_4 = 1.0\left(\frac{2.0}{1.5}\right)^{1/3} = 1.100642416298
$$

This is mathematically convenient because it preserves the same five-parameter form while absorbing the slope change into the calibrated constants. It also matches the nomogram source provided with the repository.

### 7.5 No-damage interpretation

For `Nod = 0`, the cube formula reduces to:

$$
\displaystyle \frac{H_s}{\Delta D_n} = 1.0\,s_{om}^{-0.1}
$$

for the base 1.5:1 form. Van der Meer explicitly notes that `Nod = 0` is a severe criterion and leads to larger units. In many cases, a small but nonzero damage allowance such as `Nod = 0.5` is more economical and more consistent with practical conceptual design (Van der Meer, 1999).

---

## 8. Chegini-Aghtouman concrete-unit formula for Antifer blocks

### 8.1 Role of the Antifer equations in the repository

The Antifer implementation follows the same generalized structure used for cubes:

$$
\displaystyle N_s = \left(k_1\frac{N_{od}^{k_2}}{N_z^{k_3}} + k_4\right)s_{0m}^{-k_5}
$$

but with a different coefficient set for each slope. This reflects the fact that Antifer stability was calibrated independently from cubes and that the unit geometry and hydraulic behavior are different (Chegini and Aghtouman, 2006).

### 8.2 Exact Antifer coefficient sets used by the code

| Formula ID | Description | Slope | `k1` | `k2` | `k3` | `k4` | `k5` |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3 | Chegini (2006) | 2.0:1 | 6.138 | 0.443 | 0.276 | 1.164 | 0.070 |
| 4 | Chegini (2006) | 1.5:1 | 6.951 | 0.443 | 0.291 | 1.082 | 0.082 |

### 8.3 Technical meaning of the Antifer exponents

Compared with the cube relationship, the Antifer exponents used by the repository indicate a somewhat different sensitivity to:

- accepted damage `Nod`;
- number of waves `Nz`;
- mean wave steepness `som`.

That is exactly what should be expected. Antifer blocks are not cubes with a different name; they are a distinct armor family with different shape, hydraulic interaction, and stability response. The code therefore treats them with their own coefficient sets rather than with geometric correction factors applied to the cube law.

### 8.4 Interpretation of the Antifer formula in use

The repository does not implement a separate explicit Hudson design law for Antifer units. It first sizes the Antifer block by the selected Chegini-Aghtouman equation and only then derives an equivalent `KD` value by back-substitution into Hudson form. This keeps the sizing logic consistent with the empirical tests while still reporting a `KD`-type quantity familiar to coastal engineers.

---

## 9. Exact computational chain implemented by the aligned calculators

### 9.1 Input variables

The aligned calculators expose the following principal engineering inputs:

| Symbol | Meaning | Unit |
| --- | --- | --- |
| `Hs` | significant wave height | m |
| `Tm` | mean wave period | s |
| `Nz` | number of waves in the design storm | – |
| `Nod` | relative damage number for concrete units | – |
| `Wc` | trunk concrete specific weight | kN/m³ |
| `Formula_ID` | formula selector | – |
| `UseEN13383` | standard EN 13383 underlayer selection enabled | true/false |
| `CustomFamily` | custom underlayer family | `AUTO`, `HMA`, `LMA`, `CP` |

### 9.2 Internal constants used by the code

The aligned core calculators currently use:

$$
\displaystyle g = 9.80665\ \text{m/s}^2
$$

$$
\displaystyle W_{rock} = 26.5\ \text{kN/m}^3
$$

$$
\displaystyle P_{cubes} = 0.40
$$

$$
\displaystyle P_{rock} = 0.25
$$

$$
\displaystyle \frac{K_{D,trunk}}{K_{D,head}} = 1.5
$$

The last ratio is an implementation assumption used to transfer the trunk design to the head while preserving geometry.

### 9.3 Deep-water intermediates

The code computes:

$$
\displaystyle L_0 = \frac{gT_m^2}{2\pi}
$$

$$
\displaystyle k_0 = \frac{2\pi}{L_0}
$$

$$
\displaystyle s_{0m} = \frac{H_s}{L_0}
$$

and stores `sm = s0m` in the current implementation.

### 9.4 Storm duration derived from the number of waves

The storm duration is not entered directly. Instead, the code computes:

$$
\displaystyle t_{storm} = \frac{N_z T_m}{3600}
$$

with `Tm` in seconds and the result in hours. This is consistent with the general Van der Meer emphasis on storm length as a physically relevant driver of cumulative damage, even though the concrete-unit formulae use `Nz` directly rather than duration in hours (Van der Meer, 1998; Van der Meer, 1999).

### 9.5 Trunk relative density and stability number

The trunk buoyant relative density is:

$$
\displaystyle \Delta_{trunk} = \frac{W_{c,trunk}}{W_w} - 1
$$

The trunk stability number is then evaluated by the selected coefficient set:

$$
\displaystyle N_{s,trunk} = \left(k_1\frac{N_{od}^{k_2}}{N_z^{k_3}} + k_4\right)s_{0m}^{-k_5}
$$

### 9.6 Trunk nominal dimension and weight

Once `Ns_trunk` is known, the code solves:

$$
\displaystyle D_n = \frac{H_s}{\Delta_{trunk}N_{s,trunk}}
$$

and then computes the armor-unit weight:

$$
\displaystyle W_{trunk} = W_{c,trunk}D_n^3
$$

with reported mass:

$$
\displaystyle M_{trunk} = \frac{W_{trunk}}{g}
$$

The use of `W/g` is numerically consistent with reporting mass in tonnes when `W` is expressed in kN and `g` in m/s².

### 9.7 Trunk armor thickness and packing density

The software uses the following implementation rules for the armor layer:

$$
\displaystyle r_1 = 2\times1.1\times D_n
$$

$$
\displaystyle d = \frac{100\times2\times1.1\times(1-P_{cubes})}{D_n^2}
$$

These are code-level formulas used for reporting the double-layer thickness and the notional placement density in units per 100 m². They should be interpreted as implementation conventions linked to the assumed double-layer arrangement, not as universal rules applicable to every possible construction method.

### 9.8 Equivalent trunk Hudson coefficient

The equivalent trunk Hudson coefficient reported by the code is:

$$
\displaystyle K_{D,trunk} = \frac{W_{c,trunk}H_s^3}{W_{trunk}\,\Delta_{trunk}^3\,S}
$$

where `S` in the code is the stored numerical slope ratio, i.e. `2.0` or `1.5`, corresponding to `cot α` in the Hudson formulation. Substituting `W_trunk = Wc_trunk Dn^3` gives:

$$
\displaystyle K_{D,trunk} = \frac{1}{\cot\alpha}\left(\frac{H_s}{\Delta_{trunk}D_n}\right)^3 = \frac{N_{s,trunk}^3}{\cot\alpha}
$$

This shows clearly that the reported `KD_trunk` is mathematically consistent with the solved trunk stability number.

---

## 10. Trunk-to-head transfer strategy

### 10.1 Design intent

The breakwater head is usually more vulnerable than the trunk because of three-dimensional flow concentration, stronger local turbulence, and more complex attack patterns. Manuals therefore commonly distinguish trunk and head stability. The software adopts a deliberately simple but explicit transfer model: the **head is not resized geometrically**; instead, its required density is increased so that the same `Dn` can be retained (CERC, 1984; USACE, 2006; CIRIA et al., 2007).

### 10.2 Fixed stability-ratio rule

The software sets:

$$
\displaystyle \frac{K_{D,trunk}}{K_{D,head}} = 1.5
$$

and therefore:

$$
\displaystyle K_{D,head} = \frac{K_{D,trunk}}{1.5}
$$

This is an implementation rule, not a universal closed-form law of nature. It is a practical translation of the familiar engineering principle that head units must be effectively stronger or heavier than trunk units when geometry is held unchanged.

### 10.3 Required head relative density

Holding `H`, `Dn`, and `cot α` fixed in Hudson form implies that `KD` varies inversely with `Δ^3`. The code therefore computes:

$$
\displaystyle \Delta_{head} = \Delta_{trunk}(1.5)^{1/3}
$$

This is the exact expression required for a fixed-size transfer under the adopted `KD` ratio.

### 10.4 Required head concrete specific weight

The head concrete specific weight follows directly:

$$
\displaystyle W_{c,head} = W_w(\Delta_{head}+1)
$$

### 10.5 Head unit weight at constant volume

Because the nominal size and geometry are kept constant, the head unit weight scales only with specific weight:

$$
\displaystyle W_{head} = W_{trunk}\left(\frac{W_{c,head}}{W_{c,trunk}}\right)
$$

The code then reports the corresponding head mass as:

$$
\displaystyle M_{head} = \frac{W_{head}}{g}
$$

### 10.6 Head stability number reported by the software

The head stability number is back-computed as:

$$
\displaystyle N_{s,head} = \frac{H_s}{\Delta_{head}D_n}
$$

This is not an independently fitted head formula. It is the stability number implied by the transferred density while maintaining the original `Dn`.

### 10.7 Practical interpretation

This means the repository implements the following engineering position:

- **same mold** or same geometric template at trunk and head;
- **same nominal block size**;
- **same slope**;
- **higher concrete specific weight** at the head;
- **higher block weight** at the head as a consequence of density increase.

For early-stage design this is a transparent and useful approach, especially where construction logistics favor common geometry and differential density is feasible. For final design, however, the head should still be checked against project-specific hydraulic modeling and structural-strength requirements (Van der Meer, 1995; USACE, 2006).

---

## 11. Armor-unit geometry and reported dimensions

### 11.1 Cubes

For simple cubes the code reports the armor volume as:

$$
\displaystyle V = D_n^3
$$

This is fully consistent with the nominal cube interpretation.

### 11.2 Antifer blocks

For Antifer blocks the code reports volume through:

$$
\displaystyle V = 1.0247\,H^3
$$

and derives the geometric dimensions:

$$
\displaystyle H = \left(\frac{V}{1.0247}\right)^{1/3}
$$

$$
\displaystyle A = 1.086H
$$

$$
\displaystyle B = 1.005H
$$

where `H` is reported as block height, `A` as base width, and `B` as top width. These are the exact dimension formulas currently used by the code and are therefore documented here as implementation definitions.

### 11.3 Consequence of iso-geometric transfer

Because the head design preserves geometry, the trunk and head Antifer dimensions are identical except for roundoff. The head weight changes only because the material specific weight changes.

---

## 12. Underlayer sizing logic

### 12.1 Scope of the underlayer logic

The underlayer calculations apply **only to the rock underlayers**. They do **not** alter the hydraulic trunk-sizing formula for cubes or Antifer, and they do **not** alter the trunk-to-head armor transfer logic. They only affect the selected or derived rock grading assigned to the underlayer.

### 12.2 Theoretical target weight

The software adopts the current practical rule:

$$
\displaystyle W_{target} = \frac{W_{armor}}{10}
$$

This `W/10` rule is a pragmatic underlayer-sizing convention. It is widely recognizable in preliminary breakwater practice, but it should not be confused with a full filter-design verification or a complete hydraulic-stability proof for the underlayer as an independent armor layer. When project conditions are severe, a dedicated underlayer and filter review is still required (USACE, 1995; CIRIA et al., 2007).

### 12.3 Target mass corresponding to the target weight

The code converts the target weight to an equivalent target mass:

$$
\displaystyle M_{target} = \frac{1000W_{target}}{g}
$$

This is the mass against which the grading bands are compared.

---

## 13. Standard EN 13383 grading selection

### 13.1 Grading families embedded in the code

The software contains three grading families:

- `CP`;
- `LMA`;
- `HMA`.

The embedded nominal lower and upper limits are:

#### 13.1.1 CP gradings

| Class | NLL (kg) | NUL (kg) |
| --- | ---: | ---: |
| CP32/90 | 0.868 | 19.319 |
| CP45/125 | 2.415 | 51.758 |
| CP63/180 | 6.626 | 154.548 |
| CP90/250 | 19.319 | 414.063 |
| CP45/180 | 2.415 | 154.548 |
| CP90/180 | 19.319 | 154.548 |

#### 13.1.2 LMA gradings

| Class | NLL (kg) | NUL (kg) |
| --- | ---: | ---: |
| LMA5-40 | 5 | 40 |
| LMA10-60 | 10 | 60 |
| LMA15-120 | 15 | 120 |
| LMA40-200 | 40 | 200 |
| LMA60-300 | 60 | 300 |
| LMA15-300 | 15 | 300 |

#### 13.1.3 HMA gradings

| Class | NLL (kg) | NUL (kg) |
| --- | ---: | ---: |
| HMA300-1000 | 300 | 1000 |
| HMA1000-3000 | 1000 | 3000 |
| HMA3000-6000 | 3000 | 6000 |
| HMA6000-10000 | 6000 | 10000 |
| HMA10000-15000 | 10000 | 15000 |

These values are the exact internal grading ranges used by the current code. The software documentation should therefore refer to **EN 13383**, not EN 13385, for the standard grading mode.

### 13.2 Strict containment criterion

When `UseEN13383 = true`, the code looks for classes satisfying the strict condition:

$$
\displaystyle NLL < M_{target} < NUL
$$

The comparison is intentionally **strict**. A target mass exactly equal to a boundary does not count as “inside” the class.

### 13.3 Tie-break rule when more than one class contains the target mass

If more than one class strictly contains the target mass, the code selects the **tightest** admissible class, i.e. the one with the minimum nominal bandwidth:

$$
\displaystyle \Delta M_{range} = NUL - NLL
$$

This is an important design choice. It biases the result toward the narrowest class that still contains the target, rather than toward the broadest or the first class found. In practical terms, it prevents needlessly wide gradings when a tighter class is already available.

### 13.4 Representative and extreme masses reported for standard classes

For an adopted standard class, the code reports:

$$
\displaystyle M_{50} = \frac{NLL + NUL}{2}
$$

$$
\displaystyle ELL = 0.7\,NLL
$$

$$
\displaystyle EUL = 1.5\,NUL
$$

The `ELL` and `EUL` values are therefore implementation-derived extreme bounds based on the adopted nominal limits.

---

## 14. Custom underlayer grading logic

### 14.1 When custom grading is used

A **custom grading** is used in either of the following situations:

1. `UseEN13383 = false`; or
2. `UseEN13383 = true`, but **no standard class strictly contains** `M_target`.

In that case, the code derives a custom grading inside one family: `HMA`, `LMA`, or `CP`. If the user selects `AUTO`, the family is chosen by the software.

### 14.2 Family detection in `AUTO` mode

In `AUTO` mode, the code evaluates proximity in **logarithmic mass space**. For each standard class in the embedded database, it computes a representative mass:

$$
\displaystyle M_{50,class} = \frac{NLL + NUL}{2}
$$

and a score proportional to:

$$
\displaystyle \left|\ln(M_{target}) - \ln(M_{50,class})\right|
$$

A mild additional penalty is applied to `CP` when the distances are very similar, which makes the automatic logic slightly prefer the mass-based `LMA` and `HMA` families over the CP family in ambiguous cases. This is a deliberate implementation bias, not a standard requirement.

### 14.3 Sorting within a family

Once a family is selected, the family’s classes are sorted by increasing representative mass:

$$
\displaystyle M_{50,class} = \frac{NLL + NUL}{2}
$$

This ordered sequence becomes the interpolation backbone.

### 14.4 Interpolated grading-ratio variable

The code does not interpolate `NLL` and `NUL` directly. Instead, it interpolates the **ratio**:

$$
\displaystyle R = \frac{NUL}{NLL}
$$

for the chosen family.

For each standard class in the family, the software stores the pair:

$$
\displaystyle x = \ln(M_{50,class})
$$

$$
\displaystyle y = \ln\left(\frac{NUL}{NLL}\right)
$$

and then interpolates linearly in `(x,y)` space. This means the ratio is interpolated in **log space**, which is technically appropriate when the underlying variable spans more than one order of magnitude and behaves multiplicatively rather than additively.

### 14.5 Lower-end and upper-end clamping

If the target mass lies below the smallest representative mass in the selected family, the ratio is **clamped** to the lower-end family class ratio. If it lies above the largest representative mass, it is clamped to the upper-end class ratio. This prevents uncontrolled extrapolation outside the known family range.

### 14.6 Construction of the custom nominal limits

Once the interpolated or clamped ratio `R` is known, the custom grading is built **around the target mass itself**. The code enforces:

$$
\displaystyle M_{50,custom} = M_{target}
$$

and solves for nominal limits through:

$$
\displaystyle NLL_{custom} = \frac{2M_{target}}{1+R}
$$

$$
\displaystyle NUL_{custom} = R\,NLL_{custom}
$$

This is a particularly important feature of the implementation. The custom grading is not simply “the nearest standard class.” It is a genuinely new grading centered on the target representative mass while inheriting a family-consistent nominal-width ratio.

### 14.7 Reported metadata for custom grading

When a custom grading is used, the code reports:

- `Adopted rock grading = Custom Grading`;
- `Custom family basis`;
- `Custom ratio R = NUL/NLL`.

The interpolation note string is retained internally in some code paths, but the aligned public report only exposes the family and the ratio.

### 14.8 Why the custom method is technically coherent

The custom method is analytically coherent for preliminary design because it preserves two desirable properties simultaneously:

1. the grading is centered exactly on the hydraulic underlayer target implied by the armor-unit weight;
2. the grading width remains consistent with the statistical bandwidth of the selected family.

This is more rational than simply forcing the target into a distant standard class when the standard database has gaps or when the target sits near class boundaries.

---

## 15. Underlayer representative weight, nominal diameter, thickness, and packing density

### 15.1 Representative underlayer weight

After a standard or custom grading is adopted, the code evaluates the representative rock weight:

$$
\displaystyle W_{mean,rock} = \frac{M_{50}g}{1000}
$$

### 15.2 Important implementation distinction between standard and custom mode

A subtle but important current-code detail is that the reported underlayer nominal diameter is **not** calculated the same way in standard and custom mode.

In **standard EN 13383 mode**, the code uses the representative weight derived from the adopted class:

$$
\displaystyle D_{n,rock} = \left(\frac{W_{mean,rock}}{W_{rock}}\right)^{1/3}
$$

In **custom grading mode**, however, the code deliberately uses the **theoretical target weight** rather than the representative mean weight:

$$
\displaystyle D_{n,rock} = \left(\frac{W_{target}}{W_{rock}}\right)^{1/3}
$$

This is not an accident. It ensures that a custom underlayer is dimensioned exactly from the hydraulic target, whereas a standard class is dimensioned from the representative class mass.

### 15.3 Underlayer thickness and packing density

The underlayer thickness is:

$$
\displaystyle r_2 = 2D_{n,rock}
$$

The packing density reported by the code is:

$$
\displaystyle f_2 = \frac{100\times2\times1\times(1-P_{rock})}{D_{n,rock}^2}
$$

With `P_rock = 0.25`, this becomes:

$$
\displaystyle f_2 = \frac{150}{D_{n,rock}^2}
$$

Again, this is a code-level reporting formula for a notional double layer. It should not be treated as a universal placement law independent of construction method.

---

## 16. Interpretation of the implemented outputs

### 16.1 Stability number `Ns`

For the trunk, `Ns` is the direct output of the selected empirical formula. For the head, `Ns` is a derived quantity implied by the transferred density at constant `Dn`. The head `Ns` therefore should be interpreted as a **transfer result**, not as a separately fitted head stability law.

### 16.2 Equivalent `KD`

The reported `KD_trunk` and `KD_head` are equivalent Hudson coefficients derived from the solved design. They are useful because they place the result in the language of SPM, CEM, and the Rock Manual, but they are not the primary sizing law used for the trunk armor.

### 16.3 Armor-unit mass and weight

The repository reports both **weight** in `kN` and **mass** in `t`. This is helpful because literature and practice vary between specific weights, densities, masses, and forces. The code itself is internally consistent provided the user keeps `Wc` and `Ww` as **specific weights** in `kN/m³`.

### 16.4 Underlayer grading fields

The underlayer report should be read in the following order:

1. theoretical hydraulic target (`W/10`);
2. adopted class or custom grading;
3. representative mass `M50`;
4. nominal band `NLL–NUL`;
5. extreme band `ELL–EUL`;
6. equivalent nominal diameter `Dn_rock`;
7. double-layer thickness `r2`;
8. packing density `f2`.

That sequence mirrors the actual internal logic of the code.

---

## 17. Limits of applicability and engineering cautions

### 17.1 Preliminary hydraulic sizing only

This repository is a **preliminary hydraulic sizing tool**. It does not perform:

- overtopping analysis;
- wave transformation from offshore to toe;
- diffraction or harbor agitation analysis;
- toe stability design;
- filter design from granular criteria;
- geotechnical verification;
- settlement or liquefaction assessment;
- crown-wall or parapet interaction assessment;
- structural-strength verification of large concrete units.

Those limitations are entirely consistent with the attached manuals and technical literature, which distinguish between conceptual hydraulic sizing and full design verification (Van der Meer, 1995; Van der Meer, 1998; USACE, 1995; USACE, 2006).

### 17.2 Physical-model recommendation

The attached sources are explicit that important rubble-mound structures should be validated in physical models. This is especially important for:

- breakwater heads;
- low-crested structures;
- overtopped sections;
- unusual foreshore breaking;
- structures with strong three-dimensional effects;
- very large concrete armor units where structural breakage may become relevant.

The repository should therefore be understood as an engineering calculator for **rapid, disciplined preliminary design**, not as the last design step.

### 17.3 Interpretation of the underlayer rule

The adopted underlayer rule:

$$
\displaystyle W_{target} = \frac{W_{armor}}{10}
$$

is practical and coherent for the current codebase, but it does not replace a full project-specific filter and underlayer review. In severe exposure, unusual cross-sections, or highly optimized designs, the underlayer may need its own dedicated hydraulic and geometric verification.

### 17.4 Reliability and sensitivity

The Van der Meer material emphasizes that empirical stability formulae possess scatter and should be used with sensitivity analysis. The same recommendation applies here. In practice, the user should vary at least:

- `Hs`;
- `Tm`;
- `Nod`;
- `Nz`;
- `Wc`;
- selected formula ID;
- EN versus custom underlayer mode.

That sensitivity view is usually more informative than treating one single output as an exact design truth.

---

## 18. Inputs, defaults, and constants used by the aligned core calculators

### 18.1 Exposed engineering inputs

| Symbol | Meaning | Unit |
| --- | --- | --- |
| `Hs` | significant wave height | m |
| `Tm` | mean wave period | s |
| `Nz` | number of waves in the design storm | – |
| `Nod` | damage number | – |
| `Wc` | trunk concrete specific weight | kN/m³ |
| `Formula_ID` | armor formula selector | – |
| `UseEN13383` | use standard EN 13383 underlayer grading selection | true/false |
| `CustomFamily` | custom underlayer family | `AUTO`, `HMA`, `LMA`, `CP` |

### 18.2 Internal water specific weight

| Symbol | Meaning | Unit | Default |
| --- | --- | --- | ---: |
| `Ww` | water specific weight | kN/m³ | 10.05 |

### 18.3 Default values used by the aligned core calculators

| Parameter | Default |
| --- | ---: |
| `Hs` | 11.0 m |
| `Tm` | 11.9 s |
| `Nz` | 3000 |
| `Nod` | 0.5 |
| `Wc` | 27.48 kN/m³ |
| `Ww` | 10.05 kN/m³ |
| `Formula_ID` | 1 |
| `UseEN13383` | `true` |
| `CustomFamily` | `AUTO` |

### 18.4 Internal constants currently used by the code

| Constant | Value |
| --- | ---: |
| `g` | 9.80665 m/s² |
| `W_rock_spec` | 26.5 kN/m³ |
| `P_cubes` | 0.40 |
| `P_rock` | 0.25 |
| `KD_RATIO_FIXED` | 1.5 |

---

## 19. Outputs produced

### 19.1 Trunk outputs

For the trunk the report contains:

- stability number `Ns`;
- nominal dimension `Dn`;
- block weight `W` and mass;
- equivalent `KD`;
- double-layer thickness `r1`;
- packing density `d`;
- volume and, for Antifer, the reported geometric dimensions `H`, `A`, and `B`.

### 19.2 Underlayer outputs

For each underlayer the report contains:

- theoretical `W/10` target in `kN` and `kg`;
- adopted grading or `Custom Grading`;
- representative `M50`;
- `NLL`, `NUL`, `ELL`, `EUL`;
- nominal diameter `Dn_rock`;
- double-layer thickness `r2`;
- packing density `f2`;
- and, when applicable, the custom family and ratio `R = NUL/NLL`.

### 19.3 Head outputs

For the head the report contains:

- transfer ratio `Kd_T/Kd_H`;
- nominal dimension `Dn`;
- equivalent `KD_head`;
- required head concrete specific weight `Wc_head`;
- head block weight and mass;
- packing density;
- and the head underlayer results.

---

## 20. Implementation-specific notes

### 20.1 `breakwater_calculator_cli.cpp`

This is the clearest single-file reference for:

- the formula database;
- the exact coefficient values;
- the EN 13383 and custom-grading logic;
- the report field order;
- the current append-to-file behavior.

It accepts:

```text
[Hs] [Tm] [NumberOfWaves] [Nod] [Wc] [FormulaID] [UseEN13383] [CustomFamily]
```

and appends the generated report to `output.txt`.

### 20.2 `breakwater_calculator.py`

The Python implementation mirrors the same input order, the same formula sets, the same underlayer logic, and the same report semantics as the aligned C++ CLI.

### 20.3 `breakwater_calculator_cli.f90`

The Fortran CLI version was aligned to the same reference logic. For practical purposes, the C++ CLI remains the authoritative behavioral description when there is any doubt about wording or output formatting.

### 20.4 `breakwater_calculator_gui.cpp`

The GUI uses the same hydraulic kernel but exposes the options through graphical controls. It includes the custom underlayer grading option only for the rock underlayers, leaving the armor hydraulic formulae unchanged.

### 20.5 `breakwater_calculator.xlsx`

The spreadsheet workbook reproduces the same engineering chain in worksheet form and includes controls for formula selection and underlayer-grading mode. It assembles the report inside the workbook instead of writing through executable code.

### 20.6 `breakwater_calculator.ipynb`

The notebook exposes the same main logic in cell-based form. It is suitable for transparent checking and exploratory studies, while the CLI remains the most explicit definition of the current production-style report.

### 20.7 Nomogram scripts

The nomogram scripts are companion tools:

- `nomograms_vandermeer-chegini.py` for the trunk stability equations;
- `nomogram_hudson_trunk_to_head.py` for the trunk-to-head transfer rule.

They are support tools rather than the primary reference for the full report logic.

---

## 21. Build and execution instructions

### 21.1 C++ command-line version

#### 21.1.1 Recommended build

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

#### 21.1.2 Run interactively

```bash
./breakwater_calculator_cli
```

#### 21.1.3 Run with explicit arguments

```bash
./breakwater_calculator_cli \
  11.0 \
  11.9 \
  3000 \
  0.5 \
  27.48 \
  1 \
  true \
  AUTO
```

Examples:

```bash
# Standard EN 13383 underlayer selection
./breakwater_calculator_cli 11.0 11.9 3000 0.5 27.48 1 true AUTO

# Force custom interpolated underlayer grading with automatic family choice
./breakwater_calculator_cli 11.0 11.9 3000 0.5 27.48 1 false AUTO

# Force custom interpolated underlayer grading within the CP family
./breakwater_calculator_cli 11.0 11.9 3000 0.5 27.48 1 false CP
```

### 21.2 Python command-line version

Run interactively:

```bash
python breakwater_calculator.py
```

Run with explicit arguments:

```bash
python breakwater_calculator.py \
  11.0 \
  11.9 \
  3000 \
  0.5 \
  27.48 \
  1 \
  true \
  AUTO
```

### 21.3 Fortran command-line version

#### 21.3.1 Recommended GNU Fortran build

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

#### 21.3.2 Run interactively

```bash
./breakwater_calculator_cli
```

#### 21.3.3 Run with explicit arguments

```bash
./breakwater_calculator_cli \
  11.0 \
  11.9 \
  3000 \
  0.5 \
  27.48 \
  1 \
  true \
  AUTO
```

### 21.4 GUI version

#### 21.4.1 Recommended MinGW-w64 build on Windows

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

#### 21.4.2 Running the GUI

1. Start `breakwater_calculator_gui.exe`.
2. Select the desired interface language.
3. Enter `Hs`, `Tm`, `Nz`, `Nod`, and `Wc`.
4. Select the armor formula.
5. Select whether to use standard EN 13383 underlayer grading.
6. Select the custom underlayer family if custom mode is desired.
7. Trigger the calculation.
8. Review the report in the interface and in `output.txt`.

### 21.5 If static linking fails

Some systems do not provide all static runtime libraries by default. If a fully static build fails, retry the same commands without the static-linking flags.

### 21.6 Notes on Windows execution

If compiled under `cmd.exe`, executables are typically run as:

```bat
breakwater_calculator_cli.exe
breakwater_calculator_gui.exe
```

If compiled under PowerShell, the following forms may be required:

```powershell
.\breakwater_calculator_cli.exe
.\breakwater_calculator_gui.exe
```

---

## 22. Output-file behavior

### 22.1 Append behavior

The aligned executable implementations append each new report to `output.txt`:

- `breakwater_calculator_cli.cpp`
- `breakwater_calculator.py`
- `breakwater_calculator_cli.f90`
- `breakwater_calculator_gui.cpp`

### 22.2 Workbook behavior

The Excel workbook assembles the report within the workbook instead of writing through executable code to `output.txt`.

### 22.3 Notebook behavior

Notebook behavior depends on the cells executed. For automated command-line semantics and report behavior, the CLI implementations remain the authoritative reference.

---

## 23. Compact model summary

The implemented engineering chain is:

$$
\displaystyle
(H_s, T_m, N_z, N_{od}, W_{c,trunk}, \text{Formula ID}, \text{UseEN13383}, \text{CustomFamily})
$$

$$
\displaystyle
\rightarrow L_0 \rightarrow s_{0m} \rightarrow N_{s,trunk} \rightarrow D_n \rightarrow W_{trunk} \rightarrow K_{D,trunk}
$$

$$
\displaystyle
\rightarrow \Delta_{head} \rightarrow W_{c,head} \rightarrow W_{head} \rightarrow K_{D,head}
$$

$$
\displaystyle
\rightarrow W_{target} = W_{armor}/10 \rightarrow M_{target}
$$

$$
\displaystyle
\rightarrow
\begin{cases}
\text{strict EN 13383 class selection,} & \text{if } UseEN13383 = true \text{ and } NLL < M_{target} < NUL \\
\text{custom family-based interpolated grading,} & \text{otherwise}
\end{cases}
$$

$$
\displaystyle
\rightarrow M_{50} \rightarrow D_{n,rock} \rightarrow r_2 \rightarrow f_2
$$

---

## 24. Final technical interpretation

This repository should be read as a **disciplined empirical design manual encoded in software**. Its logic is technically coherent because it respects the main empirical structure of modern rubble-mound armor design:

- wave attack is expressed through `Hs`, `Tm`, and storm length;
- stability is expressed through the dimensionless quantity `Hs/(ΔDn)`;
- damage is an explicit design variable through `Nod`;
- the trunk is sized by unit-specific empirical formulae rather than by a single universal `KD`;
- Hudson is retained as a transferable and interpretable comparative framework;
- the head is intentionally handled by an iso-geometric density transfer;
- the underlayer is treated through grading-band logic rather than as a second armor layer of identical type.

That combination makes the repository well suited to **concept design, rapid option screening, sensitivity studies, and report generation**. It does not eliminate the need for engineering judgment; it organizes that judgment in a reproducible and transparent way, which is precisely the role empirical design formulae are meant to play in coastal engineering practice (Van der Meer, 1995; Van der Meer, 1998; Van der Meer, 1999; Chegini and Aghtouman, 2006; CERC, 1984; USACE, 2006; CIRIA et al., 2007).

---

## 25. Bibliography

1. **U.S. Army Corps of Engineers, Coastal Engineering Research Center. (1984).** *Shore Protection Manual* (4th ed., Vols. I–II). U.S. Government Printing Office, Washington, DC. ([Contentdm][1])
   **Relevance:** Classical USACE design manual for rubble-mound coastal structures; foundational source for Hudson-based armour stability practice, concrete-unit stability tables, and general breakwater and revetment design guidance.
   **URL:** [Official USACE record](https://usace.contentdm.oclc.org/digital/collection/p16021coll11/id/1934/)

2. **Burcharth, H. F., & Hughes, S. A. (2003).** *Fundamentals of Design.* In *Coastal Engineering Manual* (Vol. 6, Chapter 5, pp. VI-5-i–VI-5-316). Coastal Engineering Research Center, Vicksburg, Mississippi. ([Aalborg Universitets forskningsportal][2])
   **Relevance:** Core CEM chapter on design philosophy and hydraulic response of coastal structures; directly relevant to armour stability, overtopping, wave loads, scour, and rubble-mound design methodology.
   **URL:** [Aalborg University repository entry](https://vbn.aau.dk/en/publications/fundamentals-of-design/)

3. **U.S. Army Corps of Engineers. (2002; updated through Change 3, 2011).** *Coastal Engineering Manual, Part VI: Design of Coastal Project Elements.* Washington, DC. ([Publicações do USACE][3])
   **Relevance:** Official USACE design manual consolidating coastal-structure design methods, including armour-layer stability, overtopping, crest design, toe design, and example problems used in engineering practice.
   **URL:** [USACE Part VI PDF](https://www.publications.usace.army.mil/Portals/76/Publications/EngineerManuals/EM_1110-2-1100_Part-06.pdf)

4. **U.S. Army Corps of Engineers. (1995).** *Design of Coastal Revetments, Seawalls, and Bulkheads* (EM 1110-2-1614). Washington, DC. ([Publicações do USACE][4])
   **Relevance:** Detailed USACE manual on coastal armouring systems and slope protection; useful for toe stability, underlayers, crest geometry, wave loading context, and practical cross-section detailing.
   **URL:** [Official PDF](https://www.publications.usace.army.mil/Portals/76/Publications/EngineerManuals/EM_1110-2-1614.pdf)

5. **CIRIA; CUR; CETMEF. (2007).** *The Rock Manual: The Use of Rock in Hydraulic Engineering* (2nd ed., C683). CIRIA, London. ([CIRIA][5])
   **Relevance:** Principal modern reference for rock in hydraulic engineering, including armourstone properties, grading specification, hydraulic stability, toe and filter design, damage description, and construction quality control.
   **URL:** [Official CIRIA publication page](https://www.ciria.org/ItemDetail?Category=BOOK&iProductCode=C683)

6. **Hudson, R. Y. (1958).** *Design of Quarry-Stone Cover Layers for Rubble-Mound Breakwaters: Hydraulic Laboratory Investigation* (Research Report No. 2-2). U.S. Army Engineer Waterways Experiment Station, Vicksburg, MS. ([WorldCat][6])
   **Relevance:** Original research report underlying the classic Hudson design approach for armour sizing of rubble-mound breakwaters.
   **URL:** [WorldCat bibliographic record](https://search.worldcat.org/title/design-of-quarry-stone-cover-layers-for-rubble-mound-breakwaters-hydraulic-laboratory-investigation/oclc/4578948)

7. **Hudson, R. Y. (1959).** *Laboratory Investigation of Rubble-Mound Breakwaters.* *Journal of the Waterways and Harbors Division*, 85(WW3), 93–121. ([ASCE Library][7])
   **Relevance:** Seminal paper presenting the laboratory basis of the Hudson formula and its application to rubble-mound armour stability.
   **URL:** [ASCE Library DOI record](https://ascelibrary.org/doi/10.1061/JWHEAU.0000142)

8. **Van der Meer, J. W. (1987).** *Stability of breakwater armour layers—Design formulae.* *Coastal Engineering*, 11(3), 219–239. ([Van der Meer Consulting][8])
   **Relevance:** Landmark paper introducing the Van der Meer stability formulae for rock armour under random-wave attack, including dependence on wave period, storm duration, permeability, and damage level.
   **URL:** [Author-hosted PDF](https://www.vandermeerconsulting.nl/downloads/stability_b/1987_vandermeer_armourlayers.pdf)

9. **Van der Meer, J. W. (1988).** *Stability of Cubes, Tetrapods and Accropode.* In *Design of Breakwaters* (pp. 71–80). Thomas Telford, London. ([Van der Meer Consulting][9])
   **Relevance:** Primary source for Van der Meer’s concrete armour-unit stability formulae for cubes and other units under irregular waves, including damage-number-based formulations.
   **URL:** [Author-hosted PDF](https://www.vandermeerconsulting.nl/downloads/stability_c/1988_vandermeer.pdf)

10. **Van der Meer, J. W. (1995).** *Conceptual Design of Rubble Mound Breakwaters.* In *Advances in Coastal and Ocean Engineering* (Vol. 1, pp. 221–315). World Scientific. ([World Scientific][10])
    **Relevance:** Broad conceptual-design reference covering hydraulic stability, toe stability, crest design, overtopping, reliability, and practical breakwater pre-design methodology.
    **URL:** [Author-hosted PDF](https://www.vandermeerconsulting.nl/downloads/stability_b/1995_vandermeer_conceptualdesign.pdf)

11. **Van der Meer, J. W. (1998).** *Application and stability criteria for rock and artificial units.* In K. W. Pilarczyk (ed.), *Dikes and Revetments: Design, Maintenance and Safety Assessment* (Chapter 11, pp. 191–215). A.A. Balkema, Rotterdam/Brookfield. ([Leibniz Universität Hannover][11])
    **Relevance:** Highly relevant synthesis chapter linking conceptual design to practical application criteria for both rock armour and concrete units, including Hudson and Van der Meer usage domains.
    **URL:** [Author-hosted PDF](https://www.vandermeerconsulting.nl/downloads/stability_b/1998_vandermeer_ch11.pdf)

12. **Van der Meer, J. W. (1999).** *Design of concrete armour layers.* In I. J. Losada (ed.), *Proceedings of Coastal Structures ’99* (pp. 213–221). Balkema, Rotterdam. ([UNESCO-IHE][12])
    **Relevance:** Comparative discussion of concrete armour-unit systems, including two-layer and single-layer concepts, stability behaviour, relative efficiency, and practical design implications.
    **URL:** [Author-hosted PDF](https://www.vandermeerconsulting.nl/downloads/stability_c/1999_vandermeer.pdf)

13. **Chegini, V., & Aghtouman, P. (2006).** *An investigation on the stability of rubble mound breakwaters with armour layers of Antifer cubes.* *Journal of Marine Engineering*, 2(1), 86–93. ([ResearchGate][13])
    **Relevance:** Principal reference for Antifer-cube stability relationships and physical-model evidence used in empirical design formulae for Antifer-armoured breakwaters.
    **URL:** [Journal PDF](https://marine-eng.ir/article-1-23-fa.pdf)

14. **Frens, A. B. (2007).** *The impact of placement method on Antifer-block stability.* M.Sc. thesis, Delft University of Technology. ([Core][14])
    **Relevance:** Detailed study on how Antifer placement pattern and packing density influence hydraulic stability; important for interpreting the limits of formula-only design and the role of placement quality.
    **URL:** [CORE PDF](https://files01.core.ac.uk/download/pdf/46673898.pdf)

15. **Frens, A. B., van Gent, M. R. A., & Olthof, J. (2009).** *Placement methods for Antifer armour units.* In *Coastal Engineering 2008* (pp. 3337–3345). World Scientific. ([World Scientific][15])
    **Relevance:** Conference paper extending the placement-method discussion for Antifer units and documenting the influence of orderly versus irregular placement on armour stability.
    **URL:** [World Scientific DOI page](https://www.worldscientific.com/doi/10.1142/9789814277426_0276)

16. **Latham, J.-P., Van Meulen, J. A., & Dupray, S. (2006).** *The specification of armourstone gradings and EN 13383 (2002).* *Quarterly Journal of Engineering Geology and Hydrogeology*, 39(1), 51–64. ([Lyell Collection][16])
    **Relevance:** Key paper on armourstone grading theory, statistical interpretation of grading envelopes, the (M_{50})–mean-mass relationship, and non-standard grading specification; directly relevant to custom underlayer grading logic.
    **URL:** [DOI page](https://www.lyellcollection.org/doi/10.1144/1470-9236/05-025)

17. **CEN. (2013).** *EN 13383-1:2013 Armourstone — Part 1: Specification.* European Committee for Standardization, Brussels. ([iTeh Standards][17])
    **Relevance:** Governing European standard for armourstone specification and standard gradings; directly relevant to EN-based underlayer selection and verification of quarry gradings.
    **URL:** [Preview PDF](https://webstore.ansi.org/preview-pages/bsi/preview_30195709.pdf)