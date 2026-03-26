import math
import sys
from dataclasses import dataclass, field
from typing import Dict, List


@dataclass
class FormulaParams:
    name: str
    type: str
    slope_ratio: float
    k1: float
    k2: float
    k3: float
    k4: float
    k5: float


@dataclass
class GradingDef:
    name: str
    NLL_kg: float
    NUL_kg: float


@dataclass
class Dimensions:
    H: float = 0.0
    A: float = 0.0
    B: float = 0.0


@dataclass
class UnderlayerResult:
    grading_name: str = ""
    target_W: float = 0.0
    target_M50_kg: float = 0.0
    M50_kg: float = 0.0
    ELL_kg: float = 0.0
    EUL_kg: float = 0.0
    NLL_kg: float = 0.0
    NUL_kg: float = 0.0
    W_mean_kn: float = 0.0
    Dn_rock: float = 0.0
    r2: float = 0.0
    f2: float = 0.0
    W_rock_spec: float = 0.0
    used_custom_interpolation: bool = False
    custom_family: str = ""
    custom_ratio_nul_nll: float = 0.0
    custom_ratio_note: str = ""


@dataclass
class ArmorResult:
    Ns: float = 0.0
    Dn: float = 0.0
    W: float = 0.0
    Mass_tonnes: float = 0.0
    Kd: float = 0.0
    r1: float = 0.0
    packing_density: float = 0.0
    dims: Dimensions = field(default_factory=Dimensions)
    Kd_Ratio: float = 0.0
    Delta_Required: float = 0.0
    Wc_Required: float = 0.0


@dataclass
class Intermediates:
    L0: float = 0.0
    k0: float = 0.0
    s0m: float = 0.0
    sm: float = 0.0
    Nz: float = 0.0
    Storm_Duration_hr: float = 0.0
    delta: float = 0.0
    Ns_trunk: float = 0.0


@dataclass
class Inputs:
    Hs: float
    Tm: float
    Number_of_Waves: float
    Nod: float
    Wc: float
    Ww: float
    Formula_ID: int
    use_en13383: bool
    custom_family: str


@dataclass
class FullResults:
    inputs: Inputs
    coefficients: FormulaParams
    intermediate: Intermediates
    final_trunk: ArmorResult
    underlayer_trunk: UnderlayerResult
    final_head: ArmorResult
    underlayer_head: UnderlayerResult
    P_rock: float
    P_cubes: float


class BreakwaterCalculator:
    def __init__(self) -> None:
        self.g = 9.80665
        self.W_rock_spec = 26.5
        self.P_cubes = 0.40
        self.P_rock = 0.25
        self.KD_RATIO_FIXED = 1.5

        self.formulas: Dict[int, FormulaParams] = {
            1: FormulaParams("Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)", "Cubes", 2.0, 7.374304189198, 0.4, 0.3, 1.100642416298, 0.1),
            2: FormulaParams("Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)", "Cubes", 1.5, 6.7, 0.4, 0.3, 1.0, 0.1),
            3: FormulaParams("Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)", "Antifer", 2.0, 6.138, 0.443, 0.276, 1.164, 0.07),
            4: FormulaParams("Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)", "Antifer", 1.5, 6.951, 0.443, 0.291, 1.082, 0.082),
        }

        self.standard_gradings: List[GradingDef] = [
            GradingDef("CP32/90", 0.868, 19.319),
            GradingDef("CP45/125", 2.415, 51.758),
            GradingDef("CP63/180", 6.626, 154.548),
            GradingDef("CP90/250", 19.319, 414.063),
            GradingDef("CP45/180", 2.415, 154.548),
            GradingDef("CP90/180", 19.319, 154.548),
            GradingDef("LMA5-40", 5.0, 40.0),
            GradingDef("LMA10-60", 10.0, 60.0),
            GradingDef("LMA15-120", 15.0, 120.0),
            GradingDef("LMA40-200", 40.0, 200.0),
            GradingDef("LMA60-300", 60.0, 300.0),
            GradingDef("LMA15-300", 15.0, 300.0),
            GradingDef("HMA300-1000", 300.0, 1000.0),
            GradingDef("HMA1000-3000", 1000.0, 3000.0),
            GradingDef("HMA3000-6000", 3000.0, 6000.0),
            GradingDef("HMA6000-10000", 6000.0, 10000.0),
            GradingDef("HMA10000-15000", 10000.0, 15000.0),
        ]

        self.defaults = Inputs(
            11.0,
            11.9,
            3000.0,
            0.5,
            27.48,
            10.05,
            1,
            True,
            "AUTO",
        )

    def calculate_L0(self, Tm: float) -> float:
        return (self.g * (Tm ** 2)) / (2.0 * math.pi)

    def get_formula_display_name(self, formula_id: int) -> str:
        if formula_id == 1:
            return "Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)"
        if formula_id == 2:
            return "Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)"
        if formula_id == 3:
            return "Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)"
        if formula_id == 4:
            return "Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)"
        return "Unknown formula"

    @staticmethod
    def starts_with(value: str, prefix: str) -> bool:
        return value.startswith(prefix)

    def grading_family(self, grading_name: str) -> str:
        if self.starts_with(grading_name, "HMA"):
            return "HMA"
        if self.starts_with(grading_name, "LMA"):
            return "LMA"
        if self.starts_with(grading_name, "CP"):
            return "CP"
        return "UNKNOWN"

    def get_family_gradings(self, family: str) -> List[GradingDef]:
        out = [gdef for gdef in self.standard_gradings if self.grading_family(gdef.name) == family]
        out.sort(key=lambda gdef: 0.5 * (gdef.NLL_kg + gdef.NUL_kg))
        return out

    def select_custom_family(self, target_mass: float) -> str:
        safe_mass = max(target_mass, 1e-9)
        best_score = float("inf")
        best_family = "LMA"

        for gdef in self.standard_gradings:
            fam = self.grading_family(gdef.name)
            if fam == "UNKNOWN":
                continue

            m50 = 0.5 * (gdef.NLL_kg + gdef.NUL_kg)
            score = abs(math.log(safe_mass) - math.log(max(m50, 1e-9)))
            if fam == "CP":
                score += 0.08

            if score < best_score:
                best_score = score
                best_family = fam

        return best_family

    def interpolate_family_ratio(self, target_mass: float, family: str) -> tuple[float, str]:
        family_gradings = self.get_family_gradings(family)
        if not family_gradings:
            return 3.0, "fallback ratio R=NUL/NLL = 3.0 (no family data)"

        x: List[float] = []
        y: List[float] = []
        for gdef in family_gradings:
            m50 = 0.5 * (gdef.NLL_kg + gdef.NUL_kg)
            x.append(math.log(max(m50, 1e-9)))
            y.append(math.log(max(gdef.NUL_kg / gdef.NLL_kg, 1.0 + 1e-9)))

        xt = math.log(max(target_mass, 1e-9))

        if len(family_gradings) == 1:
            ratio = family_gradings[0].NUL_kg / family_gradings[0].NLL_kg
            return ratio, f"{family} family single-point ratio used"

        if xt <= x[0]:
            ratio = math.exp(y[0])
            return ratio, f"{family} family ratio clamped to lower-end class {family_gradings[0].name}"

        if xt >= x[-1]:
            ratio = math.exp(y[-1])
            return ratio, f"{family} family ratio clamped to upper-end class {family_gradings[-1].name}"

        for i in range(len(family_gradings) - 1):
            if xt >= x[i] and xt <= x[i + 1]:
                t = (xt - x[i]) / (x[i + 1] - x[i])
                log_ratio = y[i] + t * (y[i + 1] - y[i])
                ratio = math.exp(log_ratio)
                note = (
                    f"{family} family ratio interpolated between "
                    f"{family_gradings[i].name} and {family_gradings[i + 1].name}"
                )
                return ratio, note

        ratio = math.exp(y[-1])
        return ratio, f"{family} family ratio fallback to upper-end class {family_gradings[-1].name}"

    def calculate_underlayer_params(self, W_armor: float, use_en13383: bool, requested_custom_family: str) -> UnderlayerResult:
        target_weight = W_armor / 10.0
        target_mass_kg = (target_weight * 1000.0) / self.g

        res = UnderlayerResult()
        res.target_W = target_weight
        res.target_M50_kg = target_mass_kg
        res.used_custom_interpolation = False
        res.custom_family = ""
        res.custom_ratio_nul_nll = 0.0
        res.custom_ratio_note = ""

        if use_en13383:
            selected_grading: GradingDef | None = None
            found = False
            min_range_width = float("inf")

            for grading in self.standard_gradings:
                if target_mass_kg > grading.NLL_kg and target_mass_kg < grading.NUL_kg:
                    current_range = grading.NUL_kg - grading.NLL_kg
                    if current_range < min_range_width:
                        min_range_width = current_range
                        selected_grading = grading
                        found = True

            if found and selected_grading is not None:
                res.grading_name = selected_grading.name
                res.M50_kg = 0.5 * (selected_grading.NLL_kg + selected_grading.NUL_kg)
                res.NLL_kg = selected_grading.NLL_kg
                res.NUL_kg = selected_grading.NUL_kg

        if (not use_en13383) or (res.NUL_kg <= res.NLL_kg):
            family = requested_custom_family.upper()
            if family not in {"HMA", "LMA", "CP"}:
                family = self.select_custom_family(target_mass_kg)

            ratio_nul_nll, ratio_note = self.interpolate_family_ratio(target_mass_kg, family)
            ratio_nul_nll = max(ratio_nul_nll, 1.01)

            nll_kg = (2.0 * target_mass_kg) / (1.0 + ratio_nul_nll)
            nul_kg = ratio_nul_nll * nll_kg

            res.grading_name = "Custom Grading"
            res.M50_kg = target_mass_kg
            res.NLL_kg = nll_kg
            res.NUL_kg = nul_kg
            res.used_custom_interpolation = True
            res.custom_family = family
            res.custom_ratio_nul_nll = ratio_nul_nll
            res.custom_ratio_note = ratio_note

        res.ELL_kg = 0.7 * res.NLL_kg
        res.EUL_kg = 1.5 * res.NUL_kg
        res.W_mean_kn = (res.M50_kg * self.g) / 1000.0

        reference_weight_kn = target_weight if res.used_custom_interpolation else res.W_mean_kn
        res.Dn_rock = (reference_weight_kn / self.W_rock_spec) ** (1.0 / 3.0)
        res.r2 = 2.0 * res.Dn_rock
        res.f2 = 100.0 * 2.0 * 1.0 * (1.0 - self.P_rock) / (res.Dn_rock ** 2)
        res.W_rock_spec = self.W_rock_spec
        return res

    def solve(self, formula_id: int, params: Inputs) -> FullResults:
        if formula_id not in self.formulas:
            raise RuntimeError("Invalid Formula ID")

        coeffs = self.formulas[formula_id]
        k1 = coeffs.k1
        k2 = coeffs.k2
        k3 = coeffs.k3
        k4 = coeffs.k4
        k5 = coeffs.k5

        L0 = self.calculate_L0(params.Tm)
        k0 = (2.0 * math.pi) / L0
        s0m = params.Hs / L0
        sm = s0m

        Nz = params.Number_of_Waves
        storm_duration_hr = (Nz * params.Tm) / 3600.0
        delta_trunk = (params.Wc / params.Ww) - 1.0

        term_damage = params.Nod ** k2
        term_waves = Nz ** k3
        damage_wave_ratio = term_damage / term_waves
        scaled_term = k1 * damage_wave_ratio
        inv_f = scaled_term + k4
        steepness_factor = s0m ** (-k5)
        Ns_trunk = inv_f * steepness_factor

        Dn = params.Hs / (delta_trunk * Ns_trunk)
        W_trunk = params.Wc * (Dn ** 3)
        packing_density_trunk = 100.0 * 2.0 * 1.1 * (1.0 - self.P_cubes) / (Dn ** 2)

        slope = coeffs.slope_ratio
        kd_trunk_equiv = (params.Wc * (params.Hs ** 3)) / (W_trunk * (delta_trunk ** 3) * slope)

        ul_trunk = self.calculate_underlayer_params(W_trunk, params.use_en13383, params.custom_family)

        kd_ratio = self.KD_RATIO_FIXED
        kd_head_derived = kd_trunk_equiv / kd_ratio
        delta_head = delta_trunk * (kd_ratio ** (1.0 / 3.0))
        Wc_head = params.Ww * (delta_head + 1.0)
        W_head = W_trunk * (Wc_head / params.Wc)

        Ns_head = params.Hs / (delta_head * Dn)
        packing_density_head = 100.0 * 2.0 * 1.1 * (1.0 - self.P_cubes) / (Dn ** 2)

        ul_head = self.calculate_underlayer_params(W_head, params.use_en13383, params.custom_family)

        r1 = 2.0 * 1.1 * Dn

        vol_trunk = W_trunk / params.Wc
        h_trunk = (vol_trunk / 1.0247) ** (1.0 / 3.0)
        a_trunk = 1.086 * h_trunk
        b_trunk = 1.005 * h_trunk

        vol_head = W_head / Wc_head
        h_head = (vol_head / 1.0247) ** (1.0 / 3.0)
        a_head = 1.086 * h_head
        b_head = 1.005 * h_head

        results = FullResults(
            inputs=params,
            coefficients=coeffs,
            intermediate=Intermediates(),
            final_trunk=ArmorResult(),
            underlayer_trunk=ul_trunk,
            final_head=ArmorResult(),
            underlayer_head=ul_head,
            P_rock=self.P_rock,
            P_cubes=self.P_cubes,
        )

        results.intermediate.L0 = L0
        results.intermediate.k0 = k0
        results.intermediate.s0m = s0m
        results.intermediate.sm = sm
        results.intermediate.Nz = Nz
        results.intermediate.Storm_Duration_hr = storm_duration_hr
        results.intermediate.delta = delta_trunk
        results.intermediate.Ns_trunk = Ns_trunk

        results.final_trunk.Dn = Dn
        results.final_trunk.W = W_trunk
        results.final_trunk.Mass_tonnes = W_trunk / self.g
        results.final_trunk.Kd = kd_trunk_equiv
        results.final_trunk.r1 = r1
        results.final_trunk.packing_density = packing_density_trunk
        results.final_trunk.dims = Dimensions(h_trunk, a_trunk, b_trunk)

        results.final_head.Ns = Ns_head
        results.final_head.Dn = Dn
        results.final_head.Kd = kd_head_derived
        results.final_head.Kd_Ratio = kd_ratio
        results.final_head.Delta_Required = delta_head
        results.final_head.Wc_Required = Wc_head
        results.final_head.W = W_head
        results.final_head.Mass_tonnes = W_head / self.g
        results.final_head.packing_density = packing_density_head
        results.final_head.dims = Dimensions(h_head, a_head, b_head)

        return results

    def format_report(self, results: FullResults) -> str:
        def fmt(val: float, prec: int) -> str:
            return f"{val:.{prec}f}"

        def field(label: str, value: str, width: int = 38) -> str:
            return f"   {label:<{width}} : {value}\n"

        p = results.inputs
        c = results.coefficients
        i = results.intermediate
        ft = results.final_trunk
        ut = results.underlayer_trunk
        fh = results.final_head
        uh = results.underlayer_head

        armor_vol_trunk = 1.0247 * (ft.dims.H ** 3) if c.type == "Antifer" else (ft.Dn ** 3)
        armor_vol_head = 1.0247 * (fh.dims.H ** 3) if c.type == "Antifer" else (fh.Dn ** 3)

        methodology_name = self.get_formula_display_name(p.Formula_ID)

        ss: List[str] = []
        ss.append("================================================================================\n")
        ss.append("    TECHNICAL REPORT: BREAKWATER ARMOR & UNDERLAYER DESIGN\n")
        ss.append("================================================================================\n")
        ss.append(f"Methodology: {methodology_name}\n")
        ss.append("--------------------------------------------------------------------------------\n")

        ss.append("1. INPUT PARAMETERS\n")
        ss.append(field("Hs (Significant Wave Height)", fmt(p.Hs, 2) + " m"))
        ss.append(field("Tm (Mean Wave Period)", fmt(p.Tm, 2) + " s"))
        ss.append(field("Number of waves (Nz)", fmt(p.Number_of_Waves, 0)))
        ss.append(field("Nod (Damage)", fmt(p.Nod, 2)))
        ss.append(field("Wc Trunk (Concrete Spec. Weight)", fmt(p.Wc, 2) + " kN/m3"))
        ss.append(field("Ww (Water Specific Weight)", fmt(p.Ww, 2) + " kN/m3"))
        ss.append(field("Relative Density D=(Wc/Ww)-1", fmt(i.delta, 4)))
        ss.append(field("Structure Slope (TRUNK & HEAD)", fmt(c.slope_ratio, 1) + ":1"))
        ss.append(field("Porosity (Cubes)", fmt(results.P_cubes * 100.0, 0) + "%"))
        ss.append(field("Porosity (Rock Layer)", fmt(results.P_rock * 100.0, 0) + "%"))
        ss.append("--------------------------------------------------------------------------------\n")

        ss.append("2. INTERMEDIATE PARAMETERS\n")
        ss.append(field("Wave Length (L0)", fmt(i.L0, 2) + " m"))
        ss.append(field("Wave number (k0 = 2*pi/L0)", fmt(i.k0, 4)))
        ss.append(field("Wave steepness (s0m = Hs/L0)", fmt(i.s0m, 4)))
        ss.append(field("Storm Duration (h)", fmt(i.Storm_Duration_hr, 3) + " h"))
        ss.append(field("Stability Number TRUNK (Ns)", fmt(i.Ns_trunk, 4)))
        ss.append(field("Stability Number HEAD (Ns)", fmt(fh.Ns, 4)))
        ss.append("--------------------------------------------------------------------------------\n")

        ss.append("3. ARMOR LAYER RESULTS - TRUNK\n")
        ss.append(field("BLOCK WEIGHT (W)", fmt(ft.W, 2) + " kN"))
        ss.append(field("Mass (t)", fmt(ft.Mass_tonnes, 2) + " t"))
        ss.append(field("Nominal Dimension (Dn)", fmt(ft.Dn, 3) + " m"))
        if c.type == "Cubes":
            ss.append(field("Volume = Dn^3 (V)", fmt(armor_vol_trunk, 3) + " m3"))
        else:
            ss.append(field("Volume = 1.0247 * H^3 (V)", fmt(armor_vol_trunk, 3) + " m3"))
            ss.append(field("Block Height (H)", fmt(ft.dims.H, 3) + " m"))
            ss.append(field("Block Top Width (B)", fmt(ft.dims.B, 3) + " m"))
            ss.append(field("Block Base Width (A)", fmt(ft.dims.A, 3) + " m"))
        ss.append(field("KD_TRUNK (Equivalent)", fmt(ft.Kd, 2)))
        ss.append(field("Double Layer Thickness (r1)", fmt(ft.r1, 2) + " m"))
        ss.append(field("Packing Density, d [units/100m2]", fmt(ft.packing_density, 2)))
        ss.append("\n")

        ss.append("4. UNDERLAYER RESULTS - TRUNK\n")
        ss.append(field("Theoretical Target (W/10)", fmt(ut.target_W, 2) + " kN (" + fmt(ut.target_M50_kg, 1) + " kg)"))
        ss.append(field("Adopted rock grading", ut.grading_name))
        ss.append(field("Representative M50", fmt(ut.M50_kg, 1) + " kg"))
        ss.append(field("Nominal lower limit (NLL)", fmt(ut.NLL_kg, 1) + " kg"))
        ss.append(field("Nominal upper limit (NUL)", fmt(ut.NUL_kg, 1) + " kg"))
        ss.append(field("Extreme lower limit (ELL)", fmt(ut.ELL_kg, 1) + " kg"))
        ss.append(field("Extreme upper limit (EUL)", fmt(ut.EUL_kg, 1) + " kg"))
        ss.append(field("Nominal Dimension (Dn_rock)", fmt(ut.Dn_rock, 3) + " m"))
        ss.append(field("Double Layer Thickness (r2)", fmt(ut.r2, 2) + " m"))
        ss.append(field("Packing Density, f2 [rocks/100m2]", fmt(ut.f2, 2)))
        if ut.used_custom_interpolation:
            ss.append(field("Custom family basis", ut.custom_family))
            ss.append(field("Custom ratio R=NUL/NLL", fmt(ut.custom_ratio_nul_nll, 3)))
        ss.append("-" * 80 + "\n")

        ss.append("5. ARMOR LAYER RESULTS - HEAD (High Density)\n")
        ss.append("   *Maintains same Dn and slope as Trunk*\n")
        ss.append(field("Stability Ratio (Kd_T/Kd_H)", fmt(fh.Kd_Ratio, 2)))
        ss.append(field("Nominal Dimension (Dn)", fmt(fh.Dn, 3) + " m"))
        if c.type == "Cubes":
            ss.append(field("Volume = Dn^3 (V)", fmt(armor_vol_head, 3) + " m3"))
        else:
            ss.append(field("Volume = 1.0247 * H^3 (V)", fmt(armor_vol_head, 3) + " m3"))
            ss.append(field("Block Height (H)", fmt(fh.dims.H, 3) + " m"))
            ss.append(field("Block Top Width (B)", fmt(fh.dims.B, 3) + " m"))
            ss.append(field("Block Base Width (A)", fmt(fh.dims.A, 3) + " m"))
        ss.append(field("KD_HEAD (Equivalent)", fmt(fh.Kd, 2)))
        ss.append(field("Required Concrete Density (Wc)", fmt(fh.Wc_Required, 2) + " kN/m3"))
        ss.append(field("BLOCK WEIGHT (W)", fmt(fh.W, 2) + " kN"))
        ss.append(field("Mass (t)", fmt(fh.Mass_tonnes, 2) + " t"))
        ss.append(field("Packing Density, d [units/100m2]", fmt(fh.packing_density, 2)))
        ss.append("\n")

        ss.append("6. UNDERLAYER RESULTS - HEAD\n")
        ss.append(field("Theoretical Target (W/10)", fmt(uh.target_W, 2) + " kN (" + fmt(uh.target_M50_kg, 1) + " kg)"))
        ss.append(field("Adopted rock grading", uh.grading_name))
        ss.append(field("Representative M50", fmt(uh.M50_kg, 1) + " kg"))
        ss.append(field("Nominal lower limit (NLL)", fmt(uh.NLL_kg, 1) + " kg"))
        ss.append(field("Nominal upper limit (NUL)", fmt(uh.NUL_kg, 1) + " kg"))
        ss.append(field("Extreme lower limit (ELL)", fmt(uh.ELL_kg, 1) + " kg"))
        ss.append(field("Extreme upper limit (EUL)", fmt(uh.EUL_kg, 1) + " kg"))
        ss.append(field("Nominal Dimension (Dn_rock)", fmt(uh.Dn_rock, 3) + " m"))
        ss.append(field("Double Layer Thickness (r2)", fmt(uh.r2, 2) + " m"))
        ss.append(field("Packing Density, f2 [rocks/100m2]", fmt(uh.f2, 2)))
        if uh.used_custom_interpolation:
            ss.append(field("Custom family basis", uh.custom_family))
            ss.append(field("Custom ratio R=NUL/NLL", fmt(uh.custom_ratio_nul_nll, 3)))
        ss.append("=" * 80 + "\n")

        return "".join(ss)

    def generate_report_file(self, results: FullResults, filepath: str = "output.txt") -> None:
        report_content = self.format_report(results)

        with open(filepath, "ab") as outfile:
            outfile.write(report_content.encode("utf-8"))

        sys.stdout.write(report_content)


def get_param(prompt: str, default_val: float) -> float:
    sys.stdout.write(f"{prompt} [{default_val}]: ")
    sys.stdout.flush()
    try:
        input_str = sys.stdin.readline()
    except Exception:
        return default_val

    startpos = 0
    while startpos < len(input_str) and input_str[startpos] in " \n\r\t":
        startpos += 1

    if startpos >= len(input_str):
        return default_val

    try:
        return float(input_str[startpos:])
    except Exception:
        return default_val


def trim_copy(s: str) -> str:
    start = 0
    while start < len(s) and s[start] in " \n\r\t":
        start += 1
    if start >= len(s):
        return ""
    end = len(s) - 1
    while end >= 0 and s[end] in " \n\r\t":
        end -= 1
    return s[start:end + 1]


def to_upper_copy(s: str) -> str:
    return s.upper()


def parse_bool(s: str, default_val: bool) -> bool:
    t = to_upper_copy(trim_copy(s))
    if t == "":
        return default_val
    if t in {"1", "TRUE", "T", "YES", "Y"}:
        return True
    if t in {"0", "FALSE", "F", "NO", "N"}:
        return False
    return default_val


def get_text_param(prompt: str, default_val: str) -> str:
    sys.stdout.write(f"{prompt} [{default_val}]: ")
    sys.stdout.flush()
    try:
        input_str = sys.stdin.readline()
    except Exception:
        return default_val

    startpos = 0
    while startpos < len(input_str) and input_str[startpos] in " \n\r\t":
        startpos += 1

    if startpos >= len(input_str):
        return default_val

    endpos = len(input_str) - 1
    while endpos >= 0 and input_str[endpos] in " \n\r\t":
        endpos -= 1
    return input_str[startpos:endpos + 1]


def main() -> int:
    calc = BreakwaterCalculator()
    user_inputs = Inputs(
        calc.defaults.Hs,
        calc.defaults.Tm,
        calc.defaults.Number_of_Waves,
        calc.defaults.Nod,
        calc.defaults.Wc,
        calc.defaults.Ww,
        calc.defaults.Formula_ID,
        calc.defaults.use_en13383,
        calc.defaults.custom_family,
    )
    formula_id = 1

    argc = len(sys.argv)
    argv = sys.argv

    if argc >= 7:
        try:
            user_inputs.Hs = float(argv[1])
            user_inputs.Tm = float(argv[2])
            user_inputs.Number_of_Waves = float(argv[3])
            user_inputs.Nod = float(argv[4])
            user_inputs.Wc = float(argv[5])
            formula_id = int(argv[6])

            if formula_id < 1 or formula_id > 4:
                sys.stderr.write("Error: Formula ID must be 1-4. Using default (1).\n")
                formula_id = 1
            user_inputs.Formula_ID = formula_id

            if argc >= 8:
                user_inputs.use_en13383 = parse_bool(argv[7], calc.defaults.use_en13383)
            if argc >= 9:
                user_inputs.custom_family = to_upper_copy(trim_copy(argv[8]))
                if user_inputs.custom_family not in {"AUTO", "HMA", "LMA", "CP"}:
                    user_inputs.custom_family = "AUTO"
        except Exception as e:
            sys.stderr.write(f"Error parsing command line arguments: {e}\n")
            sys.stderr.write(
                f"Usage: {argv[0]} [Hs] [Tm] [NumberOfWaves] [Nod] [Wc] [FormulaID] [UseEN13383] [CustomFamily]\n"
            )
            return 1
    else:
        sys.stdout.write("\n--- COASTAL PROTECTION BLOCK CALCULATOR (TRUNK & HEAD) ---\n")
        sys.stdout.write("1. Simple Cubes (Slope 2.0:1) - Van der Meer\n")
        sys.stdout.write("2. Simple Cubes (Slope 1.5:1) - Van der Meer\n")
        sys.stdout.write("3. Antifer (Slope 2.0:1) - Chegini\n")
        sys.stdout.write("4. Antifer (Slope 1.5:1) - Chegini\n")

        sys.stdout.write("\nOption [1-4]: ")
        sys.stdout.flush()
        selection = sys.stdin.readline()

        endpos = len(selection) - 1
        while endpos >= 0 and selection[endpos] in " \n\r\t":
            endpos -= 1
        if endpos >= 0:
            selection = selection[:endpos + 1]
        else:
            selection = ""

        if selection in {"1", "2", "3", "4"}:
            formula_id = int(selection)
        else:
            formula_id = 1
        user_inputs.Formula_ID = formula_id

        sys.stdout.write("\n--- Enter Parameters (Press ENTER for Default) ---\n")
        user_inputs.Hs = get_param("Hs (m)", calc.defaults.Hs)
        user_inputs.Tm = get_param("Tm (s)", calc.defaults.Tm)
        user_inputs.Number_of_Waves = get_param("Number of waves (Nz)", calc.defaults.Number_of_Waves)
        user_inputs.Nod = get_param("Nod (Damage)", calc.defaults.Nod)
        user_inputs.Wc = get_param("Concrete Weight Trunk (kN/m3)", calc.defaults.Wc)

        use_en_str = get_text_param(
            "Use standard EN 13383 underlayer grading? [true/false]",
            "true" if calc.defaults.use_en13383 else "false",
        )
        user_inputs.use_en13383 = parse_bool(use_en_str, calc.defaults.use_en13383)

        family = get_text_param("Custom underlayer family [AUTO/HMA/LMA/CP]", calc.defaults.custom_family)
        family = to_upper_copy(trim_copy(family))
        if family in {"AUTO", "HMA", "LMA", "CP"}:
            user_inputs.custom_family = family
        else:
            user_inputs.custom_family = "AUTO"

    try:
        results = calc.solve(formula_id, user_inputs)
        calc.generate_report_file(results, "output.txt")
    except Exception as e:
        sys.stdout.write(f"\n Calculation Error: {e}\n")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
