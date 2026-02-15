// ======================================================================================
// PROGRAM DESCRIPTION & METHODOLOGY
// ======================================================================================
//
// 1. PURPOSE:
//    This software calculates the required size and weight of artificial concrete 
//    armor units (specifically Cubes and Antifer blocks) for rubble mound breakwaters.
//    It implements an "Iso-Geometric Design Philosophy" to dimension two sections:
//    - The Trunk (the main longitudinal section).
//    - The Head (the exposed roundhead).
//
// 2. METHODOLOGY:
//    The calculator utilizes empirical power-law formulas derived from hydraulic 
//    model testing and European standards:
//    - For Simple Cubes: Uses Van der Meer (1988) stability formulas.
//    - For Antifer Blocks: Uses Chegini & Aghtouman (2006) power-law formulas.
//    - For Underlayers: Automatic sizing based on EN 13383 Standard Rock Gradings.
//
// 3. DESIGN LOGIC & STRATEGY:
//    a. Trunk Stability: 
//       Calculates the Stability Number (Ns) based on wave steepness (s0m), storm 
//       duration (Nz), and damage level (Nod). This determines the nominal 
//       diameter (Dn) and Weight (W).
//
//    b. Head Design (Iso-Geometric Transfer Strategy): 
//       The breakwater head endures higher turbulence (3D effects) than the trunk.
//       Rather than increasing block size or softening the slope—which requires 
//       different molds or crane logistics—this tool solves for the required 
//       increase in Concrete Density (rho_c).
//       
//       - Principle: Maintain constant Geometry (Dn) and Slope (alpha).
//       - Transfer Function: Delta_head = Delta_trunk * (1.5)^(1/3)
//       - Result: Head blocks use the same molds but are heavier due to higher density.
//
// 4. TECHNICAL BIBLIOGRAPHY & REFERENCES:
//
// 1. Van der Meer, J.W. (1988). "Rock Slopes and Gravel Beaches Under Wave Attack."
//    Doctoral Thesis, Delft University of Technology.
//
// 2. Van der Meer, J.W. (1988). "Stability of Cubes, Tetrapods and Accropode."
//    Proceedings of the Conference Breakwaters '88, Eastbourne, Thomas Telford.
//
// 3. Chegini, V., & Aghtouman, P. (2006). "An Investigation on the Stability of 
//    Rubble Mound Breakwaters with Armour Layers of Antifer Cubes." 
//    Journal of Marine Engineering.
//
// 4. USACE (2006). "Coastal Engineering Manual (CEM)", Chapter VI-5.
//
// 5. CEN (2002). "EN 13383-1: Armourstone - Part 1: Specification.".
//
// 5. COMPILATION INSTRUCTIONS:
//
// Compilation (MinGW on Windows):
//
// g++ -O3 -march=native -std=c++17 -municode
// breakwater_calculator_gui.cpp -o breakwater_calculator_gui
// -mwindows -static -static-libgcc -static-libstdc++
//
// 6. EXECUTION:
//
// 1. Launch the application (breakwater_calculator_gui.exe).
// 2. Enter hydraulic parameters (Hs, Tm, Duration, Damage Nod, Concrete Density).
// 3. Select the desired formula/block type from the dropdown menu.
// 4. Click "Calculate Armor" to generate the design report in the output window
//    and append it to output.txt.
//
// ======================================================================================

#define _USE_MATH_DEFINES 

#include <windows.h>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <map>
#include <algorithm>
#include <limits>
#include <fstream>
#include <codecvt> 

// --- Data Structures ---

struct FormulaParams {
    std::wstring name; 
    std::string type;
    double slope_ratio;
    double k1, k2, k3, k4, k5;
};

struct GradingDef {
    std::string name;
    double NLL_kg;
    double NUL_kg;
};

struct Dimensions {
    double H;
    double A;
    double B;
};

struct UnderlayerResult {
    std::string grading_name;
    double target_W;
    double target_M50_kg;
    double M50_kg;
    double ELL_kg;
    double EUL_kg;
    double NLL_kg;
    double NUL_kg;
    double W_mean_kn;
    double Dn_rock;
    double r2;
    double f2;
    double W_rock_spec;
};

struct ArmorResult {
    double Ns; 
    double Dn;
    double W;
    double Mass_tonnes;
    double Kd; 
    double r1;
    double packing_density;
    Dimensions dims;
    
    // Specific to Head
    double Kd_Ratio;
    double Delta_Required;
    double Wc_Required;
};

struct Intermediates {
    double L0;
    double k0;
    double s0m;
    double sm;
    double Nz;
    double delta;
    double Ns_trunk;
};

struct Inputs {
    double Hs;
    double Tm;
    double Storm_Duration_hr;
    double Nod;
    double Wc;
    double Ww;
    int Formula_ID;
};

struct FullResults {
    Inputs inputs;
    FormulaParams coefficients;
    Intermediates intermediate;
    ArmorResult final_trunk;
    UnderlayerResult underlayer_trunk;
    ArmorResult final_head;
    UnderlayerResult underlayer_head;
    double P_rock;
    double P_cubes;
};

// --- Core Calculation Logic Class ---

class BreakwaterCalculator {
private:
    const double g = 9.80665;
    const double W_rock_spec = 26.5; 
    
    const double P_cubes = 0.40;  
    const double P_rock = 0.25;   
    
    const double KD_RATIO_FIXED = 1.5;

    std::map<int, FormulaParams> formulas;
    std::vector<GradingDef> standard_gradings;

public:
    Inputs defaults;

    BreakwaterCalculator() {
        // --- FORMULA ORDER ---
        // 1. Cubes 2.0
        // 2. Cubes 1.5
        // 3. Antifer 2.0
        // 4. Antifer 1.5
        
        formulas[1] = { L"Van Der Meer (1988a) - Cubes (Slope 2.0:1)", "Cubes", 2.0, 6.7, 0.4, 0.3, 1.0, 0.1 };
        formulas[2] = { L"Van Der Meer (1988a) - Cubes (Slope 1.5:1)", "Cubes", 1.5, 6.7, 0.4, 0.3, 1.0, 0.1 };
        formulas[3] = { L"Chegini-Aghtouman (2006) - Antifer (Slope 2:1)", "Antifer", 2.0, 6.138, 0.443, 0.276, 1.164, 0.07 };
        formulas[4] = { L"Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)", "Antifer", 1.5, 6.951, 0.443, 0.291, 1.082, 0.082 };

        standard_gradings = {
            // Coarse / Light Gradings
            {"CP45/125",       0.4,   1.2},
            {"CP63/180",       1.2,   3.8},
            {"CP90/250",       3.1,   9.3},
            {"CP45/180",       0.4,   1.2},
            {"CP90/180",       2.1,   2.8},

            // Light Mass Armourstone (LMA) - Values derived from names
            {"LMA5-40",        5,     40},
            {"LMA10-60",       10,    60},
            {"LMA15-120",      15,    120},
            {"LMA40-200",      40,    200},
            {"LMA60-300",      60,    300},
            {"LMA15-300",      15,    300},

            // Heavy Mass Armourstone (HMA) - Values derived from names
            {"HMA300-1000",    300,   1000},
            {"HMA1000-3000",   1000,  3000},
            {"HMA3000-6000",   3000,  6000},
            {"HMA6000-10000",  6000,  10000},
            {"HMA10000-15000", 10000, 15000}
        };

        defaults = {
            10.0,   // Hs
            13.0,   // Tm
            12.0,   // Storm_Duration_hr
            1.0,    // Nod
            24.0,   // Wc
            10.05,  // Ww
            1       // Formula_ID
        };
    }

    double calculate_L0(double Tm) {
        return (g * std::pow(Tm, 2)) / (2 * M_PI);
    }

    UnderlayerResult calculate_underlayer_params(double W_armor) {
        double target_weight = W_armor / 10.0;
        double target_mass_kg = (target_weight * 1000.0) / g;
        
        GradingDef selected_grading;
        bool found = false;
        
        double final_M50 = 0;
        double final_NLL = 0;
        double final_NUL = 0;

        // --- CONTAINMENT & TIGHTEST RANGE LOGIC ---
        double min_range_width = std::numeric_limits<double>::max();

        for (const auto& grading : standard_gradings) {
            // Check containment: Target must be strictly inside nominal limits
            if (target_mass_kg > grading.NLL_kg && target_mass_kg < grading.NUL_kg) {
                
                double current_range = grading.NUL_kg - grading.NLL_kg;
                
                // Update if this is the first match OR if this range is tighter (smaller)
                if (current_range < min_range_width) {
                    min_range_width = current_range;
                    selected_grading = grading;
                    final_NLL = grading.NLL_kg;
                    final_NUL = grading.NUL_kg;
                    final_M50 = 0.5 * (final_NLL + final_NUL);
                    found = true;
                }
            }
        }
        
        // Fallback if no grading strictly contains the target mass
        if (!found && !standard_gradings.empty()) {
            selected_grading = standard_gradings[0];
            final_NLL = selected_grading.NLL_kg;
            final_NUL = selected_grading.NUL_kg;
            final_M50 = 0.5 * (final_NLL + final_NUL);
        }

        // --- CONTINUE WITH EXISTING LIMIT CALCULATIONS ---
        double ELL = 0.7 * final_NLL;
        double EUL = 1.5 * final_NUL;
        double W_mean_kn = (final_M50 * g) / 1000.0;
        
        double Dn_rock = std::pow(W_mean_kn / W_rock_spec, 1.0/3.0);
        double r2 = 2.0 * Dn_rock;
        double f2 = 100.0 * 2.0 * 1.0 * (1.0 - P_rock) / std::pow(Dn_rock, 2);
        
        UnderlayerResult res;
        res.grading_name = selected_grading.name;
        res.target_W = target_weight;
        res.target_M50_kg = target_mass_kg;
        res.M50_kg = final_M50;
        res.NLL_kg = final_NLL;
        res.NUL_kg = final_NUL;
        res.ELL_kg = ELL;
        res.EUL_kg = EUL;
        res.W_mean_kn = W_mean_kn;
        res.Dn_rock = Dn_rock;
        res.r2 = r2;
        res.f2 = f2;
        res.W_rock_spec = W_rock_spec;
        
        return res;
    }
	
    FullResults solve(int formula_id, Inputs params) {
        // 2. Load Coefficients
        if (formulas.find(formula_id) == formulas.end()) {
            throw std::runtime_error("Invalid Formula ID");
        }
        
        FormulaParams coeffs = formulas[formula_id];
        double k1 = coeffs.k1;
        double k2 = coeffs.k2;
        double k3 = coeffs.k3;
        double k4 = coeffs.k4;
        double k5 = coeffs.k5;

        // 3. Preliminary Hydraulic Calculations
        double L0 = calculate_L0(params.Tm);
        double k0 = (2 * M_PI) / L0;
        double s0m = params.Hs / L0;
        double sm = s0m; 

        double Nz = (params.Storm_Duration_hr * 3600.0) / params.Tm;
        double delta_trunk = (params.Wc / params.Ww) - 1.0;

        // 4. Algorithmic Core (Chegini-Aghtouman / Van der Meer) - TRUNK
        double term_damage = std::pow(params.Nod, k2);
        double term_waves = std::pow(Nz, k3);
        double damage_wave_ratio = term_damage / term_waves;
        
        double scaled_term = k1 * damage_wave_ratio;
        double inv_f = scaled_term + k4;
        
        double steepness_factor = std::pow(s0m, -k5);
        
        double Ns_trunk;
        if (formula_id == 1) {
            Ns_trunk = inv_f * steepness_factor * std::pow(2.0/1.5, 1.0/3.0);
        } else {
            Ns_trunk = inv_f * steepness_factor;
        }

        // 5. Block Sizing (Armor) - TRUNK
        double Dn = params.Hs / (delta_trunk * Ns_trunk);
        double W_trunk = params.Wc * std::pow(Dn, 3);
        double packing_density_trunk = 100.0 * 2.0 * 1.1 * (1.0 - P_cubes) / std::pow(Dn, 2);

        // 6. Hudson Comparative Calculation
        double slope = coeffs.slope_ratio;
        double kd_trunk_equiv = (params.Wc * std::pow(params.Hs, 3)) / (W_trunk * std::pow(delta_trunk, 3) * slope);

        // 7. UNDERLAYER - TRUNK
        UnderlayerResult ul_trunk = calculate_underlayer_params(W_trunk);

        // 8. HEAD CALCULATION (FIXED RATIO 1.5)
        double kd_ratio = KD_RATIO_FIXED;
        double kd_head_derived = kd_trunk_equiv / kd_ratio;
        double delta_head = delta_trunk * std::pow(kd_ratio, 1.0/3.0);
        double Wc_head = params.Ww * (delta_head + 1.0);
        double W_head = W_trunk * (Wc_head / params.Wc);

        double Ns_head = params.Hs / (delta_head * Dn);
        double packing_density_head = 100.0 * 2.0 * 1.1 * (1.0 - P_cubes) / std::pow(Dn, 2);

        // 9. UNDERLAYER - HEAD
        UnderlayerResult ul_head = calculate_underlayer_params(W_head);

        // 10. Armor Layer Details (Common)
        double r1 = 2.0 * 1.1 * Dn;

        // Dimensions
        double vol_trunk = W_trunk / params.Wc;
        double h_trunk = std::pow(vol_trunk / 1.0247, 1.0/3.0);
        double a_trunk = 1.086 * h_trunk;
        double b_trunk = 1.005 * h_trunk;

        double vol_head = W_head / Wc_head;
        double h_head = std::pow(vol_head / 1.0247, 1.0/3.0);
        double a_head = 1.086 * h_head;
        double b_head = 1.005 * h_head;

        // 11. Compile Results
        FullResults results;
        results.inputs = params;
        results.coefficients = coeffs;
        
        results.intermediate.L0 = L0;
        results.intermediate.k0 = k0;
        results.intermediate.s0m = s0m;
        results.intermediate.sm = sm;
        results.intermediate.Nz = Nz;
        results.intermediate.delta = delta_trunk;
        results.intermediate.Ns_trunk = Ns_trunk;
        
        results.final_trunk.Dn = Dn;
        results.final_trunk.W = W_trunk;
        results.final_trunk.Mass_tonnes = W_trunk / g;
        results.final_trunk.Kd = kd_trunk_equiv;
        results.final_trunk.r1 = r1;
        results.final_trunk.packing_density = packing_density_trunk;
        results.final_trunk.dims = { h_trunk, a_trunk, b_trunk };
        
        results.underlayer_trunk = ul_trunk;
        
        results.final_head.Ns = Ns_head;
        results.final_head.Kd = kd_head_derived;
        results.final_head.Kd_Ratio = kd_ratio;
        results.final_head.Delta_Required = delta_head;
        results.final_head.Wc_Required = Wc_head;
        results.final_head.W = W_head;
        results.final_head.Mass_tonnes = W_head / g;
        results.final_head.packing_density = packing_density_head;
        results.final_head.dims = { h_head, a_head, b_head };
        
        results.underlayer_head = ul_head;
        results.P_rock = P_rock;
        results.P_cubes = P_cubes;
        
        return results;
    }
};

std::wstring format_report(const FullResults& results) {
    std::wstringstream ss;
    ss << std::fixed << std::setprecision(2);
    
    auto fmt = [](double val, int prec) {
        std::wstringstream stream;
        stream << std::fixed << std::setprecision(prec) << val;
        return stream.str();
    };

    const auto& p = results.inputs;
    const auto& c = results.coefficients;
    const auto& i = results.intermediate;
    const auto& ft = results.final_trunk;
    const auto& ut = results.underlayer_trunk;
    const auto& fh = results.final_head;
    const auto& uh = results.underlayer_head;

    ss << L"================================================================================\n"
       << L"    TECHNICAL REPORT: BREAKWATER ARMOR & UNDERLAYER DESIGN                         \n"
       << L"================================================================================\n"
       << L"Methodology: " << c.name << L"\n"
       << L"--------------------------------------------------------------------------------\n";
    
    ss << L"1. INPUT PARAMETERS\n"
       << L"   Hs (Sigificant Wave Height)         : " << p.Hs << L" m\n"
       << L"   Tm (Mean Wave Period)               : " << p.Tm << L" s\n"
       << L"   Storm Duration (h)                  : " << p.Storm_Duration_hr << L" h\n"
       << L"   Nod (Damage)                        : " << p.Nod << L"\n"
       << L"   Wc Trunk (Concrete Spec. Weight)    : " << p.Wc << L" kN/m3\n"
       << L"   Ww (Water Specific Weight)          : " << p.Ww << L" kN/m3\n"
       << L"   Relative Density D=(Wc/Ww)-1        : " << fmt(i.delta, 4) << L"\n"
       << L"   Structure Slope (TRUNK & HEAD)      : " << fmt(c.slope_ratio, 1) << L":1\n"
       << L"   Porosity (Cubes)                    : " << fmt(results.P_cubes * 100, 0) << L"%\n"
       << L"   Porosity (Rock Layer)               : " << fmt(results.P_rock * 100, 0) << L"%\n"
       << L"--------------------------------------------------------------------------------\n";
    
    ss << L"2. INTERMEDIATE PARAMETERS\n"
       << L"   Wave Length (L0)                    : " << i.L0 << L" m\n"
       << L"   wave number (k0 = 2*pi/L0)          : " << fmt(i.k0, 4) << L"\n"
       << L"   wave steepness (s0m = Hs/L0)        : " << fmt(i.s0m, 4) << L"\n"
       << L"   Number of waves (Nz)                : " << fmt(i.Nz, 0) << L"\n"
       << L"   Stability Number TRUNK (Ns)         : " << fmt(i.Ns_trunk, 4) << L"\n"
       << L"   Stability Number HEAD (Ns)          : " << fmt(fh.Ns, 4) << L"\n"
       << L"--------------------------------------------------------------------------------\n";
    
    ss << L"3. ARMOR LAYER RESULTS - TRUNK\n"
       << L"   BLOCK WEIGHT (W)                    : " << ft.W << L" kN\n"
       << L"   Mass (ton)                          : " << ft.Mass_tonnes << L" t\n"
       << L"   Nominal Diameter (Dn)               : " << fmt(ft.Dn, 3) << L" m\n"
       << L"   Cube Height (H)                     : " << fmt(ft.dims.H, 3) << L" m\n"
       << L"   Cube Top Width (B)                  : " << fmt(ft.dims.B, 3) << L" m\n"
       << L"   Cube Base Width (A)                 : " << fmt(ft.dims.A, 3) << L" m\n"
       << L"   KD_TRUNK (Equivalent)               : " << ft.Kd << L"\n"
       << L"   Double Layer Thickness (r1)         : " << ft.r1 << L" m\n"
       << L"   Packing Density, d [units/100m2]    : " << ft.packing_density << L"\n"
       << L"\n";

    std::wstring ut_grading_wide(ut.grading_name.begin(), ut.grading_name.end());

        ss << L"4. UNDERLAYER RESULTS - TRUNK\n";
        ss << L"   Theoretical Target (W/10)           : " << fmt(ut.target_W, 2) << L" kN (" << fmt(ut.target_M50_kg, 1) << L" kg)\n";
        ss << L"   Adopted rock grading                : " << ut_grading_wide << L"\n";
        ss << L"   Representative M50                  : " << fmt(ut.M50_kg, 1) << L" kg\n";
        ss << L"   Nominal lower limit (NLL)           : " << fmt(ut.NLL_kg, 1) << L" kg\n";
        ss << L"   Nominal upper limit (NUL)           : " << fmt(ut.NUL_kg, 1) << L" kg\n";
        ss << L"   Extreme lower limit (ELL)           : " << fmt(ut.ELL_kg, 1) << L" kg\n";
        ss << L"   Extreme upper limit (EUL)           : " << fmt(ut.EUL_kg, 1) << L" kg\n";
        ss << L"   Nominal Diameter (Dn_rock)          : " << fmt(ut.Dn_rock, 3) << L" m\n";
        ss << L"   Double Layer Thickness (r2)         : " << fmt(ut.r2, 2) << L" m\n";
        ss << L"   Packing Density, f2 [rocks/100m2]   : " << fmt(ut.f2, 2) << L"\n";
        ss << std::wstring(80, L'-') << L"\n";

    ss << L"5. ARMOR LAYER RESULTS - HEAD (High Density)\n"
       << L"   *Maintains same Dn and Slope as Trunk*\n"
       << L"   Stability Ratio (Kd_T/Kd_H)         : " << fh.Kd_Ratio << L"\n"
       << L"   Nominal Diameter (Dn)               : " << fmt(ft.Dn, 3) << L" m\n"
       << L"   Cube Height (H)                     : " << fmt(fh.dims.H, 3) << L" m\n"
       << L"   Cube Top width (B)                  : " << fmt(fh.dims.B, 3) << L" m\n"
       << L"   Cube Base Width (A)                 : " << fmt(fh.dims.A, 3) << L" m\n"
       << L"   KD_HEAD (Equivalent)                : " << fh.Kd << L"\n"
       << L"   Required Concrete Density (Wc)      : " << fh.Wc_Required << L" kN/m3\n"
       << L"   BLOCK WEIGHT (W)                    : " << fh.W << L" kN\n"
       << L"   Mass (ton)                          : " << fh.Mass_tonnes << L" t\n"
       << L"   Packing Density, d [units/100m2]    : " << fh.packing_density << L"\n"
       << L"\n";

    std::wstring uh_grading_wide(uh.grading_name.begin(), uh.grading_name.end());

    ss << L"6. UNDERLAYER RESULTS - HEAD\n";
        ss << L"   Theoretical Target (W/10)           : " << fmt(uh.target_W, 2) << L" kN (" << fmt(uh.target_M50_kg, 1) << L" kg)\n";
        ss << L"   Adopted rock grading                : " << uh_grading_wide << L"\n";
        ss << L"   Representative M50                  : " << fmt(uh.M50_kg, 1) << L" kg\n";
        ss << L"   Nominal lower limit (NLL)           : " << fmt(uh.NLL_kg, 1) << L" kg\n";
        ss << L"   Nominal upper limit (NUL)           : " << fmt(uh.NUL_kg, 1) << L" kg\n";
        ss << L"   Extreme lower limit (ELL)           : " << fmt(uh.ELL_kg, 1) << L" kg\n";
        ss << L"   Extreme upper limit (EUL)           : " << fmt(uh.EUL_kg, 1) << L" kg\n";
        ss << L"   Nominal Diameter (Dn_rock)          : " << fmt(uh.Dn_rock, 3) << L" m\n";
        ss << L"   Double Layer Thickness (r2)         : " << fmt(uh.r2, 2) << L" m\n";
        ss << L"   Packing Density, f2 [rocks/100m2]   : " << fmt(uh.f2, 2) << L"\n";
        ss << std::wstring(80, L'=') << L"\n";

    return ss.str();
}

static std::wstring fix_newlines_for_edit_control(const std::wstring &text) {
    std::wstring out;
    out.reserve(text.size() + 100);
    for (wchar_t c : text) {
        if (c == L'\n') {
            out.push_back(L'\r');
        }
        out.push_back(c);
    }
    return out;
}

// --- FILE OUTPUT HELPER ---
void SaveReportToFile(const std::wstring& report_content) {
    // Append to "output.txt"
    std::ofstream file("output.txt", std::ios::app);
    if (file.is_open()) {
        // Convert wstring to string (narrowing cast is safe for this technical report content)
        std::string narrow_report(report_content.begin(), report_content.end());
        file << narrow_report; 
        file.close();
    }
}

// --- Win32 GUI Specifics ---

#define IDC_EDIT_HS 101
#define IDC_EDIT_TM 102
#define IDC_EDIT_DUR 103
#define IDC_EDIT_NOD 104
#define IDC_EDIT_WC 105
#define IDC_COMBO_FORMULA 106
#define IDC_BUTTON_COMPUTE 107
#define IDC_OUTPUT 108

HWND hEditHs, hEditTm, hEditDur, hEditNod, hEditWc, hComboFormula, hOutput;
BreakwaterCalculator calculator;

// Store fonts globally to delete them later
HFONT hMonoFont = NULL;
HFONT hUIFont = NULL;

void CreateControls(HWND hwnd) {
    // Font Size (20) for User Interface
    hUIFont = CreateFontW(20, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                          DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                          DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");

    int y = 20;
    int lblW = 200; // Wider labels for larger font
    int editX = 220;
    int editW = 100; // Wider edit boxes
    int step = 40;   // vertical spacing

    // Small helper to format double defaults to wstring for the inputs
    auto fmt = [](double val) {
        std::wstringstream ss;
        ss << std::fixed << std::setprecision(2) << val;
        return ss.str();
    };

    auto CreateLabel = [&](const wchar_t* text, int x, int y) {
        HWND h = CreateWindowW(L"STATIC", text, WS_CHILD | WS_VISIBLE, x, y, lblW, 25, hwnd, NULL, NULL, NULL);
        SendMessageW(h, WM_SETFONT, (WPARAM)hUIFont, TRUE);
    };

    auto CreateInput = [&](const wchar_t* def, int id, int y) -> HWND {
        // Fix for 64-bit cast warning: (HMENU)(INT_PTR)id
        HWND h = CreateWindowW(L"EDIT", def, WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL, 
                              editX, y, editW, 25, hwnd, (HMENU)(INT_PTR)id, NULL, NULL);
        SendMessageW(h, WM_SETFONT, (WPARAM)hUIFont, TRUE);
        return h;
    };

    // HS
    CreateLabel(L"Hs (m):", 10, y);
    // Dynamically use defaults
    hEditHs = CreateInput(fmt(calculator.defaults.Hs).c_str(), IDC_EDIT_HS, y);
    
    y += step;
    // TM
    CreateLabel(L"Tm (s):", 10, y);
    hEditTm = CreateInput(fmt(calculator.defaults.Tm).c_str(), IDC_EDIT_TM, y);

    y += step;
    // Duration
    CreateLabel(L"Storm Duration (h):", 10, y);
    hEditDur = CreateInput(fmt(calculator.defaults.Storm_Duration_hr).c_str(), IDC_EDIT_DUR, y);

    y += step;
    // Nod
    CreateLabel(L"Nod (Damage):", 10, y);
    hEditNod = CreateInput(fmt(calculator.defaults.Nod).c_str(), IDC_EDIT_NOD, y);

    y += step;
    // Wc
    CreateLabel(L"Wc Concrete (kN/m3):", 10, y);
    hEditWc = CreateInput(fmt(calculator.defaults.Wc).c_str(), IDC_EDIT_WC, y);

    y += step;
    // Formula Label
    CreateLabel(L"Formula / Slope:", 10, y);
    
    y += 30;
    // Formula ComboBox
    hComboFormula = CreateWindowW(L"COMBOBOX", NULL, 
        WS_CHILD | WS_VISIBLE | CBS_DROPDOWNLIST | WS_VSCROLL, 
        10, y, 350, 200, hwnd, (HMENU)IDC_COMBO_FORMULA, NULL, NULL);
    SendMessageW(hComboFormula, WM_SETFONT, (WPARAM)hUIFont, TRUE);

    // --- COMBOBOX ORDER ---
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"1. Cubes (Slope 2.0:1) - Van der Meer");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"2. Cubes (Slope 1.5:1) - Van der Meer");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"3. Antifer (Slope 2.0:1) - Chegini");
    SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)L"4. Antifer (Slope 1.5:1) - Chegini");
    SendMessageW(hComboFormula, CB_SETCURSEL, 0, 0); 

    y += 45;
    // Compute Button
    HWND hBtn = CreateWindowW(L"BUTTON", L"Calculate Armor", WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
                 10, y, 250, 40, hwnd, (HMENU)IDC_BUTTON_COMPUTE, NULL, NULL);
    SendMessageW(hBtn, WM_SETFONT, (WPARAM)hUIFont, TRUE);

    // Output Box
    hOutput = CreateWindowW(L"EDIT", L"", WS_CHILD | WS_VISIBLE | WS_BORDER |
                           ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL | ES_READONLY | WS_HSCROLL | ES_AUTOHSCROLL,
                           400, 10, 760, 720, hwnd, (HMENU)IDC_OUTPUT, NULL, NULL);

    // Monospace Font (21) for Report
    hMonoFont = CreateFontW(21, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                           DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                           DEFAULT_QUALITY, FIXED_PITCH | FF_DONTCARE, L"Courier New");
    SendMessageW(hOutput, WM_SETFONT, (WPARAM)hMonoFont, TRUE);
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_CREATE:
        CreateControls(hwnd);
        break;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDC_BUTTON_COMPUTE) {
            wchar_t buffer[64];
            Inputs inputs = calculator.defaults; 

            GetWindowTextW(hEditHs, buffer, 63); inputs.Hs = _wtof(buffer);
            GetWindowTextW(hEditTm, buffer, 63); inputs.Tm = _wtof(buffer);
            GetWindowTextW(hEditDur, buffer, 63); inputs.Storm_Duration_hr = _wtof(buffer);
            GetWindowTextW(hEditNod, buffer, 63); inputs.Nod = _wtof(buffer);
            GetWindowTextW(hEditWc, buffer, 63); inputs.Wc = _wtof(buffer);

            int selIndex = (int)SendMessageW(hComboFormula, CB_GETCURSEL, 0, 0);
            inputs.Formula_ID = selIndex + 1;

            if (inputs.Hs <= 0 || inputs.Tm <= 0 || inputs.Storm_Duration_hr <= 0 || inputs.Wc <= 0) {
                 MessageBoxW(hwnd, L"Please enter positive values for hydraulic parameters.", L"Input Error", MB_ICONERROR | MB_OK);
                 break;
            }

            try {
                FullResults res = calculator.solve(inputs.Formula_ID, inputs);
                std::wstring report = format_report(res);
                
                // --- Save to file output.txt ---
                SaveReportToFile(report); 

                std::wstring guiText = fix_newlines_for_edit_control(report);
                SetWindowTextW(hOutput, guiText.c_str());
            } catch (const std::exception& e) {
                const char* what = e.what();
                std::wstring errMsg(what, what + strlen(what));
                MessageBoxW(hwnd, errMsg.c_str(), L"Calculation Error", MB_ICONERROR | MB_OK);
            }
        }
        break;

    case WM_DESTROY:
        if (hMonoFont) DeleteObject(hMonoFont);
        if (hUIFont) DeleteObject(hUIFont);
        PostQuitMessage(0);
        break;

    default:
        return DefWindowProcW(hwnd, msg, wParam, lParam);
    }
    return 0;
}

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE, PWSTR, int nCmdShow) {
    const wchar_t CLASS_NAME[] = L"BreakwaterCalcWindow";

    WNDCLASSEXW wc = {};
    wc.cbSize = sizeof(WNDCLASSEXW);
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = CLASS_NAME;
    wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
    wc.hIconSm = LoadIcon(NULL, IDI_APPLICATION);

    if (!RegisterClassExW(&wc)) {
        return 0;
    }

    HWND hwnd = CreateWindowExW(
        0, CLASS_NAME, L"Breakwater Armor Calculator",
        WS_OVERLAPPEDWINDOW & ~WS_THICKFRAME & ~WS_MAXIMIZEBOX,
        CW_USEDEFAULT, CW_USEDEFAULT, 1200, 800,
        NULL, NULL, hInstance, NULL);

    if (!hwnd) {
        return 0;
    }

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    MSG msg = {};
    while (GetMessageW(&msg, NULL, 0, 0) > 0) {
        TranslateMessage(&msg);
        DispatchMessageW(&msg);
    }

    return static_cast<int>(msg.wParam);
}