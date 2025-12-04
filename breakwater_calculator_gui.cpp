// ======================================================================================
// PROGRAM DESCRIPTION & METHODOLOGY
// ======================================================================================
//
// 1. PURPOSE:
//    This software calculates the required size and weight of artificial concrete 
//    armor units (specifically Cubes and Antifer blocks) for rubble mound breakwaters.
//    It determines requirements for two distinct sections of the structure:
//    - The Trunk (the main longitudinal section).
//    - The Head (the exposed tip/end of the breakwater).
//
// 2. METHODOLOGY:
//    The calculator utilizes empirical stability formulas derived from hydraulic 
//    model testing:
//    - For Antifer Cubes: Uses Chegini & Aghtouman (2006) power-law formulas.
//    - For Simple Cubes: Uses Van der Meer (1988) stability formulas.
//    - For Underlayers: maps calculated requirements to EN 13383 Standard Rock Grades.
//
// 3. DESIGN LOGIC & STRATEGY:
//    a. Trunk Design: 
//       Calculates the Stability Number (Ns) based on wave height, duration, and damage 
//       allowance. This determines the nominal diameter (Dn) and Weight (W).
//
//    b. Head Design (High-Density Strategy): 
//       The breakwater head is subjected to more severe turbulence than the trunk. 
//       Standard practice usually requires larger blocks. 
//       However, this specific script uses a "Constant Geometry, Variable Density" logic:
//       - It enforces a Safety Factor (Kd Ratio = 1.5).
//       - It keeps the Block Size (Dn) identical to the Trunk (simplifying construction).
//       - It increases the Concrete Density (Wc) to achieve the required stability 
//         against the higher forces at the head.
//
// 4. TECHNICAL BIBLIOGRAPHY & REFERENCES:
//
// 1. Van der Meer, J.W. (1988). "Stability of Cubes, Tetrapods and Accropode."
//    Proceedings of the Conference Breakwaters '88, Eastbourne, Thomas Telford.
//
// 2. Van der Meer, J.W. (1998). "Application and stability criteria for rock and 
//    artificial units." Chapter 11 in: Dikes and Revetments. Balkema.
//
// 3. Chegini, V. & Aghtouman, P. (2006). "An Investigation on the Stability of Rubble 
//    Mound Breakwaters with Armour Layers of Antifer Cubes." 
//    Journal of Marine Engineering.
//
// 4. USACE (2006). "Coastal Engineering Manual" (CEM), Chapter VI-5.
//
// 5. EN 13383-1. "Armourstone - Part 1: Specification." (Standard Rock Grades).
//
// 5. COMPILATION INSTRUCTIONS:
//
// Compilation (example using MinGW/g++ on Windows):
//
// g++ -O3 -march=native -std=c++17 -Wall -Wextra -municode breakwater_calculator_gui.cpp -o breakwater_calculator_gui -mwindows -static -static-libgcc -static-libstdc++
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
    std::wstring name;
    double min_kg;
    double max_kg;
};

struct Dimensions {
    double H;
    double A;
    double B;
};

struct UnderlayerResult {
    std::wstring grading_name;
    double target_W;
    double W_mean;
    double W1;
    double W2;
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
            {L"CP45/125",       0.4,   1.2},
            {L"CP63/180",       1.2,   3.8},
            {L"CP90/250",       3.1,   9.3},
            {L"CP45/180",       0.4,   1.2},
            {L"CP90/180",       2.1,   2.8},

            // Light Rock Armour (LRA)
            {L"LRA5-40",        10,    20},
            {L"LRA10-60",       20,    35},
            {L"LRA40-200",      80,    120},
            {L"LRA60-300",      120,   190},
            {L"LRA15-300",      45,    135},

            // Heavy Rock Armour (HRA)
            {L"HRA300-1000",    540,   690},
            {L"HRA1000-3000",   1700,  2100},
            {L"HRA3000-6000",   4200,  4800},
            {L"HRA6000-10000",  7500,  8500},
            {L"HRA10000-15000", 12000, 13000}
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
        
        GradingDef selected_grading;
        double min_diff = std::numeric_limits<double>::infinity();
        
        double final_w_mean = 0;
        double final_w_min = 0;
        double final_w_max = 0;
        bool found = false;

        for (const auto& grading : standard_gradings) {
            double w_min = grading.min_kg * g / 1000.0;
            double w_max = grading.max_kg * g / 1000.0;
            double w_mean_range = (w_min + w_max) / 2.0;
            double diff = std::abs(w_mean_range - target_weight);
            
            if (diff < min_diff) {
                min_diff = diff;
                selected_grading = grading;
                final_w_mean = w_mean_range;
                final_w_min = w_min;
                final_w_max = w_max;
                found = true;
            }
        }
        
        if (!found && !standard_gradings.empty()) {
            selected_grading = standard_gradings[0];
            final_w_min = selected_grading.min_kg * g / 1000.0;
            final_w_max = selected_grading.max_kg * g / 1000.0;
            final_w_mean = (final_w_min + final_w_max) / 2.0;
        }

        double W_mean = final_w_mean;
        double Dn_rock = std::pow(W_mean / W_rock_spec, 1.0/3.0);
        double r2 = 2.0 * Dn_rock;
        double f2 = 100.0 * 2.0 * 1.0 * (1.0 - P_rock) / std::pow(Dn_rock, 2);
        
        UnderlayerResult res;
        res.grading_name = selected_grading.name;
        res.target_W = target_weight;
        res.W_mean = W_mean;
        res.W1 = final_w_min;
        res.W2 = final_w_max;
        res.Dn_rock = Dn_rock;
        res.r2 = r2;
        res.f2 = f2;
        res.W_rock_spec = W_rock_spec;
        return res;
    }

    FullResults solve(int formula_id, Inputs params) {
        if (formulas.find(formula_id) == formulas.end()) {
            throw std::runtime_error("Invalid Formula ID Selected");
        }
        
        FormulaParams coeffs = formulas[formula_id];
        double k1 = coeffs.k1;
        double k2 = coeffs.k2;
        double k3 = coeffs.k3;
        double k4 = coeffs.k4;
        double k5 = coeffs.k5;

        double L0 = calculate_L0(params.Tm);
        double k0 = (2 * M_PI) / L0;
        double s0m = params.Hs / L0;
        double sm = s0m; 

        double Nz = (params.Storm_Duration_hr * 3600.0) / params.Tm;
        double delta_trunk = (params.Wc / params.Ww) - 1.0;

        double term_damage = std::pow(params.Nod, k2);
        double term_waves = std::pow(Nz, k3);
        double damage_wave_ratio = term_damage / term_waves;
        
        double scaled_term = k1 * damage_wave_ratio;
        double inv_f = scaled_term + k4;
        
        double steepness_factor = std::pow(s0m, -k5);
        
        double Ns_trunk;
        
        // Van Der Meer (Slope 2.0) is ID 1
        if (formula_id == 1) {
            Ns_trunk = inv_f * steepness_factor * std::pow(2.0/1.5, 1.0/3.0);
        } else {
            Ns_trunk = inv_f * steepness_factor;
        }

        double Dn = params.Hs / (delta_trunk * Ns_trunk);
        double W_trunk = params.Wc * std::pow(Dn, 3);
        double packing_density_trunk = 100.0 * 2.0 * 1.1 * (1.0 - P_cubes) / std::pow(Dn, 2);

        double slope = coeffs.slope_ratio;
        double kd_trunk_equiv = (params.Wc * std::pow(params.Hs, 3)) / (W_trunk * std::pow(delta_trunk, 3) * slope);

        UnderlayerResult ul_trunk = calculate_underlayer_params(W_trunk);

        double kd_ratio = KD_RATIO_FIXED;
        double kd_head_derived = kd_trunk_equiv / kd_ratio;
        double delta_head = delta_trunk * std::pow(kd_ratio, 1.0/3.0);
        double Wc_head = params.Ww * (delta_head + 1.0);
        double W_head = W_trunk * (Wc_head / params.Wc);

        double Ns_head = params.Hs / (delta_head * Dn);
        double packing_density_head = 100.0 * 2.0 * 1.1 * (1.0 - P_cubes) / std::pow(Dn, 2);

        UnderlayerResult ul_head = calculate_underlayer_params(W_head);

        double r1 = 2.0 * 1.1 * Dn;
        double vol_trunk = W_trunk / params.Wc;
        double h_trunk = std::pow(vol_trunk / 1.0247, 1.0/3.0);
        double a_trunk = 1.086 * h_trunk;
        double b_trunk = 1.005 * h_trunk;

        double vol_head = W_head / Wc_head;
        double h_head = std::pow(vol_head / 1.0247, 1.0/3.0);
        double a_head = 1.086 * h_head;
        double b_head = 1.005 * h_head;

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

// --- Report Formatting ---

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
    
    ss << L"4. UNDERLAYER RESULTS - TRUNK\n"
       << L"   Theoretical Target (W/10)           : " << ut.target_W << L" kN\n"
       << L"   Adopted rock grading                : " << ut.grading_name << L"\n"
       << L"   Grading Min (Lower Limit)           : " << ut.W1 << L" kN\n"
       << L"   Grading Max (Upper Limit)           : " << ut.W2 << L" kN\n"
       << L"   Mean Weight (Used for thickness)    : " << ut.W_mean << L" kN\n"
       << L"   Nominal Diameter (Dn_rock)          : " << fmt(ut.Dn_rock, 3) << L" m\n"
       << L"   Double Layer Thickness (r2)         : " << ut.r2 << L" m\n"
       << L"   Packing Density, f2 [rocks/100m2]   : " << ut.f2 << L"\n"
       << L"--------------------------------------------------------------------------------\n";

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
    
    ss << L"6. UNDERLAYER RESULTS - HEAD\n"
       << L"   Theoretical Target (W/10)           : " << uh.target_W << L" kN\n"
       << L"   Adopted rock grading                : " << uh.grading_name << L"\n"
       << L"   Grading Min (Lower Limit)           : " << uh.W1 << L" kN\n"
       << L"   Grading Max (Upper Limit)           : " << uh.W2 << L" kN\n"
       << L"   Mean Weight (Used for thickness)    : " << uh.W_mean << L" kN\n"
       << L"   Nominal Diameter (Dn_rock)          : " << fmt(uh.Dn_rock, 3) << L" m\n"
       << L"   Double Layer Thickness (r2)         : " << uh.r2 << L" m\n"
       << L"   Packing Density, f2 [rocks/100m2]   : " << uh.f2 << L"\n"
       << L"================================================================================\n\n";

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