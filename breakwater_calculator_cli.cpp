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
//    - For Simple Cubes: Uses Van der Meer (1988) stability formulas.
//    - For Antifer Cubes: Uses Chegini & Aghtouman (2006) power-law formulas.
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
// Compilation (example using g++ on Windows/Linux):
//
// g++ -O3 -march=native -std=c++17 -Wall -Wextra -pedantic -Wconversion 
// -Wsign-conversion -static -static-libgcc -static-libstdc++ 
// -o breakwater_calculator_cli breakwater_calculator_cli.cpp
//
// 6. EXECUTION EXAMPLES:
//
// 1. Interactive Mode (No arguments):
// breakwater_calculator.exe
//
// 2. Command Line Arguments Mode:
// Format: breakwater_calculator.exe [Hs] [Tm] [Duration] [Nod] [Wc] [FormulaID]
//
// Example (Hs=11.5, Tm=13.0, Duration=8.5, Nod=1.0, Wc=26.6, FormulaID=1):
// breakwater_calculator.exe 11.5 13.0 8.5 1.0 26.6 1
//
// ======================================================================================

// ----------------------------------------------------------------------
// MINGW STATIC LINKING SHIM
// Fixes "undefined reference to __imp_fseeko64" errors in GCC 14/15+
// ----------------------------------------------------------------------
#if defined(_WIN32) && defined(__GNUC__)
    #include <cstdio>
    extern "C" {
        // Force the linker to use the local _fseeki64 when it asks for the import stub
        int (*__imp_fseeko64)(FILE*, long long, int) = reinterpret_cast<int(*)(FILE*, long long, int)>(&_fseeki64);
        long long (*__imp_ftello64)(FILE*) = reinterpret_cast<long long(*)(FILE*)>(&_ftelli64);
    }
#endif
// ----------------------------------------------------------------------

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <limits>

// Define M_PI if not already defined
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ----------------------------------------------------------------------
// DATA STRUCTURES
// ----------------------------------------------------------------------

struct FormulaParams {
    std::string name;
    std::string type;
    double slope_ratio;
    double k1, k2, k3, k4, k5;
};

struct GradingDef {
    std::string name;
    double min_kg;
    double max_kg;
};

struct Dimensions {
    double H;
    double A;
    double B;
};

struct UnderlayerResult {
    std::string grading_name;
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
    double Ns; // Store Ns specific to section
    double Dn;
    double W;
    double Mass_tonnes;
    double Kd; // Equivalent or Derived
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

// ----------------------------------------------------------------------
// CLASS DEFINITION
// ----------------------------------------------------------------------

class BreakwaterCalculator {
private:
    // ----------------------------------------------------------------------
    // PHYSICAL CONSTANTS
    // ----------------------------------------------------------------------
    const double g = 9.80665;  // Acceleration due to gravity (m/s^2)
    const double W_rock_spec = 26.5; // Standard Specific weight for underlayer rock (kN/m3)
    
    // ----------------------------------------------------------------------
    // LAYER CHARACTERISTICS
    // ----------------------------------------------------------------------
    const double P_cubes = 0.40;  // Porosity of armor layers (Standard for double-layer Cubes)
    const double P_rock = 0.25;   // Porosity of rock underlayers (Standard approximation)
    
    // ----------------------------------------------------------------------
    // HEAD TO TRUNK TRANSFER RATIO
    // ----------------------------------------------------------------------
    const double KD_RATIO_FIXED = 1.5;

    std::map<int, FormulaParams> formulas;
    std::vector<GradingDef> standard_gradings;

public:
    Inputs defaults;

    BreakwaterCalculator() {
        // ----------------------------------------------------------------------
        // FORMULA DATABASE (REORDERED)
        // ----------------------------------------------------------------------
        // ID 1: Van Der Meer (1988a) - Cubes (Slope 2.0:1)
        formulas[1] = { "Van Der Meer (1988a) - Cubes (Slope 2.0:1)", "Cubes", 2.0, 6.7, 0.4, 0.3, 1.0, 0.1 };
        
        // ID 2: Van Der Meer (1988a) - Cubes (Slope 1.5:1)
        formulas[2] = { "Van Der Meer (1988a) - Cubes (Slope 1.5:1)", "Cubes", 1.5, 6.7, 0.4, 0.3, 1.0, 0.1 };

        // ID 3: Chegini-Aghtouman (2006) - Antifer (Slope 2:1)
        formulas[3] = { "Chegini-Aghtouman (2006) - Antifer (Slope 2:1)", "Antifer", 2.0, 6.138, 0.443, 0.276, 1.164, 0.07 };
        
        // ID 4: Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)
        formulas[4] = { "Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)", "Antifer", 1.5, 6.951, 0.443, 0.291, 1.082, 0.082 };

        // ----------------------------------------------------------------------
        // ROCK GRADING DATABASE (EN 13383)
        // ----------------------------------------------------------------------
        standard_gradings = {
            {"EN 13383 LMA 5-40 kg",     5,    40},
            {"EN 13383 LMA 10-60 kg",    10,   60},
            {"EN 13383 LMA 40-200 kg",   40,   200},
            {"EN 13383 LMA 60-300 kg",   60,   300},
            {"EN 13383 HMA 300-1000 kg", 300,  1000},
            {"EN 13383 HMA 1-3 ton",     1000, 3000},
            {"EN 13383 HMA 3-6 ton",     3000, 6000},
            {"EN 13383 HMA 6-10 ton",    6000, 10000}
        };

        // ----------------------------------------------------------------------
        // DEFAULT INPUTS
        // ----------------------------------------------------------------------
        defaults = {
            10.0,   // Hs
            13.0,   // Tm
            12.0,   // Storm_Duration_hr
            1.0,    // Nod
            24.0,   // Wc
            10.05,  // Ww
            1       // Formula_ID (Now defaults to Cubes 2.0:1)
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
        
        // Fallback if vector empty (shouldn't happen)
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
        // 2. Load Coefficients
        if (formulas.find(formula_id) == formulas.end()) {
            std::cerr << "Invalid Formula ID: " << formula_id << std::endl;
            exit(1);
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
        
        // Updated Condition: 
        // Originally this was (formula_id == 3) for "Van Der Meer (Slope 2.0:1)"
        if (formula_id == 1) {
            Ns_trunk = inv_f * steepness_factor * std::pow(2.0/1.5, 1.0/3.0);
        } else {
            Ns_trunk = inv_f * steepness_factor;
        }

        // 5. Block Sizing (Armor) - TRUNK
        double Dn = params.Hs / (delta_trunk * Ns_trunk);
        double W_trunk = params.Wc * std::pow(Dn, 3);
        double packing_density_trunk = 100.0 * 2.0 * 1.0 * (1.0 - P_cubes) / std::pow(Dn, 2);

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

        // Calculate Stability Number for Head
        double Ns_head = params.Hs / (delta_head * Dn);
        double packing_density_head = 100.0 * 2.0 * 1.0 * (1.0 - P_cubes) / std::pow(Dn, 2);

        // 9. UNDERLAYER - HEAD
        UnderlayerResult ul_head = calculate_underlayer_params(W_head);

        // 10. Armor Layer Details (Common)
        double r1 = 2.0 * Dn;

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

    void generate_report_file(const FullResults& results, std::string filepath="output.txt") {
        std::stringstream ss;
        
        auto fmt = [](double val, int prec) {
            std::stringstream stream;
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

        ss << "================================================================================" << "\n";
        ss << "    TECHNICAL REPORT: BREAKWATER ARMOR & UNDERLAYER DESIGN                        " << "\n";
        ss << "================================================================================" << "\n";
        ss << "Methodology: " << c.name << "\n";
        ss << std::string(80, '-') << "\n";
        
        ss << "1. INPUT PARAMETERS" << "\n";
        ss << "   Hs (Sigificant Wave Height)         : " << fmt(p.Hs, 2) << " m" << "\n";
        ss << "   Tm (Mean Wave Period)               : " << fmt(p.Tm, 2) << " s" << "\n";
        ss << "   Storm Duration (h)                  : " << fmt(p.Storm_Duration_hr, 2) << " h" << "\n";
        ss << "   Nod (Damage)                        : " << fmt(p.Nod, 2) << "\n";
        ss << "   Wc Trunk (Concrete Spec. Weight)    : " << fmt(p.Wc, 2) << " kN/m3" << "\n";
        ss << "   Ww (Water Specific Weight)          : " << fmt(p.Ww, 2) << " kN/m3" << "\n";
        ss << "   Relative Density D=(Wc/Ww)-1        : " << fmt(i.delta, 5) << "\n";
        ss << "   Structure Slope (TRUNK & HEAD)      : " << fmt(c.slope_ratio, 1) << ":1" << "\n";
        ss << "   Porosity (Cubes)                    : " << fmt(results.P_cubes * 100, 0) << "%" << "\n";
        ss << "   Porosity (Rock Layer)               : " << fmt(results.P_rock * 100, 0) << "%" << "\n";
        ss << std::string(80, '-') << "\n";
        
        ss << "2. INTERMEDIATE PARAMETERS" << "\n";
        ss << "   Wave Length (L0)                    : " << fmt(i.L0, 2) << " m" << "\n";
        ss << "   wave number (k0 = 2*pi/L0)          : " << fmt(i.k0, 5) << "\n";
        ss << "   wave steepness (s0m = Hs/L0)        : " << fmt(i.s0m, 5) << "\n";
        ss << "   Number of waves (Nz)                : " << fmt(i.Nz, 0) << "\n";
        ss << "   Stability Number TRUNK (Ns)         : " << fmt(i.Ns_trunk, 5) << "\n";
        ss << "   Stability Number HEAD (Ns)          : " << fmt(fh.Ns, 5) << "\n";
        ss << std::string(80, '-') << "\n";
        
        ss << "3. ARMOR LAYER RESULTS - TRUNK" << "\n";
        ss << "   BLOCK WEIGHT (W)                    : " << fmt(ft.W, 2) << " kN" << "\n";
        ss << "   Mass (ton)                          : " << fmt(ft.Mass_tonnes, 2) << " t" << "\n";
        ss << "   Nominal Diameter (Dn)               : " << fmt(ft.Dn, 3) << " m" << "\n";
        ss << "   Cube Height (H)                     : " << fmt(ft.dims.H, 3) << " m" << "\n";
        ss << "   Cube Top Width (B)                  : " << fmt(ft.dims.B, 3) << " m" << "\n";
        ss << "   Cube Base Width (A)                 : " << fmt(ft.dims.A, 3) << " m" << "\n";
        ss << "   KD_TRUNK (Equivalent)               : " << fmt(ft.Kd, 2) << "\n";
        ss << "   Double Layer Thickness (r1)         : " << fmt(ft.r1, 2) << " m" << "\n";
        ss << "   Packing Density, d [units/100m2]    : " << fmt(ft.packing_density, 2) << "\n";
        ss << "" << "\n";
        
        ss << "4. UNDERLAYER RESULTS - TRUNK" << "\n";
        ss << "   Theoretical Target (W/10)           : " << fmt(ut.target_W, 2) << " kN" << "\n";
        ss << "   Adopted rock grading                : " << ut.grading_name << "\n";
        ss << "   Grading Min (Lower Limit)           : " << fmt(ut.W1, 2) << " kN" << "\n";
        ss << "   Grading Max (Upper Limit)           : " << fmt(ut.W2, 2) << " kN" << "\n";
        ss << "   Mean Weight (Used for thickness)    : " << fmt(ut.W_mean, 2) << " kN" << "\n";
        ss << "   Nominal Diameter (Dn_rock)          : " << fmt(ut.Dn_rock, 3) << " m" << "\n";
        ss << "   Double Layer Thickness (r2)         : " << fmt(ut.r2, 2) << " m" << "\n";
        ss << "   Packing Density, f2 [rocks/100m2]   : " << fmt(ut.f2, 2) << "\n";
        ss << std::string(80, '-') << "\n";

        ss << "5. ARMOR LAYER RESULTS - HEAD (High Density)" << "\n";
        ss << "   *Maintains same Dn and Slope as Trunk*" << "\n";
        ss << "   Stability Ratio (Kd_T/Kd_H)         : " << fmt(fh.Kd_Ratio, 2) << "\n";
        ss << "   Nominal Diameter (Dn)               : " << fmt(ft.Dn, 3) << " m" << "\n";
        ss << "   Cube Height (H)                     : " << fmt(fh.dims.H, 3) << " m" << "\n";
        ss << "   Cube Top width (B)                  : " << fmt(fh.dims.B, 3) << " m" << "\n";
        ss << "   Cube Base Width (A)                 : " << fmt(fh.dims.A, 3) << " m" << "\n";
        ss << "   KD_HEAD (Equivalent)                : " << fmt(fh.Kd, 2) << "\n";
        ss << "   Required Concrete Density (Wc)      : " << fmt(fh.Wc_Required, 2) << " kN/m3" << "\n";
        ss << "   BLOCK WEIGHT (W)                    : " << fmt(fh.W, 2) << " kN" << "\n";
        ss << "   Mass (ton)                          : " << fmt(fh.Mass_tonnes, 2) << " t" << "\n";
        ss << "   Packing Density, d [units/100m2]    : " << fmt(fh.packing_density, 2) << "\n";
        ss << "" << "\n";
        
        ss << "6. UNDERLAYER RESULTS - HEAD" << "\n";
        ss << "   Theoretical Target (W/10)           : " << fmt(uh.target_W, 2) << " kN" << "\n";
        ss << "   Adopted rock grading                : " << uh.grading_name << "\n";
        ss << "   Grading Min (Lower Limit)           : " << fmt(uh.W1, 2) << " kN" << "\n";
        ss << "   Grading Max (Upper Limit)           : " << fmt(uh.W2, 2) << " kN" << "\n";
        ss << "   Mean Weight (Used for thickness)    : " << fmt(uh.W_mean, 2) << " kN" << "\n";
        ss << "   Nominal Diameter (Dn_rock)          : " << fmt(uh.Dn_rock, 3) << " m" << "\n";
        ss << "   Double Layer Thickness (r2)         : " << fmt(uh.r2, 2) << " m" << "\n";
        ss << "   Packing Density, f2 [rocks/100m2]   : " << fmt(uh.f2, 2) << "\n";
        ss << "================================================================================" << "\n";

        std::string report_content = ss.str();

        std::ofstream outfile(filepath);
        if (outfile.is_open()) {
            outfile << report_content;
            outfile.close();
            std::cout << "\n Report generated successfully: " << filepath << std::endl;
            std::cout << "Report content:\n\n";
            std::cout << report_content << std::endl;
        } else {
            std::cerr << "\n Error saving file." << std::endl;
        }
    }
};

// ==============================================================================
// MAIN EXECUTION BLOCK
// ==============================================================================
int main(int argc, char* argv[]) {
    BreakwaterCalculator calc;

    // Helper to get param with default
    auto get_param = [](std::string prompt, double default_val) -> double {
        std::cout << prompt << " [" << default_val << "]: ";
        std::string input_str;
        std::getline(std::cin, input_str);
        if (input_str.empty()) {
            return default_val;
        }
        try {
            return std::stod(input_str);
        } catch (...) {
            return default_val;
        }
    };

    Inputs user_inputs = calc.defaults;
    int formula_id = 1;

    // Check if command line arguments are provided
    if (argc >= 7) {
        // Parse CLI arguments
        try {
            user_inputs.Hs = std::stod(argv[1]);
            user_inputs.Tm = std::stod(argv[2]);
            user_inputs.Storm_Duration_hr = std::stod(argv[3]);
            user_inputs.Nod = std::stod(argv[4]);
            user_inputs.Wc = std::stod(argv[5]);
            formula_id = std::stoi(argv[6]);
            
            if (formula_id < 1 || formula_id > 4) {
                 std::cerr << "Error: Formula ID must be 1-4. Using default (1)." << std::endl;
                 formula_id = 1;
            }
            user_inputs.Formula_ID = formula_id;
        } catch (const std::exception& e) {
            std::cerr << "Error parsing command line arguments: " << e.what() << std::endl;
            std::cerr << "Usage: " << argv[0] << " [Hs] [Tm] [Duration] [Nod] [Wc] [FormulaID]" << std::endl;
            return 1;
        }
    } else {
        // Interactive Mode
        std::cout << "\n--- COASTAL PROTECTION BLOCK CALCULATOR (TRUNK & HEAD) ---" << std::endl;
        std::cout << "1. Van Der Meer (1988) - Cubes (Slope 2.0:1)" << std::endl;
        std::cout << "2. Van Der Meer (1988) - Cubes (Slope 1.5:1)" << std::endl;
        std::cout << "3. Chegini-Aghtouman (2006) - Antifer (Slope 2:1)" << std::endl;
        std::cout << "4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)" << std::endl;

        std::cout << "\nOption [1-4]: ";
        std::string selection;
        std::getline(std::cin, selection);
        selection.erase(selection.find_last_not_of(" \n\r\t") + 1);

        if (selection == "1" || selection == "2" || selection == "3" || selection == "4") {
            formula_id = std::stoi(selection);
            user_inputs.Formula_ID = formula_id;

            std::cout << "\n--- Enter Parameters (Press ENTER for Default) ---" << std::endl;
            user_inputs.Hs = get_param("Hs (m)", calc.defaults.Hs);
            user_inputs.Tm = get_param("Tm (s)", calc.defaults.Tm);
            user_inputs.Nod = get_param("Nod (Damage)", calc.defaults.Nod);
            user_inputs.Storm_Duration_hr = get_param("Duration (h)", calc.defaults.Storm_Duration_hr);
            user_inputs.Wc = get_param("Concrete Weight Trunk (kN/m3)", calc.defaults.Wc);
        }
    }

    try {
        FullResults results = calc.solve(formula_id, user_inputs);
        calc.generate_report_file(results, "output.txt");
    } catch (const std::exception& e) {
        std::cout << "\n Calculation Error: " << e.what() << std::endl;
    }

    return 0;
}