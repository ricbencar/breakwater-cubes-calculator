// ======================================================================================
// PROGRAM DESCRIPTION & METHODOLOGY
// ======================================================================================
//
// 1. PURPOSE:
//    This software performs preliminary hydraulic sizing of rubble-mound
//    breakwater armor made with artificial concrete units and the associated
//    rock underlayers.
//
//    The currently implemented armor-unit families are:
//    - Simple Cubes
//    - Antifer Blocks
//
//    The program dimensions two structural zones:
//    - The Trunk (main longitudinal section)
//    - The Head (exposed roundhead / end section)
//
//    The implemented design philosophy is iso-geometric, but not iso-weight:
//    - the nominal diameter is kept the same at Trunk and Head;
//    - the armor-unit geometry is kept the same at Trunk and Head;
//    - the slope is kept the same at Trunk and Head;
//    - the Head unit weight increases through higher required concrete specific
//      weight rather than through larger unit size.
//
// 2. METHODOLOGY:
//    The calculator uses empirical hydraulic-stability formulas derived from
//    physical model testing together with an internal EN 13383 rock-grading
//    database.
//
//    Implemented armor formulas:
//    - Simple Cubes: Van der Meer (1988)
//    - Antifer Blocks: Chegini & Aghtouman (2006)
//
//    Underlayers:
//    - standard EN 13383 grading selection from an embedded mass-based
//      grading table, with optional custom interpolated grading by family
//      (AUTO / HMA / LMA / CP) for the rock underlayers only.
//
// 3. DESIGN LOGIC & COMPUTATIONAL STRATEGY:
//
//    a. Exposed storm input:
//       The storm variable entered by the user is the number of waves, Nz.
//       Storm duration is not entered directly; it is computed internally as:
//
//           Storm_Duration_hr = Nz * Tm / 3600
//
//    b. Trunk stability:
//       For the selected formula, the program evaluates the Stability Number:
//
//           Ns = (k1 * (Nod^k2 / Nz^k3) + k4) * s0m^(-k5)
//
//       where:
//
//           L0  = g * Tm^2 / (2 * pi)
//           s0m = Hs / L0
//
//       The trunk buoyant relative density is:
//
//           Delta_trunk = (Wc_trunk / Ww) - 1
//
//       The nominal diameter is then:
//
//           Dn = Hs / (Delta_trunk * Ns)
//
//       and the trunk armor-unit weight is:
//
//           W_trunk = Wc_trunk * Dn^3
//
//       The program also derives an equivalent Hudson stability coefficient
//       for the trunk.
//
//    c. Head design - iso-geometric transfer:
//       The breakwater head is subjected to stronger three-dimensional flow
//       effects than the trunk. Instead of increasing block size or flattening
//       the slope, the program preserves armor geometry and transfers the trunk
//       design to the head through a fixed Hudson-coefficient ratio:
//
//           Kd_trunk / Kd_head = 1.5
//
//       With geometry kept fixed, the required head relative density becomes:
//
//           Delta_head = Delta_trunk * (1.5)^(1/3)
//
//       and the required concrete specific weight at the head is:
//
//           Wc_head = Ww * (Delta_head + 1)
//
//       Since block volume is kept constant, the head block weight scales only
//       with concrete specific weight:
//
//           W_head = W_trunk * (Wc_head / Wc_trunk)
//
//       Practical consequence:
//       - the same molds can be used at Trunk and Head;
//       - the same nominal diameter is maintained;
//       - Head units are heavier because concrete specific weight increases.
//
//    d. Underlayer selection:
//       The theoretical underlayer target is one tenth of the armor-unit weight:
//
//           W_target = W_armor / 10
//
//       The target mass is compared against the embedded EN 13383 grading
//       limits. The program selects a grading only if the target mass is
//       strictly inside the nominal lower and upper limits.
//
//       If more than one grading contains the target mass, the program selects
//       the tightest range, i.e. the grading with the smallest:
//
//           NUL - NLL
//
//       If no grading strictly contains the target mass, or if the user
//       explicitly disables EN 13383 mode, the program derives a custom
//       interpolated underlayer grading by family (AUTO / HMA / LMA / CP).
//
// 4. DEFAULT INPUTS, CONSTANTS, AND CURRENT FORMULA SETS USED IN THE CODE:
//
//    a. Default user inputs currently used by the CLI code:
//       - Hs              = 11.0 m
//       - Tm              = 11.9 s
//       - Number_of_Waves = 3000
//       - Nod             = 0.5
//       - Wc              = 27.48 kN/m3
//       - Ww              = 10.05 kN/m3
//       - Formula_ID      = 1
//       - UseEN13383    = true
//       - CustomFamily  = AUTO
//
//       Default Formula_ID = 1 corresponds to:
//       Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)
//
//    b. Internal constants currently used by the code:
//       - g              = 9.80665 m/s2
//       - W_rock_spec    = 26.5 kN/m3
//       - P_cubes        = 0.40
//       - P_rock         = 0.25
//       - KD_RATIO_FIXED = 1.5
//
//    c. Formula sets currently embedded in the code:
//       1. Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)
//       2. Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)
//       3. Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)
//       4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)
//
// 5. OUTPUTS PRODUCED:
//    The generated report contains, for both Trunk and Head:
//    - stability number;
//    - nominal diameter;
//    - unit weight and unit mass;
//    - equivalent Kd value;
//    - armor layer thickness;
//    - packing density;
//    - geometric dimensions reported by the code;
//    - adopted underlayer grading and derived rock parameters.
//
//    Volume reporting follows the current code logic:
//
//    - for Simple Cubes:
//
//          V = Dn^3
//
//    - for Antifer blocks:
//
//          V = 1.0247 * H^3
//
//    The CLI version writes the technical report to output.txt and overwrites
//    that file at each execution.
//
// 6. TECHNICAL BIBLIOGRAPHY & REFERENCES:
//
//    1. Van der Meer, J.W. (1988). "Rock Slopes and Gravel Beaches Under Wave
//       Attack." Doctoral Thesis, Delft University of Technology.
//
//    2. Van der Meer, J.W. (1988). "Stability of Cubes, Tetrapods and
//       Accropode." Proceedings of the Conference Breakwaters '88, Eastbourne,
//       Thomas Telford.
//
//    3. Chegini, V., & Aghtouman, P. (2006). "An Investigation on the
//       Stability of Rubble Mound Breakwaters with Armour Layers of Antifer
//       Cubes."
//
//    4. USACE (2006). "Coastal Engineering Manual (CEM)", Chapter VI-5.
//
//    5. CEN (2002). "EN 13383-1: Armourstone - Part 1: Specification."
//
// 7. COMPILATION INSTRUCTIONS:
//
//    Recommended GCC / MinGW-w64 build with static runtime linking:
//
//    g++ -O3 -march=native -std=c++17 -Wall -Wextra -static -static-libgcc
//    -static-libstdc++ -o breakwater_calculator_cli breakwater_calculator_cli.cpp
//
//    Meaning of the main compilation flags:
//    - -O3               -> high optimization level;
//    - -march=native     -> optimize for the local CPU;
//    - -std=c++17        -> compile as C++17;
//    - -Wall -Wextra     -> enable common compiler warnings;
//    - -static           -> prefer static linking;
//    - -static-libgcc    -> statically link GCC runtime;
//    - -static-libstdc++ -> statically link C++ runtime.
//
//    If full static linking fails on the local toolchain, a non-static fallback
//    build can be used:
//
//    g++ -O3 -march=native -std=c++17 -Wall -Wextra -o
//    breakwater_calculator_cli breakwater_calculator_cli.cpp
//
//    Note:
//    This source file already contains a MinGW/GCC Windows compatibility shim
//    intended to help avoid "__imp_fseeko64" / "__imp_ftello64" static-linking
//    errors on newer GCC toolchains.
//
// 8. EXECUTION MODES AND EXAMPLES:
//
//    a. Interactive mode:
//       Run the executable with no command-line arguments:
//
//           ./breakwater_calculator_cli
//
//       The program will:
//       - display the 4 available formula options;
//       - ask for the formula selection;
//       - ask for Hs, Tm, Nz, Nod, Wc, and underlayer grading mode;
//       - accept ENTER to use each default value;
//       - calculate results and write output.txt.
//
//    b. Command-line argument mode:
//       Format:
//
//           ./breakwater_calculator_cli [Hs] [Tm] [NumberOfWaves] [Nod] [Wc] [FormulaID] [UseEN13383] [CustomFamily]
//
//       Example using the current code defaults:
//
//           ./breakwater_calculator_cli 11.0 11.9 3000 0.5 27.48 1 true AUTO
//
//       FormulaID must be one of:
//       - 1 -> Simple Cubes, slope 2.0:1
//       - 2 -> Simple Cubes, slope 1.5:1
//       - 3 -> Antifer, slope 2.0:1
//       - 4 -> Antifer, slope 1.5:1
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

    // Custom interpolated grading metadata
    bool used_custom_interpolation;
    std::string custom_family;
    double custom_ratio_nul_nll;
    std::string custom_ratio_note;
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
    double Storm_Duration_hr;
    double delta;
    double Ns_trunk;
};

struct Inputs {
    double Hs;
    double Tm;
    double Number_of_Waves;
    double Nod;
    double Wc;
    double Ww;
    int Formula_ID;
    bool use_en13383;          // true = standard EN 13383 grading; false = custom interpolated grading
    std::string custom_family; // AUTO, HMA, LMA, or CP when using custom grading
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
        // ID 1: Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)
        formulas[1] = { "Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)", "Cubes", 2.0, 7.374304189198, 0.4, 0.3, 1.100642416298, 0.1 };
        
        // ID 2: Van Der Meer (1988a) - Simple Cubes (Slope 1.5:1)
        formulas[2] = { "Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)", "Cubes", 1.5, 6.7, 0.4, 0.3, 1.0, 0.1 };

        // ID 3: Chegini-Aghtouman (2006) - Antifer (Slope 2:1)
        formulas[3] = { "Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)", "Antifer", 2.0, 6.138, 0.443, 0.276, 1.164, 0.07 };
        
        // ID 4: Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)
        formulas[4] = { "Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)", "Antifer", 1.5, 6.951, 0.443, 0.291, 1.082, 0.082 };

        // ----------------------------------------------------------------------
        // ROCK GRADING DATABASE (EN 13383)
        // ----------------------------------------------------------------------
        standard_gradings = {
            // Coarse / Light Gradings
            {"CP32/90",        0.868,   19.319},
            {"CP45/125",       2.415,	51.758},
            {"CP63/180",       6.626,	154.548},
            {"CP90/250",       19.319,	414.063},
            {"CP45/180",       2.415,	154.548},
            {"CP90/180",       19.319,	154.548},

            // Light Mass Armourstone (LMA)
            {"LMA5-40",        5,     40},
            {"LMA10-60",       10,    60},
            {"LMA15-120",      15,    120},
            {"LMA40-200",      40,    200},
            {"LMA60-300",      60,    300},
            {"LMA15-300",      15,    300},

            // Heavy Mass Armourstone (HMA)
            {"HMA300-1000",    300,   1000},
            {"HMA1000-3000",   1000,  3000},
            {"HMA3000-6000",   3000,  6000},
            {"HMA6000-10000",  6000,  10000},
            {"HMA10000-15000", 10000, 15000}
        };
		
        // ----------------------------------------------------------------------
        // DEFAULT INPUTS
        // ----------------------------------------------------------------------
        defaults = {
            11.0,   // Hs
            11.9,   // Tm
            3000.0, // Number_of_Waves
            0.5,    // Nod
            27.48,  // Wc
            10.05,  // Ww
            1,      // Formula_ID (Now defaults to Cubes 2.0:1)
            true,   // use_en13383
            "AUTO"  // custom_family for custom grading mode
        };
    }

    double calculate_L0(double Tm) {
        return (g * std::pow(Tm, 2)) / (2 * M_PI);
    }

    std::string get_formula_display_name(int formula_id) const {
        switch (formula_id) {
            case 1: return "Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)";
            case 2: return "Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)";
            case 3: return "Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)";
            case 4: return "Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)";
            default: return "Unknown formula";
        }
    }

    static bool starts_with(const std::string& value, const std::string& prefix) {
        return value.rfind(prefix, 0) == 0;
    }

    std::string grading_family(const std::string& grading_name) const {
        if (starts_with(grading_name, "HMA")) return "HMA";
        if (starts_with(grading_name, "LMA")) return "LMA";
        if (starts_with(grading_name, "CP"))  return "CP";
        return "UNKNOWN";
    }

    std::vector<GradingDef> get_family_gradings(const std::string& family) const {
        std::vector<GradingDef> out;
        for (const auto& gdef : standard_gradings) {
            if (grading_family(gdef.name) == family) out.push_back(gdef);
        }
        std::sort(out.begin(), out.end(), [](const GradingDef& a, const GradingDef& b) {
            const double ma = 0.5 * (a.NLL_kg + a.NUL_kg);
            const double mb = 0.5 * (b.NLL_kg + b.NUL_kg);
            return ma < mb;
        });
        return out;
    }

    std::string select_custom_family(double target_mass) const {
        double safe_mass = std::max(target_mass, 1e-9);
        double best_score = std::numeric_limits<double>::max();
        std::string best_family = "LMA";

        for (const auto& gdef : standard_gradings) {
            const std::string fam = grading_family(gdef.name);
            if (fam == "UNKNOWN") continue;

            const double m50 = 0.5 * (gdef.NLL_kg + gdef.NUL_kg);
            double score = std::abs(std::log(safe_mass) - std::log(std::max(m50, 1e-9)));

            // Mild engineering preference: favour mass-based armourstone families
            // over CP when distances are very similar.
            if (fam == "CP") score += 0.08;

            if (score < best_score) {
                best_score = score;
                best_family = fam;
            }
        }
        return best_family;
    }

    double interpolate_family_ratio(double target_mass, const std::string& family, std::string& note) const {
        std::vector<GradingDef> family_gradings = get_family_gradings(family);
        if (family_gradings.empty()) {
            note = "fallback ratio R=NUL/NLL = 3.0 (no family data)";
            return 3.0;
        }

        std::vector<double> x;
        std::vector<double> y;
        x.reserve(family_gradings.size());
        y.reserve(family_gradings.size());

        for (const auto& gdef : family_gradings) {
            const double m50 = 0.5 * (gdef.NLL_kg + gdef.NUL_kg);
            x.push_back(std::log(std::max(m50, 1e-9)));
            y.push_back(std::log(std::max(gdef.NUL_kg / gdef.NLL_kg, 1.0 + 1e-9)));
        }

        double xt = std::log(std::max(target_mass, 1e-9));
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3);

        if (family_gradings.size() == 1) {
            double ratio = family_gradings.front().NUL_kg / family_gradings.front().NLL_kg;
            oss << family << " family single-point ratio used";
            note = oss.str();
            return ratio;
        }

        if (xt <= x.front()) {
            double ratio = std::exp(y.front());
            oss << family << " family ratio clamped to lower-end class " << family_gradings.front().name;
            note = oss.str();
            return ratio;
        }

        if (xt >= x.back()) {
            double ratio = std::exp(y.back());
            oss << family << " family ratio clamped to upper-end class " << family_gradings.back().name;
            note = oss.str();
            return ratio;
        }

        for (size_t i = 0; i + 1 < family_gradings.size(); ++i) {
            if (xt >= x[i] && xt <= x[i + 1]) {
                double t = (xt - x[i]) / (x[i + 1] - x[i]);
                double log_ratio = y[i] + t * (y[i + 1] - y[i]);
                double ratio = std::exp(log_ratio);
                oss << family << " family ratio interpolated between "
                    << family_gradings[i].name << " and " << family_gradings[i + 1].name;
                note = oss.str();
                return ratio;
            }
        }

        double ratio = std::exp(y.back());
        oss << family << " family ratio fallback to upper-end class " << family_gradings.back().name;
        note = oss.str();
        return ratio;
    }

    UnderlayerResult calculate_underlayer_params(double W_armor, bool use_en13383, const std::string& requested_custom_family) {
        double target_weight = W_armor / 10.0;
        double target_mass_kg = (target_weight * 1000.0) / g;

        UnderlayerResult res{};
        res.target_W = target_weight;
        res.target_M50_kg = target_mass_kg;
        res.used_custom_interpolation = false;
        res.custom_family = "";
        res.custom_ratio_nul_nll = 0.0;
        res.custom_ratio_note = "";

        if (use_en13383) {
            GradingDef selected_grading{};
            bool found = false;
            double min_range_width = std::numeric_limits<double>::max();

            for (const auto& grading : standard_gradings) {
                if (target_mass_kg > grading.NLL_kg && target_mass_kg < grading.NUL_kg) {
                    double current_range = grading.NUL_kg - grading.NLL_kg;
                    if (current_range < min_range_width) {
                        min_range_width = current_range;
                        selected_grading = grading;
                        found = true;
                    }
                }
            }

            if (found) {
                res.grading_name = selected_grading.name;
                res.M50_kg = 0.5 * (selected_grading.NLL_kg + selected_grading.NUL_kg);
                res.NLL_kg = selected_grading.NLL_kg;
                res.NUL_kg = selected_grading.NUL_kg;
            }
        }

        if (!use_en13383 || res.NUL_kg <= res.NLL_kg) {
            std::string family = requested_custom_family;
            std::transform(family.begin(), family.end(), family.begin(),
                           [](unsigned char c) { return static_cast<char>(std::toupper(c)); });

            if (family != "HMA" && family != "LMA" && family != "CP") {
                family = select_custom_family(target_mass_kg);
            }

            std::string ratio_note;
            double ratio_nul_nll = interpolate_family_ratio(target_mass_kg, family, ratio_note);
            ratio_nul_nll = std::max(ratio_nul_nll, 1.01);

            double nll_kg = (2.0 * target_mass_kg) / (1.0 + ratio_nul_nll);
            double nul_kg = ratio_nul_nll * nll_kg;

            res.grading_name = "Custom Grading";
            res.M50_kg = target_mass_kg;
            res.NLL_kg = nll_kg;
            res.NUL_kg = nul_kg;
            res.used_custom_interpolation = true;
            res.custom_family = family;
            res.custom_ratio_nul_nll = ratio_nul_nll;
            res.custom_ratio_note = ratio_note;
        }

        res.ELL_kg = 0.7 * res.NLL_kg;
        res.EUL_kg = 1.5 * res.NUL_kg;
        res.W_mean_kn = (res.M50_kg * g) / 1000.0;

        double reference_weight_kn = res.used_custom_interpolation ? target_weight : res.W_mean_kn;
        res.Dn_rock = std::pow(reference_weight_kn / W_rock_spec, 1.0 / 3.0);
        res.r2 = 2.0 * res.Dn_rock;
        res.f2 = 100.0 * 2.0 * 1.0 * (1.0 - P_rock) / std::pow(res.Dn_rock, 2);
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

        double Nz = params.Number_of_Waves;
        double storm_duration_hr = (Nz * params.Tm) / 3600.0;
        double delta_trunk = (params.Wc / params.Ww) - 1.0;

        // 4. Algorithmic Core (Chegini-Aghtouman / Van der Meer) - TRUNK
        double term_damage = std::pow(params.Nod, k2);
        double term_waves = std::pow(Nz, k3);
        double damage_wave_ratio = term_damage / term_waves;
        
        double scaled_term = k1 * damage_wave_ratio;
        double inv_f = scaled_term + k4;
        
        double steepness_factor = std::pow(s0m, -k5);
        
        double Ns_trunk = inv_f * steepness_factor;

        // 5. Block Sizing (Armor) - TRUNK
        double Dn = params.Hs / (delta_trunk * Ns_trunk);
        double W_trunk = params.Wc * std::pow(Dn, 3);
        double packing_density_trunk = 100.0 * 2.0 * 1.1 * (1.0 - P_cubes) / std::pow(Dn, 2);

        // 6. Hudson Comparative Calculation
        double slope = coeffs.slope_ratio;
        double kd_trunk_equiv = (params.Wc * std::pow(params.Hs, 3)) / (W_trunk * std::pow(delta_trunk, 3) * slope);

        // 7. UNDERLAYER - TRUNK
        UnderlayerResult ul_trunk = calculate_underlayer_params(W_trunk, params.use_en13383, params.custom_family);

        // 8. HEAD CALCULATION (FIXED RATIO 1.5)
        double kd_ratio = KD_RATIO_FIXED;
        double kd_head_derived = kd_trunk_equiv / kd_ratio;
        double delta_head = delta_trunk * std::pow(kd_ratio, 1.0/3.0);
        double Wc_head = params.Ww * (delta_head + 1.0);
        double W_head = W_trunk * (Wc_head / params.Wc);

        // Calculate Stability Number for Head
        double Ns_head = params.Hs / (delta_head * Dn);
        double packing_density_head = 100.0 * 2.0 * 1.1 * (1.0 - P_cubes) / std::pow(Dn, 2);

        // 9. UNDERLAYER - HEAD
        UnderlayerResult ul_head = calculate_underlayer_params(W_head, params.use_en13383, params.custom_family);

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
        results.intermediate.Storm_Duration_hr = storm_duration_hr;
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
        results.final_head.Dn = Dn;
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

    std::string format_report(const FullResults& results) const {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(2);

        auto fmt = [](double val, int prec) {
            std::stringstream stream;
            stream << std::fixed << std::setprecision(prec) << val;
            return stream.str();
        };

        auto field = [&](const std::string& label, const std::string& value, int width = 38) {
            ss << "   " << std::left << std::setw(width) << label << " : " << value << "\n";
        };

        const auto& p = results.inputs;
        const auto& c = results.coefficients;
        const auto& i = results.intermediate;
        const auto& ft = results.final_trunk;
        const auto& ut = results.underlayer_trunk;
        const auto& fh = results.final_head;
        const auto& uh = results.underlayer_head;

        double armor_vol_trunk = (c.type == "Antifer")
            ? 1.0247 * std::pow(ft.dims.H, 3)
            : std::pow(ft.Dn, 3);
        double armor_vol_head = (c.type == "Antifer")
            ? 1.0247 * std::pow(fh.dims.H, 3)
            : std::pow(fh.Dn, 3);

        std::string methodology_name = get_formula_display_name(p.Formula_ID);

        ss << "================================================================================\n"
           << "    TECHNICAL REPORT: BREAKWATER ARMOR & UNDERLAYER DESIGN\n"
           << "================================================================================\n"
           << "Methodology: " << methodology_name << "\n"
           << "--------------------------------------------------------------------------------\n";

        ss << "1. INPUT PARAMETERS\n";
        field("Hs (Significant Wave Height)", fmt(p.Hs, 2) + " m");
        field("Tm (Mean Wave Period)", fmt(p.Tm, 2) + " s");
        field("Number of waves (Nz)", fmt(p.Number_of_Waves, 0));
        field("Nod (Damage)", fmt(p.Nod, 2));
        field("Wc Trunk (Concrete Spec. Weight)", fmt(p.Wc, 2) + " kN/m3");
        field("Ww (Water Specific Weight)", fmt(p.Ww, 2) + " kN/m3");
        field("Relative Density D=(Wc/Ww)-1", fmt(i.delta, 4));
        field("Structure Slope (TRUNK & HEAD)", fmt(c.slope_ratio, 1) + ":1");
        field("Porosity (Cubes)", fmt(results.P_cubes * 100, 0) + "%");
        field("Porosity (Rock Layer)", fmt(results.P_rock * 100, 0) + "%");
        ss << "--------------------------------------------------------------------------------\n";

        ss << "2. INTERMEDIATE PARAMETERS\n";
        field("Wave Length (L0)", fmt(i.L0, 2) + " m");
        field("Wave number (k0 = 2*pi/L0)", fmt(i.k0, 4));
        field("Wave steepness (s0m = Hs/L0)", fmt(i.s0m, 4));
        field("Storm Duration (h)", fmt(i.Storm_Duration_hr, 3) + " h");
        field("Stability Number TRUNK (Ns)", fmt(i.Ns_trunk, 4));
        field("Stability Number HEAD (Ns)", fmt(fh.Ns, 4));
        ss << "--------------------------------------------------------------------------------\n";

        ss << "3. ARMOR LAYER RESULTS - TRUNK\n";
        field("BLOCK WEIGHT (W)", fmt(ft.W, 2) + " kN");
        field("Mass (t)", fmt(ft.Mass_tonnes, 2) + " t");
        field("Nominal Dimension (Dn)", fmt(ft.Dn, 3) + " m");
        if (c.type == "Cubes") {
            field("Volume = Dn^3 (V)", fmt(armor_vol_trunk, 3) + " m3");
        } else {
            field("Volume = 1.0247 * H^3 (V)", fmt(armor_vol_trunk, 3) + " m3");
            field("Block Height (H)", fmt(ft.dims.H, 3) + " m");
            field("Block Top Width (B)", fmt(ft.dims.B, 3) + " m");
            field("Block Base Width (A)", fmt(ft.dims.A, 3) + " m");
        }
        field("KD_TRUNK (Equivalent)", fmt(ft.Kd, 2));
        field("Double Layer Thickness (r1)", fmt(ft.r1, 2) + " m");
        field("Packing Density, d [units/100m2]", fmt(ft.packing_density, 2));
        ss << "\n";

        ss << "4. UNDERLAYER RESULTS - TRUNK\n";
        field("Theoretical Target (W/10)", fmt(ut.target_W, 2) + " kN (" + fmt(ut.target_M50_kg, 1) + " kg)");
        field("Adopted rock grading", ut.grading_name);
        field("Representative M50", fmt(ut.M50_kg, 1) + " kg");
        field("Nominal lower limit (NLL)", fmt(ut.NLL_kg, 1) + " kg");
        field("Nominal upper limit (NUL)", fmt(ut.NUL_kg, 1) + " kg");
        field("Extreme lower limit (ELL)", fmt(ut.ELL_kg, 1) + " kg");
        field("Extreme upper limit (EUL)", fmt(ut.EUL_kg, 1) + " kg");
        field("Nominal Dimension (Dn_rock)", fmt(ut.Dn_rock, 3) + " m");
        field("Double Layer Thickness (r2)", fmt(ut.r2, 2) + " m");
        field("Packing Density, f2 [rocks/100m2]", fmt(ut.f2, 2));
        if (ut.used_custom_interpolation) {
            field("Custom family basis", ut.custom_family);
            field("Custom ratio R=NUL/NLL", fmt(ut.custom_ratio_nul_nll, 3));
        }
        ss << std::string(80, '-') << "\n";

        ss << "5. ARMOR LAYER RESULTS - HEAD (High Density)\n"
           << "   *Maintains same Dn and slope as Trunk*\n";
        field("Stability Ratio (Kd_T/Kd_H)", fmt(fh.Kd_Ratio, 2));
        field("Nominal Dimension (Dn)", fmt(fh.Dn, 3) + " m");
        if (c.type == "Cubes") {
            field("Volume = Dn^3 (V)", fmt(armor_vol_head, 3) + " m3");
        } else {
            field("Volume = 1.0247 * H^3 (V)", fmt(armor_vol_head, 3) + " m3");
            field("Block Height (H)", fmt(fh.dims.H, 3) + " m");
            field("Block Top Width (B)", fmt(fh.dims.B, 3) + " m");
            field("Block Base Width (A)", fmt(fh.dims.A, 3) + " m");
        }
        field("KD_HEAD (Equivalent)", fmt(fh.Kd, 2));
        field("Required Concrete Density (Wc)", fmt(fh.Wc_Required, 2) + " kN/m3");
        field("BLOCK WEIGHT (W)", fmt(fh.W, 2) + " kN");
        field("Mass (t)", fmt(fh.Mass_tonnes, 2) + " t");
        field("Packing Density, d [units/100m2]", fmt(fh.packing_density, 2));
        ss << "\n";

        ss << "6. UNDERLAYER RESULTS - HEAD\n";
        field("Theoretical Target (W/10)", fmt(uh.target_W, 2) + " kN (" + fmt(uh.target_M50_kg, 1) + " kg)");
        field("Adopted rock grading", uh.grading_name);
        field("Representative M50", fmt(uh.M50_kg, 1) + " kg");
        field("Nominal lower limit (NLL)", fmt(uh.NLL_kg, 1) + " kg");
        field("Nominal upper limit (NUL)", fmt(uh.NUL_kg, 1) + " kg");
        field("Extreme lower limit (ELL)", fmt(uh.ELL_kg, 1) + " kg");
        field("Extreme upper limit (EUL)", fmt(uh.EUL_kg, 1) + " kg");
        field("Nominal Dimension (Dn_rock)", fmt(uh.Dn_rock, 3) + " m");
        field("Double Layer Thickness (r2)", fmt(uh.r2, 2) + " m");
        field("Packing Density, f2 [rocks/100m2]", fmt(uh.f2, 2));
        if (uh.used_custom_interpolation) {
            field("Custom family basis", uh.custom_family);
            field("Custom ratio R=NUL/NLL", fmt(uh.custom_ratio_nul_nll, 3));
        }
        ss << std::string(80, '=') << "\n";

        return ss.str();
    }

    void generate_report_file(const FullResults& results, std::string filepath="output.txt") {
        std::string report_content = format_report(results);

        std::ofstream outfile(filepath, std::ios::app | std::ios::binary);
        if (outfile.is_open()) {
            outfile << report_content;
            outfile.close();
        }

        std::cout << report_content;
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

        size_t startpos = input_str.find_first_not_of(" \n\r\t");
        if (startpos == std::string::npos) {
            return default_val; // Empty or whitespace only
        }

        try {
            return std::stod(input_str.substr(startpos));
        } catch (...) {
            return default_val;
        }
    };

    auto trim_copy = [](std::string s) -> std::string {
        size_t start = s.find_first_not_of(" \n\r\t");
        if (start == std::string::npos) return "";
        size_t end = s.find_last_not_of(" \n\r\t");
        return s.substr(start, end - start + 1);
    };

    auto to_upper_copy = [](std::string s) -> std::string {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::toupper(c); });
        return s;
    };

    auto parse_bool = [&](const std::string& s, bool default_val) -> bool {
        std::string t = to_upper_copy(trim_copy(s));
        if (t.empty()) return default_val;
        if (t == "1" || t == "TRUE" || t == "T" || t == "YES" || t == "Y") return true;
        if (t == "0" || t == "FALSE" || t == "F" || t == "NO" || t == "N") return false;
        return default_val;
    };

    auto get_text_param = [](const std::string& prompt, const std::string& default_val) -> std::string {
        std::cout << prompt << " [" << default_val << "]: ";
        std::string input_str;
        std::getline(std::cin, input_str);

        size_t startpos = input_str.find_first_not_of(" \n\r\t");
        if (startpos == std::string::npos) {
            return default_val;
        }

        size_t endpos = input_str.find_last_not_of(" \n\r\t");
        return input_str.substr(startpos, endpos - startpos + 1);
    };

    Inputs user_inputs = calc.defaults;
    int formula_id = 1;

    // Check if command line arguments are provided
    if (argc >= 7) {
        // Parse CLI arguments
        try {
            user_inputs.Hs = std::stod(argv[1]);
            user_inputs.Tm = std::stod(argv[2]);
            user_inputs.Number_of_Waves = std::stod(argv[3]);
            user_inputs.Nod = std::stod(argv[4]);
            user_inputs.Wc = std::stod(argv[5]);
            formula_id = std::stoi(argv[6]);

            if (formula_id < 1 || formula_id > 4) {
                 std::cerr << "Error: Formula ID must be 1-4. Using default (1)." << std::endl;
                 formula_id = 1;
            }
            user_inputs.Formula_ID = formula_id;

            if (argc >= 8) {
                user_inputs.use_en13383 = parse_bool(argv[7], calc.defaults.use_en13383);
            }
            if (argc >= 9) {
                user_inputs.custom_family = to_upper_copy(trim_copy(argv[8]));
                if (user_inputs.custom_family != "AUTO" &&
                    user_inputs.custom_family != "HMA" &&
                    user_inputs.custom_family != "LMA" &&
                    user_inputs.custom_family != "CP") {
                    user_inputs.custom_family = "AUTO";
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "Error parsing command line arguments: " << e.what() << std::endl;
            std::cerr << "Usage: " << argv[0] << " [Hs] [Tm] [NumberOfWaves] [Nod] [Wc] [FormulaID] [UseEN13383] [CustomFamily]" << std::endl;
            return 1;
        }
    } else {
        // Interactive Mode
        std::cout << "\n--- COASTAL PROTECTION BLOCK CALCULATOR (TRUNK & HEAD) ---" << std::endl;
        std::cout << "1. Simple Cubes (Slope 2.0:1) - Van der Meer" << std::endl;
        std::cout << "2. Simple Cubes (Slope 1.5:1) - Van der Meer" << std::endl;
        std::cout << "3. Antifer (Slope 2.0:1) - Chegini" << std::endl;
        std::cout << "4. Antifer (Slope 1.5:1) - Chegini" << std::endl;

        std::cout << "\nOption [1-4]: ";
        std::string selection;
        std::getline(std::cin, selection);
        
        // Safely trim whitespace
        size_t endpos = selection.find_last_not_of(" \n\r\t");
        if (endpos != std::string::npos) {
            selection.erase(endpos + 1);
        } else {
            selection.clear();
        }

        // Set formula ID, default to 1 if user just pressed ENTER or typed nonsense
        if (selection == "1" || selection == "2" || selection == "3" || selection == "4") {
            formula_id = std::stoi(selection);
        } else {
            formula_id = 1; 
        }
        user_inputs.Formula_ID = formula_id;

        // Unconditionally ask for parameters (matching Fortran behavior)
        std::cout << "\n--- Enter Parameters (Press ENTER for Default) ---" << std::endl;
        user_inputs.Hs = get_param("Hs (m)", calc.defaults.Hs);
        user_inputs.Tm = get_param("Tm (s)", calc.defaults.Tm);
        user_inputs.Number_of_Waves = get_param("Number of waves (Nz)", calc.defaults.Number_of_Waves);
        user_inputs.Nod = get_param("Nod (Damage)", calc.defaults.Nod);
        user_inputs.Wc = get_param("Concrete Weight Trunk (kN/m3)", calc.defaults.Wc);

        std::string use_en_str = get_text_param("Use standard EN 13383 underlayer grading? [true/false]", calc.defaults.use_en13383 ? "true" : "false");
        user_inputs.use_en13383 = parse_bool(use_en_str, calc.defaults.use_en13383);

        std::string family = get_text_param("Custom underlayer family [AUTO/HMA/LMA/CP]", calc.defaults.custom_family);
        family = to_upper_copy(trim_copy(family));
        if (family == "AUTO" || family == "HMA" || family == "LMA" || family == "CP") {
            user_inputs.custom_family = family;
        } else {
            user_inputs.custom_family = "AUTO";
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
