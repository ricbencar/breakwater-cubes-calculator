// ======================================================================================
// PROGRAM DESCRIPTION, MULTILINGUAL CAPABILITY & METHODOLOGY
// ======================================================================================
//
// 1. PURPOSE:
//    This software performs preliminary hydraulic sizing of rubble-mound
//    breakwater armor made with artificial concrete units and the associated
//    rock underlayers.
//
//    The armor-unit families currently implemented in the code are:
//    - Simple Cubes
//    - Antifer Blocks
//
//    The program dimensions two structural zones:
//    - The Trunk (main longitudinal section)
//    - The Head (the exposed roundhead / end section)
//
//    The implemented design philosophy is iso-geometric, but not iso-weight:
//    - the nominal dimension is kept the same at Trunk and Head;
//    - the armor-unit geometry is kept the same at Trunk and Head;
//    - the slope is kept the same at Trunk and Head;
//    - the Head unit weight increases through higher required concrete specific
//      weight rather than through larger unit size.
//
//    This GUI version is multilingual and supports three interface languages:
//    - English
//    - Portuguese
//    - French
//
//    The selected language changes the GUI labels, dropdown texts, messages,
//    and generated technical report. The hydraulic calculations, constants,
//    formulas, and sizing logic remain exactly the same regardless of the
//    selected language.
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
//    - automatic grading selection from an embedded EN 13383 mass-based
//      grading table.
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
//       The nominal dimension is then:
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
//       - the same nominal dimension is maintained;
//       - Head units are heavier because concrete specific weight increases.
//
//    d. Underlayer selection:
//       The theoretical underlayer target is one tenth of the armor-unit weight:
//
//           W_target = W_armor / 10
//
//       The target mass is compared against the embedded EN 13383 grading
//       limits. A grading is selected only if the target mass is strictly
//       contained within the nominal lower and upper limits.
//
//       If more than one grading contains the target mass, the program selects
//       the tightest range, i.e. the grading with the smallest:
//
//           NUL - NLL
//
//       If no grading strictly contains the target mass, the program falls back
//       to the first grading in the internal grading database.
//
// 4. DEFAULT INPUTS, CONSTANTS, AND CURRENT FORMULA SETS USED IN THE CODE:
//
//    a. Default user inputs currently used by the GUI:
//       - Hs         = 11.0 m
//       - Tm         = 11.9 s
//       - Nz         = 3000
//       - Nod        = 0.5
//       - Wc         = 27.48 kN/m3
//       - Ww         = 10.05 kN/m3
//       - Formula_ID = 1
//       - Language   = English (default startup language unless changed in code)
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
//       2. Van Der Meer (1988a) - Simple Cubes (Slope 1.5:1)
//       3. Chegini-Aghtouman (2006) - Antifer (Slope 2:1)
//       4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)
//
// 5. OUTPUTS PRODUCED:
//    The generated report contains, for both Trunk and Head:
//    - stability number;
//    - nominal dimension;
//    - unit weight and unit mass;
//    - equivalent Kd value;
//    - armor layer thickness;
//    - packing density;
//    - geometric dimensions reported by the code;
//    - adopted underlayer grading and derived rock parameters.
//
//    Volume reporting follows the current GUI code logic:
//
//    - for Simple Cubes:
//
//          V = Dn^3
//
//    - for Antifer blocks:
//
//          V = 1.0247 * H^3
//
//    The GUI displays the report in the output window and appends the same
//    report to output.txt in the execution folder.
//
//    The language used in the report matches the language currently selected
//    in the GUI. Numerical values and engineering results are identical in all
//    languages; only the displayed text changes.
//
// 6. MULTILINGUAL USER INTERFACE CAPABILITY:
//
//    The code contains an internal language dictionary used to centralize all
//    user-visible text strings. This includes, among others:
//    - input labels;
//    - button captions;
//    - formula names shown in the dropdown list;
//    - warnings and error messages;
//    - report headings and engineering output descriptions.
//
//    The multilingual implementation is interface-level only. It does not alter
//    the computational kernel, formula coefficients, default constants, or
//    engineering methodology.
//
//    The purpose of the dictionary-based approach is to:
//    - avoid maintaining separate source files for each language;
//    - ensure calculation consistency across all languages;
//    - simplify future translation updates or addition of new languages.
//
// 7. TECHNICAL BIBLIOGRAPHY & REFERENCES:
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
// 8. COMPILATION INSTRUCTIONS:
//
//    Recommended MinGW-w64 build on Windows:
//
//    g++ -O3 -march=native -std=c++17 -municode -mwindows -static
//    -static-libgcc -static-libstdc++ -o breakwater_calculator_gui.exe
//    breakwater_calculator_gui.cpp
//
//    Meaning of the main compilation flags:
//    - -O3               -> high optimization level;
//    - -march=native     -> optimize for the local CPU;
//    - -std=c++17        -> compile as C++17;
//    - -municode         -> enable Win32 Unicode entry point;
//    - -mwindows         -> build a Windows GUI application without console;
//    - -static           -> prefer static linking;
//    - -static-libgcc    -> statically link GCC runtime;
//    - -static-libstdc++ -> statically link C++ runtime.
//
//    If full static linking fails on the local toolchain, a fallback build is:
//
//    g++ -O3 -march=native -std=c++17 -municode -mwindows -o
//    breakwater_calculator_gui.exe breakwater_calculator_gui.cpp
//
// 9. EXECUTION:
//
//    1. Launch breakwater_calculator_gui.exe.
//    2. Select the desired interface language:
//       - English
//       - Portuguese
//       - French
//    3. Enter Hs, Tm, Number of Waves (Nz), Nod, and Wc.
//    4. Select the required formula / armor type from the dropdown list.
//    5. Click "Calculate Armor".
//    6. The technical report is displayed in the output window.
//    7. The same report is appended to output.txt in the execution folder,
//       using the currently selected GUI language.
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
#include <cstring>
#include <stdexcept>

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
        
        formulas[1] = { L"Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)", "Cubes", 2.0, 7.374304189198, 0.4, 0.3, 1.100642416298, 0.1 };
        formulas[2] = { L"Van Der Meer (1988a) - Simple Cubes (Slope 1.5:1)", "Cubes", 1.5, 6.7, 0.4, 0.3, 1.0, 0.1 };
        formulas[3] = { L"Chegini-Aghtouman (2006) - Antifer (Slope 2:1)", "Antifer", 2.0, 6.138, 0.443, 0.276, 1.164, 0.07 };
        formulas[4] = { L"Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)", "Antifer", 1.5, 6.951, 0.443, 0.291, 1.082, 0.082 };

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

        defaults = {
            11.0,   // Hs
            11.9,   // Tm
            3000.0, // Number_of_Waves
            0.5,    // Nod
            27.48,  // Wc
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
};


enum class Language {
    English = 0,
    Portuguese = 1,
    French = 2
};

struct LangPack {
    std::wstring window_title;
    std::wstring input_error_title;
    std::wstring input_error_positive;
    std::wstring calc_error_title;

    std::wstring language_label;
    std::wstring hs_label;
    std::wstring tm_label;
    std::wstring nz_label;
    std::wstring nod_label;
    std::wstring wc_label;
    std::wstring formula_label;
    std::wstring button_calculate;

    std::wstring report_title;
    std::wstring methodology;
    std::wstring section_input;
    std::wstring section_intermediate;
    std::wstring section_armor_trunk;
    std::wstring section_underlayer_trunk;
    std::wstring section_armor_head;
    std::wstring section_underlayer_head;

    std::wstring input_hs;
    std::wstring input_tm;
    std::wstring input_nz;
    std::wstring input_nod;
    std::wstring input_wc;
    std::wstring input_ww;
    std::wstring input_delta;
    std::wstring input_slope;
    std::wstring input_porosity_cubes;
    std::wstring input_porosity_rock;

    std::wstring inter_l0;
    std::wstring inter_k0;
    std::wstring inter_s0m;
    std::wstring inter_storm_duration;
    std::wstring inter_ns_trunk;
    std::wstring inter_ns_head;

    std::wstring armor_weight;
    std::wstring armor_mass;
    std::wstring armor_dn;
    std::wstring armor_vol_dn3;
    std::wstring armor_vol_antifer;
    std::wstring armor_h;
    std::wstring armor_b;
    std::wstring armor_a;
    std::wstring armor_kd_trunk;
    std::wstring armor_kd_head;
    std::wstring armor_r1;
    std::wstring armor_packing;

    std::wstring under_target;
    std::wstring under_grading;
    std::wstring under_m50;
    std::wstring under_nll;
    std::wstring under_nul;
    std::wstring under_ell;
    std::wstring under_eul;
    std::wstring under_dn;
    std::wstring under_r2;
    std::wstring under_f2;

    std::wstring head_note;
    std::wstring head_ratio;
    std::wstring head_required_wc;
};

static Language g_currentLanguage = Language::English;
static bool g_hasLastResults = false;
static FullResults g_lastResults = {};

const LangPack& GetLangPack(Language lang) {
    static const LangPack en = {
        L"Breakwater Armor Calculator",
        L"Input Error",
        L"Please enter positive values for hydraulic parameters.",
        L"Calculation Error",

        L"Language:",
        L"Wave Height, Hs (m)",
        L"Mean Period, Tm (m)",
        L"Number of Waves (Nz):",
        L"Nod (Damage):",
        L"Wc Concrete (kN/m3):",
        L"Formula / Slope:",
        L"Calculate Armor",

        L"TECHNICAL REPORT: BREAKWATER ARMOR & UNDERLAYER DESIGN",
        L"Methodology:",
        L"1. INPUT PARAMETERS",
        L"2. INTERMEDIATE PARAMETERS",
        L"3. ARMOR LAYER RESULTS - TRUNK",
        L"4. UNDERLAYER RESULTS - TRUNK",
        L"5. ARMOR LAYER RESULTS - HEAD (High Density)",
        L"6. UNDERLAYER RESULTS - HEAD",

        L"Hs (Significant Wave Height)",
        L"Tm (Mean Wave Period)",
        L"Number of waves (Nz)",
        L"Nod (Damage)",
        L"Wc Trunk (Concrete Spec. Weight)",
        L"Ww (Water Specific Weight)",
        L"Relative Density D=(Wc/Ww)-1",
        L"Structure Slope (TRUNK & HEAD)",
        L"Porosity (Cubes)",
        L"Porosity (Rock Layer)",

        L"Wave Length (L0)",
        L"Wave number (k0 = 2*pi/L0)",
        L"Wave steepness (s0m = Hs/L0)",
        L"Storm Duration (h)",
        L"Stability Number TRUNK (Ns)",
        L"Stability Number HEAD (Ns)",

        L"BLOCK WEIGHT (W)",
        L"Mass (t)",
        L"Nominal Dimension (Dn)",
        L"Volume = Dn^3 (V)",
        L"Volume = 1.0247 * H^3 (V)",
        L"Block Height (H)",
        L"Block Top Width (B)",
        L"Block Base Width (A)",
        L"KD_TRUNK (Equivalent)",
        L"KD_HEAD (Equivalent)",
        L"Double Layer Thickness (r1)",
        L"Packing Density, d [units/100m2]",

        L"Theoretical Target (W/10)",
        L"Adopted rock grading",
        L"Representative M50",
        L"Nominal lower limit (NLL)",
        L"Nominal upper limit (NUL)",
        L"Extreme lower limit (ELL)",
        L"Extreme upper limit (EUL)",
        L"Nominal Dimension (Dn_rock)",
        L"Double Layer Thickness (r2)",
        L"Packing Density, f2 [rocks/100m2]",

        L"*Maintains same Dn and slope as Trunk*",
        L"Stability Ratio (Kd_T/Kd_H)",
        L"Required Concrete Density (Wc)"
    };

    static const LangPack pt = {
        L"Calculadora do Manto de Quebra-Mares",
        L"Erro de entrada",
        L"Introduza valores positivos para os parâmetros hidráulicos.",
        L"Erro de cálculo",

        L"Idioma:",
        L"Altura de Onda, Hs (m)",
        L"Período Médio, Tm (m)",
        L"Número de ondas (Nz):",
        L"Nod (Dano):",
        L"Wc do betão (kN/m3):",
        L"Fórmula / Talude:",
        L"Calcular Manto",

        L"RELATÓRIO TÉCNICO: DIMENSIONAMENTO DO MANTO E DA SUBCAMADA DO QUEBRA-MAR",
        L"Metodologia:",
        L"1. PARÂMETROS DE ENTRADA",
        L"2. PARÂMETROS INTERMÉDIOS",
        L"3. RESULTADOS DA CAMADA DO MANTO - TRONCO",
        L"4. RESULTADOS DA SUBCAMADA - TRONCO",
        L"5. RESULTADOS DA CAMADA DO MANTO - CABEÇA (ALTA DENSIDADE)",
        L"6. RESULTADOS DA SUBCAMADA - CABEÇA",

        L"Hs (Altura significativa da onda)",
        L"Tm (Período médio da onda)",
        L"Número de ondas (Nz)",
        L"Nod (Dano)",
        L"Wc do tronco (peso esp. do betão)",
        L"Ww (peso específico da água)",
        L"Densidade relativa D=(Wc/Ww)-1",
        L"Talude da estrutura (TRONCO/CABEÇA)",
        L"Porosidade (Cubos)",
        L"Porosidade (Camada de enrocamento)",

        L"Comprimento de onda (L0)",
        L"Número de onda (k0 = 2*pi/L0)",
        L"Declividade da onda (s0m = Hs/L0)",
        L"Duração da tempestade (h)",
        L"Número de estabilidade TRONCO (Ns)",
        L"Número de estabilidade CABEÇA (Ns)",

        L"PESO DO BLOCO (W)",
        L"Massa (t)",
        L"Dimensão nominal (Dn)",
        L"Volume = Dn^3 (V)",
        L"Volume = 1.0247 * H^3 (V)",
        L"Altura do bloco (H)",
        L"Largura do bloco no topo (B)",
        L"Largura do bloco na base (A)",
        L"KD_TRONCO (equivalente)",
        L"KD_CABEÇA (equivalente)",
        L"Espessura da camada dupla (r1)",
        L"Densidade de colocação, d [un/100m2]",

        L"Alvo teórico (W/10)",
        L"Classe de enrocamento adoptada",
        L"M50 representativo",
        L"Limite nominal inferior (NLL)",
        L"Limite nominal superior (NUL)",
        L"Limite extremo inferior (ELL)",
        L"Limite extremo superior (EUL)",
        L"Diâmetro nominal (Dn_enrocam)",
        L"Espessura da camada dupla (r2)",
        L"Densidade, f2 [enrocam/100m2]",

        L"*Mantém o mesmo Dn e talude do tronco*",
        L"Razão de estabilidade (Kd_T/Kd_H)",
        L"Densidade de betão requerida (Wc)"
    };

    static const LangPack fr = {
        L"Calculateur de carapace de brise-lames",
        L"Erreur de saisie",
        L"Veuillez saisir des valeurs positives pour les paramètres hydrauliques.",
        L"Erreur de calcul",

        L"Langue :",
        L"Hauteur de Vague, Hs (m)",
        L"Période Moyenne, Tm (m)",
        L"Nombre de vagues (Nz) :",
        L"Nod (Endommagement) :",
        L"Wc béton (kN/m3) :",
        L"Formule / Talus :",
        L"Calculer la carapace",

        L"RAPPORT TECHNIQUE : DIMENSIONNEMENT DE LA CARAPACE ET DE LA SOUS-COUCHE DU BRISE-LAMES",
        L"Méthodologie :",
        L"1. PARAMÈTRES D'ENTRÉE",
        L"2. PARAMÈTRES INTERMÉDIAIRES",
        L"3. RÉSULTATS DE LA CARAPACE - TRONC",
        L"4. RÉSULTATS DE LA SOUS-COUCHE - TRONC",
        L"5. RÉSULTATS DE LA CARAPACE - TÊTE (HAUTE DENSITÉ)",
        L"6. RÉSULTATS DE LA SOUS-COUCHE - TÊTE",

        L"Hs (Hauteur significative de vague)",
        L"Tm (Période moyenne de vague)",
        L"Nombre de vagues (Nz)",
        L"Nod (Endommagement)",
        L"Wc du tronc (poids volumique du béton)",
        L"Ww (poids volumique de l'eau)",
        L"Densité relative D=(Wc/Ww)-1",
        L"Talus de l'ouvrage (TRONC/TÊTE)",
        L"Porosité (Cubes)",
        L"Porosité (Couche d'enrochement)",

        L"Longueur d'onde (L0)",
        L"Nombre d'onde (k0 = 2*pi/L0)",
        L"Raideur de vague (s0m = Hs/L0)",
        L"Durée de tempête (h)",
        L"Nombre de stabilité TRONC (Ns)",
        L"Nombre de stabilité TÊTE (Ns)",

        L"POIDS DU BLOC (W)",
        L"Masse (t)",
        L"Dimension nominal (Dn)",
        L"Volume = Dn^3 (V)",
        L"Volume = 1.0247 * H^3 (V)",
        L"Hauteur du bloc (H)",
        L"Largeur supérieure du bloc (B)",
        L"Largeur de base du bloc (A)",
        L"KD_TRONC (équivalent)",
        L"KD_TÊTE (équivalent)",
        L"Épaisseur de la double couche (r1)",
        L"Densité de pose, d [u/100m2]",

        L"Cible théorique (W/10)",
        L"Classe d'enrochement adoptée",
        L"M50 représentatif",
        L"Limite nominale inférieure (NLL)",
        L"Limite nominale supérieure (NUL)",
        L"Limite extrême inférieure (ELL)",
        L"Limite extrême supérieure (EUL)",
        L"Diamètre nominal (Dn_enroch)",
        L"Épaisseur de la double couche (r2)",
        L"Densité, f2 [enroch/100m2]",

        L"*Conserve le même Dn et le même talus que le tronc*",
        L"Rapport de stabilité (Kd_T/Kd_H)",
        L"Poids volumique de béton requis (Wc)"
    };

    switch (lang) {
        case Language::Portuguese: return pt;
        case Language::French:     return fr;
        case Language::English:
        default:                   return en;
    }
}

std::wstring GetFormulaDisplayName(int formula_id, Language lang) {
    switch (lang) {
        case Language::Portuguese:
            switch (formula_id) {
                case 1: return L"Van der Meer (1988a) - Cubos Simples (Talude 2.0:1)";
                case 2: return L"Van der Meer (1988a) - Cubos Simples (Talude 1.5:1)";
                case 3: return L"Chegini-Aghtouman (2006) - Antifer (Talude 2.0:1)";
                case 4: return L"Chegini-Aghtouman (2006) - Antifer (Talude 1.5:1)";
                default: return L"Fórmula desconhecida";
            }
        case Language::French:
            switch (formula_id) {
                case 1: return L"Van der Meer (1988a) - Cubes Simples (Talus 2.0:1)";
                case 2: return L"Van der Meer (1988a) - Cubes Simples (Talus 1.5:1)";
                case 3: return L"Chegini-Aghtouman (2006) - Antifer (Talus 2.0:1)";
                case 4: return L"Chegini-Aghtouman (2006) - Antifer (Talus 1.5:1)";
                default: return L"Formule inconnue";
            }
        case Language::English:
        default:
            switch (formula_id) {
                case 1: return L"Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)";
                case 2: return L"Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)";
                case 3: return L"Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)";
                case 4: return L"Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)";
                default: return L"Unknown formula";
            }
    }
}

std::wstring GetFormulaComboLabel(int formula_id, Language lang) {
    switch (lang) {
        case Language::Portuguese:
            switch (formula_id) {
                case 1: return L"1. Cubos (Talude 2.0:1) - Van der Meer";
                case 2: return L"2. Cubos (Talude 1.5:1) - Van der Meer";
                case 3: return L"3. Antifer (Talude 2.0:1) - Chegini";
                case 4: return L"4. Antifer (Talude 1.5:1) - Chegini";
                default: return L"";
            }
        case Language::French:
            switch (formula_id) {
                case 1: return L"1. Cubes (Talus 2.0:1) - Van der Meer";
                case 2: return L"2. Cubes (Talus 1.5:1) - Van der Meer";
                case 3: return L"3. Antifer (Talus 2.0:1) - Chegini";
                case 4: return L"4. Antifer (Talus 1.5:1) - Chegini";
                default: return L"";
            }
        case Language::English:
        default:
            switch (formula_id) {
                case 1: return L"1. Cubes (Slope 2.0:1) - Van der Meer";
                case 2: return L"2. Cubes (Slope 1.5:1) - Van der Meer";
                case 3: return L"3. Antifer (Slope 2.0:1) - Chegini";
                case 4: return L"4. Antifer (Slope 1.5:1) - Chegini";
                default: return L"";
            }
    }
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

static std::string wide_to_utf8(const std::wstring& text) {
    if (text.empty()) {
        return std::string();
    }

    int size_needed = WideCharToMultiByte(
        CP_UTF8,
        0,
        text.c_str(),
        static_cast<int>(text.size()),
        nullptr,
        0,
        nullptr,
        nullptr
    );

    if (size_needed <= 0) {
        throw std::runtime_error("UTF-8 conversion failed");
    }

    std::string utf8(static_cast<size_t>(size_needed), '\0');
    int converted = WideCharToMultiByte(
        CP_UTF8,
        0,
        text.c_str(),
        static_cast<int>(text.size()),
        utf8.data(),
        size_needed,
        nullptr,
        nullptr
    );

    if (converted != size_needed) {
        throw std::runtime_error("UTF-8 conversion failed");
    }

    return utf8;
}

std::wstring format_report(const FullResults& results, Language lang) {
    const LangPack& t = GetLangPack(lang);

    std::wstringstream ss;
    ss << std::fixed << std::setprecision(2);

    auto fmt = [](double val, int prec) {
        std::wstringstream stream;
        stream << std::fixed << std::setprecision(prec) << val;
        return stream.str();
    };

    auto field = [&](const std::wstring& label, const std::wstring& value, int width = 38) {
        ss << L"   " << std::left << std::setw(width) << label << L" : " << value << L"\n";
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

    std::wstring methodology_name = GetFormulaDisplayName(p.Formula_ID, lang);
    std::wstring ut_grading_wide(ut.grading_name.begin(), ut.grading_name.end());
    std::wstring uh_grading_wide(uh.grading_name.begin(), uh.grading_name.end());

    ss << L"================================================================================\n"
       << L"    " << t.report_title << L"\n"
       << L"================================================================================\n"
       << t.methodology << L" " << methodology_name << L"\n"
       << L"--------------------------------------------------------------------------------\n";

    ss << t.section_input << L"\n";
    field(t.input_hs, fmt(p.Hs, 2) + L" m");
    field(t.input_tm, fmt(p.Tm, 2) + L" s");
    field(t.input_nz, fmt(p.Number_of_Waves, 0));
    field(t.input_nod, fmt(p.Nod, 2));
    field(t.input_wc, fmt(p.Wc, 2) + L" kN/m3");
    field(t.input_ww, fmt(p.Ww, 2) + L" kN/m3");
    field(t.input_delta, fmt(i.delta, 4));
    field(t.input_slope, fmt(c.slope_ratio, 1) + L":1");
    field(t.input_porosity_cubes, fmt(results.P_cubes * 100, 0) + L"%");
    field(t.input_porosity_rock, fmt(results.P_rock * 100, 0) + L"%");
    ss << L"--------------------------------------------------------------------------------\n";

    ss << t.section_intermediate << L"\n";
    field(t.inter_l0, fmt(i.L0, 2) + L" m");
    field(t.inter_k0, fmt(i.k0, 4));
    field(t.inter_s0m, fmt(i.s0m, 4));
    field(t.inter_storm_duration, fmt(i.Storm_Duration_hr, 3) + L" h");
    field(t.inter_ns_trunk, fmt(i.Ns_trunk, 4));
    field(t.inter_ns_head, fmt(fh.Ns, 4));
    ss << L"--------------------------------------------------------------------------------\n";

    ss << t.section_armor_trunk << L"\n";
    field(t.armor_weight, fmt(ft.W, 2) + L" kN");
    field(t.armor_mass, fmt(ft.Mass_tonnes, 2) + L" t");
    field(t.armor_dn, fmt(ft.Dn, 3) + L" m");
    if (c.type == "Cubes") {
        field(t.armor_vol_dn3, fmt(armor_vol_trunk, 3) + L" m3");
    } else {
        field(t.armor_vol_antifer, fmt(armor_vol_trunk, 3) + L" m3");
        field(t.armor_h, fmt(ft.dims.H, 3) + L" m");
        field(t.armor_b, fmt(ft.dims.B, 3) + L" m");
        field(t.armor_a, fmt(ft.dims.A, 3) + L" m");
    }
    field(t.armor_kd_trunk, fmt(ft.Kd, 2));
    field(t.armor_r1, fmt(ft.r1, 2) + L" m");
    field(t.armor_packing, fmt(ft.packing_density, 2));
    ss << L"\n";

    ss << t.section_underlayer_trunk << L"\n";
    field(t.under_target, fmt(ut.target_W, 2) + L" kN (" + fmt(ut.target_M50_kg, 1) + L" kg)");
    field(t.under_grading, ut_grading_wide);
    field(t.under_m50, fmt(ut.M50_kg, 1) + L" kg");
    field(t.under_nll, fmt(ut.NLL_kg, 1) + L" kg");
    field(t.under_nul, fmt(ut.NUL_kg, 1) + L" kg");
    field(t.under_ell, fmt(ut.ELL_kg, 1) + L" kg");
    field(t.under_eul, fmt(ut.EUL_kg, 1) + L" kg");
    field(t.under_dn, fmt(ut.Dn_rock, 3) + L" m");
    field(t.under_r2, fmt(ut.r2, 2) + L" m");
    field(t.under_f2, fmt(ut.f2, 2));
    ss << std::wstring(80, L'-') << L"\n";

    ss << t.section_armor_head << L"\n"
       << L"   " << t.head_note << L"\n";
    field(t.head_ratio, fmt(fh.Kd_Ratio, 2));
    field(t.armor_dn, fmt(fh.Dn, 3) + L" m");
    if (c.type == "Cubes") {
        field(t.armor_vol_dn3, fmt(armor_vol_head, 3) + L" m3");
    } else {
        field(t.armor_vol_antifer, fmt(armor_vol_head, 3) + L" m3");
        field(t.armor_h, fmt(fh.dims.H, 3) + L" m");
        field(t.armor_b, fmt(fh.dims.B, 3) + L" m");
        field(t.armor_a, fmt(fh.dims.A, 3) + L" m");
    }
    field(t.armor_kd_head, fmt(fh.Kd, 2));
    field(t.head_required_wc, fmt(fh.Wc_Required, 2) + L" kN/m3");
    field(t.armor_weight, fmt(fh.W, 2) + L" kN");
    field(t.armor_mass, fmt(fh.Mass_tonnes, 2) + L" t");
    field(t.armor_packing, fmt(fh.packing_density, 2));
    ss << L"\n";

    ss << t.section_underlayer_head << L"\n";
    field(t.under_target, fmt(uh.target_W, 2) + L" kN (" + fmt(uh.target_M50_kg, 1) + L" kg)");
    field(t.under_grading, uh_grading_wide);
    field(t.under_m50, fmt(uh.M50_kg, 1) + L" kg");
    field(t.under_nll, fmt(uh.NLL_kg, 1) + L" kg");
    field(t.under_nul, fmt(uh.NUL_kg, 1) + L" kg");
    field(t.under_ell, fmt(uh.ELL_kg, 1) + L" kg");
    field(t.under_eul, fmt(uh.EUL_kg, 1) + L" kg");
    field(t.under_dn, fmt(uh.Dn_rock, 3) + L" m");
    field(t.under_r2, fmt(uh.r2, 2) + L" m");
    field(t.under_f2, fmt(uh.f2, 2));
    ss << std::wstring(80, L'=') << L"\n";

    return ss.str();
}

// --- FILE OUTPUT HELPER ---
void SaveReportToFile(const std::wstring& report_content) {
    std::ofstream file("output.txt", std::ios::app | std::ios::binary);
    if (!file.is_open()) {
        return;
    }

    std::string utf8 = wide_to_utf8(report_content);
    file.write(utf8.data(), static_cast<std::streamsize>(utf8.size()));
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
#define IDC_COMBO_LANGUAGE 109

HWND hLabelLanguage, hLabelHs, hLabelTm, hLabelDur, hLabelNod, hLabelWc, hLabelFormula;
HWND hEditHs, hEditTm, hEditDur, hEditNod, hEditWc, hComboFormula, hComboLanguage, hButtonCompute, hOutput;
BreakwaterCalculator calculator;

// Store fonts globally to delete them later
HFONT hMonoFont = NULL;
HFONT hUIFont = NULL;

void PopulateFormulaCombo(Language lang) {
    int sel = (int)SendMessageW(hComboFormula, CB_GETCURSEL, 0, 0);
    if (sel < 0 || sel > 3) {
        sel = calculator.defaults.Formula_ID - 1;
    }

    SendMessageW(hComboFormula, CB_RESETCONTENT, 0, 0);
    for (int i = 1; i <= 4; ++i) {
        std::wstring label = GetFormulaComboLabel(i, lang);
        SendMessageW(hComboFormula, CB_ADDSTRING, 0, (LPARAM)label.c_str());
    }
    SendMessageW(hComboFormula, CB_SETCURSEL, (WPARAM)sel, 0);
}

void ApplyLanguage(HWND hwnd) {
    const LangPack& t = GetLangPack(g_currentLanguage);

    SetWindowTextW(hwnd, t.window_title.c_str());
    SetWindowTextW(hLabelLanguage, t.language_label.c_str());
    SetWindowTextW(hLabelHs, t.hs_label.c_str());
    SetWindowTextW(hLabelTm, t.tm_label.c_str());
    SetWindowTextW(hLabelDur, t.nz_label.c_str());
    SetWindowTextW(hLabelNod, t.nod_label.c_str());
    SetWindowTextW(hLabelWc, t.wc_label.c_str());
    SetWindowTextW(hLabelFormula, t.formula_label.c_str());
    SetWindowTextW(hButtonCompute, t.button_calculate.c_str());

    PopulateFormulaCombo(g_currentLanguage);
    SendMessageW(hComboLanguage, CB_SETCURSEL, (WPARAM)g_currentLanguage, 0);

    if (g_hasLastResults) {
        std::wstring report = format_report(g_lastResults, g_currentLanguage);
        std::wstring guiText = fix_newlines_for_edit_control(report);
        SetWindowTextW(hOutput, guiText.c_str());
    }
}

void CreateControls(HWND hwnd) {
    hUIFont = CreateFontW(20, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                          DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                          DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");

    int y = 20;
    int lblW = 230;
    int editX = 245;
    int editW = 110;
    int step = 40;

    auto fmt = [](double val) {
        std::wstringstream ss;
        ss << std::fixed << std::setprecision(2) << val;
        return ss.str();
    };

    auto CreateLabel = [&](const wchar_t* text, int x, int y, int width = -1) -> HWND {
        if (width < 0) width = lblW;
        HWND h = CreateWindowW(L"STATIC", text, WS_CHILD | WS_VISIBLE, x, y, width, 25, hwnd, NULL, NULL, NULL);
        SendMessageW(h, WM_SETFONT, (WPARAM)hUIFont, TRUE);
        return h;
    };

    auto CreateInput = [&](const wchar_t* def, int id, int y) -> HWND {
        HWND h = CreateWindowW(L"EDIT", def, WS_CHILD | WS_VISIBLE | WS_BORDER | ES_AUTOHSCROLL,
                               editX, y, editW, 25, hwnd, (HMENU)(INT_PTR)id, NULL, NULL);
        SendMessageW(h, WM_SETFONT, (WPARAM)hUIFont, TRUE);
        return h;
    };

    hLabelLanguage = CreateLabel(L"Language:", 10, y);
    hComboLanguage = CreateWindowW(L"COMBOBOX", NULL,
        WS_CHILD | WS_VISIBLE | CBS_DROPDOWNLIST | WS_VSCROLL,
        editX, y, editW, 200, hwnd, (HMENU)(INT_PTR)IDC_COMBO_LANGUAGE, NULL, NULL);
    SendMessageW(hComboLanguage, WM_SETFONT, (WPARAM)hUIFont, TRUE);
    SendMessageW(hComboLanguage, CB_ADDSTRING, 0, (LPARAM)L"English");
    SendMessageW(hComboLanguage, CB_ADDSTRING, 0, (LPARAM)L"Português");
    SendMessageW(hComboLanguage, CB_ADDSTRING, 0, (LPARAM)L"Français");
    SendMessageW(hComboLanguage, CB_SETCURSEL, (WPARAM)g_currentLanguage, 0);

    y += step;
    hLabelHs = CreateLabel(L"Wave Height, Hs (m)", 10, y);
    hEditHs = CreateInput(fmt(calculator.defaults.Hs).c_str(), IDC_EDIT_HS, y);

    y += step;
    hLabelTm = CreateLabel(L"Mean Period, Tm (m)", 10, y);
    hEditTm = CreateInput(fmt(calculator.defaults.Tm).c_str(), IDC_EDIT_TM, y);

    y += step;
    hLabelDur = CreateLabel(L"Number of Waves (Nz):", 10, y);
    hEditDur = CreateInput(fmt(calculator.defaults.Number_of_Waves).c_str(), IDC_EDIT_DUR, y);

    y += step;
    hLabelNod = CreateLabel(L"Nod (Damage):", 10, y);
    hEditNod = CreateInput(fmt(calculator.defaults.Nod).c_str(), IDC_EDIT_NOD, y);

    y += step;
    hLabelWc = CreateLabel(L"Wc Concrete (kN/m3):", 10, y);
    hEditWc = CreateInput(fmt(calculator.defaults.Wc).c_str(), IDC_EDIT_WC, y);

    y += step;
    hLabelFormula = CreateLabel(L"Formula / Slope:", 10, y);

    y += 30;
    hComboFormula = CreateWindowW(L"COMBOBOX", NULL,
        WS_CHILD | WS_VISIBLE | CBS_DROPDOWNLIST | WS_VSCROLL,
        10, y, 350, 200, hwnd, (HMENU)(INT_PTR)IDC_COMBO_FORMULA, NULL, NULL);
    SendMessageW(hComboFormula, WM_SETFONT, (WPARAM)hUIFont, TRUE);

    y += 45;
    hButtonCompute = CreateWindowW(L"BUTTON", L"Calculate Armor",
        WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON,
        10, y, 260, 40, hwnd, (HMENU)(INT_PTR)IDC_BUTTON_COMPUTE, NULL, NULL);
    SendMessageW(hButtonCompute, WM_SETFONT, (WPARAM)hUIFont, TRUE);

    hOutput = CreateWindowW(L"EDIT", L"", WS_CHILD | WS_VISIBLE | WS_BORDER |
                           ES_MULTILINE | ES_AUTOVSCROLL | WS_VSCROLL | ES_READONLY | WS_HSCROLL | ES_AUTOHSCROLL,
                           400, 10, 790, 760, hwnd, (HMENU)IDC_OUTPUT, NULL, NULL);

    hMonoFont = CreateFontW(21, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE,
                           DEFAULT_CHARSET, OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS,
                           DEFAULT_QUALITY, FIXED_PITCH | FF_DONTCARE, L"Courier New");
    SendMessageW(hOutput, WM_SETFONT, (WPARAM)hMonoFont, TRUE);

    ApplyLanguage(hwnd);
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_CREATE:
        CreateControls(hwnd);
        break;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDC_COMBO_LANGUAGE && HIWORD(wParam) == CBN_SELCHANGE) {
            int sel = (int)SendMessageW(hComboLanguage, CB_GETCURSEL, 0, 0);
            if (sel < 0) sel = 0;
            g_currentLanguage = static_cast<Language>(sel);
            ApplyLanguage(hwnd);
            break;
        }

        if (LOWORD(wParam) == IDC_BUTTON_COMPUTE) {
            const LangPack& t = GetLangPack(g_currentLanguage);

            wchar_t buffer[64];
            Inputs inputs = calculator.defaults;

            GetWindowTextW(hEditHs, buffer, 63);  inputs.Hs = _wtof(buffer);
            GetWindowTextW(hEditTm, buffer, 63);  inputs.Tm = _wtof(buffer);
            GetWindowTextW(hEditDur, buffer, 63); inputs.Number_of_Waves = _wtof(buffer);
            GetWindowTextW(hEditNod, buffer, 63); inputs.Nod = _wtof(buffer);
            GetWindowTextW(hEditWc, buffer, 63);  inputs.Wc = _wtof(buffer);

            int selIndex = (int)SendMessageW(hComboFormula, CB_GETCURSEL, 0, 0);
            if (selIndex < 0 || selIndex > 3) {
                selIndex = calculator.defaults.Formula_ID - 1;
            }
            inputs.Formula_ID = selIndex + 1;

            if (inputs.Hs <= 0 || inputs.Tm <= 0 || inputs.Number_of_Waves <= 0 || inputs.Wc <= 0) {
                MessageBoxW(hwnd, t.input_error_positive.c_str(), t.input_error_title.c_str(), MB_ICONERROR | MB_OK);
                break;
            }

            try {
                FullResults res = calculator.solve(inputs.Formula_ID, inputs);
                g_lastResults = res;
                g_hasLastResults = true;

                std::wstring report = format_report(res, g_currentLanguage);
                SaveReportToFile(report);

                std::wstring guiText = fix_newlines_for_edit_control(report);
                SetWindowTextW(hOutput, guiText.c_str());
            } catch (const std::exception& e) {
                const char* what = e.what();
                std::wstring errMsg(what, what + strlen(what));
                MessageBoxW(hwnd, errMsg.c_str(), t.calc_error_title.c_str(), MB_ICONERROR | MB_OK);
            } catch (...) {
                MessageBoxW(hwnd, L"Unexpected internal error.", t.calc_error_title.c_str(), MB_ICONERROR | MB_OK);
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
        CW_USEDEFAULT, CW_USEDEFAULT, 1240, 840,
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