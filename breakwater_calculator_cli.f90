! ======================================================================================
! PROGRAM DESCRIPTION & METHODOLOGY
! ======================================================================================
!
! 1. PURPOSE:
!    This software performs preliminary hydraulic sizing of rubble-mound
!    breakwater armor made with artificial concrete units and the associated
!    rock underlayers.
!
!    The armor-unit families currently implemented in the code are:
!    - Simple Cubes
!    - Antifer Blocks
!
!    The program dimensions two structural zones:
!    - The Trunk (main longitudinal section)
!    - The Head (the exposed roundhead / end section)
!
!    The implemented design philosophy is iso-geometric, but not iso-weight:
!    - the nominal diameter is kept the same at Trunk and Head;
!    - the armor-unit geometry is kept the same at Trunk and Head;
!    - the slope is kept the same at Trunk and Head;
!    - the Head unit weight increases through higher required concrete specific
!      weight rather than through larger unit size.
!
! 2. METHODOLOGY:
!    The calculator uses empirical hydraulic-stability formulas derived from
!    physical model testing together with an internal EN 13383 rock-grading
!    database.
!
!    Implemented armor formulas:
!    - Simple Cubes: Van der Meer (1988)
!    - Antifer Blocks: Chegini & Aghtouman (2006)
!
!    Underlayers:
!    - automatic grading selection from an embedded EN 13383 mass-based
!      grading table.
!
! 3. DESIGN LOGIC & COMPUTATIONAL STRATEGY:
!
!    a. Exposed storm input:
!       The storm variable entered by the user is the number of waves, Nz.
!       Storm duration is not entered directly; it is computed internally as:
!
!           Storm_Duration_hr = Nz * Tm / 3600
!
!    b. Trunk stability:
!       For the selected formula, the program evaluates the Stability Number:
!
!           Ns = (k1 * (Nod^k2 / Nz^k3) + k4) * s0m^(-k5)
!
!       where:
!
!           L0  = g * Tm^2 / (2 * pi)
!           s0m = Hs / L0
!
!       The trunk buoyant relative density is:
!
!           Delta_trunk = (Wc_trunk / Ww) - 1
!
!       The nominal diameter is then:
!
!           Dn = Hs / (Delta_trunk * Ns)
!
!       and the trunk armor-unit weight is:
!
!           W_trunk = Wc_trunk * Dn^3
!
!       The program also derives an equivalent Hudson stability coefficient
!       for the trunk.
!
!    c. Head design - iso-geometric transfer:
!       The breakwater head is subjected to stronger three-dimensional flow
!       effects than the trunk. Instead of increasing block size or flattening
!       the slope, the program preserves armor geometry and transfers the trunk
!       design to the head through a fixed Hudson-coefficient ratio:
!
!           Kd_trunk / Kd_head = 1.5
!
!       With geometry kept fixed, the required head relative density becomes:
!
!           Delta_head = Delta_trunk * (1.5)^(1/3)
!
!       and the required concrete specific weight at the head is:
!
!           Wc_head = Ww * (Delta_head + 1)
!
!       Since block volume is kept constant, the head block weight scales only
!       with concrete specific weight:
!
!           W_head = W_trunk * (Wc_head / Wc_trunk)
!
!       Practical consequence:
!       - the same molds can be used at Trunk and Head;
!       - the same nominal diameter is maintained;
!       - Head units are heavier because concrete specific weight increases.
!
!    d. Underlayer selection:
!       The theoretical underlayer target is one tenth of the armor-unit weight:
!
!           W_target = W_armor / 10
!
!       The target mass is compared against the embedded EN 13383 grading
!       limits. A grading is selected only if the target mass is strictly
!       contained within the nominal lower and upper limits.
!
!       If more than one grading contains the target mass, the program selects
!       the tightest range, i.e. the grading with the smallest:
!
!           NUL - NLL
!
!       If no grading strictly contains the target mass, the program falls back
!       to the first grading in the internal grading database.
!
! 4. DEFAULT INPUTS, CONSTANTS, AND CURRENT FORMULA SETS USED IN THE CODE:
!
!    a. Default user inputs currently used by the Fortran CLI:
!       - Hs              = 11.0 m
!       - Tm              = 11.9 s
!       - Number_of_Waves = 3000
!       - Nod             = 0.5
!       - Wc              = 27.48 kN/m3
!       - Ww              = 10.05 kN/m3
!       - Formula_ID      = 1
!
!       Default Formula_ID = 1 corresponds to:
!       Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)
!
!    b. Internal constants currently used by the code:
!       - g              = 9.80665 m/s2
!       - W_rock_spec    = 26.5 kN/m3
!       - P_cubes        = 0.40
!       - P_rock         = 0.25
!       - KD_RATIO_FIXED = 1.5
!
!    c. Formula sets currently embedded in the code:
!       1. Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)
!       2. Van Der Meer (1988a) - Simple Cubes (Slope 1.5:1)
!       3. Chegini-Aghtouman (2006) - Antifer (Slope 2:1)
!       4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)
!
! 5. OUTPUTS PRODUCED:
!    The generated report contains, for both Trunk and Head:
!    - stability number;
!    - nominal diameter;
!    - unit weight and unit mass;
!    - equivalent Kd value;
!    - armor layer thickness;
!    - packing density;
!    - geometric dimensions reported by the code;
!    - adopted underlayer grading and derived rock parameters.
!
!    Volume reporting follows the current Fortran code logic:
!
!    - for Simple Cubes:
!
!          V = Dn^3
!
!    - for Antifer blocks:
!
!          V = 1.0247 * H^3
!
!    The Fortran CLI appends the technical report to output.txt at each
!    execution.
!
! 6. TECHNICAL BIBLIOGRAPHY & REFERENCES:
!
!    1. Van der Meer, J.W. (1988). "Rock Slopes and Gravel Beaches Under Wave
!       Attack." Doctoral Thesis, Delft University of Technology.
!
!    2. Van der Meer, J.W. (1988). "Stability of Cubes, Tetrapods and
!       Accropode." Proceedings of the Conference Breakwaters '88, Eastbourne,
!       Thomas Telford.
!
!    3. Chegini, V., & Aghtouman, P. (2006). "An Investigation on the
!       Stability of Rubble Mound Breakwaters with Armour Layers of Antifer
!       Cubes."
!
!    4. USACE (2006). "Coastal Engineering Manual (CEM)", Chapter VI-5.
!
!    5. CEN (2002). "EN 13383-1: Armourstone - Part 1: Specification."
!
! 7. COMPILATION INSTRUCTIONS:
!
!    Recommended GNU Fortran build:
!
!    gfortran -O3 -march=native -std=f2008 -Wall -Wextra -pedantic
!    -Wconversion -static -static-libgfortran -static-libgcc
!    -o breakwater_calculator_cli breakwater_calculator_cli.f90
!
!    Meaning of the main compilation flags:
!    - -O3                  -> high optimization level;
!    - -march=native        -> optimize for the local CPU;
!    - -std=f2008           -> compile as Fortran 2008;
!    - -Wall -Wextra        -> enable common warnings;
!    - -pedantic            -> request stricter standards checking;
!    - -Wconversion         -> warn about implicit numeric conversions;
!    - -static              -> prefer static linking;
!    - -static-libgfortran  -> statically link Fortran runtime;
!    - -static-libgcc       -> statically link GCC runtime.
!
!    If full static linking fails on the local toolchain, a fallback build is:
!
!    gfortran -O3 -march=native -std=f2008 -Wall -Wextra -pedantic
!    -Wconversion -o breakwater_calculator_cli breakwater_calculator_cli.f90
!
! 8. EXECUTION MODES AND EXAMPLES:
!
!    a. Interactive mode:
!       Run the executable with no command-line arguments:
!
!           ./breakwater_calculator_cli
!
!       The program will:
!       - display the 4 available formula options;
!       - ask for the formula selection;
!       - ask for Hs, Tm, Nz, Nod, and Wc;
!       - accept ENTER to use each default value;
!       - calculate results and append output.txt.
!
!    b. Command-line argument mode:
!       Format:
!
!           ./breakwater_calculator_cli [Hs] [Tm] [NumberOfWaves] [Nod] [Wc] [FormulaID]
!
!       Example using the current code defaults:
!
!           ./breakwater_calculator_cli 11.0 11.9 3000 0.5 27.48 1
!
!       FormulaID must be one of:
!       - 1 -> Simple Cubes, slope 2.0:1
!       - 2 -> Simple Cubes, slope 1.5:1
!       - 3 -> Antifer, slope 2:1
!       - 4 -> Antifer, slope 1.5:1
!
! ======================================================================================

PROGRAM BreakwaterCalculator
    IMPLICIT NONE

    ! ----------------------------------------------------------------------
    ! DATA STRUCTURES
    ! ----------------------------------------------------------------------

    ! Define M_PI if not already defined
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307)
    REAL(dp), PARAMETER :: PI = 3.14159265358979323846_dp

    ! ----------------------------------------------------------------------
    ! PHYSICAL CONSTANTS
    ! ----------------------------------------------------------------------
    REAL(dp), PARAMETER :: g = 9.80665_dp  ! Acceleration due to gravity (m/s^2)
    REAL(dp), PARAMETER :: W_rock_spec_val = 26.5_dp ! Standard Specific weight for underlayer rock (kN/m3)

    ! ----------------------------------------------------------------------
    ! LAYER CHARACTERISTICS
    ! ----------------------------------------------------------------------
    REAL(dp), PARAMETER :: P_cubes = 0.40_dp ! Porosity of armor layers (Standard for double-layer Cubes)
    REAL(dp), PARAMETER :: P_rock = 0.25_dp  ! Porosity of rock underlayers (Standard approximation)

    ! ----------------------------------------------------------------------
    ! HEAD TO TRUNK TRANSFER RATIO
    ! ----------------------------------------------------------------------
    REAL(dp), PARAMETER :: KD_RATIO_FIXED = 1.5_dp

    ! Struct definitions mirroring C++
    
    TYPE :: FormulaParams
        CHARACTER(LEN=100) :: name
        CHARACTER(LEN=64) :: type_name
        REAL(dp) :: slope_ratio
        REAL(dp) :: k1, k2, k3, k4, k5
    END TYPE FormulaParams

    TYPE :: GradingDef
        CHARACTER(LEN=64) :: name
        REAL(dp) :: NLL_kg
        REAL(dp) :: NUL_kg
    END TYPE GradingDef

    TYPE :: Dimensions
        REAL(dp) :: H
        REAL(dp) :: A
        REAL(dp) :: B
    END TYPE Dimensions

    TYPE :: UnderlayerResult
        CHARACTER(LEN=64) :: grading_name
        REAL(dp) :: target_W
        REAL(dp) :: target_M50_kg
        REAL(dp) :: M50_kg
        REAL(dp) :: ELL_kg
        REAL(dp) :: EUL_kg
        REAL(dp) :: NLL_kg
        REAL(dp) :: NUL_kg
        REAL(dp) :: W_mean_kn
        REAL(dp) :: Dn_rock
        REAL(dp) :: r2
        REAL(dp) :: f2
        REAL(dp) :: W_rock_spec
    END TYPE UnderlayerResult

    TYPE :: ArmorResult
        REAL(dp) :: Ns
        REAL(dp) :: Dn
        REAL(dp) :: W
        REAL(dp) :: Mass_tonnes
        REAL(dp) :: Kd
        REAL(dp) :: r1
        REAL(dp) :: packing_density
        TYPE(Dimensions) :: dims
        
        ! Specific to Head
        REAL(dp) :: Kd_Ratio
        REAL(dp) :: Delta_Required
        REAL(dp) :: Wc_Required
    END TYPE ArmorResult

    TYPE :: Intermediates
        REAL(dp) :: L0
        REAL(dp) :: k0
        REAL(dp) :: s0m
        REAL(dp) :: sm
        REAL(dp) :: Nz
        REAL(dp) :: Storm_Duration_hr
        REAL(dp) :: delta
        REAL(dp) :: Ns_trunk
    END TYPE Intermediates

    TYPE :: Inputs
        REAL(dp) :: Hs
        REAL(dp) :: Tm
        REAL(dp) :: Number_Of_Waves
        REAL(dp) :: Nod
        REAL(dp) :: Wc
        REAL(dp) :: Ww
        INTEGER :: Formula_ID
    END TYPE Inputs

    TYPE :: FullResults
        TYPE(Inputs) :: inputs
        TYPE(FormulaParams) :: coefficients
        TYPE(Intermediates) :: intermediate
        TYPE(ArmorResult) :: final_trunk
        TYPE(UnderlayerResult) :: underlayer_trunk
        TYPE(ArmorResult) :: final_head
        TYPE(UnderlayerResult) :: underlayer_head
        REAL(dp) :: P_rock_val
        REAL(dp) :: P_cubes_val
    END TYPE FullResults

    ! Data containers
    TYPE(FormulaParams) :: formulas(4)
    TYPE(GradingDef) :: standard_gradings(17)
    
    ! Local Variables for Main
    TYPE(Inputs) :: user_inputs
    TYPE(FullResults) :: results
    INTEGER :: formula_id_arg
    INTEGER :: num_args
    CHARACTER(LEN=100) :: arg_val
    CHARACTER(LEN=100) :: buffer
    INTEGER :: iostat

    ! ----------------------------------------------------------------------
    ! DATA INITIALIZATION (Runtime Assignment)
    ! ----------------------------------------------------------------------
    ! FORMULA DATABASE
    
    ! ID 1: Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)
    formulas(1)%name = "Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)"
    formulas(1)%type_name = "Cubes"
    formulas(1)%slope_ratio = 2.0_dp
    formulas(1)%k1 = 7.374304189198_dp
    formulas(1)%k2 = 0.4_dp
    formulas(1)%k3 = 0.3_dp
    formulas(1)%k4 = 1.100642416298_dp
    formulas(1)%k5 = 0.1_dp
    
    ! ID 2: Van Der Meer (1988a) - Simple Cubes (Slope 1.5:1)
    formulas(2)%name = "Van Der Meer (1988a) - Simple Cubes (Slope 1.5:1)"
    formulas(2)%type_name = "Cubes"
    formulas(2)%slope_ratio = 1.5_dp
    formulas(2)%k1 = 6.7_dp
    formulas(2)%k2 = 0.4_dp
    formulas(2)%k3 = 0.3_dp
    formulas(2)%k4 = 1.0_dp
    formulas(2)%k5 = 0.1_dp

    ! ID 3: Chegini-Aghtouman (2006) - Antifer (Slope 2:1)
    formulas(3)%name = "Chegini-Aghtouman (2006) - Antifer (Slope 2:1)"
    formulas(3)%type_name = "Antifer"
    formulas(3)%slope_ratio = 2.0_dp
    formulas(3)%k1 = 6.138_dp
    formulas(3)%k2 = 0.443_dp
    formulas(3)%k3 = 0.276_dp
    formulas(3)%k4 = 1.164_dp
    formulas(3)%k5 = 0.07_dp
    
    ! ID 4: Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)
    formulas(4)%name = "Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)"
    formulas(4)%type_name = "Antifer"
    formulas(4)%slope_ratio = 1.5_dp
    formulas(4)%k1 = 6.951_dp
    formulas(4)%k2 = 0.443_dp
    formulas(4)%k3 = 0.291_dp
    formulas(4)%k4 = 1.082_dp
    formulas(4)%k5 = 0.082_dp

    ! ROCK GRADING DATABASE (EN 13383)
    ! Coarse / Light Gradings
    standard_gradings(1)%name="CP32/90";       standard_gradings(1)%NLL_kg=0.868_dp;    standard_gradings(1)%NUL_kg=19.319_dp
    standard_gradings(2)%name="CP45/125";      standard_gradings(2)%NLL_kg=2.415_dp;    standard_gradings(2)%NUL_kg=51.758_dp
    standard_gradings(3)%name="CP63/180";      standard_gradings(3)%NLL_kg=6.626_dp;    standard_gradings(3)%NUL_kg=154.548_dp
    standard_gradings(4)%name="CP90/250";      standard_gradings(4)%NLL_kg=19.319_dp;   standard_gradings(4)%NUL_kg=414.063_dp
    standard_gradings(5)%name="CP45/180";      standard_gradings(5)%NLL_kg=2.415_dp;    standard_gradings(5)%NUL_kg=154.548_dp
    standard_gradings(6)%name="CP90/180";      standard_gradings(6)%NLL_kg=19.319_dp;   standard_gradings(6)%NUL_kg=154.548_dp
    ! Light Mass Armourstone (LMA)
    standard_gradings(7)%name="LMA5-40";       standard_gradings(7)%NLL_kg=5.0_dp;      standard_gradings(7)%NUL_kg=40.0_dp
    standard_gradings(8)%name="LMA10-60";      standard_gradings(8)%NLL_kg=10.0_dp;     standard_gradings(8)%NUL_kg=60.0_dp
    standard_gradings(9)%name="LMA15-120";     standard_gradings(9)%NLL_kg=15.0_dp;     standard_gradings(9)%NUL_kg=120.0_dp
    standard_gradings(10)%name="LMA40-200";    standard_gradings(10)%NLL_kg=40.0_dp;    standard_gradings(10)%NUL_kg=200.0_dp
    standard_gradings(11)%name="LMA60-300";    standard_gradings(11)%NLL_kg=60.0_dp;    standard_gradings(11)%NUL_kg=300.0_dp
    standard_gradings(12)%name="LMA15-300";    standard_gradings(12)%NLL_kg=15.0_dp;    standard_gradings(12)%NUL_kg=300.0_dp
    ! Heavy Mass Armourstone (HMA)
    standard_gradings(13)%name="HMA300-1000";  standard_gradings(13)%NLL_kg=300.0_dp;   standard_gradings(13)%NUL_kg=1000.0_dp
    standard_gradings(14)%name="HMA1000-3000"; standard_gradings(14)%NLL_kg=1000.0_dp;  standard_gradings(14)%NUL_kg=3000.0_dp
    standard_gradings(15)%name="HMA3000-6000"; standard_gradings(15)%NLL_kg=3000.0_dp;  standard_gradings(15)%NUL_kg=6000.0_dp
    standard_gradings(16)%name="HMA6000-10000";standard_gradings(16)%NLL_kg=6000.0_dp;  standard_gradings(16)%NUL_kg=10000.0_dp
    standard_gradings(17)%name="HMA10000-15000";standard_gradings(17)%NLL_kg=10000.0_dp;standard_gradings(17)%NUL_kg=15000.0_dp

    ! ----------------------------------------------------------------------
    ! DEFAULT INPUTS
    ! ----------------------------------------------------------------------
    user_inputs%Hs = 11.0_dp
    user_inputs%Tm = 11.9_dp
    user_inputs%Number_Of_Waves = 3000.0_dp
    user_inputs%Nod = 0.5_dp
    user_inputs%Wc = 27.48_dp
    user_inputs%Ww = 10.05_dp
    user_inputs%Formula_ID = 1

    ! ==============================================================================
    ! MAIN EXECUTION BLOCK
    ! ==============================================================================
    
    num_args = COMMAND_ARGUMENT_COUNT()

    ! Check if command line arguments are provided
    IF (num_args >= 6) THEN
        ! Parse CLI arguments
        CALL GET_COMMAND_ARGUMENT(1, arg_val); READ(arg_val, *) user_inputs%Hs
        CALL GET_COMMAND_ARGUMENT(2, arg_val); READ(arg_val, *) user_inputs%Tm
        CALL GET_COMMAND_ARGUMENT(3, arg_val); READ(arg_val, *) user_inputs%Number_Of_Waves
        CALL GET_COMMAND_ARGUMENT(4, arg_val); READ(arg_val, *) user_inputs%Nod
        CALL GET_COMMAND_ARGUMENT(5, arg_val); READ(arg_val, *) user_inputs%Wc
        CALL GET_COMMAND_ARGUMENT(6, arg_val); READ(arg_val, *) formula_id_arg
        
        IF (formula_id_arg < 1 .OR. formula_id_arg > 4) THEN
             PRINT *, "Error: Formula ID must be 1-4. Using default (1)."
             formula_id_arg = 1
        END IF
        user_inputs%Formula_ID = formula_id_arg
    ELSE
        ! Interactive Mode
        PRINT *, ""
        PRINT *, "--- COASTAL PROTECTION BLOCK CALCULATOR (TRUNK & HEAD) ---"
        PRINT *, "1. Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)"
        PRINT *, "2. Van Der Meer (1988a) - Cubes (Slope 1.5:1)"
        PRINT *, "3. Chegini-Aghtouman (2006) - Antifer (Slope 2:1)"
        PRINT *, "4. Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)"
        
        PRINT *, ""
        WRITE(*, '(A)', ADVANCE='NO') "Option [1-4]: "
        READ(*, '(A)') buffer
        
        IF (LEN_TRIM(buffer) > 0) THEN
            READ(buffer, *, IOSTAT=iostat) formula_id_arg
            IF (iostat /= 0 .OR. formula_id_arg < 1 .OR. formula_id_arg > 4) formula_id_arg = 1
        ELSE
            formula_id_arg = 1
        END IF
        user_inputs%Formula_ID = formula_id_arg

        PRINT *, ""
        PRINT *, "--- Enter Parameters (Press ENTER for Default) ---"
        
        CALL get_param("Hs (m)", user_inputs%Hs)
        CALL get_param("Tm (s)", user_inputs%Tm)
        CALL get_param("Number of waves (Nz)", user_inputs%Number_Of_Waves)
        CALL get_param("Nod (Damage)", user_inputs%Nod)
        CALL get_param("Concrete Weight Trunk (kN/m3)", user_inputs%Wc)
    END IF

    CALL solve(user_inputs%Formula_ID, user_inputs, results)
    CALL generate_report_file(results, "output.txt")

CONTAINS

    FUNCTION calculate_L0(Tm) RESULT(L0)
        REAL(dp), INTENT(IN) :: Tm
        REAL(dp) :: L0
        L0 = (g * Tm**2) / (2.0_dp * PI)
    END FUNCTION calculate_L0

    FUNCTION calculate_underlayer_params(W_armor) RESULT(res)
        REAL(dp), INTENT(IN) :: W_armor
        TYPE(UnderlayerResult) :: res
        
        REAL(dp) :: target_weight
        REAL(dp) :: target_mass_kg
        REAL(dp) :: min_range_width, current_range
        REAL(dp) :: final_M50, final_NLL, final_NUL
        REAL(dp) :: ELL, EUL, W_mean_kn
        INTEGER :: i
        CHARACTER(LEN=64) :: best_name
        LOGICAL :: found

        target_weight = W_armor / 10.0_dp
        target_mass_kg = (target_weight * 1000.0_dp) / g
        
        min_range_width = HUGE(1.0_dp)
        found = .FALSE.
        
        final_M50 = 0.0_dp
        final_NLL = 0.0_dp
        final_NUL = 0.0_dp

        ! --- CONTAINMENT & TIGHTEST RANGE LOGIC ---
        DO i = 1, SIZE(standard_gradings)
            ! Check containment: Target must be strictly inside nominal limits
            IF (target_mass_kg > standard_gradings(i)%NLL_kg .AND. target_mass_kg < standard_gradings(i)%NUL_kg) THEN
                current_range = standard_gradings(i)%NUL_kg - standard_gradings(i)%NLL_kg
                
                ! Update if this is the first match OR if this range is tighter (smaller)
                IF (current_range < min_range_width) THEN
                    min_range_width = current_range
                    best_name = standard_gradings(i)%name
                    final_NLL = standard_gradings(i)%NLL_kg
                    final_NUL = standard_gradings(i)%NUL_kg
                    final_M50 = 0.5_dp * (final_NLL + final_NUL)
                    found = .TRUE.
                END IF
            END IF
        END DO
        
        ! Fallback if no grading strictly contains the target mass
        IF (.NOT. found .AND. SIZE(standard_gradings) > 0) THEN
            best_name = standard_gradings(1)%name
            final_NLL = standard_gradings(1)%NLL_kg
            final_NUL = standard_gradings(1)%NUL_kg
            final_M50 = 0.5_dp * (final_NLL + final_NUL)
        END IF

        ! --- EN 13383 LIMIT CALCULATIONS ---
        ELL = 0.7_dp * final_NLL
        EUL = 1.5_dp * final_NUL
        W_mean_kn = (final_M50 * g) / 1000.0_dp
        
        res%grading_name = best_name
        res%target_W = target_weight
        res%target_M50_kg = target_mass_kg
        res%M50_kg = final_M50
        res%NLL_kg = final_NLL
        res%NUL_kg = final_NUL
        res%ELL_kg = ELL
        res%EUL_kg = EUL
        res%W_mean_kn = W_mean_kn
        res%Dn_rock = (W_mean_kn / W_rock_spec_val)**(1.0_dp/3.0_dp)
        res%r2 = 2.0_dp * res%Dn_rock
        res%f2 = 100.0_dp * 2.0_dp * 1.0_dp * (1.0_dp - P_rock) / (res%Dn_rock**2)
        res%W_rock_spec = W_rock_spec_val
    END FUNCTION calculate_underlayer_params

    SUBROUTINE solve(formula_id, params, results_out)
        INTEGER, INTENT(IN) :: formula_id
        TYPE(Inputs), INTENT(IN) :: params
        TYPE(FullResults), INTENT(OUT) :: results_out
        
        TYPE(FormulaParams) :: coeffs
        
        REAL(dp) :: L0, k0, s0m, sm, Nz, storm_duration_hr, delta_trunk
        REAL(dp) :: term_damage, term_waves, damage_wave_ratio, scaled_term, inv_f
        REAL(dp) :: steepness_factor, Ns_trunk
        REAL(dp) :: Dn, W_trunk, packing_density_trunk, kd_trunk_equiv
        REAL(dp) :: h_trunk, a_trunk, b_trunk
        
        ! Head vars
        REAL(dp) :: kd_ratio, kd_head_derived, delta_head, Wc_head, W_head
        REAL(dp) :: Ns_head, packing_density_head, h_head, a_head, b_head
        
        ! 2. Load Coefficients
        IF (formula_id < 1 .OR. formula_id > 4) THEN
            PRINT *, "Invalid Formula ID"
            STOP 1
        END IF
        
        ! Copy to local struct for easy access
        coeffs = formulas(formula_id)
        
        ! 3. Preliminary Hydraulic Calculations
        L0 = calculate_L0(params%Tm)
        k0 = (2.0_dp * PI) / L0
        s0m = params%Hs / L0
        sm = s0m
        
        Nz = params%Number_Of_Waves
        storm_duration_hr = (Nz * params%Tm) / 3600.0_dp
        delta_trunk = (params%Wc / params%Ww) - 1.0_dp
        
        ! 4. Algorithmic Core (Chegini-Aghtouman / Van der Meer) - TRUNK
        term_damage = params%Nod**coeffs%k2
        term_waves = Nz**coeffs%k3
        damage_wave_ratio = term_damage / term_waves
        
        scaled_term = coeffs%k1 * damage_wave_ratio
        inv_f = scaled_term + coeffs%k4
        
        steepness_factor = s0m**(-coeffs%k5)
        
        Ns_trunk = inv_f * steepness_factor
        
        ! 5. Block Sizing (Armor) - TRUNK
        Dn = params%Hs / (delta_trunk * Ns_trunk)
        W_trunk = params%Wc * (Dn**3)
        packing_density_trunk = 100.0_dp * 2.0_dp * 1.1_dp * (1.0_dp - P_cubes) / (Dn**2)
        
        ! 6. Hudson Comparative Calculation
        kd_trunk_equiv = (params%Wc * (params%Hs**3)) / (W_trunk * (delta_trunk**3) * coeffs%slope_ratio)
        
        ! 7. UNDERLAYER - TRUNK
        results_out%underlayer_trunk = calculate_underlayer_params(W_trunk)
        
        ! 8. HEAD CALCULATION (FIXED RATIO 1.5)
        kd_ratio = KD_RATIO_FIXED
        kd_head_derived = kd_trunk_equiv / kd_ratio
        delta_head = delta_trunk * (kd_ratio**(1.0_dp/3.0_dp))
        Wc_head = params%Ww * (delta_head + 1.0_dp)
        W_head = W_trunk * (Wc_head / params%Wc)
        
        ! Calculate Stability Number for Head
        Ns_head = params%Hs / (delta_head * Dn)
        packing_density_head = 100.0_dp * 2.0_dp * 1.1_dp * (1.0_dp - P_cubes) / (Dn**2)
        
        ! 9. UNDERLAYER - HEAD
        results_out%underlayer_head = calculate_underlayer_params(W_head)
        
        ! 10. Armor Layer Details (Common)
        results_out%final_trunk%r1 = 2.0_dp * 1.1_dp * Dn
        
        ! Dimensions
        h_trunk = ((W_trunk / params%Wc) / 1.0247_dp)**(1.0_dp/3.0_dp)
        a_trunk = 1.086_dp * h_trunk
        b_trunk = 1.005_dp * h_trunk
        
        h_head = ((W_head / Wc_head) / 1.0247_dp)**(1.0_dp/3.0_dp)
        a_head = 1.086_dp * h_head
        b_head = 1.005_dp * h_head
        
        ! 11. Compile Results
        results_out%inputs = params
        results_out%coefficients = coeffs
        
        results_out%intermediate%L0 = L0
        results_out%intermediate%k0 = k0
        results_out%intermediate%s0m = s0m
        results_out%intermediate%sm = sm
        results_out%intermediate%Nz = Nz
        results_out%intermediate%Storm_Duration_hr = storm_duration_hr
        results_out%intermediate%delta = delta_trunk
        results_out%intermediate%Ns_trunk = Ns_trunk
        
        results_out%final_trunk%Ns = Ns_trunk
        results_out%final_trunk%Dn = Dn
        results_out%final_trunk%W = W_trunk
        results_out%final_trunk%Mass_tonnes = W_trunk / g
        results_out%final_trunk%Kd = kd_trunk_equiv
        results_out%final_trunk%packing_density = packing_density_trunk
        results_out%final_trunk%dims%H = h_trunk
        results_out%final_trunk%dims%A = a_trunk
        results_out%final_trunk%dims%B = b_trunk
        
        results_out%final_head%Ns = Ns_head
        results_out%final_head%Dn = Dn
        results_out%final_head%Kd = kd_head_derived
        results_out%final_head%Kd_Ratio = kd_ratio
        results_out%final_head%Delta_Required = delta_head
        results_out%final_head%Wc_Required = Wc_head
        results_out%final_head%W = W_head
        results_out%final_head%Mass_tonnes = W_head / g
        results_out%final_head%packing_density = packing_density_head
        results_out%final_head%dims%H = h_head
        results_out%final_head%dims%A = a_head
        results_out%final_head%dims%B = b_head
        
        results_out%P_rock_val = P_rock
        results_out%P_cubes_val = P_cubes
        
    END SUBROUTINE solve

    SUBROUTINE generate_report_file(results, filepath)
        TYPE(FullResults), INTENT(IN) :: results
        CHARACTER(LEN=*), INTENT(IN) :: filepath
        
        INTEGER :: u, iostatus_val
        CHARACTER(LEN=200) :: line
        CHARACTER(LEN=80) :: separator
        REAL(dp) :: armor_vol_trunk, armor_vol_head
        
        separator = "--------------------------------------------------------------------------------"

        IF (TRIM(results%coefficients%type_name) == "Antifer") THEN
            armor_vol_trunk = 1.0247_dp * (results%final_trunk%dims%H**3)
            armor_vol_head = 1.0247_dp * (results%final_head%dims%H**3)
        ELSE
            armor_vol_trunk = results%final_trunk%Dn**3
            armor_vol_head = results%final_head%Dn**3
        END IF
        
        OPEN(NEWUNIT=u, FILE=filepath, STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND', IOSTAT=iostatus_val)
        IF (iostatus_val /= 0) THEN
            PRINT *, "Error saving file."
            RETURN
        END IF
        
        WRITE(u, '(A)') "================================================================================"
        WRITE(u, '(A)') "    TECHNICAL REPORT: BREAKWATER ARMOR & UNDERLAYER DESIGN                        "
        WRITE(u, '(A)') "================================================================================"
        WRITE(u, '(A)') "Methodology: " // TRIM(results%coefficients%name)
        WRITE(u, '(A)') TRIM(separator)
        
        WRITE(u, '(A)') "1. INPUT PARAMETERS"
        WRITE(u, '(A, F0.2, A)') "   Hs (Sigificant Wave Height)         : ", results%inputs%Hs, " m"
        WRITE(u, '(A, F0.2, A)') "   Tm (Mean Wave Period)               : ", results%inputs%Tm, " s"
        WRITE(u, '(A, I0)')      "   Number of waves (Nz)                : ", NINT(results%inputs%Number_Of_Waves)
        WRITE(u, '(A, F0.2)')    "   Nod (Damage)                        : ", results%inputs%Nod
        WRITE(u, '(A, F0.2, A)') "   Wc Trunk (Concrete Spec. Weight)    : ", results%inputs%Wc, " kN/m3"
        WRITE(u, '(A, F0.2, A)') "   Ww (Water Specific Weight)          : ", results%inputs%Ww, " kN/m3"
        WRITE(u, '(A, F0.4)')    "   Relative Density D=(Wc/Ww)-1        : ", results%intermediate%delta
        WRITE(u, '(A, F0.1, A)') "   Structure Slope (TRUNK & HEAD)      : ", results%coefficients%slope_ratio, ":1"
        WRITE(u, '(A, I0, A)')   "   Porosity (Cubes)                    : ", INT(results%P_cubes_val * 100.0_dp), "%"
        WRITE(u, '(A, I0, A)')   "   Porosity (Rock Layer)               : ", INT(results%P_rock_val * 100.0_dp), "%"
        WRITE(u, '(A)') TRIM(separator)
        
        WRITE(u, '(A)') "2. INTERMEDIATE PARAMETERS"
        WRITE(u, '(A, F0.2, A)') "   Wave Length (L0)                    : ", results%intermediate%L0, " m"
        WRITE(u, '(A, F0.4)')    "   wave number (k0 = 2*pi/L0)          : ", results%intermediate%k0
        WRITE(u, '(A, F0.4)')    "   wave steepness (s0m = Hs/L0)        : ", results%intermediate%s0m
        WRITE(u, '(A, F0.3, A)') "   Storm Duration (h)                  : ", results%intermediate%Storm_Duration_hr, " h"
        WRITE(u, '(A, F0.4)')    "   Stability Number TRUNK (Ns)         : ", results%intermediate%Ns_trunk
        WRITE(u, '(A, F0.4)')    "   Stability Number HEAD (Ns)          : ", results%final_head%Ns
        WRITE(u, '(A)') TRIM(separator)
        
        WRITE(u, '(A)') "3. ARMOR LAYER RESULTS - TRUNK"
        WRITE(u, '(A, F0.2, A)') "   BLOCK WEIGHT (W)                    : ", results%final_trunk%W, " kN"
        WRITE(u, '(A, F0.2, A)') "   Mass (ton)                          : ", results%final_trunk%Mass_tonnes, " t"
        WRITE(u, '(A, F0.3, A)') "   Nominal Diameter (Dn)               : ", results%final_trunk%Dn, " m"
        IF (TRIM(results%coefficients%type_name) == "Antifer") THEN
            WRITE(u, '(A, F0.3, A)') "   Volume = 1.0247 * H^3 (V)           : ", armor_vol_trunk, " m3"
        ELSE
            WRITE(u, '(A, F0.3, A)') "   Volume = Dn^3 (V)                   : ", armor_vol_trunk, " m3"
            WRITE(u, '(A, F0.3, A)') "   Cube Height (H)                     : ", results%final_trunk%dims%H, " m"
            WRITE(u, '(A, F0.3, A)') "   Cube Top Width (B)                  : ", results%final_trunk%dims%B, " m"
            WRITE(u, '(A, F0.3, A)') "   Cube Base Width (A)                 : ", results%final_trunk%dims%A, " m"
        END IF
        WRITE(u, '(A, F0.2)')    "   KD_TRUNK (Equivalent)               : ", results%final_trunk%Kd
        WRITE(u, '(A, F0.2, A)') "   Double Layer Thickness (r1)         : ", results%final_trunk%r1, " m"
        WRITE(u, '(A, F0.2)')    "   Packing Density, d [units/100m2]    : ", results%final_trunk%packing_density
        WRITE(u, '(A)') ""
        
        WRITE(u, '(A)') "4. UNDERLAYER RESULTS - TRUNK"
        WRITE(u, '(A, F0.2, A, F0.1, A)') "   Theoretical Target (W/10)           : ", &
            results%underlayer_trunk%target_W, " kN (", &
            results%underlayer_trunk%target_M50_kg, " kg)"
        WRITE(u, '(A, A)')                "   Adopted rock grading                : ", TRIM(results%underlayer_trunk%grading_name)
        WRITE(u, '(A, F0.1, A)')          "   Representative M50                  : ", results%underlayer_trunk%M50_kg, " kg"
        WRITE(u, '(A, F0.1, A)')          "   Nominal lower limit (NLL)           : ", results%underlayer_trunk%NLL_kg, " kg"
        WRITE(u, '(A, F0.1, A)')          "   Nominal upper limit (NUL)           : ", results%underlayer_trunk%NUL_kg, " kg"
        WRITE(u, '(A, F0.1, A)')          "   Extreme lower limit (ELL)           : ", results%underlayer_trunk%ELL_kg, " kg"
        WRITE(u, '(A, F0.1, A)')          "   Extreme upper limit (EUL)           : ", results%underlayer_trunk%EUL_kg, " kg"
        WRITE(u, '(A, F0.3, A)')          "   Nominal Diameter (Dn_rock)          : ", results%underlayer_trunk%Dn_rock, " m"
        WRITE(u, '(A, F0.2, A)')          "   Double Layer Thickness (r2)         : ", results%underlayer_trunk%r2, " m"
        WRITE(u, '(A, F0.2)')             "   Packing Density, f2 [rocks/100m2]   : ", results%underlayer_trunk%f2
        WRITE(u, '(A)') TRIM(separator)
        
        WRITE(u, '(A)') "5. ARMOR LAYER RESULTS - HEAD (High Density)"
        WRITE(u, '(A)') "   *Maintains same Dn and Slope as Trunk*"
        WRITE(u, '(A, F0.2)')    "   Stability Ratio (Kd_T/Kd_H)         : ", results%final_head%Kd_Ratio
        WRITE(u, '(A, F0.3, A)') "   Nominal Diameter (Dn)               : ", results%final_head%Dn, " m"
        IF (TRIM(results%coefficients%type_name) == "Antifer") THEN
            WRITE(u, '(A, F0.3, A)') "   Volume = 1.0247 * H^3 (V)           : ", armor_vol_head, " m3"
        ELSE
            WRITE(u, '(A, F0.3, A)') "   Volume = Dn^3 (V)                   : ", armor_vol_head, " m3"
            WRITE(u, '(A, F0.3, A)') "   Cube Height (H)                     : ", results%final_head%dims%H, " m"
           WRITE(u, '(A, F0.3, A)') "   Cube Top width (B)                  : ", results%final_head%dims%B, " m"
           WRITE(u, '(A, F0.3, A)') "   Cube Base Width (A)                 : ", results%final_head%dims%A, " m"
        END IF
        WRITE(u, '(A, F0.2)')    "   KD_HEAD (Equivalent)                : ", results%final_head%Kd
        WRITE(u, '(A, F0.2, A)') "   Required Concrete Density (Wc)      : ", results%final_head%Wc_Required, " kN/m3"
        WRITE(u, '(A, F0.2, A)') "   BLOCK WEIGHT (W)                    : ", results%final_head%W, " kN"
        WRITE(u, '(A, F0.2, A)') "   Mass (ton)                          : ", results%final_head%Mass_tonnes, " t"
        WRITE(u, '(A, F0.2)')    "   Packing Density, d [units/100m2]    : ", results%final_head%packing_density
        WRITE(u, '(A)') ""
        
        WRITE(u, '(A)') "6. UNDERLAYER RESULTS - HEAD"
        WRITE(u, '(A, F0.2, A, F0.1, A)') "   Theoretical Target (W/10)           : ", &
            results%underlayer_head%target_W, " kN (", &
            results%underlayer_head%target_M50_kg, " kg)"
        WRITE(u, '(A, A)')                "   Adopted rock grading                : ", TRIM(results%underlayer_head%grading_name)
        WRITE(u, '(A, F0.1, A)')          "   Representative M50                  : ", results%underlayer_head%M50_kg, " kg"
        WRITE(u, '(A, F0.1, A)')          "   Nominal lower limit (NLL)           : ", results%underlayer_head%NLL_kg, " kg"
        WRITE(u, '(A, F0.1, A)')          "   Nominal upper limit (NUL)           : ", results%underlayer_head%NUL_kg, " kg"
        WRITE(u, '(A, F0.1, A)')          "   Extreme lower limit (ELL)           : ", results%underlayer_head%ELL_kg, " kg"
        WRITE(u, '(A, F0.1, A)')          "   Extreme upper limit (EUL)           : ", results%underlayer_head%EUL_kg, " kg"
        WRITE(u, '(A, F0.3, A)')          "   Nominal Diameter (Dn_rock)          : ", results%underlayer_head%Dn_rock, " m"
        WRITE(u, '(A, F0.2, A)')          "   Double Layer Thickness (r2)         : ", results%underlayer_head%r2, " m"
        WRITE(u, '(A, F0.2)')             "   Packing Density, f2 [rocks/100m2]   : ", results%underlayer_head%f2
        WRITE(u, '(A)') "================================================================================"
        
        CLOSE(u)
        
        PRINT *, ""
        PRINT *, " Report generated successfully: ", TRIM(filepath)
        PRINT *, "Report content:"
        PRINT *, ""
        
        OPEN(NEWUNIT=u, FILE=filepath, STATUS='OLD', ACTION='READ')
        DO
            READ(u, '(A)', IOSTAT=iostatus_val) line
            IF (iostatus_val /= 0) EXIT
            PRINT '(A)', TRIM(line)
        END DO
        CLOSE(u)
    END SUBROUTINE generate_report_file

    ! Helper to get param with default
    SUBROUTINE get_param(prompt, val)
        CHARACTER(LEN=*), INTENT(IN) :: prompt
        REAL(dp), INTENT(INOUT) :: val
        CHARACTER(LEN=100) :: buffer
        INTEGER :: iostatus_val
        
        WRITE(*, '(A, " [", F0.2, "]: ")', ADVANCE='NO') prompt, val
        READ(*, '(A)') buffer
        
        IF (LEN_TRIM(buffer) > 0) THEN
            READ(buffer, *, IOSTAT=iostatus_val) val
            IF (iostatus_val /= 0) THEN
                ! If read fails, keep default
            END IF
        END IF
    END SUBROUTINE get_param

END PROGRAM BreakwaterCalculator