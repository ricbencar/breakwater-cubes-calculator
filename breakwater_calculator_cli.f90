! ======================================================================================
! PROGRAM DESCRIPTION & METHODOLOGY
! ======================================================================================
!
! 1. PURPOSE:
!    This software calculates the required size and weight of artificial concrete 
!    armor units (specifically Cubes and Antifer blocks) for rubble mound breakwaters.
!    It implements an "Iso-Geometric Design Philosophy" to dimension two sections:
!    - The Trunk (the main longitudinal section).
!    - The Head (the exposed roundhead).
!
! 2. METHODOLOGY:
!    The calculator utilizes empirical power-law formulas derived from hydraulic 
!    model testing and European standards:
!    - For Simple Cubes: Uses Van der Meer (1988) stability formulas.
!    - For Antifer Blocks: Uses Chegini & Aghtouman (2006) power-law formulas.
!    - For Underlayers: Automatic sizing based on EN 13383 Standard Rock Gradings.
!
! 3. DESIGN LOGIC & STRATEGY:
!    a. Trunk Stability: 
!       Calculates the Stability Number (Ns) based on wave steepness (s0m), storm 
!       duration (Nz), and damage level (Nod). This determines the nominal 
!       diameter (Dn) and Weight (W).
!
!    b. Head Design (Iso-Geometric Transfer Strategy): 
!       The breakwater head endures higher turbulence (3D effects) than the trunk.
!       Rather than increasing block size or softening the slope—which requires 
!       different molds or crane logistics—this tool solves for the required 
!       increase in Concrete Density (rho_c).
!       
!       - Principle: Maintain constant Geometry (Dn) and Slope (alpha).
!       - Transfer Function: Delta_head = Delta_trunk * (1.5)^(1/3)
!       - Result: Head blocks use the same molds but are heavier due to higher density.
!
! 4. TECHNICAL BIBLIOGRAPHY & REFERENCES:
!
! 1. Van der Meer, J.W. (1988). "Rock Slopes and Gravel Beaches Under Wave Attack."
!    Doctoral Thesis, Delft University of Technology.
!
! 2. Van der Meer, J.W. (1988). "Stability of Cubes, Tetrapods and Accropode."
!    Proceedings of the Conference Breakwaters '88, Eastbourne, Thomas Telford.
!
! 3. Chegini, V., & Aghtouman, P. (2006). "An Investigation on the Stability of 
!    Rubble Mound Breakwaters with Armour Layers of Antifer Cubes." 
!    Journal of Marine Engineering.
!
! 4. USACE (2006). "Coastal Engineering Manual (CEM)", Chapter VI-5.
!
! 5. CEN (2002). "EN 13383-1: Armourstone - Part 1: Specification."
!
! 5. COMPILATION INSTRUCTIONS:
!
! Compilation (example using gfortran on Windows/Linux):
!
! gfortran -O3 -march=native -std=f2008 -Wall -Wextra -pedantic -Wconversion -static
! -static-libgfortran -static-libgcc -o breakwater_calculator_cli breakwater_calculator_cli.f90
!
! 6. EXECUTION EXAMPLES:
!
! 1. Interactive Mode (No arguments):
! ./breakwater_calculator_cli
!
! 2. Command Line Arguments Mode:
! Format: ./breakwater_calculator_cli [Hs] [Tm] [Duration] [Nod] [Wc] [FormulaID]
!
! Example (Hs=10.0, Tm=13.0, Duration=12.0, Nod=1.0, Wc=24.0, FormulaID=1):
! ./breakwater_calculator_cli 10.0 13.0 12.0 1.0 24.0 1
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
        REAL(dp) :: min_kg
        REAL(dp) :: max_kg
    END TYPE GradingDef

    TYPE :: Dimensions
        REAL(dp) :: H
        REAL(dp) :: A
        REAL(dp) :: B
    END TYPE Dimensions

    TYPE :: UnderlayerResult
        CHARACTER(LEN=64) :: grading_name
        REAL(dp) :: target_W
        REAL(dp) :: W_mean
        REAL(dp) :: W1
        REAL(dp) :: W2
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
        REAL(dp) :: delta
        REAL(dp) :: Ns_trunk
    END TYPE Intermediates

    TYPE :: Inputs
        REAL(dp) :: Hs
        REAL(dp) :: Tm
        REAL(dp) :: Storm_Duration_hr
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
    TYPE(GradingDef) :: standard_gradings(16)
    
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
    
    ! ID 1: Van Der Meer (1988a) - Cubes (Slope 2.0:1)
    formulas(1)%name = "Van Der Meer (1988a) - Cubes (Slope 2.0:1)"
    formulas(1)%type_name = "Cubes"
    formulas(1)%slope_ratio = 2.0_dp
    formulas(1)%k1 = 6.7_dp
    formulas(1)%k2 = 0.4_dp
    formulas(1)%k3 = 0.3_dp
    formulas(1)%k4 = 1.0_dp
    formulas(1)%k5 = 0.1_dp
    
    ! ID 2: Van Der Meer (1988a) - Cubes (Slope 1.5:1)
    formulas(2)%name = "Van Der Meer (1988a) - Cubes (Slope 1.5:1)"
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
    standard_gradings(1)%name="CP45/125";      standard_gradings(1)%min_kg=0.4_dp;    standard_gradings(1)%max_kg=1.2_dp
    standard_gradings(2)%name="CP63/180";      standard_gradings(2)%min_kg=1.2_dp;    standard_gradings(2)%max_kg=3.8_dp
    standard_gradings(3)%name="CP90/250";      standard_gradings(3)%min_kg=3.1_dp;    standard_gradings(3)%max_kg=9.3_dp
    standard_gradings(4)%name="CP45/180";      standard_gradings(4)%min_kg=0.4_dp;    standard_gradings(4)%max_kg=1.2_dp
    standard_gradings(5)%name="CP90/180";      standard_gradings(5)%min_kg=2.1_dp;    standard_gradings(5)%max_kg=2.8_dp
    ! Light Mass Armourstone (LMA)
    standard_gradings(6)%name="LMA5-40";       standard_gradings(6)%min_kg=10.0_dp;   standard_gradings(6)%max_kg=20.0_dp
    standard_gradings(7)%name="LMA10-60";      standard_gradings(7)%min_kg=20.0_dp;   standard_gradings(7)%max_kg=35.0_dp
    standard_gradings(8)%name="LMA15-120";     standard_gradings(7)%min_kg=35.0_dp;   standard_gradings(7)%max_kg=60.0_dp
    standard_gradings(9)%name="LMA40-200";     standard_gradings(8)%min_kg=80.0_dp;   standard_gradings(8)%max_kg=120.0_dp
    standard_gradings(10)%name="LMA60-300";     standard_gradings(9)%min_kg=120.0_dp;  standard_gradings(9)%max_kg=190.0_dp
    standard_gradings(11)%name="LMA15-300";    standard_gradings(10)%min_kg=45.0_dp;  standard_gradings(10)%max_kg=135.0_dp
    ! Heavy Masss Armourstone (HMA)
    standard_gradings(12)%name="HMA300-1000";  standard_gradings(11)%min_kg=540.0_dp; standard_gradings(11)%max_kg=690.0_dp
    standard_gradings(13)%name="HMA1000-3000"; standard_gradings(12)%min_kg=1700.0_dp;standard_gradings(12)%max_kg=2100.0_dp
    standard_gradings(14)%name="HMA3000-6000"; standard_gradings(13)%min_kg=4200.0_dp;standard_gradings(13)%max_kg=4800.0_dp
    standard_gradings(15)%name="HMA6000-10000";standard_gradings(14)%min_kg=7500.0_dp;standard_gradings(14)%max_kg=8500.0_dp
    standard_gradings(16)%name="HMA10000-15000";standard_gradings(15)%min_kg=12000.0_dp;standard_gradings(15)%max_kg=13000.0_dp

    ! ----------------------------------------------------------------------
    ! DEFAULT INPUTS
    ! ----------------------------------------------------------------------
    user_inputs%Hs = 10.0_dp
    user_inputs%Tm = 13.0_dp
    user_inputs%Storm_Duration_hr = 12.0_dp
    user_inputs%Nod = 1.0_dp
    user_inputs%Wc = 24.0_dp
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
        CALL GET_COMMAND_ARGUMENT(3, arg_val); READ(arg_val, *) user_inputs%Storm_Duration_hr
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
        PRINT *, "1. Van Der Meer (1988) - Cubes (Slope 2.0:1)"
        PRINT *, "2. Van Der Meer (1988) - Cubes (Slope 1.5:1)"
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
        CALL get_param("Nod (Damage)", user_inputs%Nod)
        CALL get_param("Duration (h)", user_inputs%Storm_Duration_hr)
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
        REAL(dp) :: min_diff, diff_val
        REAL(dp) :: w_min, w_max, w_mean_range
        REAL(dp) :: final_w_mean, final_w_min, final_w_max
        INTEGER :: i

        CHARACTER(LEN=64) :: best_name
        LOGICAL :: found

        target_weight = W_armor / 10.0_dp
        min_diff = HUGE(1.0_dp)
        found = .FALSE.
        
        ! Initialize defaults
        best_name = standard_gradings(1)%name
        final_w_mean = 0.0_dp
        final_w_min = 0.0_dp
        final_w_max = 0.0_dp
        
        DO i = 1, SIZE(standard_gradings)
            w_min = standard_gradings(i)%min_kg * g / 1000.0_dp
            w_max = standard_gradings(i)%max_kg * g / 1000.0_dp
            w_mean_range = (w_min + w_max) / 2.0_dp
            diff_val = ABS(w_mean_range - target_weight)
            
            IF (diff_val < min_diff) THEN
                min_diff = diff_val
                best_name = standard_gradings(i)%name
                final_w_mean = w_mean_range
                final_w_min = w_min
                final_w_max = w_max
                found = .TRUE.
            END IF
        END DO
        
        IF (.NOT. found) THEN
            final_w_min = standard_gradings(1)%min_kg * g / 1000.0_dp
            final_w_max = standard_gradings(1)%max_kg * g / 1000.0_dp
            final_w_mean = (final_w_min + final_w_max) / 2.0_dp
        END IF

        res%grading_name = best_name
        res%target_W = target_weight
        res%W_mean = final_w_mean
        res%W1 = final_w_min
        res%W2 = final_w_max
        res%Dn_rock = (final_w_mean / W_rock_spec_val)**(1.0_dp/3.0_dp)
        res%r2 = 2.0_dp * res%Dn_rock
        res%f2 = 100.0_dp * 2.0_dp * 1.0_dp * (1.0_dp - P_rock) / (res%Dn_rock**2)
        res%W_rock_spec = W_rock_spec_val
    END FUNCTION calculate_underlayer_params

    SUBROUTINE solve(formula_id, params, results_out)
        INTEGER, INTENT(IN) :: formula_id
        TYPE(Inputs), INTENT(IN) :: params
        TYPE(FullResults), INTENT(OUT) :: results_out
        
        TYPE(FormulaParams) :: coeffs
        
        REAL(dp) :: L0, k0, s0m, sm, Nz, delta_trunk
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
        
        Nz = (params%Storm_Duration_hr * 3600.0_dp) / params%Tm
        delta_trunk = (params%Wc / params%Ww) - 1.0_dp
        
        ! 4. Algorithmic Core (Chegini-Aghtouman / Van der Meer) - TRUNK
        term_damage = params%Nod**coeffs%k2
        term_waves = Nz**coeffs%k3
        damage_wave_ratio = term_damage / term_waves
        
        scaled_term = coeffs%k1 * damage_wave_ratio
        inv_f = scaled_term + coeffs%k4
        
        steepness_factor = s0m**(-coeffs%k5)
        
        ! Updated Condition: 
        IF (formula_id == 1) THEN
            Ns_trunk = inv_f * steepness_factor * ((2.0_dp/1.5_dp)**(1.0_dp/3.0_dp))
        ELSE
            Ns_trunk = inv_f * steepness_factor
        END IF
        
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
        results_out%intermediate%delta = delta_trunk
        results_out%intermediate%Ns_trunk = Ns_trunk
        
        results_out%final_trunk%Dn = Dn
        results_out%final_trunk%W = W_trunk
        results_out%final_trunk%Mass_tonnes = W_trunk / g
        results_out%final_trunk%Kd = kd_trunk_equiv
        results_out%final_trunk%packing_density = packing_density_trunk
        results_out%final_trunk%dims%H = h_trunk
        results_out%final_trunk%dims%A = a_trunk
        results_out%final_trunk%dims%B = b_trunk
        
        results_out%final_head%Ns = Ns_head
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
        
        separator = "--------------------------------------------------------------------------------"
        
        OPEN(NEWUNIT=u, FILE=filepath, STATUS='REPLACE', ACTION='WRITE', IOSTAT=iostatus_val)
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
        WRITE(u, '(A, F0.2, A)') "   Storm Duration (h)                  : ", results%inputs%Storm_Duration_hr, " h"
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
        WRITE(u, '(A, F6.4)')    "   wave number (k0 = 2*pi/L0)          : ", results%intermediate%k0
        WRITE(u, '(A, F6.4)')    "   wave steepness (s0m = Hs/L0)        : ", results%intermediate%s0m
        WRITE(u, '(A, I0)')      "   Number of waves (Nz)                : ", INT(results%intermediate%Nz)
        WRITE(u, '(A, F0.4)')    "   Stability Number TRUNK (Ns)         : ", results%intermediate%Ns_trunk
        WRITE(u, '(A, F0.4)')    "   Stability Number HEAD (Ns)          : ", results%final_head%Ns
        WRITE(u, '(A)') TRIM(separator)
        
        WRITE(u, '(A)') "3. ARMOR LAYER RESULTS - TRUNK"
        WRITE(u, '(A, F0.2, A)') "   BLOCK WEIGHT (W)                    : ", results%final_trunk%W, " kN"
        WRITE(u, '(A, F0.2, A)') "   Mass (ton)                          : ", results%final_trunk%Mass_tonnes, " t"
        WRITE(u, '(A, F0.3, A)') "   Nominal Diameter (Dn)               : ", results%final_trunk%Dn, " m"
        WRITE(u, '(A, F0.3, A)') "   Cube Height (H)                     : ", results%final_trunk%dims%H, " m"
        WRITE(u, '(A, F0.3, A)') "   Cube Top Width (B)                  : ", results%final_trunk%dims%B, " m"
        WRITE(u, '(A, F0.3, A)') "   Cube Base Width (A)                 : ", results%final_trunk%dims%A, " m"
        WRITE(u, '(A, F0.2)')    "   KD_TRUNK (Equivalent)               : ", results%final_trunk%Kd
        WRITE(u, '(A, F0.2, A)') "   Double Layer Thickness (r1)         : ", results%final_trunk%r1, " m"
        WRITE(u, '(A, F0.2)')    "   Packing Density, d [units/100m2]    : ", results%final_trunk%packing_density
        WRITE(u, '(A)') ""
        
        WRITE(u, '(A)') "4. UNDERLAYER RESULTS - TRUNK"
        WRITE(u, '(A, F0.2, A)') "   Theoretical Target (W/10)           : ", results%underlayer_trunk%target_W, " kN"
        WRITE(u, '(A, A)')       "   Adopted rock grading                : ", TRIM(results%underlayer_trunk%grading_name)
        WRITE(u, '(A, F0.2, A)') "   Grading Min (Lower Limit)           : ", results%underlayer_trunk%W1, " kN"
        WRITE(u, '(A, F0.2, A)') "   Grading Max (Upper Limit)           : ", results%underlayer_trunk%W2, " kN"
        WRITE(u, '(A, F0.2, A)') "   Mean Weight (Used for thickness)    : ", results%underlayer_trunk%W_mean, " kN"
        WRITE(u, '(A, F0.3, A)') "   Nominal Diameter (Dn_rock)          : ", results%underlayer_trunk%Dn_rock, " m"
        WRITE(u, '(A, F0.2, A)') "   Double Layer Thickness (r2)         : ", results%underlayer_trunk%r2, " m"
        WRITE(u, '(A, F0.2)')    "   Packing Density, f2 [rocks/100m2]   : ", results%underlayer_trunk%f2
        WRITE(u, '(A)') TRIM(separator)
        
        WRITE(u, '(A)') "5. ARMOR LAYER RESULTS - HEAD (High Density)"
        WRITE(u, '(A)') "   *Maintains same Dn and Slope as Trunk*"
        WRITE(u, '(A, F0.2)')    "   Stability Ratio (Kd_T/Kd_H)         : ", results%final_head%Kd_Ratio
        WRITE(u, '(A, F0.3, A)') "   Nominal Diameter (Dn)               : ", results%final_trunk%Dn, " m"
        WRITE(u, '(A, F0.3, A)') "   Cube Height (H)                     : ", results%final_head%dims%H, " m"
        WRITE(u, '(A, F0.3, A)') "   Cube Top width (B)                  : ", results%final_head%dims%B, " m"
        WRITE(u, '(A, F0.3, A)') "   Cube Base Width (A)                 : ", results%final_head%dims%A, " m"
        WRITE(u, '(A, F0.2)')    "   KD_HEAD (Equivalent)                : ", results%final_head%Kd
        WRITE(u, '(A, F0.2, A)') "   Required Concrete Density (Wc)      : ", results%final_head%Wc_Required, " kN/m3"
        WRITE(u, '(A, F0.2, A)') "   BLOCK WEIGHT (W)                    : ", results%final_head%W, " kN"
        WRITE(u, '(A, F0.2, A)') "   Mass (ton)                          : ", results%final_head%Mass_tonnes, " t"
        WRITE(u, '(A, F0.2)')    "   Packing Density, d [units/100m2]    : ", results%final_head%packing_density
        WRITE(u, '(A)') ""
        
        WRITE(u, '(A)') "6. UNDERLAYER RESULTS - HEAD"
        WRITE(u, '(A, F0.2, A)') "   Theoretical Target (W/10)           : ", results%underlayer_head%target_W, " kN"
        WRITE(u, '(A, A)')       "   Adopted rock grading                : ", TRIM(results%underlayer_head%grading_name)
        WRITE(u, '(A, F0.2, A)') "   Grading Min (Lower Limit)           : ", results%underlayer_head%W1, " kN"
        WRITE(u, '(A, F0.2, A)') "   Grading Max (Upper Limit)           : ", results%underlayer_head%W2, " kN"
        WRITE(u, '(A, F0.2, A)') "   Mean Weight (Used for thickness)    : ", results%underlayer_head%W_mean, " kN"
        WRITE(u, '(A, F0.3, A)') "   Nominal Diameter (Dn_rock)          : ", results%underlayer_head%Dn_rock, " m"
        WRITE(u, '(A, F0.2, A)') "   Double Layer Thickness (r2)         : ", results%underlayer_head%r2, " m"
        WRITE(u, '(A, F0.2)')    "   Packing Density, f2 [rocks/100m2]   : ", results%underlayer_head%f2
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