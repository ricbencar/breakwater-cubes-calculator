! ======================================================================================
! PROGRAM DESCRIPTION & METHODOLOGY
! ======================================================================================
!
! 1. PURPOSE:
!    This software performs preliminary hydraulic sizing of rubble-mound
!    breakwater armor made with artificial concrete units and the associated
!    rock underlayers.
!
!    The currently implemented armor-unit families are:
!    - Simple Cubes
!    - Antifer Blocks
!
!    The program dimensions two structural zones:
!    - The Trunk (main longitudinal section)
!    - The Head (exposed roundhead / end section)
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
!    - standard EN 13383 grading selection from an embedded mass-based
!      grading table, with optional custom interpolated grading by family
!      (AUTO / HMA / LMA / CP) for the rock underlayers only.
!
! ======================================================================================

program BreakwaterCalculator
    implicit none

    integer, parameter :: dp = selected_real_kind(15, 307)
    real(dp), parameter :: PI = 3.14159265358979323846_dp
    real(dp), parameter :: g = 9.80665_dp
    real(dp), parameter :: W_rock_spec_val = 26.5_dp
    real(dp), parameter :: P_cubes = 0.40_dp
    real(dp), parameter :: P_rock = 0.25_dp
    real(dp), parameter :: KD_RATIO_FIXED = 1.5_dp

    type :: FormulaParams
        character(len=100) :: name = ""
        character(len=64)  :: type_name = ""
        real(dp) :: slope_ratio = 0.0_dp
        real(dp) :: k1 = 0.0_dp
        real(dp) :: k2 = 0.0_dp
        real(dp) :: k3 = 0.0_dp
        real(dp) :: k4 = 0.0_dp
        real(dp) :: k5 = 0.0_dp
    end type FormulaParams

    type :: GradingDef
        character(len=64) :: name = ""
        real(dp) :: NLL_kg = 0.0_dp
        real(dp) :: NUL_kg = 0.0_dp
    end type GradingDef

    type :: Dimensions
        real(dp) :: H = 0.0_dp
        real(dp) :: A = 0.0_dp
        real(dp) :: B = 0.0_dp
    end type Dimensions

    type :: UnderlayerResult
        character(len=64) :: grading_name = ""
        real(dp) :: target_W = 0.0_dp
        real(dp) :: target_M50_kg = 0.0_dp
        real(dp) :: M50_kg = 0.0_dp
        real(dp) :: ELL_kg = 0.0_dp
        real(dp) :: EUL_kg = 0.0_dp
        real(dp) :: NLL_kg = 0.0_dp
        real(dp) :: NUL_kg = 0.0_dp
        real(dp) :: W_mean_kn = 0.0_dp
        real(dp) :: Dn_rock = 0.0_dp
        real(dp) :: r2 = 0.0_dp
        real(dp) :: f2 = 0.0_dp
        real(dp) :: W_rock_spec = 0.0_dp
        logical  :: used_custom_interpolation = .false.
        character(len=16) :: custom_family = ""
        real(dp) :: custom_ratio_nul_nll = 0.0_dp
        character(len=128) :: custom_ratio_note = ""
    end type UnderlayerResult

    type :: ArmorResult
        real(dp) :: Ns = 0.0_dp
        real(dp) :: Dn = 0.0_dp
        real(dp) :: W = 0.0_dp
        real(dp) :: Mass_tonnes = 0.0_dp
        real(dp) :: Kd = 0.0_dp
        real(dp) :: r1 = 0.0_dp
        real(dp) :: packing_density = 0.0_dp
        type(Dimensions) :: dims
        real(dp) :: Kd_Ratio = 0.0_dp
        real(dp) :: Delta_Required = 0.0_dp
        real(dp) :: Wc_Required = 0.0_dp
    end type ArmorResult

    type :: Intermediates
        real(dp) :: L0 = 0.0_dp
        real(dp) :: k0 = 0.0_dp
        real(dp) :: s0m = 0.0_dp
        real(dp) :: sm = 0.0_dp
        real(dp) :: Nz = 0.0_dp
        real(dp) :: Storm_Duration_hr = 0.0_dp
        real(dp) :: delta = 0.0_dp
        real(dp) :: Ns_trunk = 0.0_dp
    end type Intermediates

    type :: Inputs
        real(dp) :: Hs = 0.0_dp
        real(dp) :: Tm = 0.0_dp
        real(dp) :: Number_of_Waves = 0.0_dp
        real(dp) :: Nod = 0.0_dp
        real(dp) :: Wc = 0.0_dp
        real(dp) :: Ww = 0.0_dp
        integer  :: Formula_ID = 1
        logical  :: use_en13383 = .true.
        character(len=16) :: custom_family = "AUTO"
    end type Inputs

    type :: FullResults
        type(Inputs) :: inputs
        type(FormulaParams) :: coefficients
        type(Intermediates) :: intermediate
        type(ArmorResult) :: final_trunk
        type(UnderlayerResult) :: underlayer_trunk
        type(ArmorResult) :: final_head
        type(UnderlayerResult) :: underlayer_head
        real(dp) :: P_rock_val = 0.0_dp
        real(dp) :: P_cubes_val = 0.0_dp
    end type FullResults

    type(FormulaParams) :: formulas(4)
    type(GradingDef)    :: standard_gradings(17)
    type(Inputs)        :: defaults, user_inputs
    type(FullResults)   :: results
    integer :: num_args, formula_id
    character(len=256) :: arg

    call init_data(formulas, standard_gradings, defaults)
    user_inputs = defaults
    formula_id = defaults%Formula_ID

    num_args = command_argument_count()

    if (num_args >= 6) then
        if (.not. read_real_arg(1, user_inputs%Hs)) call cli_parse_fail()
        if (.not. read_real_arg(2, user_inputs%Tm)) call cli_parse_fail()
        if (.not. read_real_arg(3, user_inputs%Number_of_Waves)) call cli_parse_fail()
        if (.not. read_real_arg(4, user_inputs%Nod)) call cli_parse_fail()
        if (.not. read_real_arg(5, user_inputs%Wc)) call cli_parse_fail()
        if (.not. read_int_arg(6, formula_id)) call cli_parse_fail()

        if (formula_id < 1 .or. formula_id > 4) then
            write(*,'(A)') 'Error: Formula ID must be 1-4. Using default (1).'
            formula_id = 1
        end if
        user_inputs%Formula_ID = formula_id

        if (num_args >= 7) then
            call get_command_argument(7, arg)
            user_inputs%use_en13383 = parse_bool(trim(arg), defaults%use_en13383)
        end if
        if (num_args >= 8) then
            call get_command_argument(8, arg)
            user_inputs%custom_family = normalize_family(trim(arg))
        end if
    else
        call run_interactive(defaults, user_inputs, formula_id)
    end if

    call solve(formulas, standard_gradings, formula_id, user_inputs, results)
    call generate_report_file(results, 'output.txt')

contains

    subroutine init_data(formulas, gradings, defaults)
        type(FormulaParams), intent(out) :: formulas(4)
        type(GradingDef), intent(out)    :: gradings(17)
        type(Inputs), intent(out)        :: defaults

        formulas(1)%name = 'Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)'
        formulas(1)%type_name = 'Cubes'
        formulas(1)%slope_ratio = 2.0_dp
        formulas(1)%k1 = 7.374304189198_dp
        formulas(1)%k2 = 0.4_dp
        formulas(1)%k3 = 0.3_dp
        formulas(1)%k4 = 1.100642416298_dp
        formulas(1)%k5 = 0.1_dp

        formulas(2)%name = 'Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)'
        formulas(2)%type_name = 'Cubes'
        formulas(2)%slope_ratio = 1.5_dp
        formulas(2)%k1 = 6.7_dp
        formulas(2)%k2 = 0.4_dp
        formulas(2)%k3 = 0.3_dp
        formulas(2)%k4 = 1.0_dp
        formulas(2)%k5 = 0.1_dp

        formulas(3)%name = 'Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)'
        formulas(3)%type_name = 'Antifer'
        formulas(3)%slope_ratio = 2.0_dp
        formulas(3)%k1 = 6.138_dp
        formulas(3)%k2 = 0.443_dp
        formulas(3)%k3 = 0.276_dp
        formulas(3)%k4 = 1.164_dp
        formulas(3)%k5 = 0.07_dp

        formulas(4)%name = 'Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)'
        formulas(4)%type_name = 'Antifer'
        formulas(4)%slope_ratio = 1.5_dp
        formulas(4)%k1 = 6.951_dp
        formulas(4)%k2 = 0.443_dp
        formulas(4)%k3 = 0.291_dp
        formulas(4)%k4 = 1.082_dp
        formulas(4)%k5 = 0.082_dp

        gradings(1)  = GradingDef('CP32/90',          0.868_dp,   19.319_dp)
        gradings(2)  = GradingDef('CP45/125',         2.415_dp,   51.758_dp)
        gradings(3)  = GradingDef('CP63/180',         6.626_dp,  154.548_dp)
        gradings(4)  = GradingDef('CP90/250',        19.319_dp,  414.063_dp)
        gradings(5)  = GradingDef('CP45/180',         2.415_dp,  154.548_dp)
        gradings(6)  = GradingDef('CP90/180',        19.319_dp,  154.548_dp)
        gradings(7)  = GradingDef('LMA5-40',          5.0_dp,     40.0_dp)
        gradings(8)  = GradingDef('LMA10-60',        10.0_dp,     60.0_dp)
        gradings(9)  = GradingDef('LMA15-120',       15.0_dp,    120.0_dp)
        gradings(10) = GradingDef('LMA40-200',       40.0_dp,    200.0_dp)
        gradings(11) = GradingDef('LMA60-300',       60.0_dp,    300.0_dp)
        gradings(12) = GradingDef('LMA15-300',       15.0_dp,    300.0_dp)
        gradings(13) = GradingDef('HMA300-1000',    300.0_dp,   1000.0_dp)
        gradings(14) = GradingDef('HMA1000-3000',  1000.0_dp,   3000.0_dp)
        gradings(15) = GradingDef('HMA3000-6000',  3000.0_dp,   6000.0_dp)
        gradings(16) = GradingDef('HMA6000-10000', 6000.0_dp,  10000.0_dp)
        gradings(17) = GradingDef('HMA10000-15000',10000.0_dp, 15000.0_dp)

        defaults%Hs = 11.0_dp
        defaults%Tm = 11.9_dp
        defaults%Number_of_Waves = 3000.0_dp
        defaults%Nod = 0.5_dp
        defaults%Wc = 27.48_dp
        defaults%Ww = 10.05_dp
        defaults%Formula_ID = 1
        defaults%use_en13383 = .true.
        defaults%custom_family = 'AUTO'
    end subroutine init_data

    logical function read_real_arg(idx, val)
        integer, intent(in) :: idx
        real(dp), intent(out) :: val
        character(len=256) :: s
        integer :: ios_local
        call get_command_argument(idx, s)
        read(s, *, iostat=ios_local) val
        read_real_arg = (ios_local == 0)
    end function read_real_arg

    logical function read_int_arg(idx, val)
        integer, intent(in) :: idx
        integer, intent(out) :: val
        character(len=256) :: s
        integer :: ios_local
        call get_command_argument(idx, s)
        read(s, *, iostat=ios_local) val
        read_int_arg = (ios_local == 0)
    end function read_int_arg

    subroutine cli_parse_fail()
        write(*,'(A)') 'Error parsing command line arguments.'
        write(*,'(A)') 'Usage: ./breakwater_calculator_cli [Hs] [Tm] [NumberOfWaves] [Nod] [Wc] [FormulaID] '// &
            '[UseEN13383] [CustomFamily]'
        stop 1
    end subroutine cli_parse_fail

    subroutine run_interactive(defaults, user_inputs, formula_id)
        type(Inputs), intent(in) :: defaults
        type(Inputs), intent(inout) :: user_inputs
        integer, intent(out) :: formula_id
        character(len=256) :: selection, s

        write(*,'(A)') ''
        write(*,'(A)') '--- COASTAL PROTECTION BLOCK CALCULATOR (TRUNK & HEAD) ---'
        write(*,'(A)') '1. Simple Cubes (Slope 2.0:1) - Van der Meer'
        write(*,'(A)') '2. Simple Cubes (Slope 1.5:1) - Van der Meer'
        write(*,'(A)') '3. Antifer (Slope 2.0:1) - Chegini'
        write(*,'(A)') '4. Antifer (Slope 1.5:1) - Chegini'
        write(*,'(A)') ''
        write(*,'(A)', advance='no') 'Option [1-4]: '
        read(*,'(A)') selection
        selection = trim(adjustl(selection))
        select case (selection)
            case ('1','2','3','4')
                read(selection,*) formula_id
            case default
                formula_id = 1
        end select
        user_inputs%Formula_ID = formula_id

        write(*,'(A)') ''
        write(*,'(A)') '--- Enter Parameters (Press ENTER for Default) ---'
        call get_param('Hs (m)', defaults%Hs, user_inputs%Hs)
        call get_param('Tm (s)', defaults%Tm, user_inputs%Tm)
        call get_param('Number of waves (Nz)', defaults%Number_of_Waves, user_inputs%Number_of_Waves)
        call get_param('Nod (Damage)', defaults%Nod, user_inputs%Nod)
        call get_param('Concrete Weight Trunk (kN/m3)', defaults%Wc, user_inputs%Wc)
        if (defaults%use_en13383) then
            call get_text_param('Use standard EN 13383 underlayer grading? [true/false]', 'true', s)
        else
            call get_text_param('Use standard EN 13383 underlayer grading? [true/false]', 'false', s)
        end if
        user_inputs%use_en13383 = parse_bool(s, defaults%use_en13383)
        call get_text_param('Custom underlayer family [AUTO/HMA/LMA/CP]', trim(defaults%custom_family), s)
        user_inputs%custom_family = normalize_family(s)
        user_inputs%Ww = defaults%Ww
    end subroutine run_interactive

    subroutine get_param(prompt, default_val, out_val)
        character(len=*), intent(in) :: prompt
        real(dp), intent(in) :: default_val
        real(dp), intent(out) :: out_val
        character(len=256) :: s
        integer :: ios_local
        out_val = default_val
        write(*,'(A," [",F0.2,"]: ")', advance='no') trim(prompt), default_val
        read(*,'(A)') s
        s = trim(adjustl(s))
        if (len_trim(s) > 0) then
            read(s, *, iostat=ios_local) out_val
            if (ios_local /= 0) out_val = default_val
        end if
    end subroutine get_param

    subroutine get_text_param(prompt, default_val, out_val)
        character(len=*), intent(in) :: prompt, default_val
        character(len=*), intent(out) :: out_val
        character(len=256) :: s
        write(*,'(A," [",A,"]: ")', advance='no') trim(prompt), trim(default_val)
        read(*,'(A)') s
        s = trim(adjustl(s))
        if (len_trim(s) == 0) then
            out_val = trim(default_val)
        else
            out_val = s
        end if
    end subroutine get_text_param

    pure function upper_copy(s) result(out)
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, c
        out = s
        do i = 1, len(s)
            c = iachar(out(i:i))
            if (c >= iachar('a') .and. c <= iachar('z')) out(i:i) = achar(c - 32)
        end do
    end function upper_copy

    pure logical function starts_with(value, prefix)
        character(len=*), intent(in) :: value, prefix
        integer :: lp
        lp = len_trim(prefix)
        if (len_trim(value) < lp) then
            starts_with = .false.
        else
            starts_with = (value(1:lp) == prefix(1:lp))
        end if
    end function starts_with

    pure function normalize_family(s) result(out)
        character(len=*), intent(in) :: s
        character(len=16) :: out
        character(len=256) :: t
        t = upper_copy(trim(adjustl(s)))
        select case (trim(t))
            case ('AUTO','HMA','LMA','CP')
                out = trim(t)
            case default
                out = 'AUTO'
        end select
    end function normalize_family

    logical function parse_bool(s, default_val)
        character(len=*), intent(in) :: s
        logical, intent(in) :: default_val
        character(len=256) :: t
        t = upper_copy(trim(adjustl(s)))
        if (len_trim(t) == 0) then
            parse_bool = default_val
        else if (trim(t) == '1' .or. trim(t) == 'TRUE' .or. trim(t) == 'T' .or. trim(t) == 'YES' .or. trim(t) == 'Y') then
            parse_bool = .true.
        else if (trim(t) == '0' .or. trim(t) == 'FALSE' .or. trim(t) == 'F' .or. trim(t) == 'NO' .or. trim(t) == 'N') then
            parse_bool = .false.
        else
            parse_bool = default_val
        end if
    end function parse_bool

    pure real(dp) function calculate_L0(Tm)
        real(dp), intent(in) :: Tm
        calculate_L0 = (g * Tm**2) / (2.0_dp * PI)
    end function calculate_L0

    pure function get_formula_display_name(formula_id) result(name)
        integer, intent(in) :: formula_id
        character(len=100) :: name
        select case (formula_id)
            case (1)
                name = 'Van der Meer (1988a) - Simple Cubes (Slope 2.0:1)'
            case (2)
                name = 'Van der Meer (1988a) - Simple Cubes (Slope 1.5:1)'
            case (3)
                name = 'Chegini-Aghtouman (2006) - Antifer (Slope 2.0:1)'
            case (4)
                name = 'Chegini-Aghtouman (2006) - Antifer (Slope 1.5:1)'
            case default
                name = 'Unknown formula'
        end select
    end function get_formula_display_name

    pure function grading_family(grading_name) result(family)
        character(len=*), intent(in) :: grading_name
        character(len=8) :: family
        if (starts_with(trim(grading_name), 'HMA')) then
            family = 'HMA'
        else if (starts_with(trim(grading_name), 'LMA')) then
            family = 'LMA'
        else if (starts_with(trim(grading_name), 'CP')) then
            family = 'CP'
        else
            family = 'UNKNOWN'
        end if
    end function grading_family

    subroutine get_family_gradings(gradings, family, out, nout)
        type(GradingDef), intent(in) :: gradings(:)
        character(len=*), intent(in) :: family
        type(GradingDef), intent(out) :: out(size(gradings))
        integer, intent(out) :: nout
        integer :: i, j
        type(GradingDef) :: tmp
        real(dp) :: ma, mb

        nout = 0
        do i = 1, size(gradings)
            if (trim(grading_family(gradings(i)%name)) == trim(family)) then
                nout = nout + 1
                out(nout) = gradings(i)
            end if
        end do
        if (nout <= 1) return
        do i = 1, nout - 1
            do j = i + 1, nout
                ma = 0.5_dp * (out(i)%NLL_kg + out(i)%NUL_kg)
                mb = 0.5_dp * (out(j)%NLL_kg + out(j)%NUL_kg)
                if (mb < ma) then
                    tmp = out(i)
                    out(i) = out(j)
                    out(j) = tmp
                end if
            end do
        end do
    end subroutine get_family_gradings

    function select_custom_family(gradings, target_mass) result(best_family)
        type(GradingDef), intent(in) :: gradings(:)
        real(dp), intent(in) :: target_mass
        character(len=8) :: best_family
        integer :: i
        real(dp) :: safe_mass, best_score, score, m50
        character(len=8) :: fam

        safe_mass = max(target_mass, 1.0e-9_dp)
        best_score = huge(1.0_dp)
        best_family = 'LMA'
        do i = 1, size(gradings)
            fam = grading_family(gradings(i)%name)
            if (trim(fam) == 'UNKNOWN') cycle
            m50 = 0.5_dp * (gradings(i)%NLL_kg + gradings(i)%NUL_kg)
            score = abs(log(safe_mass) - log(max(m50, 1.0e-9_dp)))
            if (trim(fam) == 'CP') score = score + 0.08_dp
            if (score < best_score) then
                best_score = score
                best_family = fam
            end if
        end do
    end function select_custom_family

    function interpolate_family_ratio(gradings, target_mass, family, note) result(ratio)
        type(GradingDef), intent(in) :: gradings(:)
        real(dp), intent(in) :: target_mass
        character(len=*), intent(in) :: family
        character(len=*), intent(out) :: note
        real(dp) :: ratio
        type(GradingDef) :: fam(size(gradings))
        integer :: n, i
        real(dp) :: x(size(gradings)), y(size(gradings)), xt, t, m50, log_ratio

        call get_family_gradings(gradings, family, fam, n)
        if (n == 0) then
            note = 'fallback ratio R=NUL/NLL = 3.0 (no family data)'
            ratio = 3.0_dp
            return
        end if

        do i = 1, n
            m50 = 0.5_dp * (fam(i)%NLL_kg + fam(i)%NUL_kg)
            x(i) = log(max(m50, 1.0e-9_dp))
            y(i) = log(max(fam(i)%NUL_kg / fam(i)%NLL_kg, 1.0_dp + 1.0e-9_dp))
        end do
        xt = log(max(target_mass, 1.0e-9_dp))

        if (n == 1) then
            ratio = fam(1)%NUL_kg / fam(1)%NLL_kg
            note = trim(family)//' family single-point ratio used'
            return
        end if

        if (xt <= x(1)) then
            ratio = exp(y(1))
            note = trim(family)//' family ratio clamped to lower-end class '//trim(fam(1)%name)
            return
        end if

        if (xt >= x(n)) then
            ratio = exp(y(n))
            note = trim(family)//' family ratio clamped to upper-end class '//trim(fam(n)%name)
            return
        end if

        do i = 1, n - 1
            if (xt >= x(i) .and. xt <= x(i+1)) then
                t = (xt - x(i)) / (x(i+1) - x(i))
                log_ratio = y(i) + t * (y(i+1) - y(i))
                ratio = exp(log_ratio)
                note = trim(family)//' family ratio interpolated between '//trim(fam(i)%name)//' and '//trim(fam(i+1)%name)
                return
            end if
        end do

        ratio = exp(y(n))
        note = trim(family)//' family ratio fallback to upper-end class '//trim(fam(n)%name)
    end function interpolate_family_ratio

    function calculate_underlayer_params(gradings, W_armor, use_en13383, requested_custom_family) result(res)
        type(GradingDef), intent(in) :: gradings(:)
        real(dp), intent(in) :: W_armor
        logical, intent(in) :: use_en13383
        character(len=*), intent(in) :: requested_custom_family
        type(UnderlayerResult) :: res
        integer :: i
        logical :: found
        type(GradingDef) :: selected
        real(dp) :: target_weight, target_mass_kg, min_range_width, current_range
        character(len=16) :: family
        character(len=128) :: ratio_note
        real(dp) :: ratio_nul_nll, nll_kg, nul_kg, reference_weight_kn

        target_weight = W_armor / 10.0_dp
        target_mass_kg = (target_weight * 1000.0_dp) / g

        res%target_W = target_weight
        res%target_M50_kg = target_mass_kg
        res%used_custom_interpolation = .false.
        res%custom_family = ''
        res%custom_ratio_nul_nll = 0.0_dp
        res%custom_ratio_note = ''

        if (use_en13383) then
            found = .false.
            min_range_width = huge(1.0_dp)
            do i = 1, size(gradings)
                if (target_mass_kg > gradings(i)%NLL_kg .and. target_mass_kg < gradings(i)%NUL_kg) then
                    current_range = gradings(i)%NUL_kg - gradings(i)%NLL_kg
                    if (current_range < min_range_width) then
                        min_range_width = current_range
                        selected = gradings(i)
                        found = .true.
                    end if
                end if
            end do
            if (found) then
                res%grading_name = selected%name
                res%M50_kg = 0.5_dp * (selected%NLL_kg + selected%NUL_kg)
                res%NLL_kg = selected%NLL_kg
                res%NUL_kg = selected%NUL_kg
            end if
        end if

        if ((.not. use_en13383) .or. (res%NUL_kg <= res%NLL_kg)) then
            family = normalize_family(requested_custom_family)
            if (trim(family) /= 'HMA' .and. trim(family) /= 'LMA' .and. trim(family) /= 'CP') then
                family = select_custom_family(gradings, target_mass_kg)
            end if

            ratio_nul_nll = interpolate_family_ratio(gradings, target_mass_kg, trim(family), ratio_note)
            ratio_nul_nll = max(ratio_nul_nll, 1.01_dp)

            nll_kg = (2.0_dp * target_mass_kg) / (1.0_dp + ratio_nul_nll)
            nul_kg = ratio_nul_nll * nll_kg

            res%grading_name = 'Custom Grading'
            res%M50_kg = target_mass_kg
            res%NLL_kg = nll_kg
            res%NUL_kg = nul_kg
            res%used_custom_interpolation = .true.
            res%custom_family = family
            res%custom_ratio_nul_nll = ratio_nul_nll
            res%custom_ratio_note = ratio_note
        end if

        res%ELL_kg = 0.7_dp * res%NLL_kg
        res%EUL_kg = 1.5_dp * res%NUL_kg
        res%W_mean_kn = (res%M50_kg * g) / 1000.0_dp

        if (res%used_custom_interpolation) then
            reference_weight_kn = target_weight
        else
            reference_weight_kn = res%W_mean_kn
        end if
        res%Dn_rock = (reference_weight_kn / W_rock_spec_val)**(1.0_dp / 3.0_dp)
        res%r2 = 2.0_dp * res%Dn_rock
        res%f2 = 100.0_dp * 2.0_dp * 1.0_dp * (1.0_dp - P_rock) / (res%Dn_rock**2)
        res%W_rock_spec = W_rock_spec_val
    end function calculate_underlayer_params

    subroutine solve(formulas, gradings, formula_id, params, results)
        type(FormulaParams), intent(in) :: formulas(4)
        type(GradingDef), intent(in) :: gradings(:)
        integer, intent(in) :: formula_id
        type(Inputs), intent(in) :: params
        type(FullResults), intent(out) :: results
        type(FormulaParams) :: coeffs
        real(dp) :: k1, k2, k3, k4, k5
        real(dp) :: L0, k0, s0m, sm, Nz, storm_duration_hr, delta_trunk
        real(dp) :: term_damage, term_waves, damage_wave_ratio, scaled_term, inv_f, steepness_factor
        real(dp) :: Ns_trunk, Dn, W_trunk, packing_density_trunk
        real(dp) :: slope, kd_trunk_equiv
        type(UnderlayerResult) :: ul_trunk, ul_head
        real(dp) :: kd_ratio, kd_head_derived, delta_head, Wc_head, W_head, Ns_head, packing_density_head
        real(dp) :: r1, vol_trunk, h_trunk, a_trunk, b_trunk, vol_head, h_head, a_head, b_head

        if (formula_id < 1 .or. formula_id > 4) stop 'Invalid Formula ID'

        coeffs = formulas(formula_id)
        k1 = coeffs%k1
        k2 = coeffs%k2
        k3 = coeffs%k3
        k4 = coeffs%k4
        k5 = coeffs%k5

        L0 = calculate_L0(params%Tm)
        k0 = (2.0_dp * PI) / L0
        s0m = params%Hs / L0
        sm = s0m
        Nz = params%Number_of_Waves
        storm_duration_hr = (Nz * params%Tm) / 3600.0_dp
        delta_trunk = (params%Wc / params%Ww) - 1.0_dp

        term_damage = params%Nod**k2
        term_waves = Nz**k3
        damage_wave_ratio = term_damage / term_waves
        scaled_term = k1 * damage_wave_ratio
        inv_f = scaled_term + k4
        steepness_factor = s0m**(-k5)
        Ns_trunk = inv_f * steepness_factor

        Dn = params%Hs / (delta_trunk * Ns_trunk)
        W_trunk = params%Wc * Dn**3
        packing_density_trunk = 100.0_dp * 2.0_dp * 1.1_dp * (1.0_dp - P_cubes) / Dn**2

        slope = coeffs%slope_ratio
        kd_trunk_equiv = (params%Wc * params%Hs**3) / (W_trunk * delta_trunk**3 * slope)

        ul_trunk = calculate_underlayer_params(gradings, W_trunk, params%use_en13383, params%custom_family)

        kd_ratio = KD_RATIO_FIXED
        kd_head_derived = kd_trunk_equiv / kd_ratio
        delta_head = delta_trunk * kd_ratio**(1.0_dp / 3.0_dp)
        Wc_head = params%Ww * (delta_head + 1.0_dp)
        W_head = W_trunk * (Wc_head / params%Wc)
        Ns_head = params%Hs / (delta_head * Dn)
        packing_density_head = 100.0_dp * 2.0_dp * 1.1_dp * (1.0_dp - P_cubes) / Dn**2

        ul_head = calculate_underlayer_params(gradings, W_head, params%use_en13383, params%custom_family)

        r1 = 2.0_dp * 1.1_dp * Dn

        vol_trunk = W_trunk / params%Wc
        h_trunk = (vol_trunk / 1.0247_dp)**(1.0_dp / 3.0_dp)
        a_trunk = 1.086_dp * h_trunk
        b_trunk = 1.005_dp * h_trunk

        vol_head = W_head / Wc_head
        h_head = (vol_head / 1.0247_dp)**(1.0_dp / 3.0_dp)
        a_head = 1.086_dp * h_head
        b_head = 1.005_dp * h_head

        results%inputs = params
        results%coefficients = coeffs
        results%intermediate%L0 = L0
        results%intermediate%k0 = k0
        results%intermediate%s0m = s0m
        results%intermediate%sm = sm
        results%intermediate%Nz = Nz
        results%intermediate%Storm_Duration_hr = storm_duration_hr
        results%intermediate%delta = delta_trunk
        results%intermediate%Ns_trunk = Ns_trunk

        results%final_trunk%Ns = Ns_trunk
        results%final_trunk%Dn = Dn
        results%final_trunk%W = W_trunk
        results%final_trunk%Mass_tonnes = W_trunk / g
        results%final_trunk%Kd = kd_trunk_equiv
        results%final_trunk%r1 = r1
        results%final_trunk%packing_density = packing_density_trunk
        results%final_trunk%dims = Dimensions(h_trunk, a_trunk, b_trunk)

        results%underlayer_trunk = ul_trunk

        results%final_head%Ns = Ns_head
        results%final_head%Dn = Dn
        results%final_head%Kd = kd_head_derived
        results%final_head%Kd_Ratio = kd_ratio
        results%final_head%Delta_Required = delta_head
        results%final_head%Wc_Required = Wc_head
        results%final_head%W = W_head
        results%final_head%Mass_tonnes = W_head / g
        results%final_head%packing_density = packing_density_head
        results%final_head%dims = Dimensions(h_head, a_head, b_head)

        results%underlayer_head = ul_head
        results%P_rock_val = P_rock
        results%P_cubes_val = P_cubes
    end subroutine solve

    subroutine write_line_pair(u, line)
        integer, intent(in) :: u
        character(len=*), intent(in) :: line
        write(u,'(A)') trim(line)
        write(*,'(A)') trim(line)
    end subroutine write_line_pair

    function fmt_real(val, prec) result(out)
        real(dp), intent(in) :: val
        integer, intent(in) :: prec
        character(len=64) :: out
        character(len=32) :: f
        integer :: width
        if (prec == 0) then
            write(out,'(I0)') nint(val)
        else
            width = max(prec + 8, 8)
            write(f,'("(F",I0,".",I0,")")') width, prec
            write(out, f) val
            out = adjustl(out)
        end if
    end function fmt_real

    function fmt_int(val) result(out)
        integer, intent(in) :: val
        character(len=64) :: out
        write(out,'(I0)') val
        out = adjustl(out)
    end function fmt_int

    function make_field(label, value, width) result(line)
        character(len=*), intent(in) :: label, value
        integer, intent(in) :: width
        character(len=256) :: line
        line = '   ' // pad_right(trim(label), width) // ' : ' // trim(value)
    end function make_field

    pure function pad_right(s, width) result(out)
        character(len=*), intent(in) :: s
        integer, intent(in) :: width
        character(len=:), allocatable :: out
        integer :: ls
        ls = len_trim(s)
        if (ls >= width) then
            out = s(1:ls)
        else
            out = s(1:ls) // repeat(' ', width - ls)
        end if
    end function pad_right


    subroutine generate_report_file(results, filepath)
        type(FullResults), intent(in) :: results
        character(len=*), intent(in) :: filepath
        integer :: u, iostat_val
        real(dp) :: armor_vol_trunk, armor_vol_head
        character(len=100) :: methodology_name

        if (trim(results%coefficients%type_name) == 'Antifer') then
            armor_vol_trunk = 1.0247_dp * results%final_trunk%dims%H**3
            armor_vol_head = 1.0247_dp * results%final_head%dims%H**3
        else
            armor_vol_trunk = results%final_trunk%Dn**3
            armor_vol_head = results%final_head%Dn**3
        end if

        methodology_name = get_formula_display_name(results%inputs%Formula_ID)

        open(newunit=u, file=filepath, status='unknown', action='write', &
             position='append', iostat=iostat_val)
        if (iostat_val /= 0) stop 'Failed to open output.txt'

        call write_line_pair(u, '================================================================================')
        call write_line_pair(u, '    TECHNICAL REPORT: BREAKWATER ARMOR & UNDERLAYER DESIGN')
        call write_line_pair(u, '================================================================================')
        call write_line_pair(u, 'Methodology: ' // trim(methodology_name))
        call write_line_pair(u, '--------------------------------------------------------------------------------')

        call write_line_pair(u, '1. INPUT PARAMETERS')
        call write_line_pair(u, make_field('Hs (Significant Wave Height)', &
            trim(fmt_real(results%inputs%Hs,2))//' m', 38))
        call write_line_pair(u, make_field('Tm (Mean Wave Period)', &
            trim(fmt_real(results%inputs%Tm,2))//' s', 38))
        call write_line_pair(u, make_field('Number of waves (Nz)', &
            trim(fmt_int(nint(results%inputs%Number_of_Waves))), 38))
        call write_line_pair(u, make_field('Nod (Damage)', &
            trim(fmt_real(results%inputs%Nod,2)), 38))
        call write_line_pair(u, make_field('Wc Trunk (Concrete Spec. Weight)', &
            trim(fmt_real(results%inputs%Wc,2))//' kN/m3', 38))
        call write_line_pair(u, make_field('Ww (Water Specific Weight)', &
            trim(fmt_real(results%inputs%Ww,2))//' kN/m3', 38))
        call write_line_pair(u, make_field('Relative Density D=(Wc/Ww)-1', &
            trim(fmt_real(results%intermediate%delta,4)), 38))
        call write_line_pair(u, make_field('Structure Slope (TRUNK & HEAD)', &
            trim(fmt_real(results%coefficients%slope_ratio,1))//':1', 38))
        call write_line_pair(u, make_field('Porosity (Cubes)', &
            trim(fmt_real(results%P_cubes_val*100.0_dp,0))//'%', 38))
        call write_line_pair(u, make_field('Porosity (Rock Layer)', &
            trim(fmt_real(results%P_rock_val*100.0_dp,0))//'%', 38))
        call write_line_pair(u, '--------------------------------------------------------------------------------')

        call write_line_pair(u, '2. INTERMEDIATE PARAMETERS')
        call write_line_pair(u, make_field('Wave Length (L0)', &
            trim(fmt_real(results%intermediate%L0,2))//' m', 38))
        call write_line_pair(u, make_field('Wave number (k0 = 2*pi/L0)', &
            trim(fmt_real(results%intermediate%k0,4)), 38))
        call write_line_pair(u, make_field('Wave steepness (s0m = Hs/L0)', &
            trim(fmt_real(results%intermediate%s0m,4)), 38))
        call write_line_pair(u, make_field('Storm Duration (h)', &
            trim(fmt_real(results%intermediate%Storm_Duration_hr,3))//' h', 38))
        call write_line_pair(u, make_field('Stability Number TRUNK (Ns)', &
            trim(fmt_real(results%intermediate%Ns_trunk,4)), 38))
        call write_line_pair(u, make_field('Stability Number HEAD (Ns)', &
            trim(fmt_real(results%final_head%Ns,4)), 38))
        call write_line_pair(u, '--------------------------------------------------------------------------------')

        call write_line_pair(u, '3. ARMOR LAYER RESULTS - TRUNK')
        call write_line_pair(u, make_field('BLOCK WEIGHT (W)', &
            trim(fmt_real(results%final_trunk%W,2))//' kN', 38))
        call write_line_pair(u, make_field('Mass (t)', &
            trim(fmt_real(results%final_trunk%Mass_tonnes,2))//' t', 38))
        call write_line_pair(u, make_field('Nominal Dimension (Dn)', &
            trim(fmt_real(results%final_trunk%Dn,3))//' m', 38))
        if (trim(results%coefficients%type_name) == 'Cubes') then
            call write_line_pair(u, make_field('Volume = Dn^3 (V)', &
                trim(fmt_real(armor_vol_trunk,3))//' m3', 38))
        else
            call write_line_pair(u, make_field('Volume = 1.0247 * H^3 (V)', &
                trim(fmt_real(armor_vol_trunk,3))//' m3', 38))
            call write_line_pair(u, make_field('Block Height (H)', &
                trim(fmt_real(results%final_trunk%dims%H,3))//' m', 38))
            call write_line_pair(u, make_field('Block Top Width (B)', &
                trim(fmt_real(results%final_trunk%dims%B,3))//' m', 38))
            call write_line_pair(u, make_field('Block Base Width (A)', &
                trim(fmt_real(results%final_trunk%dims%A,3))//' m', 38))
        end if
        call write_line_pair(u, make_field('KD_TRUNK (Equivalent)', &
            trim(fmt_real(results%final_trunk%Kd,2)), 38))
        call write_line_pair(u, make_field('Double Layer Thickness (r1)', &
            trim(fmt_real(results%final_trunk%r1,2))//' m', 38))
        call write_line_pair(u, make_field('Packing Density, d [units/100m2]', &
            trim(fmt_real(results%final_trunk%packing_density,2)), 38))
        call write_line_pair(u, '')

        call write_line_pair(u, '4. UNDERLAYER RESULTS - TRUNK')
        call write_line_pair(u, make_field('Theoretical Target (W/10)', &
            trim(fmt_real(results%underlayer_trunk%target_W,2))//' kN ('// &
            trim(fmt_real(results%underlayer_trunk%target_M50_kg,1))//' kg)', 38))
        call write_line_pair(u, make_field('Adopted rock grading', &
            trim(results%underlayer_trunk%grading_name), 38))
        call write_line_pair(u, make_field('Representative M50', &
            trim(fmt_real(results%underlayer_trunk%M50_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Nominal lower limit (NLL)', &
            trim(fmt_real(results%underlayer_trunk%NLL_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Nominal upper limit (NUL)', &
            trim(fmt_real(results%underlayer_trunk%NUL_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Extreme lower limit (ELL)', &
            trim(fmt_real(results%underlayer_trunk%ELL_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Extreme upper limit (EUL)', &
            trim(fmt_real(results%underlayer_trunk%EUL_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Nominal Dimension (Dn_rock)', &
            trim(fmt_real(results%underlayer_trunk%Dn_rock,3))//' m', 38))
        call write_line_pair(u, make_field('Double Layer Thickness (r2)', &
            trim(fmt_real(results%underlayer_trunk%r2,2))//' m', 38))
        call write_line_pair(u, make_field('Packing Density, f2 [rocks/100m2]', &
            trim(fmt_real(results%underlayer_trunk%f2,2)), 38))
        if (results%underlayer_trunk%used_custom_interpolation) then
            call write_line_pair(u, make_field('Custom family basis', &
                trim(results%underlayer_trunk%custom_family), 38))
            call write_line_pair(u, make_field('Custom ratio R=NUL/NLL', &
                trim(fmt_real(results%underlayer_trunk%custom_ratio_nul_nll,3)), 38))
        end if
        call write_line_pair(u, repeat('-',80))

        call write_line_pair(u, '5. ARMOR LAYER RESULTS - HEAD (High Density)')
        call write_line_pair(u, '   *Maintains same Dn and slope as Trunk*')
        call write_line_pair(u, make_field('Stability Ratio (Kd_T/Kd_H)', &
            trim(fmt_real(results%final_head%Kd_Ratio,2)), 38))
        call write_line_pair(u, make_field('Nominal Dimension (Dn)', &
            trim(fmt_real(results%final_head%Dn,3))//' m', 38))
        if (trim(results%coefficients%type_name) == 'Cubes') then
            call write_line_pair(u, make_field('Volume = Dn^3 (V)', &
                trim(fmt_real(armor_vol_head,3))//' m3', 38))
        else
            call write_line_pair(u, make_field('Volume = 1.0247 * H^3 (V)', &
                trim(fmt_real(armor_vol_head,3))//' m3', 38))
            call write_line_pair(u, make_field('Block Height (H)', &
                trim(fmt_real(results%final_head%dims%H,3))//' m', 38))
            call write_line_pair(u, make_field('Block Top Width (B)', &
                trim(fmt_real(results%final_head%dims%B,3))//' m', 38))
            call write_line_pair(u, make_field('Block Base Width (A)', &
                trim(fmt_real(results%final_head%dims%A,3))//' m', 38))
        end if
        call write_line_pair(u, make_field('KD_HEAD (Equivalent)', &
            trim(fmt_real(results%final_head%Kd,2)), 38))
        call write_line_pair(u, make_field('Required Concrete Density (Wc)', &
            trim(fmt_real(results%final_head%Wc_Required,2))//' kN/m3', 38))
        call write_line_pair(u, make_field('BLOCK WEIGHT (W)', &
            trim(fmt_real(results%final_head%W,2))//' kN', 38))
        call write_line_pair(u, make_field('Mass (t)', &
            trim(fmt_real(results%final_head%Mass_tonnes,2))//' t', 38))
        call write_line_pair(u, make_field('Packing Density, d [units/100m2]', &
            trim(fmt_real(results%final_head%packing_density,2)), 38))
        call write_line_pair(u, '')

        call write_line_pair(u, '6. UNDERLAYER RESULTS - HEAD')
        call write_line_pair(u, make_field('Theoretical Target (W/10)', &
            trim(fmt_real(results%underlayer_head%target_W,2))//' kN ('// &
            trim(fmt_real(results%underlayer_head%target_M50_kg,1))//' kg)', 38))
        call write_line_pair(u, make_field('Adopted rock grading', &
            trim(results%underlayer_head%grading_name), 38))
        call write_line_pair(u, make_field('Representative M50', &
            trim(fmt_real(results%underlayer_head%M50_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Nominal lower limit (NLL)', &
            trim(fmt_real(results%underlayer_head%NLL_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Nominal upper limit (NUL)', &
            trim(fmt_real(results%underlayer_head%NUL_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Extreme lower limit (ELL)', &
            trim(fmt_real(results%underlayer_head%ELL_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Extreme upper limit (EUL)', &
            trim(fmt_real(results%underlayer_head%EUL_kg,1))//' kg', 38))
        call write_line_pair(u, make_field('Nominal Dimension (Dn_rock)', &
            trim(fmt_real(results%underlayer_head%Dn_rock,3))//' m', 38))
        call write_line_pair(u, make_field('Double Layer Thickness (r2)', &
            trim(fmt_real(results%underlayer_head%r2,2))//' m', 38))
        call write_line_pair(u, make_field('Packing Density, f2 [rocks/100m2]', &
            trim(fmt_real(results%underlayer_head%f2,2)), 38))
        if (results%underlayer_head%used_custom_interpolation) then
            call write_line_pair(u, make_field('Custom family basis', &
                trim(results%underlayer_head%custom_family), 38))
            call write_line_pair(u, make_field('Custom ratio R=NUL/NLL', &
                trim(fmt_real(results%underlayer_head%custom_ratio_nul_nll,3)), 38))
        end if
        call write_line_pair(u, repeat('=',80))

        close(u)
    end subroutine generate_report_file

end program BreakwaterCalculator
