# Cirrus-cloud-studies
Model of cirrus cloud dynamics 
Here is the full fixed and complete code with the long WRITE format repaired by moving it into a named format string parameter.

program cirrus_parcel_two_bin
  implicit none

  !=============================================================
  ! DATA DICTIONARY
  !
  ! PURPOSE
  !   Compact 0-D cirrus parcel demo with:
  !   - sinusoidal vertical velocity forcing
  !   - diagnosed RH wrt ice
  !   - continuous nucleation into small-crystal bin
  !   - two ice bins:
  !       bin 1 = small crystals
  !       bin 2 = large crystals
  !   - per-bin diffusional/depositional growth
  !   - transfer from bin 1 to bin 2
  !   - latent heat feedback
  !   - separate sedimentation for each bin
  !   - CSV output
  !   - summary statistics
  !   - ASCII RH plot
  !
  ! IMPORTANT
  !   This is a teaching/demo parcel model, not a full research
  !   microphysics model.
  !=============================================================

  integer, parameter :: outu = 20
  integer, parameter :: max_store = 5000

  character(len=*), parameter :: header_line = &
       'time_s,temp_k,w_ms,rh_i_percent,s_i_percent,qv_kgkg,' // &
       'q1_kgkg,n1_perkg,r1_um,q2_kgkg,n2_perkg,r2_um,' // &
       'j_nuc_perkgps,dq1dt_dep,dq2dt_dep,dtemp_latent_Kps'

  character(len=*), parameter :: csv_fmt = &
       '(F10.1,A,F10.3,A,F9.4,A,F10.3,A,F10.3,A,' // &
       'ES12.4,A,ES12.4,A,ES12.4,A,F10.3,A,' // &
       'ES12.4,A,ES12.4,A,F10.3,A,ES12.4,A,' // &
       'ES12.4,A,ES12.4,A,ES12.4)'

  integer :: i, nsteps
  real(8) :: dt, tmax, t, pi

  ! Parcel state
  real(8) :: p, temp_k, qv
  real(8) :: rh_i, s_i
  real(8) :: w, w_amp, period, gamma_m

  ! Two ice bins
  real(8) :: n1, q1, r1
  real(8) :: n2, q2, r2
  real(8) :: qi, ni_ice

  ! Nucleation
  real(8) :: rh_nuc, c_nuc, alpha_nuc, m_seed
  real(8) :: j_nuc, dn1dt_nuc, dq1dt_nuc

  ! Growth / deposition
  real(8) :: alpha_dep, dep_base
  real(8) :: dm1_dt, dm2_dt
  real(8) :: dq1dt_dep, dq2dt_dep
  real(8) :: gcoef1, gcoef2

  ! Transfer between bins
  real(8) :: r_transfer, tau_transfer
  real(8) :: dn12dt, dq12dt

  ! Sedimentation
  real(8) :: sed_tau_q1, sed_tau_q2
  real(8) :: sed_tau_n1, sed_tau_n2

  ! Thermodynamics
  real(8) :: es_i, qsi
  real(8) :: dqvdt
  real(8) :: dtemp_dyn, dtemp_latent, dtempdt
  real(8) :: l_sub, cp_air

  ! Tendency containers
  real(8) :: dq1dt, dq2dt, dn1dt, dn2dt

  ! Physical bounds
  real(8) :: rho_ice
  real(8) :: q_min, n_min
  real(8) :: r1_min, r1_max, r2_min, r2_max

  ! Diagnostics / stats
  real(8) :: rh_max, rh_sum, rh_mean
  real(8) :: time_supersat, cloud_time
  real(8) :: max_q1, max_q2, max_qi
  real(8) :: max_n1, max_n2, max_ni
  real(8) :: max_r1, max_r2
  real(8) :: first_nuc_time
  real(8) :: max_j_nuc, total_nucleated
  real(8) :: max_dtemp_latent, latent_temp_integral
  real(8) :: max_dep_total
  integer :: supersat_episodes
  logical :: was_supersat, is_supersat
  logical :: first_nuc_recorded

  ! Plot storage
  real(8), dimension(max_store) :: rh_store
  integer :: store_count, plot_stride

  pi = 4.0d0 * atan(1.0d0)

  !----------------------------
  ! Numerical controls
  !----------------------------
  dt     = 1.0d0
  tmax   = 4.0d0 * 3600.0d0
  nsteps = int(tmax / dt)

  !----------------------------
  ! Initial state
  !----------------------------
  p      = 25000.0d0
  temp_k = 220.0d0
  qv     = 4.5d-5

  n1     = 0.0d0
  q1     = 0.0d0
  r1     = 5.0d-6

  n2     = 0.0d0
  q2     = 0.0d0
  r2     = 30.0d-6

  !----------------------------
  ! Forcing
  !----------------------------
  w_amp   = 0.25d0
  period  = 900.0d0
  gamma_m = 9.8d-3

  !----------------------------
  ! Nucleation parameters
  !----------------------------
  rh_nuc    = 1.35d0
  c_nuc     = 5.0d5
  alpha_nuc = 2.0d0
  m_seed    = 1.0d-18

  !----------------------------
  ! Deposition parameters
  !----------------------------
  alpha_dep = 0.30d0
  dep_base  = 1.0d-10

  !----------------------------
  ! Bin transfer parameters
  !----------------------------
  r_transfer   = 20.0d-6
  tau_transfer = 600.0d0

  !----------------------------
  ! Sedimentation parameters
  !----------------------------
  sed_tau_q1 = 12000.0d0
  sed_tau_q2 = 3600.0d0
  sed_tau_n1 = 14400.0d0
  sed_tau_n2 = 5400.0d0

  !----------------------------
  ! Latent heat constants
  !----------------------------
  l_sub  = 2.834d6
  cp_air = 1004.0d0

  !----------------------------
  ! Physical bounds
  !----------------------------
  rho_ice = 917.0d0
  q_min   = 1.0d-20
  n_min   = 1.0d-12
  r1_min  = 1.0d-6
  r1_max  = 50.0d-6
  r2_min  = 10.0d-6
  r2_max  = 300.0d-6

  !----------------------------
  ! Statistics initialization
  !----------------------------
  rh_max               = -1.0d30
  rh_sum               = 0.0d0
  rh_mean              = 0.0d0
  time_supersat        = 0.0d0
  cloud_time           = 0.0d0
  max_q1               = 0.0d0
  max_q2               = 0.0d0
  max_qi               = 0.0d0
  max_n1               = 0.0d0
  max_n2               = 0.0d0
  max_ni               = 0.0d0
  max_r1               = 0.0d0
  max_r2               = 0.0d0
  first_nuc_time       = -1.0d0
  max_j_nuc            = 0.0d0
  total_nucleated      = 0.0d0
  max_dtemp_latent     = 0.0d0
  latent_temp_integral = 0.0d0
  max_dep_total        = 0.0d0
  supersat_episodes    = 0
  was_supersat         = .false.
  first_nuc_recorded   = .false.

  !----------------------------
  ! Plot controls
  !----------------------------
  plot_stride = max(1, nsteps / max_store)
  store_count = 0

  t = 0.0d0

  open(unit=outu, file='cirrus_two_bin_output.csv', status='replace', &
       action='write', form='formatted')

  write(outu,'(A)') header_line

  do i = 0, nsteps

     !-------------------------------------------
     ! 1) Vertical velocity forcing
     !-------------------------------------------
     w = w_amp * sin(2.0d0 * pi * t / period)

     !-------------------------------------------
     ! 2) Diagnose RH wrt ice from current state
     !-------------------------------------------
     es_i = sat_vapor_ice(temp_k)
     qsi  = eps_ratio() * es_i / max(1.0d-6, p - es_i)
     rh_i = qv / max(qsi, 1.0d-20)
     s_i  = rh_i - 1.0d0

     !-------------------------------------------
     ! 3) Continuous nucleation into small bin
     !-------------------------------------------
     if (rh_i > rh_nuc) then
        j_nuc = c_nuc * (rh_i - rh_nuc) ** alpha_nuc
     else
        j_nuc = 0.0d0
     end if

     dn1dt_nuc = j_nuc
     dq1dt_nuc = j_nuc * m_seed

     if (j_nuc > 0.0d0 .and. .not. first_nuc_recorded) then
        first_nuc_time = t
        first_nuc_recorded = .true.
     end if

     !-------------------------------------------
     ! 4) Diagnose bin radii
     !-------------------------------------------
     if (n1 > n_min .and. q1 > q_min) then
        r1 = ((3.0d0 * q1) / (4.0d0 * pi * rho_ice * n1)) ** (1.0d0 / 3.0d0)
        r1 = min(max(r1, r1_min), r1_max)
     else if (n1 > n_min .and. q1 <= q_min) then
        r1 = r1_min
     else
        r1 = 5.0d-6
     end if

     if (n2 > n_min .and. q2 > q_min) then
        r2 = ((3.0d0 * q2) / (4.0d0 * pi * rho_ice * n2)) ** (1.0d0 / 3.0d0)
        r2 = min(max(r2, r2_min), r2_max)
     else if (n2 > n_min .and. q2 <= q_min) then
        r2 = r2_min
     else
        r2 = 30.0d-6
     end if

     !-------------------------------------------
     ! 5) Per-bin diffusional growth
     !-------------------------------------------
     if (n1 > n_min .and. r1 > 0.0d0) then
        gcoef1 = growth_coeff(temp_k, p, alpha_dep, dep_base)
        dm1_dt = gcoef1 * r1 * s_i
     else
        dm1_dt = 0.0d0
     end if

     if (n2 > n_min .and. r2 > 0.0d0) then
        gcoef2 = growth_coeff(temp_k, p, alpha_dep, dep_base)
        dm2_dt = gcoef2 * r2 * s_i
     else
        dm2_dt = 0.0d0
     end if

     dq1dt_dep = n1 * dm1_dt
     dq2dt_dep = n2 * dm2_dt

     if (q1 + dt * dq1dt_dep < 0.0d0) dq1dt_dep = -q1 / dt
     if (q2 + dt * dq2dt_dep < 0.0d0) dq2dt_dep = -q2 / dt

     !-------------------------------------------
     ! 6) Transfer from small bin to large bin
     !-------------------------------------------
     if (n1 > n_min .and. q1 > q_min .and. r1 >= r_transfer) then
        dn12dt = n1 / tau_transfer
        dq12dt = q1 / tau_transfer
     else
        dn12dt = 0.0d0
        dq12dt = 0.0d0
     end if

     if (n1 + dt * (-dn12dt) < 0.0d0) dn12dt = n1 / dt
     if (q1 + dt * (-dq12dt) < 0.0d0) dq12dt = q1 / dt

     !-------------------------------------------
     ! 7) Bin tendencies
     !-------------------------------------------
     dq1dt = dq1dt_dep + dq1dt_nuc - dq12dt
     dq2dt = dq2dt_dep + dq12dt

     if (q1 > q_min) dq1dt = dq1dt - q1 / sed_tau_q1
     if (q2 > q_min) dq2dt = dq2dt - q2 / sed_tau_q2

     dn1dt = dn1dt_nuc - dn12dt
     dn2dt = dn12dt

     if (n1 > n_min) dn1dt = dn1dt - n1 / sed_tau_n1
     if (n2 > n_min) dn2dt = dn2dt - n2 / sed_tau_n2

     ! Vapor loss only from deposition
     dqvdt = -(dq1dt_dep + dq2dt_dep)

     !-------------------------------------------
     ! 8) Temperature tendencies
     !-------------------------------------------
     dtemp_dyn    = -gamma_m * w
     dtemp_latent = (l_sub / cp_air) * (dq1dt_dep + dq2dt_dep)
     dtempdt      = dtemp_dyn + dtemp_latent

     !-------------------------------------------
     ! 9) State update
     !-------------------------------------------
     temp_k = temp_k + dt * dtempdt
     qv     = qv + dt * dqvdt

     q1 = q1 + dt * dq1dt
     q2 = q2 + dt * dq2dt
     n1 = n1 + dt * dn1dt
     n2 = n2 + dt * dn2dt

     !-------------------------------------------
     ! 10) Floors / cleanup
     !-------------------------------------------
     if (qv < 1.0d-12) qv = 1.0d-12

     if (q1 < 0.0d0) q1 = 0.0d0
     if (q2 < 0.0d0) q2 = 0.0d0
     if (n1 < 0.0d0) n1 = 0.0d0
     if (n2 < 0.0d0) n2 = 0.0d0

     if (q1 <= q_min) q1 = 0.0d0
     if (q2 <= q_min) q2 = 0.0d0
     if (n1 <= n_min) n1 = 0.0d0
     if (n2 <= n_min) n2 = 0.0d0

     !-------------------------------------------
     ! 11) Post-update diagnostics
     !-------------------------------------------
     qi     = q1 + q2
     ni_ice = n1 + n2

     es_i = sat_vapor_ice(temp_k)
     qsi  = eps_ratio() * es_i / max(1.0d-6, p - es_i)
     rh_i = qv / max(qsi, 1.0d-20)
     s_i  = rh_i - 1.0d0

     if (n1 > n_min .and. q1 > q_min) then
        r1 = ((3.0d0 * q1) / (4.0d0 * pi * rho_ice * n1)) ** (1.0d0 / 3.0d0)
        r1 = min(max(r1, r1_min), r1_max)
     else
        r1 = 0.0d0
     end if

     if (n2 > n_min .and. q2 > q_min) then
        r2 = ((3.0d0 * q2) / (4.0d0 * pi * rho_ice * n2)) ** (1.0d0 / 3.0d0)
        r2 = min(max(r2, r2_min), r2_max)
     else
        r2 = 0.0d0
     end if

     !-------------------------------------------
     ! 12) CSV output
     !-------------------------------------------
     write(outu, csv_fmt) &
          t, ',', temp_k, ',', w, ',', 100.0d0 * rh_i, ',', 100.0d0 * s_i, ',', &
          qv, ',', q1, ',', n1, ',', 1.0d6 * r1, ',', &
          q2, ',', n2, ',', 1.0d6 * r2, ',', &
          j_nuc, ',', dq1dt_dep, ',', dq2dt_dep, ',', dtemp_latent

     !-------------------------------------------
     ! 13) Statistics
     !-------------------------------------------
     rh_max = max(rh_max, rh_i)
     rh_sum = rh_sum + rh_i

     is_supersat = (rh_i > 1.0d0)
     if (is_supersat) time_supersat = time_supersat + dt
     if (qi > q_min) cloud_time = cloud_time + dt

     if (is_supersat .and. .not. was_supersat) supersat_episodes = supersat_episodes + 1
     was_supersat = is_supersat

     max_q1 = max(max_q1, q1)
     max_q2 = max(max_q2, q2)
     max_qi = max(max_qi, qi)

     max_n1 = max(max_n1, n1)
     max_n2 = max(max_n2, n2)
     max_ni = max(max_ni, ni_ice)

     max_r1 = max(max_r1, r1)
     max_r2 = max(max_r2, r2)

     max_j_nuc = max(max_j_nuc, j_nuc)
     total_nucleated = total_nucleated + dt * j_nuc
     max_dtemp_latent = max(max_dtemp_latent, abs(dtemp_latent))
     latent_temp_integral = latent_temp_integral + dt * dtemp_latent
     max_dep_total = max(max_dep_total, abs(dq1dt_dep + dq2dt_dep))

     !-------------------------------------------
     ! 14) Plot storage
     !-------------------------------------------
     if (mod(i, plot_stride) == 0) then
        if (store_count < max_store) then
           store_count = store_count + 1
           rh_store(store_count) = rh_i
        end if
     end if

     t = t + dt
  end do

  close(outu)

  rh_mean = rh_sum / dble(nsteps + 1)

  print *, '==================================================='
  print *, 'CIRRUS PARCEL SUMMARY WITH TWO ICE BINS'
  print *, '==================================================='
  print '(A,F10.4)',  'Max RH_ice:                      ', rh_max
  print '(A,F10.4)',  'Mean RH_ice:                     ', rh_mean
  print '(A,F10.1)',  'Time supersaturated (s):         ', time_supersat
  print '(A,F10.1)',  'Cloud lifetime (s):              ', cloud_time
  if (first_nuc_time >= 0.0d0) then
     print '(A,F10.1)', 'First nucleation time (s):       ', first_nuc_time
  else
     print '(A)',       'First nucleation time (s):       none'
  end if

  print '(A,ES12.4)', 'Max q1 small-bin mass:           ', max_q1
  print '(A,ES12.4)', 'Max q2 large-bin mass:           ', max_q2
  print '(A,ES12.4)', 'Max total ice mass:              ', max_qi

  print '(A,ES12.4)', 'Max n1 small-bin number:         ', max_n1
  print '(A,ES12.4)', 'Max n2 large-bin number:         ', max_n2
  print '(A,ES12.4)', 'Max total ice number:            ', max_ni

  print '(A,F10.3)',  'Max r1 small-bin radius (um):    ', 1.0d6 * max_r1
  print '(A,F10.3)',  'Max r2 large-bin radius (um):    ', 1.0d6 * max_r2

  print '(A,I10)',    'Supersaturation episodes:        ', supersat_episodes
  print '(A,ES12.4)', 'Peak nucleation rate:            ', max_j_nuc
  print '(A,ES12.4)', 'Total crystals nucleated:        ', total_nucleated
  print '(A,ES12.4)', 'Max |total deposition dq_i/dt|:  ', max_dep_total
  print '(A,ES12.4)', 'Max |latent dT/dt| (K/s):        ', max_dtemp_latent
  print '(A,F10.4)',  'Integral latent dT (K):          ', latent_temp_integral
  print *, 'CSV file: cirrus_two_bin_output.csv'
  print *, ' '

  print *, 'ASCII PLOT OF RH_ICE'
  print *, 'Marker | indicates RH_ice = 1.0'
  print *, 'Scale roughly spans RH 0.80 to 1.60'
  print *, '---------------------------------------------------'

  do i = 1, store_count
     call ascii_bar(rh_store(i))
  end do

contains

  subroutine ascii_bar(rh)
    implicit none
    real(8), intent(in) :: rh
    integer :: pos, sat_pos, width
    character(len=80) :: line

    width   = 60
    sat_pos = 15
    line    = ' '

    pos = 1 + int(((rh - 0.80d0) / 0.80d0) * dble(width - 1))
    pos = max(1, min(width, pos))

    if (sat_pos >= 1 .and. sat_pos <= width) then
       line(sat_pos:sat_pos) = '|'
    end if

    line(pos:pos) = '*'

    print '(F7.3,2X,A)', rh, line(1:width)
  end subroutine ascii_bar

  real(8) function sat_vapor_ice(tk)
    implicit none
    real(8), intent(in) :: tk

    ! Simple saturation vapor pressure over ice approximation, Pa
    sat_vapor_ice = 6.112d0 * exp(22.46d0 * (tk - 273.15d0) / (tk - 0.55d0))
    sat_vapor_ice = sat_vapor_ice * 100.0d0
  end function sat_vapor_ice

  real(8) function eps_ratio()
    implicit none
    eps_ratio = 0.622d0
  end function eps_ratio

  real(8) function growth_coeff(tk, pres_pa, alpha_dep_in, dep_base_in)
    implicit none
    real(8), intent(in) :: tk, pres_pa, alpha_dep_in, dep_base_in
    real(8) :: temp_factor, press_factor, bounded_tk

    bounded_tk = min(max(tk, 180.0d0), 273.15d0)

    temp_factor  = (bounded_tk / 220.0d0) ** 1.75d0
    press_factor = (25000.0d0 / max(pres_pa, 1000.0d0)) ** 0.30d0

    growth_coeff = dep_base_in * alpha_dep_in * temp_factor * press_factor
  end function growth_coeff

end program cirrus_parcel_two_bin

Compile with:

gfortran -O2 -Wall -Wextra -std=f2008 cirrusmodel3.f90 -o cirrusmodel3
./cirrusmodel3

This version fixes the line-truncation problem by removing the overlong inline format string.
