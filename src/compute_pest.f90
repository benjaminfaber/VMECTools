! This code is based off of the 
! vmec2sfl.f90
! Written by Matt Landreman, University of Maryland
! Initial code written August 2017.


! This code has been further modified by B.J. Faber to generalize to different
! coordinate systems.  As such, some of the variable names have bene modified.
! Skip down ~25 lines for detailed description of the input and output parameters.

module compute_pest 

  use types, only: dp, pi, mu_0
  use vmec_object, only: VMEC_Obj
  use pest_object, only: PEST_Obj
  implicit none

  public :: compute_pest_surface, compute_pest_sfl
  
  private

  ! Module wide to variables
  integer :: sign_toroidal_flux
  real(dp) :: theta_pest_target, zeta0, edge_toroidal_flux_over_2pi
  real(dp) :: ds, d_pressure_d_s, d_iota_d_s
  real(dp), dimension(2) :: vmec_radial_weight_full, vmec_radial_weight_half
  integer, dimension(2) :: vmec_radial_index_full, vmec_radial_index_half

contains

  subroutine compute_pest_surface(desired_normalized_toroidal_flux,&
    &vmec_surface_option, pest)

    !use read_wout_mod,pest%iz2_vmec => pest%nzeta ! VMEC has a variable ngzrid which conflicts with our pest%nzeta, so rename vmec's version.

    implicit none
    !*********************************************************************
    ! Input parameters
    !*********************************************************************
    ! The parameter desired_normalized_toroidal_flux determines which flux surface from the VMEC file will be used
    ! for the computation. This parameter should lie in the interval [0,1].
    real(dp), intent(in) :: desired_normalized_toroidal_flux


    ! If vmec_surface_option = 0, the magnetic surface specified by desired_normalized_toroidal_flux will be used,
    ! by intedpolating between the surfaces available in the vmec file.
    ! If vmec_surface_option = 1, the magnetic surface on vmec's HALF radial mesh will be used that is closest to desired_normalized_toroidal_flux.
    ! If vmec_surface_option = 2, the magnetic surface on vmec's FULL radial mesh will be used that is closest to desired_normalized_toroidal_flux.    
    ! Other values of vmec_surface_option will cause the program to abort with an error.
    integer, intent(in) :: vmec_surface_option

    type(PEST_Obj), intent(inout) :: pest


    !***************************************************************************
    ! Output parameters
    !***************************************************************************
    !! Pass back number of field periods to the interface routine
    !real(dp), intent(out) :: local_nfp

    integer :: j, index
    !integer :: ierr, iopen 
    real(dp) :: dphi, min_dr2
    real(dp), dimension(:), allocatable :: dr2, normalized_toroidal_flux_full_grid, normalized_toroidal_flux_half_grid
    real(dp), dimension(:), allocatable :: d_pressure_d_s_on_half_grid, d_iota_d_s_on_half_grid

    integer :: ns
    logical :: verbose, test
    !*********************************************************************
    ! VMEC variables of interest:
    ! ns = number of flux surfaces used by VMEC
    ! nfp = number of field periods, e.g. 5 for W7-X, 4 for HSX
    ! pest%iotas = rotational transform (1/q) on the half grid.
    ! pest%iotaf = rotational transform on the full grid.
    ! presf = pressure on the full grid.
    !
    ! All VMEC quantities (B, pressure, etc) are in SI units.
    ! 
    ! In VMEC, quantities on the half grid have the same number of array elements (ns) as quantities on the full grid,
    ! but the first array element is 0.
    !
    !*********************************************************************

    !*********************************************************************
    ! Beginning of executable statements.
    !*********************************************************************
    ns = pest%vmec%ns
    verbose = .false.
    test = .false.

    if (verbose) print *,"Entering subroutine compute_pest_surface."
    !*********************************************************************
    ! Do some validation.
    !*********************************************************************
    if (desired_normalized_toroidal_flux <= 0) then
       print *,"Error! desired_normalized_toroidal_flux must be >0. Instead it is",desired_normalized_toroidal_flux
       stop
    end if

    if (desired_normalized_toroidal_flux > 1) then
       print *,"Error! desired_normalized_toroidal_flux must be <= 1. Instead it is",desired_normalized_toroidal_flux
       stop
    end if

    edge_toroidal_flux_over_2pi = pest%vmec%phi(ns) / (2*pi) * pest%vmec%isigng ! isigns is called signgs in the wout*.nc file. Why is this signgs here?

    ! this gives the sign of the edge toroidal flux
    sign_toroidal_flux = int(sign(1.1,edge_toroidal_flux_over_2pi))
    write (*,*) 'sign_toroidal_flux', sign_toroidal_flux


    allocate(normalized_toroidal_flux_full_grid(ns))
    normalized_toroidal_flux_full_grid = [( real(j-1)/(ns-1), j=1,ns )]

    ! Build an array of the half grid points:
    allocate(normalized_toroidal_flux_half_grid(ns-1))
    do j = 1,ns-1
       normalized_toroidal_flux_half_grid(j) = (normalized_toroidal_flux_full_grid(j) + normalized_toroidal_flux_full_grid(j+1))*(0.5d+0)
    end do

    !*********************************************************************
    ! Determine which flux surface to use, based on 
    ! desired_normalized_toroidal_flux and vmec_surface_option.
    !*********************************************************************

    ! Possible values of vmec_surface_option:
    ! 0 = Use the exact radius requested.
    ! 1 = Use the nearest value of the VMEC half grid.
    ! 2 = Use the nearest value of the VMEC full grid.

    select case (vmec_surface_option)
    case (0)
       ! Use exact radius requested.
       pest%s0 = desired_normalized_toroidal_flux

    case (1)
       ! Use nearest value of the VMEC half grid

       ! Compute differences
       allocate(dr2(ns-1))
       dr2 = (normalized_toroidal_flux_half_grid - desired_normalized_toroidal_flux) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,ns-1
          if (dr2(j)<min_dr2) then
             index = j
             min_dr2 = dr2(j)
          end if
       end do

       pest%s0 = normalized_toroidal_flux_half_grid(index)
       deallocate(dr2)

    case (2)
       ! Use nearest value of the VMEC full grid

       ! Compute differences
       allocate(dr2(ns))
       dr2 = (normalized_toroidal_flux_full_grid - desired_normalized_toroidal_flux) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,ns
          if (dr2(j)<min_dr2) then
             index = j
             min_dr2 = dr2(j)
          end if
       end do

       pest%s0 = normalized_toroidal_flux_full_grid(index)
       deallocate(dr2)

    case default
       print *,"Error! vmec_surface_option must be 0, 1, or 2. It is instead ",vmec_surface_option
       stop
    end select

    ! --------------------------------------------------------------------------------
    ! Done choosing the actual radius to use.
    ! --------------------------------------------------------------------------------

    ! In general, we get quantities for gs2 by linear intedpolation, taking a weighted average of the quantity from
    ! 2 surfaces in the VMEC file. Sometimes the weights are 0 and 1, i.e. no intedpolation is needed.

    ! For any VMEC quantity Q on the full grid, the value used in GS2 will be
    !  Q_gs2 = Q(vmec_radial_index_full(1))*vmec_radial_weight_full(1) + Q(vmec_radial_index_full(2))*vmec_radial_weight_full(2)

    ! For any VMEC quantity Q on the half grid, the value used in GS2 will be
    !  Q_gs2 = Q(vmec_radial_index_half(1))*vmec_radial_weight_half(1) + Q(vmec_radial_index_half(2))*vmec_radial_weight_half(2)


    ! Handle quantities for the full grid
    if (pest%s0>1) then
       stop "Error! pest%s0 cannot be >1"
    elseif (pest%s0<0) then
       stop "Error! pest%s0 cannot be <0"
    elseif (pest%s0==1) then
       vmec_radial_index_full(1) = ns-1
       vmec_radial_index_full(2) = ns
       vmec_radial_weight_full(1) = 0.0d0
    else
       ! pest%s0 is >= 0 and <1
       ! This is the most common case.
       vmec_radial_index_full(1) = floor(pest%s0*(ns-1))+1
       vmec_radial_index_full(2) = vmec_radial_index_full(1) + 1
       vmec_radial_weight_full(1) = vmec_radial_index_full(1) - pest%s0*(ns-1.0d0)
    end if
    vmec_radial_weight_full(2) = 1.0d0 - vmec_radial_weight_full(1)

    ! Handle quantities for the half grid
    if (pest%s0 < normalized_toroidal_flux_half_grid(1)) then
       print *,"Warning: extrapolating beyond the end of VMEC's half grid."
       print *,"(Extrapolating towards the magnetic axis.) Results are likely to be inaccurate."

       ! We start at element 2 since element 1 is always 0 for quantities on the half grid.
       vmec_radial_index_half(1) = 2
       vmec_radial_index_half(2) = 3
       vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(2) - pest%s0) / (normalized_toroidal_flux_half_grid(2) - normalized_toroidal_flux_half_grid(1))

    elseif (pest%s0 > normalized_toroidal_flux_half_grid(ns-1)) then
       print *,"Warning: extrapolating beyond the end of VMEC's half grid."
       print *,"(Extrapolating towards the last closed flux surface.) Results may be inaccurate."
       vmec_radial_index_half(1) = ns-1
       vmec_radial_index_half(2) = ns
       vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(ns-1) - pest%s0) &
            / (normalized_toroidal_flux_half_grid(ns-1) - normalized_toroidal_flux_half_grid(ns-2))

    elseif (pest%s0 == normalized_toroidal_flux_half_grid(ns-1)) then
       ! We are exactly at the last point of the half grid
       vmec_radial_index_half(1) = ns-1
       vmec_radial_index_half(2) = ns
       vmec_radial_weight_half(1) = 0.0d0
    else
       ! pest%s0 is inside the half grid.
       ! This is the most common case.
       vmec_radial_index_half(1) = floor(pest%s0*(ns-1) + 0.5d+0)+1
       if (vmec_radial_index_half(1) < 2) then
          ! This can occur sometimes due to roundoff error.
          vmec_radial_index_half(1) = 2
       end if
       vmec_radial_index_half(2) = vmec_radial_index_half(1) + 1
       vmec_radial_weight_half(1) = vmec_radial_index_half(1) - pest%s0*(ns-1.0d0) - (0.5d+0)
    end if
    vmec_radial_weight_half(2) = 1.0d0-vmec_radial_weight_half(1)

    if (verbose) then
       if (abs(vmec_radial_weight_half(1)) < 1e-14) then
          print "(a,i4,a,i4,a)","   Using radial index ",vmec_radial_index_half(2)," of ",ns," from vmec's half mesh."
       elseif (abs(vmec_radial_weight_half(2)) < 1e-14) then
          print "(a,i4,a,i4,a)","   Using radial index ",vmec_radial_index_half(1)," of ",ns," from vmec's half mesh."
       else
          print "(a,i4,a,i4,a,i4,a)", "   Intedpolating using radial indices ",vmec_radial_index_half(1)," and ",vmec_radial_index_half(2),&
               " of ",ns," from vmec's half mesh."
          print "(a,f17.14,a,f17.14)", "   Weights for half mesh = ",vmec_radial_weight_half(1)," and ",vmec_radial_weight_half(2)
          print "(a,i4,a,i4,a,i4,a)", "   Intedpolating using radial indices ",vmec_radial_index_full(1)," and ",vmec_radial_index_full(2),&
               " of ",ns," from vmec's full mesh."
          print "(a,f17.14,a,f17.14)", "   Weights for full mesh = ",vmec_radial_weight_full(1)," and ",vmec_radial_weight_full(2)
       end if
    end if

    !*********************************************************************
    ! Evaluate several radial-profile functions at the flux surface
    ! we ended up choosing.
    !*********************************************************************

    pest%iota = pest%vmec%iotas(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + pest%vmec%iotas(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    if (verbose) print *,"  pest%iota =",pest%iota
    pest%safety_factor_q = 1/pest%iota

    allocate(d_iota_d_s_on_half_grid(ns))
    d_iota_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    if (verbose) print *,"  ds =",ds
    d_iota_d_s_on_half_grid(2:ns) = (pest%vmec%iotaf(2:ns) - pest%vmec%iotaf(1:ns-1)) / ds
    d_iota_d_s =  &
         d_iota_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + d_iota_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    deallocate(d_iota_d_s_on_half_grid)
    if (verbose) print *,"  d pest%iota / d s =",d_iota_d_s
    ! shat = (r/q)(dq/dr) where r = a sqrt(s).
    !      = - (r/pest%iota) (d pest%iota / d r) = -2 (s/pest%iota) (d pest%iota / d s)
    pest%shat = (-2 * pest%s0 / pest%iota) * d_iota_d_s

    allocate(d_pressure_d_s_on_half_grid(ns))
    d_pressure_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    d_pressure_d_s_on_half_grid(2:ns) = (pest%vmec%presf(2:ns) - pest%vmec%presf(1:ns-1)) / ds
    d_pressure_d_s =  &
         d_pressure_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + d_pressure_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    deallocate(d_pressure_d_s_on_half_grid)
    if (verbose) print *,"  d pressure / d s =",d_pressure_d_s

    deallocate(normalized_toroidal_flux_full_grid)
    deallocate(normalized_toroidal_flux_half_grid)
  end subroutine compute_pest_surface


!*******************************************************************************
! Subroutine compute_pest_sfl
!*******************************************************************************
  subroutine compute_pest_sfl(zeta_center, &
    & number_of_field_periods_to_include,pest)

    implicit none
    !***************************************************************************
    ! Input parameters
    !***************************************************************************

    ! The pest%zeta domain is centered at zeta_center. Setting zeta_center = 2*pi*N/nfp for any integer N should
    ! yield identical results to setting zeta_center = 0, where nfp is the number of field periods (as in VMEC).
    real(dp), intent(in) :: zeta_center

    ! If number_of_field_periods_to_include is > 0, then this parameter does what you think:
    ! the extent of the toroidal in pest%zeta will be 2*pi*number_of_field_periods_to_include/nfp.
    ! If number_of_field_periods_to_include is <= 0, the entire 2*pi toroidal domain will be included.
    real(dp), intent(in) :: number_of_field_periods_to_include
    type(PEST_Obj), intent(inout) :: pest


    !*********************************************************************
    ! Variables used internally by this subroutine
    !*********************************************************************

    integer :: j, iz, iztemp, ia, which_surface, isurf, m, n, imn, imn_nyq, ns
    real(dp) :: angle, sin_angle, cos_angle, temp 
    real(dp), dimension(:,:), allocatable :: theta_vmec
    !integer :: ierr, iopen, fzero_flag
    integer :: fzero_flag
    real(dp) :: number_of_field_periods_to_include_final
    real(dp) :: d_pressure_d_s, scale_factor
    real(dp) :: theta_vmec_min, theta_vmec_max, sqrt_s
    real(dp) :: root_solve_absolute_tolerance, root_solve_relative_tolerance
    logical :: non_Nyquist_mode_available, found_imn
    real(dp), dimension(:,:), allocatable :: B, sqrt_g, R, Z, B_dot_grad_theta_pest_over_B_dot_grad_zeta, temp2D
    real(dp), dimension(:,:), allocatable :: d_B_d_theta_vmec, d_B_d_zeta, d_B_d_s
    real(dp), dimension(:,:), allocatable :: d_R_d_theta_vmec, d_R_d_zeta, d_R_d_s
    real(dp), dimension(:,:), allocatable :: d_Z_d_theta_vmec, d_Z_d_zeta, d_Z_d_s
    real(dp), dimension(:,:), allocatable :: d_X_d_theta_vmec, d_X_d_zeta, d_X_d_s
    real(dp), dimension(:,:), allocatable :: d_Y_d_theta_vmec, d_Y_d_zeta, d_Y_d_s
    real(dp), dimension(:,:), allocatable :: d_Lambda_d_theta_vmec, d_Lambda_d_zeta, d_Lambda_d_s
    real(dp), dimension(:,:), allocatable :: B_sub_s, B_sub_theta_vmec, B_sub_zeta
    real(dp), dimension(:,:), allocatable :: B_sup_theta_vmec, B_sup_zeta
    real(dp), dimension(:), allocatable :: d_B_d_s_mnc, d_B_d_s_mns
    real(dp), dimension(:), allocatable :: d_R_d_s_mnc, d_R_d_s_mns
    real(dp), dimension(:), allocatable :: d_Z_d_s_mnc, d_Z_d_s_mns
    real(dp), dimension(:), allocatable :: d_Lambda_d_s_mnc, d_Lambda_d_s_mns
    real(dp), dimension(:,:), allocatable :: grad_s_X, grad_s_Y, grad_s_Z
    real(dp), dimension(:,:), allocatable :: grad_theta_vmec_X, grad_theta_vmec_Y, grad_theta_vmec_Z
    real(dp), dimension(:,:), allocatable :: grad_zeta_X, grad_zeta_Y, grad_zeta_Z
    real(dp), dimension(:,:), allocatable :: grad_psi_X, grad_psi_Y, grad_psi_Z
    real(dp), dimension(:,:), allocatable :: grad_alpha_X, grad_alpha_Y, grad_alpha_Z
    real(dp), dimension(:,:), allocatable :: B_cross_grad_B_dot_grad_alpha, B_cross_grad_B_dot_grad_alpha_alternate
    real(dp), dimension(:,:), allocatable :: B_cross_grad_s_dot_grad_alpha, B_cross_grad_s_dot_grad_alpha_alternate
    real(dp), dimension(:,:), allocatable :: grad_B_X, grad_B_Y, grad_B_Z
    real(dp), dimension(:,:), allocatable :: B_X, B_Y, B_Z
    logical :: verbose, test
    verbose = .false.
    test = .false.
    !*********************************************************************
    ! Read in everything from the vmec wout file using libstell.
    !*********************************************************************

    !if (verbose) print *,"  About to read VMEC wout file ",trim(vmec_filename)
    !call read_wout_file(vmec_filename, ierr, iopen)
    !if (iopen .ne. 0) stop 'error opening wout file'
    !if (ierr .ne. 0) stop 'error reading wout file'
    !if (verbose) print *,"  Successfully read VMEC data from ",trim(vmec_filename)
!
    if (verbose) print *,"  Number of field periods (nfp):",pest%vmec%nfp

    !*********************************************************************
    ! Do some validation.
    !*********************************************************************

    if (pest%nalpha<1) then
       print *,"Error! pest%nalpha must be >= 1. Instead it is",pest%nalpha
       stop
    end if

    if (pest%nzeta<1) then
       print *,"Error! pest%nzeta must be >= 1. Instead it is",pest%nzeta
       stop
    end if

    ns = pest%vmec%ns
    !*********************************************************************
    ! Set up the coordinate grids.
    !*********************************************************************

    pest%alpha = [( (j*2*pi) / pest%nalpha, j=pest%ia1, pest%ia2 )]

!!$    if (number_of_field_periods_to_include > nfp) then
!!$       print *,"Error! number_of_field_periods_to_include > nfp"
!!$       print *,"  number_of_field_periods_to_include =",number_of_field_periods_to_include
!!$       print *,"  nfp =",nfp
!!$       stop
!!$    end if
    number_of_field_periods_to_include_final = number_of_field_periods_to_include
    if (number_of_field_periods_to_include <= 0) then
       number_of_field_periods_to_include_final = pest%vmec%nfp
       if (verbose) print *,"  Since number_of_field_periods_to_include was <= 0, it is being reset to nfp =",pest%vmec%nfp
    end if

    pest%zeta = [( zeta_center + (pi*j*number_of_field_periods_to_include_final)/(pest%vmec%nfp*pest%nzeta), j=pest%iz1,pest%iz2 )]
print *, pest%zeta

    !*********************************************************************
    ! We know theta_pest = alpha + pest%iota * pest%zeta, but we need to determine
    ! theta_vmec = theta_pest - Lambda.
    !*********************************************************************

    allocate(theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))

    if (verbose) print *,"  Beginning root solves to determine theta_vmec."
    root_solve_absolute_tolerance = 1.0d-10
    root_solve_relative_tolerance = 1.0d-10
    do iz = pest%iz1,pest%iz2
       zeta0 = pest%zeta(iz)
       do ia = pest%ia1,pest%ia2
          theta_pest_target = pest%alpha(ia) + pest%iota * zeta0
          ! Guess that theta_vmec will be within 0.3 radians of theta_pest:
          theta_vmec_min = theta_pest_target - 0.3
          theta_vmec_max = theta_pest_target + 0.3

          ! In the 4th argument, we are telling the root-finder (fzero) to use theta_pest as the initial guess for theta_vmec.
          call fzero(fzero_residual, theta_vmec_min, theta_vmec_max, theta_pest_target, &
               root_solve_relative_tolerance, root_solve_absolute_tolerance, fzero_flag)
          ! Note: fzero returns its answer in theta_vmec_min.
          theta_vmec(ia,iz) = theta_vmec_min
          if (fzero_flag == 4) then
             stop "ERROR: fzero returned error 4: no sign change in residual"
          else if (fzero_flag > 2) then
             print *,"WARNING: fzero returned an error code:",fzero_flag
          end if
       end do
    end do
    if (verbose) then
       print *,"  Done with root solves. Here comes theta_vmec:"
       do j = pest%ia1,pest%ia2
          print *,theta_vmec(j,:)
       end do
    end if

    !*********************************************************************
    ! Initialize geometry arrays
    !*********************************************************************

!    bmag = 0
!    gradpar = 0
!    gds2 = 0
!    gds21 = 0
!    gds22 = 0
!    gbdrift = 0
!    gbdrift0 = 0
!    cvdrift = 0
!    cvdrift0 = 0
!    jac_gist_inv = 0
!    d_B_d_par = 0


    allocate(B(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(temp2D(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(sqrt_g(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(R(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(Z(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_B_d_theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_B_d_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_B_d_s(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_R_d_theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_R_d_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_R_d_s(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Z_d_theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Z_d_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Z_d_s(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Lambda_d_theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Lambda_d_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Lambda_d_s(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_sub_s(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_sub_theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_sub_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_sup_theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_sup_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))

    allocate(d_B_d_s_mnc(ns))
    allocate(d_B_d_s_mns(ns))
    allocate(d_R_d_s_mnc(ns))
    allocate(d_R_d_s_mns(ns))
    allocate(d_Z_d_s_mnc(ns))
    allocate(d_Z_d_s_mns(ns))
    allocate(d_Lambda_d_s_mnc(ns))
    allocate(d_Lambda_d_s_mns(ns))

    allocate(d_X_d_s(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_X_d_theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_X_d_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Y_d_s(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Y_d_theta_vmec(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(d_Y_d_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))

    allocate(grad_s_X(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_s_Y(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_s_Z(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_theta_vmec_X(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_theta_vmec_Y(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_theta_vmec_Z(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_zeta_X(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_zeta_Y(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_zeta_Z(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_psi_X(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_psi_Y(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_psi_Z(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_alpha_X(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_alpha_Y(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_alpha_Z(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    
    allocate(B_X(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_Y(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_Z(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_B_X(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_B_Y(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(grad_B_Z(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_cross_grad_B_dot_grad_alpha(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_cross_grad_B_dot_grad_alpha_alternate(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_cross_grad_s_dot_grad_alpha(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
    allocate(B_cross_grad_s_dot_grad_alpha_alternate(pest%ia1:pest%ia2,pest%iz1:pest%iz2))

    B = 0
    sqrt_g = 0
    R = 0
    Z = 0
    d_B_d_theta_vmec = 0
    d_B_d_zeta = 0
    d_B_d_s = 0
    d_R_d_theta_vmec = 0
    d_R_d_zeta = 0
    d_R_d_s = 0
    d_Z_d_theta_vmec = 0
    d_Z_d_zeta = 0
    d_Z_d_s = 0
    d_Lambda_d_theta_vmec = 0
    d_Lambda_d_zeta = 0
    d_Lambda_d_s = 0
    B_sub_s = 0
    B_sub_theta_vmec = 0
    B_sub_zeta = 0
    B_sup_theta_vmec = 0
    B_sup_zeta = 0

    !*********************************************************************
    ! Now that we know the grid points in theta_vmec, we can evaluate
    ! all the geometric quantities on the grid points.
    !*********************************************************************
    
    do imn_nyq = 1, pest%vmec%mnmax_nyq ! All the quantities we need except R, Z, and Lambda use the _nyq mode numbers.
       m = pest%vmec%xm_nyq(imn_nyq)
       n = pest%vmec%xn_nyq(imn_nyq)/pest%vmec%nfp

print *, m,' ',pest%vmec%mpol,' ',n,' ',pest%vmec%ntor
       if (abs(m) >= pest%vmec%mpol .or. abs(n) > pest%vmec%ntor) then
          non_Nyquist_mode_available = .false.
       else
          non_Nyquist_mode_available = .true.
          ! Find the imn in the non-Nyquist arrays that corresponds to the same m and n.
          found_imn = .false.
          do imn = 1,pest%vmec%mnmax
             if (pest%vmec%xm(imn)==m .and. pest%vmec%xn(imn)==n*pest%vmec%nfp) then
                found_imn = .true.
                exit
             end if
          end do
          if ((pest%vmec%xm(imn) .ne. m) .or. (pest%vmec%xn(imn) .ne. n*pest%vmec%nfp)) stop "Something went wrong!"
          if (.not. found_imn) stop "Error! imn could not be found matching the given imn_nyq."
       end if

       ! All quantities are multiplied by a variable scale_factor which can in principle depend on m and n.
       ! For now we just set scale_factor = 1. In the future, scale_factor could be used to lower the
       ! symmetry-breaking Fourier components, or filter out certain Fourier components in some way.
       scale_factor = 1

       ! -----------------------------------------------------
       ! First, consider just the stellarator-symmetric terms:
       ! -----------------------------------------------------

       ! Evaluate the radial derivatives we will need:

       ! B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
       ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.

       d_B_d_s_mnc(2:ns-1) = (pest%vmec%bmnc(imn_nyq,3:ns) - pest%vmec%bmnc(imn_nyq,2:ns-1)) / ds
       ! Simplistic extrapolation at the endpoints:
       d_B_d_s_mnc(1) = d_B_d_s_mnc(2)
       d_B_d_s_mnc(ns) = d_B_d_s_mnc(ns-1)

       if (non_Nyquist_mode_available) then
          ! R is on the full mesh:
          d_R_d_s_mnc(2:ns) = (pest%vmec%rmnc(imn,2:ns) - pest%vmec%rmnc(imn,1:ns-1)) / ds
          d_R_d_s_mnc(1) = 0

          ! Z is on the full mesh:
          d_Z_d_s_mns(2:ns) = (pest%vmec%zmns(imn,2:ns) - pest%vmec%zmns(imn,1:ns-1)) / ds
          d_Z_d_s_mns(1) = 0

          ! Lambda is on the half mesh:
          d_Lambda_d_s_mns(2:ns-1) = (pest%vmec%lmns(imn,3:ns) - pest%vmec%lmns(imn,2:ns-1)) / ds
          ! Simplistic extrapolation at the endpoints:
          d_Lambda_d_s_mns(1) = d_Lambda_d_s_mns(2)
          d_Lambda_d_s_mns(ns) = d_Lambda_d_s_mns(ns-1)
       else
          d_R_d_s_mnc = 0
          d_Z_d_s_mns = 0
          d_Lambda_d_s_mns = 0
       end if

       ! End of evaluating radial derivatives.

       do iz = pest%iz1,pest%iz2
          do ia = pest%ia1,pest%ia2
             angle = m * theta_vmec(ia,iz) - n * pest%vmec%nfp * pest%zeta(iz)
             cos_angle = cos(angle)
             sin_angle = sin(angle)

             do isurf = 1,2

                ! Handle |B|:
                temp = pest%vmec%bmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B(ia,iz) = B(ia,iz) + temp * cos_angle
                d_B_d_theta_vmec(ia,iz) = d_B_d_theta_vmec(ia,iz) - m * temp * sin_angle
                d_B_d_zeta(ia,iz) = d_B_d_zeta(ia,iz) + n * pest%vmec%nfp * temp * sin_angle       

                ! Handle Jacobian:
                temp = pest%vmec%gmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                sqrt_g(ia,iz) = sqrt_g(ia,iz) + temp * cos_angle

                ! Handle B sup theta:
                temp = pest%vmec%bsupumnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sup_theta_vmec(ia,iz) = B_sup_theta_vmec(ia,iz) + temp * cos_angle

                ! Handle B sup pest%zeta:
                temp = pest%vmec%bsupvmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sup_zeta(ia,iz) = B_sup_zeta(ia,iz) + temp * cos_angle

                ! Handle B sub theta:
                temp = pest%vmec%bsubumnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sub_theta_vmec(ia,iz) = B_sub_theta_vmec(ia,iz) + temp * cos_angle

                ! Handle B sub pest%zeta:
                temp = pest%vmec%bsubvmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sub_zeta(ia,iz) = B_sub_zeta(ia,iz) + temp * cos_angle

                ! Handle B sub psi.
                ! Unlike the other components of B, this one is on the full mesh.
                temp = pest%vmec%bsubsmns(imn_nyq,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                temp = temp*scale_factor
                B_sub_s(ia,iz) = B_sub_s(ia,iz) + temp * sin_angle

                ! Handle d B / d s
                ! Since bmnc is on the half mesh, its radial derivative is on the full mesh.
                temp = d_B_d_s_mnc(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                temp = temp*scale_factor
                d_B_d_s(ia,iz) = d_B_d_s(ia,iz) + temp * cos_angle

                ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                if (non_Nyquist_mode_available) then

                   ! Handle R, which is on the full mesh
                   temp = pest%vmec%rmnc(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   R(ia,iz) = R(ia,iz) + temp * cos_angle
                   d_R_d_theta_vmec(ia,iz) = d_R_d_theta_vmec(ia,iz) - temp * m * sin_angle
                   d_R_d_zeta(ia,iz)  = d_R_d_zeta(ia,iz)  + temp * n * pest%vmec%nfp * sin_angle

                   ! Handle Z, which is on the full mesh
                   temp = pest%vmec%zmns(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   Z(ia,iz) = Z(ia,iz) + temp * sin_angle  ! We don't actually need Z itself, only derivatives of Z.
                   d_Z_d_theta_vmec(ia,iz) = d_Z_d_theta_vmec(ia,iz) + temp * m * cos_angle
                   d_Z_d_zeta(ia,iz)  = d_Z_d_zeta(ia,iz)  - temp * n * pest%vmec%nfp * cos_angle

                   ! Handle Lambda:
                   temp = pest%vmec%lmns(imn,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   ! We don't need Lambda itself, just its derivatives.
                   d_Lambda_d_theta_vmec(ia,iz) = d_Lambda_d_theta_vmec(ia,iz) + m * temp * cos_angle
                   d_Lambda_d_zeta(ia,iz) = d_Lambda_d_zeta(ia,iz) - n * pest%vmec%nfp * temp * cos_angle       

                   ! Handle d R / d s
                   ! Since R is on the full mesh, its radial derivative is on the half mesh.
                   temp = d_R_d_s_mnc(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   d_R_d_s(ia,iz) = d_R_d_s(ia,iz) + temp * cos_angle

                   ! Handle d Z / d s
                   ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                   temp = d_Z_d_s_mns(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   d_Z_d_s(ia,iz) = d_Z_d_s(ia,iz) + temp * sin_angle

                   ! Handle d Lambda / d s
                   ! Since Lambda is on the half mesh, its radial derivative is on the full mesh.
                   temp = d_Lambda_d_s_mns(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   d_Lambda_d_s(ia,iz) = d_Lambda_d_s(ia,iz) + temp * sin_angle

                end if
             end do ! End surface loop
          end do ! End x3 loop
       end do ! End x2 loop

       ! -----------------------------------------------------
       ! Now consider the stellarator-asymmetric terms.
       ! -----------------------------------------------------

       if (pest%vmec%lasym) then

          ! Evaluate the radial derivatives we will need:

          ! B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
          ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.

          d_B_d_s_mns(2:ns-1) = (pest%vmec%bmns(imn_nyq,3:ns) - pest%vmec%bmns(imn_nyq,2:ns-1)) / ds
          ! Simplistic extrapolation at the endpoints:
          d_B_d_s_mns(1) = d_B_d_s_mns(2)
          d_B_d_s_mns(ns) = d_B_d_s_mns(ns-1)

          if (non_Nyquist_mode_available) then
             ! R is on the full mesh:
             d_R_d_s_mns(2:ns) = (pest%vmec%rmns(imn,2:ns) - pest%vmec%rmns(imn,1:ns-1)) / ds
             d_R_d_s_mns(1) = 0

             ! Z is on the full mesh:
             d_Z_d_s_mnc(2:ns) = (pest%vmec%zmnc(imn,2:ns) - pest%vmec%zmnc(imn,1:ns-1)) / ds
             d_Z_d_s_mnc(1) = 0

             ! Lambda is on the half mesh:
             d_Lambda_d_s_mnc(2:ns-1) = (pest%vmec%lmnc(imn_nyq,3:ns) - pest%vmec%lmnc(imn_nyq,2:ns-1)) / ds
             ! Simplistic extrapolation at the endpoints:
             d_Lambda_d_s_mnc(1) = d_Lambda_d_s_mnc(2)
             d_Lambda_d_s_mnc(ns) = d_Lambda_d_s_mnc(ns-1)
          else
             d_R_d_s_mns = 0
             d_Z_d_s_mnc = 0
             d_Lambda_d_s_mnc = 0
          end if

          ! End of evaluating radial derivatives.

          do iz = pest%iz1,pest%iz2
             do ia = pest%ia1,pest%ia2
                angle = m * theta_vmec(ia,iz) - n * pest%vmec%nfp * pest%zeta(iz)
                cos_angle = cos(angle)
                sin_angle = sin(angle)

                do isurf = 1,2

                   ! Handle |B|:
                   temp = pest%vmec%bmns(imn_nyq,vmec_radial_index_half(1)) * vmec_radial_weight_half(1)
                   temp = temp * scale_factor
                   B(ia,iz) = B(ia,iz) + temp * sin_angle
                   d_B_d_theta_vmec(ia,iz) = d_B_d_theta_vmec(ia,iz) + m * temp * cos_angle
                   d_B_d_zeta(ia,iz) = d_B_d_zeta(ia,iz) - n * pest%vmec%nfp * temp * cos_angle

                   ! Handle Jacobian:
                   temp = pest%vmec%gmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   sqrt_g(ia,iz) = sqrt_g(ia,iz) + temp * sin_angle

                   ! Handle B sup theta:
                   temp = pest%vmec%bsupumns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sup_theta_vmec(ia,iz) = B_sup_theta_vmec(ia,iz) + temp * sin_angle

                   ! Handle B sup pest%zeta:
                   temp = pest%vmec%bsupvmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sup_zeta(ia,iz) = B_sup_zeta(ia,iz) + temp * sin_angle

                   ! Handle B sub theta:
                   temp = pest%vmec%bsubumns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sub_theta_vmec(ia,iz) = B_sub_theta_vmec(ia,iz) + temp * sin_angle

                   ! Handle B sub pest%zeta:
                   temp = pest%vmec%bsubvmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sub_zeta(ia,iz) = B_sub_zeta(ia,iz) + temp * sin_angle

                   ! Handle B sub psi.
                   ! Unlike the other components of B, this one is on the full mesh.
                   temp = pest%vmec%bsubsmnc(imn_nyq,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   B_sub_s(ia,iz) = B_sub_s(ia,iz) + temp * cos_angle

                   ! Handle d B / d s.
                   ! Since bmns is on the half mesh, its radial derivative is on the full mesh.
                   temp = d_B_d_s_mns(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   d_B_d_s(ia,iz) = d_B_d_s(ia,iz) + temp * sin_angle

                   ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                   if (non_Nyquist_mode_available) then

                      ! Handle R, which is on the full mesh
                      temp = pest%vmec%rmns(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                      temp = temp*scale_factor
                      R(ia,iz) = R(ia,iz) + temp * sin_angle
                      d_R_d_theta_vmec(ia,iz) = d_R_d_theta_vmec(ia,iz) + temp * m * cos_angle
                      d_R_d_zeta(ia,iz)  = d_R_d_zeta(ia,iz)  - temp * n * pest%vmec%nfp * cos_angle

                      ! Handle Z, which is on the full mesh
                      temp = pest%vmec%zmnc(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                      temp = temp*scale_factor
                      Z(ia,iz) = Z(ia,iz) + temp * cos_angle   ! We don't actually need Z itself, only derivatives of Z.
                      d_Z_d_theta_vmec(ia,iz) = d_Z_d_theta_vmec(ia,iz) - temp * m * sin_angle
                      d_Z_d_zeta(ia,iz)  = d_Z_d_zeta(ia,iz)  + temp * n * pest%vmec%nfp * sin_angle

                      ! Handle Lambda, which is on the half mesh
                      temp = pest%vmec%lmnc(imn,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                      temp = temp*scale_factor
                      ! We don't actually need Lambda itself, only derivatives of Lambda.
                      d_Lambda_d_theta_vmec(ia,iz) = d_Lambda_d_theta_vmec(ia,iz) - temp * m * sin_angle
                      d_Lambda_d_zeta(ia,iz)  = d_Lambda_d_zeta(ia,iz)  + temp * n * pest%vmec%nfp * sin_angle

                      ! Handle d R / d s.
                      ! Since R is on the full mesh, its radial derivative is on the half mesh.
                      temp = d_R_d_s_mns(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                      temp = temp*scale_factor
                      d_R_d_s(ia,iz) = d_R_d_s(ia,iz) + temp * sin_angle

                      ! Handle d Z / d s.
                      ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                      temp = d_Z_d_s_mnc(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                      temp = temp*scale_factor
                      d_Z_d_s(ia,iz) = d_Z_d_s(ia,iz) + temp * cos_angle

                      ! Handle d Lambda / d s.
                      ! Since Lambda is on the half mesh, its radial derivative is on the full mesh.
                      temp = d_Lambda_d_s_mnc(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                      temp = temp*scale_factor
                      d_Lambda_d_s(ia,iz) = d_Lambda_d_s(ia,iz) + temp * cos_angle
                   end if
                end do ! End surface loop
             end do ! End x3 loop
          end do ! End x2 loop
       end if ! End lasym
    end do ! End imn_nyq loop 

    !*********************************************************************
    ! Sanity check: If the conversion to theta_pest has been done 
    ! correctly, we should find that 
    ! (B dot grad theta_pest) / (B dot grad pest%zeta) = pest%iota.
    ! Let's verify this:
    !*********************************************************************
    
    if (test) then
      allocate(B_dot_grad_theta_pest_over_B_dot_grad_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2))
      ! Compute (B dot grad theta_pest) / (B dot grad pest%zeta):
      B_dot_grad_theta_pest_over_B_dot_grad_zeta = &
        & (B_sup_theta_vmec * (1 + d_Lambda_d_theta_vmec) + B_sup_zeta * d_Lambda_d_zeta) / B_sup_zeta 
      temp2D = pest%iota
      call test_arrays(B_dot_grad_theta_pest_over_B_dot_grad_zeta, temp2D, .false., 0.01, 'pest%iota')
      deallocate(B_dot_grad_theta_pest_over_B_dot_grad_zeta)
    end if

    !*********************************************************************
    ! Using R(theta,pest%zeta) and Z(theta,pest%zeta), compute the Cartesian
    ! components of the gradient basis vectors using the dual relations:
    !*********************************************************************

    sqrt_s = sqrt(pest%s0)

    do iz = pest%iz1,pest%iz2
       cos_angle = cos(pest%zeta(iz))
       sin_angle = sin(pest%zeta(iz))

       ! X = R * cos(pest%zeta)
       d_X_d_theta_vmec(:,iz) = d_R_d_theta_vmec(:,iz) * cos_angle
       d_X_d_zeta(:,iz) = d_R_d_zeta(:,iz) * cos_angle - R(:,iz) * sin_angle
       d_X_d_s(:,iz) = d_R_d_s(:,iz) * cos_angle

       ! Y = R * sin(pest%zeta)
       d_Y_d_theta_vmec(:,iz) = d_R_d_theta_vmec(:,iz) * sin_angle
       d_Y_d_zeta(:,iz) = d_R_d_zeta(:,iz) * sin_angle + R(:,iz) * cos_angle
       d_Y_d_s(:,iz) = d_R_d_s(:,iz) * sin_angle

    end do

    ! Use the dual relations to get the Cartesian components of grad s, grad theta_vmec, and grad pest%zeta:
    grad_s_X = (d_Y_d_theta_vmec * d_Z_d_zeta - d_Z_d_theta_vmec * d_Y_d_zeta) / sqrt_g
    grad_s_Y = (d_Z_d_theta_vmec * d_X_d_zeta - d_X_d_theta_vmec * d_Z_d_zeta) / sqrt_g
    grad_s_Z = (d_X_d_theta_vmec * d_Y_d_zeta - d_Y_d_theta_vmec * d_X_d_zeta) / sqrt_g

    grad_theta_vmec_X = (d_Y_d_zeta * d_Z_d_s - d_Z_d_zeta * d_Y_d_s) / sqrt_g
    grad_theta_vmec_Y = (d_Z_d_zeta * d_X_d_s - d_X_d_zeta * d_Z_d_s) / sqrt_g
    grad_theta_vmec_Z = (d_X_d_zeta * d_Y_d_s - d_Y_d_zeta * d_X_d_s) / sqrt_g

    grad_zeta_X = (d_Y_d_s * d_Z_d_theta_vmec - d_Z_d_s * d_Y_d_theta_vmec) / sqrt_g
    grad_zeta_Y = (d_Z_d_s * d_X_d_theta_vmec - d_X_d_s * d_Z_d_theta_vmec) / sqrt_g
    grad_zeta_Z = (d_X_d_s * d_Y_d_theta_vmec - d_Y_d_s * d_X_d_theta_vmec) / sqrt_g
    ! End of the dual relations.

    if (test) then
      ! Sanity check: grad_zeta_X should be -sin(pest%zeta) / R:
      do iz = pest%iz1,pest%iz2
         temp2D(:,iz) = -sin(pest%zeta(iz)) / R(:,iz)
      end do
      call test_arrays(grad_zeta_X, temp2D, .false., 1.0e-2, 'grad_zeta_X')
      grad_zeta_X = temp2D ! We might as well use the exact value, which is in temp2D.

      ! Sanity check: grad_zeta_Y should be cos(pest%zeta) / R:
      do iz = pest%iz1,pest%iz2
         temp2D(:,iz) = cos(pest%zeta(iz)) / R(:,iz)
      end do
      call test_arrays(grad_zeta_Y, temp2D, .false., 1.0e-2, 'grad_zeta_Y')
      grad_zeta_Y = temp2D ! We might as well use the exact value, which is in temp2D.

      ! grad_zeta_Z should be 0:
      call test_arrays(grad_zeta_Z, temp2D, .true., 1.0e-14, 'grad_zeta_Z')
    end if
    grad_zeta_Z = 0
    
    !*********************************************************************
    ! Compute the Cartesian components of other quantities we need:
    !*********************************************************************

    grad_psi_X = grad_s_X * edge_toroidal_flux_over_2pi
    grad_psi_Y = grad_s_Y * edge_toroidal_flux_over_2pi
    grad_psi_Z = grad_s_Z * edge_toroidal_flux_over_2pi

    ! Form grad alpha = grad (theta_vmec + Lambda - pest%iota * pest%zeta)
    do iz = pest%iz1,pest%iz2
       grad_alpha_X(:,iz) = (d_Lambda_d_s(:,iz) - pest%zeta(iz) * d_iota_d_s) * grad_s_X(:,iz)
       grad_alpha_Y(:,iz) = (d_Lambda_d_s(:,iz) - pest%zeta(iz) * d_iota_d_s) * grad_s_Y(:,iz)
       grad_alpha_Z(:,iz) = (d_Lambda_d_s(:,iz) - pest%zeta(iz) * d_iota_d_s) * grad_s_Z(:,iz)
    end do
    grad_alpha_X = grad_alpha_X + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_X + (-pest%iota + d_Lambda_d_zeta) * grad_zeta_X
    grad_alpha_Y = grad_alpha_Y + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Y + (-pest%iota + d_Lambda_d_zeta) * grad_zeta_Y
    grad_alpha_Z = grad_alpha_Z + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Z + (-pest%iota + d_Lambda_d_zeta) * grad_zeta_Z

    grad_B_X = d_B_d_s * grad_s_X + d_B_d_theta_vmec * grad_theta_vmec_X + d_B_d_zeta * grad_zeta_X
    grad_B_Y = d_B_d_s * grad_s_Y + d_B_d_theta_vmec * grad_theta_vmec_Y + d_B_d_zeta * grad_zeta_Y
    grad_B_Z = d_B_d_s * grad_s_Z + d_B_d_theta_vmec * grad_theta_vmec_Z + d_B_d_zeta * grad_zeta_Z

    B_X = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_X_d_zeta + (pest%iota - d_Lambda_d_zeta) * d_X_d_theta_vmec) / sqrt_g
    B_Y = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_Y_d_zeta + (pest%iota - d_Lambda_d_zeta) * d_Y_d_theta_vmec) / sqrt_g
    B_Z = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_Z_d_zeta + (pest%iota - d_Lambda_d_zeta) * d_Z_d_theta_vmec) / sqrt_g

    if (test) then
    !*********************************************************************
    ! Sanity tests: Verify that the Jacobian equals the appropriate
    ! cross product of the basis vectors.
    !*********************************************************************

      temp2D = 0 &
           + d_X_d_s * d_Y_d_theta_vmec * d_Z_d_zeta &
           + d_Y_d_s * d_Z_d_theta_vmec * d_X_d_zeta &
           + d_Z_d_s * d_X_d_theta_vmec * d_Y_d_zeta &
           - d_Z_d_s * d_Y_d_theta_vmec * d_X_d_zeta &
           - d_X_d_s * d_Z_d_theta_vmec * d_Y_d_zeta &
           - d_Y_d_s * d_X_d_theta_vmec * d_Z_d_zeta
      call test_arrays(sqrt_g, temp2D, .false., 1.0e-2, 'sqrt_g')

      temp2D = 0 &
           + grad_s_X * grad_theta_vmec_Y * grad_zeta_Z &
           + grad_s_Y * grad_theta_vmec_Z * grad_zeta_X &
           + grad_s_Z * grad_theta_vmec_X * grad_zeta_Y &
           - grad_s_Z * grad_theta_vmec_Y * grad_zeta_X &
           - grad_s_X * grad_theta_vmec_Z * grad_zeta_Y &
           - grad_s_Y * grad_theta_vmec_X * grad_zeta_Z
      call test_arrays(1/sqrt_g, temp2D, .false., 1.0e-2, '1/sqrt_g')
    end if

    ! Things for GIST
    !jac_gist_inv = (L_ref**3)/(2*safety_factor_q)*(1 + d_Lambda_d_theta_vmec)/sqrt_g

    !d_B_d_par = -safety_factor_q/L_ref *jac_gist_inv / B * d_B_d_zeta 

        
    if (test) then
    !*********************************************************************
    ! Sanity tests: Verify that 
    ! \vec{B} dot (each of the covariant and contravariant basis vectors)
    ! matches the corresponding term from VMEC.
    !*********************************************************************

      call test_arrays(B_X * d_X_d_theta_vmec + B_Y * d_Y_d_theta_vmec + B_Z * d_Z_d_theta_vmec, B_sub_theta_vmec, .false., 1.0e-2, 'B_sub_theta_vmec')
      call test_arrays(B_X * d_X_d_s          + B_Y * d_Y_d_s          + B_Z * d_Z_d_s,          B_sub_s,          .false., 1.0e-2, 'B_sub_s')
      call test_arrays(B_X * d_X_d_zeta       + B_Y * d_Y_d_zeta       + B_Z * d_Z_d_zeta,       B_sub_zeta,       .false., 1.0e-2, 'B_sub_zeta')

      call test_arrays(B_X *          grad_s_X + B_Y *          grad_s_Y + B_Z *          grad_s_Z,           temp2D,  .true., 1.0e-2, 'B_sup_s')
      call test_arrays(B_X *       grad_zeta_X + B_Y *       grad_zeta_Y + B_Z *       grad_zeta_Z,       B_sup_zeta, .false., 1.0e-2, 'B_sup_zeta')
      call test_arrays(B_X * grad_theta_vmec_X + B_Y * grad_theta_vmec_Y + B_Z * grad_theta_vmec_Z, B_sup_theta_vmec, .false., 1.0e-2, 'B_sup_theta_vmec')
    end if

    !*********************************************************************
    ! For gbdrift, we need \vect{B} cross grad |B| dot grad alpha.
    ! For cvdrift, we also need \vect{B} cross grad s dot grad alpha.
    ! Let us compute both of these quantities 2 ways, and make sure the two
    ! approaches give the same answer (within some tolerance).
    !*********************************************************************

    B_cross_grad_s_dot_grad_alpha = (B_sub_zeta * (1 + d_Lambda_d_theta_vmec) &
         - B_sub_theta_vmec * (d_Lambda_d_zeta - pest%iota) ) / sqrt_g

    if (test) then
    B_cross_grad_s_dot_grad_alpha_alternate = 0 &
         + B_X * grad_s_Y * grad_alpha_Z &
         + B_Y * grad_s_Z * grad_alpha_X &
         + B_Z * grad_s_X * grad_alpha_Y &
         - B_Z * grad_s_Y * grad_alpha_X &
         - B_X * grad_s_Z * grad_alpha_Y &
         - B_Y * grad_s_X * grad_alpha_Z 

      call test_arrays(B_cross_grad_s_dot_grad_alpha, B_cross_grad_s_dot_grad_alpha_alternate, &
         .false., 1.0e-2, 'B_cross_grad_s_dot_grad_alpha')
    end if

    do iz = pest%iz1,pest%iz2
       B_cross_grad_B_dot_grad_alpha(:,iz) = 0 &
            + (B_sub_s(:,iz) * d_B_d_theta_vmec(:,iz) * (d_Lambda_d_zeta(:,iz) - pest%iota) &
            + B_sub_theta_vmec(:,iz) * d_B_d_zeta(:,iz) * (d_Lambda_d_s(:,iz) - pest%zeta(iz) * d_iota_d_s) &
            + B_sub_zeta(:,iz) * d_B_d_s(:,iz) * (1 + d_Lambda_d_theta_vmec(:,iz)) &
            - B_sub_zeta(:,iz) * d_B_d_theta_vmec(:,iz) * (d_Lambda_d_s(:,iz) - pest%zeta(iz) * d_iota_d_s) &
            - B_sub_theta_vmec(:,iz) * d_B_d_s(:,iz) * (d_Lambda_d_zeta(:,iz) - pest%iota) &
            - B_sub_s(:,iz) * d_B_d_zeta(:,iz) * (1 + d_Lambda_d_theta_vmec(:,iz))) / sqrt_g(:,iz)
    end do

    if (test) then
    B_cross_grad_B_dot_grad_alpha_alternate = 0 &
         + B_X * grad_B_Y * grad_alpha_Z &
         + B_Y * grad_B_Z * grad_alpha_X &
         + B_Z * grad_B_X * grad_alpha_Y &
         - B_Z * grad_B_Y * grad_alpha_X &
         - B_X * grad_B_Z * grad_alpha_Y &
         - B_Y * grad_B_X * grad_alpha_Z 

    call test_arrays(B_cross_grad_B_dot_grad_alpha, B_cross_grad_B_dot_grad_alpha_alternate, &
         .false., 1.0e-2, 'B_cross_grad_B_dot_grad_alpha')
    end if

    pest%bmag(:,:,1) = B    
    pest%jac(:,:,1) = sqrt_g
    pest%gss(:,:,1) = 4.0/((pest%B_ref**2)*(pest%L_ref**4))*(grad_psi_X * grad_psi_X + grad_psi_Y * grad_psi_Y + grad_psi_Z * grad_psi_Z)
    pest%gsa(:,:,1) = 2.0/(pest%B_ref*(pest%L_ref**2))*(grad_psi_X * grad_alpha_X + grad_psi_Y * grad_alpha_Y + grad_psi_Z * grad_alpha_Z)
    pest%gaa(:,:,1) = grad_alpha_X * grad_alpha_X + grad_alpha_Y * grad_alpha_Y + grad_alpha_Z * grad_alpha_Z
    pest%gsz(:,:,1) = 2.0/(pest%B_ref*(pest%L_ref**2))*(grad_psi_X * grad_zeta_X + grad_psi_Y * grad_zeta_Y + grad_psi_Z * grad_zeta_Z)
    pest%gaz(:,:,1) = grad_alpha_X * grad_zeta_X + grad_alpha_Y * grad_zeta_Y + grad_alpha_Z * grad_zeta_Z
    pest%gzz(:,:,1) = grad_zeta_X * grad_zeta_X + grad_zeta_Y * grad_zeta_Y + grad_zeta_Z * grad_zeta_Z
    pest%d_L_d_theta_v(:,:,1) = d_Lambda_d_theta_vmec
    print *, d_Lambda_d_theta_vmec
    !*********************************************************************
    ! Finally, assemble the quantities needed for gs2.
    !*********************************************************************

    ! See the latex note gs2_full_surface_stellarator_geometry in the "doc" directory for a derivation of the formulae that follow.

!    bmag = B / pest%B_ref
!
!    gradpar = pest%L_ref * B_sup_zeta / B
!
!    gds2 = (grad_alpha_X * grad_alpha_X + grad_alpha_Y * grad_alpha_Y + grad_alpha_Z * grad_alpha_Z) &
!         * pest%L_ref * pest%L_ref * pest%s0
!
!    gds21 = (grad_alpha_X * grad_psi_X + grad_alpha_Y * grad_psi_Y + grad_alpha_Z * grad_psi_Z) &
!         * sign_toroidal_flux * shat / pest%B_ref
!
!    gds22 = (grad_psi_X * grad_psi_X + grad_psi_Y * grad_psi_Y + grad_psi_Z * grad_psi_Z) &
!         * shat * shat / (pest%L_ref * pest%L_ref * pest%B_ref * pest%B_ref * pest%s0)
!
!    gbdrift = sign_toroidal_flux * 2 * pest%B_ref * pest%L_ref * pest%L_ref * sqrt_s * B_cross_grad_B_dot_grad_alpha &
!         / (B * B * B)
!
!    gbdrift0 = abs(edge_toroidal_flux_over_2pi) &
!      & * (B_sub_theta_vmec * d_B_d_zeta - B_sub_zeta * d_B_d_theta_vmec) / sqrt_g &
!      & * 2 * shat / (B * B * B * sqrt_s)
!    ! In the above 2-line expression for gbdrift0, the first line is \vec{B} \times \nabla B \cdot \nabla \psi.
!
!    cvdrift = gbdrift + sign_toroidal_flux * 2 * pest%B_ref * pest%L_ref * pest%L_ref * sqrt_s * mu_0 * d_pressure_d_s &
!         * B_cross_grad_s_dot_grad_alpha / (B * B * B * B)
!
!    cvdrift0 = gbdrift0

    !*********************************************************************
    ! Copy the surface quantities
    !*********************************************************************

    !Rsurf = R
    !Zsurf = Z

    !*********************************************************************
    ! Free all arrays that were allocated.
    !*********************************************************************

    deallocate(B)
    deallocate(temp2D)
    deallocate(sqrt_g)
    deallocate(R)
    deallocate(Z)
    deallocate(d_B_d_theta_vmec)
    deallocate(d_B_d_zeta)
    deallocate(d_B_d_s)
    deallocate(d_R_d_theta_vmec)
    deallocate(d_R_d_zeta)
    deallocate(d_R_d_s)
    deallocate(d_Z_d_theta_vmec)
    deallocate(d_Z_d_zeta)
    deallocate(d_Z_d_s)
    deallocate(d_Lambda_d_theta_vmec)
    deallocate(d_Lambda_d_zeta)
    deallocate(d_Lambda_d_s)
    deallocate(B_sub_s)
    deallocate(B_sub_theta_vmec)
    deallocate(B_sub_zeta)
    deallocate(B_sup_theta_vmec)
    deallocate(B_sup_zeta)

    deallocate(d_B_d_s_mnc)
    deallocate(d_B_d_s_mns)
    deallocate(d_R_d_s_mnc)
    deallocate(d_R_d_s_mns)
    deallocate(d_Z_d_s_mnc)
    deallocate(d_Z_d_s_mns)
    deallocate(d_Lambda_d_s_mnc)
    deallocate(d_Lambda_d_s_mns)

    deallocate(d_X_d_s)
    deallocate(d_X_d_theta_vmec)
    deallocate(d_X_d_zeta)
    deallocate(d_Y_d_s)
    deallocate(d_Y_d_theta_vmec)
    deallocate(d_Y_d_zeta)

    deallocate(grad_s_X)
    deallocate(grad_s_Y)
    deallocate(grad_s_Z)
    deallocate(grad_theta_vmec_X)
    deallocate(grad_theta_vmec_Y)
    deallocate(grad_theta_vmec_Z)
    deallocate(grad_zeta_X)
    deallocate(grad_zeta_Y)
    deallocate(grad_zeta_Z)
    deallocate(grad_psi_X)
    deallocate(grad_psi_Y)
    deallocate(grad_psi_Z)
    deallocate(grad_alpha_X)
    deallocate(grad_alpha_Y)
    deallocate(grad_alpha_Z)

    deallocate(B_X)
    deallocate(B_Y)
    deallocate(B_Z)
    deallocate(grad_B_X)
    deallocate(grad_B_Y)
    deallocate(grad_B_Z)
    deallocate(B_cross_grad_B_dot_grad_alpha)
    deallocate(B_cross_grad_B_dot_grad_alpha_alternate)
    deallocate(B_cross_grad_s_dot_grad_alpha)
    deallocate(B_cross_grad_s_dot_grad_alpha_alternate)

    deallocate(theta_vmec)

    if (verbose) print *,"Leaving compute_pest_sfl."


  contains

    subroutine test_arrays(array1, array2, should_be_0, tolerance, name)
      ! This subroutine is used for verifying the geometry arrays.
      ! When should_be_0 = .true., the subroutine verifies that |array1| = 0 to within 
      !     an absolute tolerance specified by 'tolerance'. array2 is ignored in this case.
      ! When should_be_0 = .false., the subroutine verifies that array1 = array2 
      !     to within a relative tolerance specified by 'tolerance'.

      implicit none

      real(dp), dimension(pest%ia1:pest%ia2,pest%iz1:pest%iz2) :: array1, array2
      real(dp) :: tolerance
      character(len=*) :: name
      logical :: should_be_0
      real(dp) :: max_value, max_difference

      if (should_be_0) then
         max_value = maxval(abs(array1))
         if (verbose) print *,"  maxval(abs(",trim(name),")):",max_value,"(should be << 1.)"
         if (max_value > tolerance) then
            print *,"Error! ",trim(name)," should be 0, but instead it is:"
            do ia = pest%ia1,pest%ia2
               print *,array1(ia,:)
            end do
            stop
         end if
      else
         max_difference = maxval(abs(array1 - array2)) / maxval(abs(array1) + abs(array2))
         if (verbose) print *,"  Relative difference between two methods for computing ",trim(name),":",max_difference,"(should be << 1.)"
         if (max_difference > tolerance) then
            print *,"Error! Two methods for computing ",trim(name)," disagree. Here comes method 1:"
            do ia = pest%ia1,pest%ia2
               print *,array1(ia,:)
            end do
            print *,"Here comes method 2:"
            do ia = pest%ia1,pest%ia2
               print *,array2(ia,:)
            end do
            print *,"Here comes the difference:"
            do ia = pest%ia1,pest%ia2
               print *,array1(ia,:) - array2(ia,:)
            end do
            stop
         end if
      end if

    end subroutine test_arrays

  end subroutine compute_pest_sfl

  ! --------------------------------------------------------------------------



  real(dp) function fzero_residual(theta_vmec_try, vmec) result(residual)
    ! Note that lmns and lmnc use the non-Nyquist xm, xn, and pest%vmec%mnmax.
    ! Also note that lmns and lmnc are on the HALF grid.

    implicit none
    
    real(dp), intent(in) :: theta_vmec_try
    real(dp) :: angle, sinangle, cosangle
    integer :: imn, which_surface
    type(VMEC_Obj), intent(in) :: vmec

    ! residual = (theta_pest based on theta_vmec_try) - theta_pest_target = theta_vmec_try + Lambda - theta_pest_target
    residual = theta_vmec_try - theta_pest_target

    do imn = 1, vmec%mnmax
       angle = vmec%xm(imn)*theta_vmec_try - vmec%xn(imn)*zeta0
       sinangle = sin(angle)
       cosangle = cos(angle)
       do which_surface = 1,2
          residual = residual + vmec_radial_weight_half(which_surface) * vmec%lmns(imn,vmec_radial_index_half(which_surface)) * sinangle
          if (vmec%lasym) then
             residual = residual + vmec_radial_weight_half(which_surface) * vmec%lmnc(imn,vmec_radial_index_half(which_surface)) * cosangle
          end if
       end do
    end do

  end function

end module
