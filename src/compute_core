! This code is based off of the 
! vmec2sfl.f90
! Written by Matt Landreman, University of Maryland
! Initial code written August 2017.


! This code has been further modified by B.J. Faber to generalize to different
! coordinate systems.  As such, some of the variable names have bene modified.
! Skip down ~25 lines for detailed description of the input and output parameters.

module compute_core 

  use types: dp, pi, mu_0
  use eq_object: Eq_Object
  implicit none

  public :: compute_surface_quantities, compute_sfl_quantities
  
  private

  ! Module wide to variables
  integer :: sign_toroidal_flux
  real(dp) :: theta_pest_target, zeta0, edge_toroidal_flux_over_2pi
  real(dp) :: iota, ds, d_pressure_d_s, d_iota_d_s
  real(dp), dimension(2) :: vmec_radial_weight_full, vmec_radial_weight_half
  integer, dimension(2) :: vmec_radial_index_full, vmec_radial_index_half

contains

  subroutine compute_surface_quantities(desired_normalized_toroidal_flux,&
    &vmec_surface_option, Eq_Obj)

    !use read_wout_mod, nzgrid_vmec => Eq_Obj%nx3 ! VMEC has a variable ngzrid which conflicts with our Eq_Obj%nx3, so rename vmec's version.

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

    type(Eq_Object), intent(inout) :: Eq_Obj


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
    !*********************************************************************
    ! VMEC variables of interest:
    ! ns = number of flux surfaces used by VMEC
    ! nfp = number of field periods, e.g. 5 for W7-X, 4 for HSX
    ! iotas = rotational transform (1/q) on the half grid.
    ! iotaf = rotational transform on the full grid.
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

    if (verbose) print *,"Entering subroutine compute_surface_quantities."

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

    !*********************************************************************
    ! Read in everything from the vmec wout file using libstell.
    !*********************************************************************
    
    !nfp = nfp_vmec
    !local_nfp = nfp

    ! There is a bug in libstell read_wout_file for ASCII-format wout files, in which the xm_nyq and xn_nyq arrays are sometimes
    ! not populated. The next few lines here provide a workaround:
    if (maxval(abs(xm_nyq)) < 1 .and. maxval(abs(xn_nyq)) < 1) then
       if (mnmax_nyq == mnmax) then
          if (verbose) print *,"xm_nyq and xn_nyq arrays are not populated in the wout file. Using xm and xn instead."
          xm_nyq = xm
          xn_nyq = xn
       else
          print *,"Error! xm_nyq and xn_nyq arrays are not populated in the wout file, and mnmax_nyq != mnmax."
          stop
       end if
    end if

    edge_toroidal_flux_over_2pi = Eq_Obj%VMEC_Obj%phi(ns) / (2*pi) * Eq_Obj%VMEC_Obj%isigng ! isigns is called signgs in the wout*.nc file. Why is this signgs here?

    ! this gives the sign of the edge toroidal flux
    sign_toroidal_flux = int(sign(1.1,edge_toroidal_flux_over_2pi))
    write (*,*) 'sign_toroidal_flux', sign_toroidal_flux

    ! Set reference length and magnetic field for GS2's normalization, 
    ! using the choices made by Pavlos Xanthopoulos in GIST:
    Eq_Obj%major_R = Rmajor
    Eq_Obj%minor_r = Aminor ! Note that 'Aminor' in read_wout_mod is called 'Aminor_p' in the wout*.nc file.

    if (norm_type == 'minor_r') then
      Eq_Obj%L_reference => Eq_Obj%minor_r
    else
      Eq_Obj%L_reference => Eq_Obj%major_R
    end if

    Eq_Obj%B_reference = 2 * abs(edge_toroidal_flux_over_2pi) / (Eq_Obj%L_reference * Eq_Obj%L_reference)
    if (verbose) then
       print *,"  Reference length for normalization:",Eq_Obj%L_reference," meters."
       print *,"  Reference magnetic field strength normalization:",Eq_Obj%B_reference," Tesla."
    end if

    ! --------------------------------------------------------------------------------
    ! Do some sanity checking to ensure the VMEC arrays have some expected properties.
    ! --------------------------------------------------------------------------------

    ! 'phi' is vmec's array of the toroidal flux (not divided by 2pi!) on vmec's radial grid.
    if (abs(phi(1)) > 1d-14) then
       print *,"Error! VMEC phi array does not begin with 0."
       print *,"phi:",phi
       stop
    end if

    dphi = phi(2) - phi(1)
    do j=3,ns
       if (abs(phi(j)-phi(j-1)-dphi) > 1d-11) then
          print *,"Error! VMEC phi array is not uniformly spaced."
          print *,"phi:",phi
          stop
       end if
    end do

    ! The variable called 'phips' in the wout file is called just 'phip' in read_wout_mod.F.
    ! phips is on the half-mesh, so skip first point.
    do j=2,ns
       if (abs(phip(j)+phi(ns)/(2*pi)) > 1d-11) then
          print *,"Error! VMEC phips array is not constant and equal to -phi(ns)/(2*pi)."
          print *,"phip(s):",phip
          stop
       end if
    end do

    ! The first mode in the m and n arrays should be m=n=0:
    if (xm(1) .ne. 0) stop "First element of xm in the wout file should be 0."
    if (xn(1) .ne. 0) stop "First element of xn in the wout file should be 0."
    if (xm_nyq(1) .ne. 0) stop "First element of xm_nyq in the wout file should be 0."
    if (xn_nyq(1) .ne. 0) stop "First element of xn_nyq in the wout file should be 0."

    ! Lambda should be on the half mesh, so its value at radial index 1 should be 0 for all (m,n)
    if (maxval(abs(lmns(:,1))) > 0) then
       print *,"Error! Expected lmns to be on the half mesh, but its value at radial index 1 is nonzero."
       print *,"Here comes lmns(:,1):", lmns(:,1)
       stop
    end if
    if (lasym) then
       if (maxval(abs(lmnc(:,1))) > 0) then
          print *,"Error! Expected lmnc to be on the half mesh, but its value at radial index 1 is nonzero."
          print *,"Here comes lmnc(:,1):", lmnc(:,1)
          stop
       end if
    end if

    ! --------------------------------------------------------------------------------
    ! End of sanity checks.
    ! --------------------------------------------------------------------------------

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
       normalized_toroidal_flux_used = desired_normalized_toroidal_flux

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

       normalized_toroidal_flux_used = normalized_toroidal_flux_half_grid(index)
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

       normalized_toroidal_flux_used = normalized_toroidal_flux_full_grid(index)
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
    if (normalized_toroidal_flux_used>1) then
       stop "Error! normalized_toroidal_flux_used cannot be >1"
    elseif (normalized_toroidal_flux_used<0) then
       stop "Error! normalized_toroidal_flux_used cannot be <0"
    elseif (normalized_toroidal_flux_used==1) then
       vmec_radial_index_full(1) = ns-1
       vmec_radial_index_full(2) = ns
       vmec_radial_weight_full(1) = 0.0d0
    else
       ! normalized_toroidal_flux_used is >= 0 and <1
       ! This is the most common case.
       vmec_radial_index_full(1) = floor(normalized_toroidal_flux_used*(ns-1))+1
       vmec_radial_index_full(2) = vmec_radial_index_full(1) + 1
       vmec_radial_weight_full(1) = vmec_radial_index_full(1) - normalized_toroidal_flux_used*(ns-1.0d0)
    end if
    vmec_radial_weight_full(2) = 1.0d0 - vmec_radial_weight_full(1)

    ! Handle quantities for the half grid
    if (normalized_toroidal_flux_used < normalized_toroidal_flux_half_grid(1)) then
       print *,"Warning: extrapolating beyond the end of VMEC's half grid."
       print *,"(Extrapolating towards the magnetic axis.) Results are likely to be inaccurate."

       ! We start at element 2 since element 1 is always 0 for quantities on the half grid.
       vmec_radial_index_half(1) = 2
       vmec_radial_index_half(2) = 3
       vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(2) - normalized_toroidal_flux_used) / (normalized_toroidal_flux_half_grid(2) - normalized_toroidal_flux_half_grid(1))

    elseif (normalized_toroidal_flux_used > normalized_toroidal_flux_half_grid(ns-1)) then
       print *,"Warning: extrapolating beyond the end of VMEC's half grid."
       print *,"(Extrapolating towards the last closed flux surface.) Results may be inaccurate."
       vmec_radial_index_half(1) = ns-1
       vmec_radial_index_half(2) = ns
       vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(ns-1) - normalized_toroidal_flux_used) &
            / (normalized_toroidal_flux_half_grid(ns-1) - normalized_toroidal_flux_half_grid(ns-2))

    elseif (normalized_toroidal_flux_used == normalized_toroidal_flux_half_grid(ns-1)) then
       ! We are exactly at the last point of the half grid
       vmec_radial_index_half(1) = ns-1
       vmec_radial_index_half(2) = ns
       vmec_radial_weight_half(1) = 0.0d0
    else
       ! normalized_toroidal_flux_used is inside the half grid.
       ! This is the most common case.
       vmec_radial_index_half(1) = floor(normalized_toroidal_flux_used*(ns-1) + 0.5d+0)+1
       if (vmec_radial_index_half(1) < 2) then
          ! This can occur sometimes due to roundoff error.
          vmec_radial_index_half(1) = 2
       end if
       vmec_radial_index_half(2) = vmec_radial_index_half(1) + 1
       vmec_radial_weight_half(1) = vmec_radial_index_half(1) - normalized_toroidal_flux_used*(ns-1.0d0) - (0.5d+0)
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

    Eq_Obj%iota = iotas(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + iotas(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    if (verbose) print *,"  iota =",Eq_Obj%iota
    Eq_Obj%safety_factor_q = 1/Eq_Obj%iota

    allocate(d_iota_d_s_on_half_grid(ns))
    d_iota_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    if (verbose) print *,"  ds =",ds
    d_iota_d_s_on_half_grid(2:ns) = (iotaf(2:ns) - iotaf(1:ns-1)) / ds
    d_iota_d_s =  &
         d_iota_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + d_iota_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    deallocate(d_iota_d_s_on_half_grid)
    if (verbose) print *,"  d iota / d s =",d_iota_d_s
    ! shat = (r/q)(dq/dr) where r = a sqrt(s).
    !      = - (r/iota) (d iota / d r) = -2 (s/iota) (d iota / d s)
    Eq_Obj%shat = (-2 * normalized_toroidal_flux_used / iota) * d_iota_d_s

    allocate(d_pressure_d_s_on_half_grid(ns))
    d_pressure_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    d_pressure_d_s_on_half_grid(2:ns) = (presf(2:ns) - presf(1:ns-1)) / ds
    d_pressure_d_s =  &
         d_pressure_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + d_pressure_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    deallocate(d_pressure_d_s_on_half_grid)
    if (verbose) print *,"  d pressure / d s =",d_pressure_d_s

    deallocate(normalized_toroidal_flux_full_grid)
    deallocate(normalized_toroidal_flux_half_grid)
  end subroutine compute_surface_quantities


!*******************************************************************************
! Subroutine compute_sfl_quantities
!*******************************************************************************
  subroutine compute_sfl_quantities(Eq_Obj%nx2, Eq_Obj%nx3, zeta_center, &
    & number_of_field_periods_to_include,Eq_Obj)

    use read_wout_mod, Eq_Obj%nx3_vmec => Eq_Obj%nx3  ! VMEC has a variable Eq_Obj%nx3 which conflicts with our Eq_Obj%nx3, so rename vmec's version.

    implicit none
    !***************************************************************************
    ! Input parameters
    !***************************************************************************

    ! Eq_Obj%nx2 is the number of grid points in the alpha coordinate:
    integer, intent(in) :: Eq_Obj%nx2

    ! The zeta grid has Eq_Obj%nx3*2+1 points, including the "repeated" point at index -Eq_Obj%nx3 and +Eq_Obj%nx3.
    integer, intent(in) :: Eq_Obj%nx3

    ! The zeta domain is centered at zeta_center. Setting zeta_center = 2*pi*N/nfp for any integer N should
    ! yield identical results to setting zeta_center = 0, where nfp is the number of field periods (as in VMEC).
    real(dp) :: zeta_center

    ! If number_of_field_periods_to_include is > 0, then this parameter does what you think:
    ! the extent of the toroidal in zeta will be 2*pi*number_of_field_periods_to_include/nfp.
    ! If number_of_field_periods_to_include is <= 0, the entire 2*pi toroidal domain will be included.
    real(dp) :: number_of_field_periods_to_include


    !*********************************************************************
    ! Variables used internally by this subroutine
    !*********************************************************************

    integer :: j, iz, ia, which_surface, isurf, m, n, imn, imn_nyq
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
    !*********************************************************************
    ! Read in everything from the vmec wout file using libstell.
    !*********************************************************************

    !if (verbose) print *,"  About to read VMEC wout file ",trim(vmec_filename)
    !call read_wout_file(vmec_filename, ierr, iopen)
    !if (iopen .ne. 0) stop 'error opening wout file'
    !if (ierr .ne. 0) stop 'error reading wout file'
    !if (verbose) print *,"  Successfully read VMEC data from ",trim(vmec_filename)
!
    if (verbose) print *,"  Number of field periods (nfp):",nfp

    !*********************************************************************
    ! Do some validation.
    !*********************************************************************

    if (Eq_Obj%nx2<1) then
       print *,"Error! Eq_Obj%nx2 must be >= 1. Instead it is",Eq_Obj%nx2
       stop
    end if

    if (Eq_Obj%nx3<1) then
       print *,"Error! Eq_Obj%nx3 must be >= 1. Instead it is",Eq_Obj%nx3
       stop
    end if

!    ! There is a bug in libstell read_wout_file for ASCII-format wout files, in which the xm_nyq and xn_nyq arrays are sometimes
!    ! not populated. The next few lines here provide a workaround:
!    if (maxval(abs(xm_nyq)) < 1 .and. maxval(abs(xn_nyq)) < 1) then
!       if (mnmax_nyq == mnmax) then
!          if (verbose) print *,"xm_nyq and xn_nyq arrays are not populated in the wout file. Using xm and xn instead."
!          xm_nyq = xm
!          xn_nyq = xn
!       else
!          print *,"Error! xm_nyq and xn_nyq arrays are not populated in the wout file, and mnmax_nyq != mnmax."
!          stop
!       end if
!    end if


    !*********************************************************************
    ! Set up the coordinate grids.
    !*********************************************************************

    alpha = [( ((j-1)*2*pi) / Eq_Obj%nx2, j=1, Eq_Obj%nx2 )]

!!$    if (number_of_field_periods_to_include > nfp) then
!!$       print *,"Error! number_of_field_periods_to_include > nfp"
!!$       print *,"  number_of_field_periods_to_include =",number_of_field_periods_to_include
!!$       print *,"  nfp =",nfp
!!$       stop
!!$    end if
    number_of_field_periods_to_include_final = number_of_field_periods_to_include
    if (number_of_field_periods_to_include <= 0) then
       number_of_field_periods_to_include_final = Eq_Obj%VMEC_Obj%nfp
       if (verbose) print *,"  Since number_of_field_periods_to_include was <= 0, it is being reset to nfp =",nfp
    end if

    zeta = [( zeta_center + (pi*j*number_of_field_periods_to_include_final)/(Eq_Obj%VMEC_Obj%nfp*Eq_Obj%nx3), j=-Eq_Obj%nx3, Eq_Obj%nx3 )]

    !*********************************************************************
    ! We know theta_pest = alpha + iota * zeta, but we need to determine
    ! theta_vmec = theta_pest - Lambda.
    !*********************************************************************

    allocate(theta_vmec(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))

    if (verbose) print *,"  Beginning root solves to determine theta_vmec."
    root_solve_absolute_tolerance = 1.0d-10
    root_solve_relative_tolerance = 1.0d-10
    do iz = -Eq_Obj%nx3, Eq_Obj%nx3
       zeta0 = zeta(iz)
       do ia = 1,Eq_Obj%nx2
          theta_pest_target = alpha(ia) + iota * zeta0
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
       do j = 1, Eq_Obj%nx2
          print *,theta_vmec(j,:)
       end do
    end if

    !*********************************************************************
    ! Initialize geometry arrays
    !*********************************************************************

    bmag = 0
    gradpar = 0
    gds2 = 0
    gds21 = 0
    gds22 = 0
    gbdrift = 0
    gbdrift0 = 0
    cvdrift = 0
    cvdrift0 = 0
    jac_gist_inv = 0
    d_B_d_par = 0


    allocate(B(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(temp2D(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(sqrt_g(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(R(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(Z(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_B_d_theta_vmec(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_B_d_zeta(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_B_d_s(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_R_d_theta_vmec(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_R_d_zeta(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_R_d_s(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Z_d_theta_vmec(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Z_d_zeta(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Z_d_s(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Lambda_d_theta_vmec(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Lambda_d_zeta(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Lambda_d_s(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_sub_s(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_sub_theta_vmec(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_sub_zeta(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_sup_theta_vmec(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_sup_zeta(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3))

    allocate(d_B_d_s_mnc(ns))
    allocate(d_B_d_s_mns(ns))
    allocate(d_R_d_s_mnc(ns))
    allocate(d_R_d_s_mns(ns))
    allocate(d_Z_d_s_mnc(ns))
    allocate(d_Z_d_s_mns(ns))
    allocate(d_Lambda_d_s_mnc(ns))
    allocate(d_Lambda_d_s_mns(ns))

    allocate(d_X_d_s(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_X_d_theta_vmec(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_X_d_zeta(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Y_d_s(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Y_d_theta_vmec(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(d_Y_d_zeta(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))

    allocate(grad_s_X(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_s_Y(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_s_Z(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_theta_vmec_X(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_theta_vmec_Y(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_theta_vmec_Z(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_zeta_X(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_zeta_Y(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_zeta_Z(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_psi_X(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_psi_Y(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_psi_Z(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_alpha_X(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_alpha_Y(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_alpha_Z(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    
    allocate(B_X(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_Y(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_Z(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_B_X(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_B_Y(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(grad_B_Z(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_cross_grad_B_dot_grad_alpha(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_cross_grad_B_dot_grad_alpha_alternate(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_cross_grad_s_dot_grad_alpha(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
    allocate(B_cross_grad_s_dot_grad_alpha_alternate(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))

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
    
    do imn_nyq = 1, mnmax_nyq ! All the quantities we need except R, Z, and Lambda use the _nyq mode numbers.
       m = xm_nyq(imn_nyq)
       n = xn_nyq(imn_nyq)/nfp

       if (abs(m) >= mpol .or. abs(n) > ntor) then
          non_Nyquist_mode_available = .false.
       else
          non_Nyquist_mode_available = .true.
          ! Find the imn in the non-Nyquist arrays that corresponds to the same m and n.
          found_imn = .false.
          do imn = 1,mnmax
             if (xm(imn)==m .and. xn(imn)==n*nfp) then
                found_imn = .true.
                exit
             end if
          end do
          if ((xm(imn) .ne. m) .or. (xn(imn) .ne. n*nfp)) stop "Something went wrong!"
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

       d_B_d_s_mnc(2:ns-1) = (bmnc(imn_nyq,3:ns) - bmnc(imn_nyq,2:ns-1)) / ds
       ! Simplistic extrapolation at the endpoints:
       d_B_d_s_mnc(1) = d_B_d_s_mnc(2)
       d_B_d_s_mnc(ns) = d_B_d_s_mnc(ns-1)

       if (non_Nyquist_mode_available) then
          ! R is on the full mesh:
          d_R_d_s_mnc(2:ns) = (rmnc(imn,2:ns) - rmnc(imn,1:ns-1)) / ds
          d_R_d_s_mnc(1) = 0

          ! Z is on the full mesh:
          d_Z_d_s_mns(2:ns) = (zmns(imn,2:ns) - zmns(imn,1:ns-1)) / ds
          d_Z_d_s_mns(1) = 0

          ! Lambda is on the half mesh:
          d_Lambda_d_s_mns(2:ns-1) = (lmns(imn,3:ns) - lmns(imn,2:ns-1)) / ds
          ! Simplistic extrapolation at the endpoints:
          d_Lambda_d_s_mns(1) = d_Lambda_d_s_mns(2)
          d_Lambda_d_s_mns(ns) = d_Lambda_d_s_mns(ns-1)
       else
          d_R_d_s_mnc = 0
          d_Z_d_s_mns = 0
          d_Lambda_d_s_mns = 0
       end if

       ! End of evaluating radial derivatives.

       do ia = 1,Eq_Obj%nx2
          do iz = -Eq_Obj%nx3, Eq_Obj%nx3
             angle = m * theta_vmec(ia,iz) - n * nfp * zeta(iz)
             cos_angle = cos(angle)
             sin_angle = sin(angle)

             do isurf = 1,2

                ! Handle |B|:
                temp = bmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B(ia,iz) = B(ia,iz) + temp * cos_angle
                d_B_d_theta_vmec(ia,iz) = d_B_d_theta_vmec(ia,iz) - m * temp * sin_angle
                d_B_d_zeta(ia,iz) = d_B_d_zeta(ia,iz) + n * nfp * temp * sin_angle       

                ! Handle Jacobian:
                temp = gmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                sqrt_g(ia,iz) = sqrt_g(ia,iz) + temp * cos_angle

                ! Handle B sup theta:
                temp = bsupumnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sup_theta_vmec(ia,iz) = B_sup_theta_vmec(ia,iz) + temp * cos_angle

                ! Handle B sup zeta:
                temp = bsupvmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sup_zeta(ia,iz) = B_sup_zeta(ia,iz) + temp * cos_angle

                ! Handle B sub theta:
                temp = bsubumnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sub_theta_vmec(ia,iz) = B_sub_theta_vmec(ia,iz) + temp * cos_angle

                ! Handle B sub zeta:
                temp = bsubvmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                temp = temp*scale_factor
                B_sub_zeta(ia,iz) = B_sub_zeta(ia,iz) + temp * cos_angle

                ! Handle B sub psi.
                ! Unlike the other components of B, this one is on the full mesh.
                temp = bsubsmns(imn_nyq,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
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
                   temp = rmnc(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   R(ia,iz) = R(ia,iz) + temp * cos_angle
                   d_R_d_theta_vmec(ia,iz) = d_R_d_theta_vmec(ia,iz) - temp * m * sin_angle
                   d_R_d_zeta(ia,iz)  = d_R_d_zeta(ia,iz)  + temp * n * nfp * sin_angle

                   ! Handle Z, which is on the full mesh
                   temp = zmns(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                   temp = temp*scale_factor
                   Z(ia,iz) = Z(ia,iz) + temp * sin_angle  ! We don't actually need Z itself, only derivatives of Z.
                   d_Z_d_theta_vmec(ia,iz) = d_Z_d_theta_vmec(ia,iz) + temp * m * cos_angle
                   d_Z_d_zeta(ia,iz)  = d_Z_d_zeta(ia,iz)  - temp * n * nfp * cos_angle

                   ! Handle Lambda:
                   temp = lmns(imn,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   ! We don't need Lambda itself, just its derivatives.
                   d_Lambda_d_theta_vmec(ia,iz) = d_Lambda_d_theta_vmec(ia,iz) + m * temp * cos_angle
                   d_Lambda_d_zeta(ia,iz) = d_Lambda_d_zeta(ia,iz) - n * nfp * temp * cos_angle       

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
             end do
          end do
       end do

       ! -----------------------------------------------------
       ! Now consider the stellarator-asymmetric terms.
       ! -----------------------------------------------------

       if (lasym) then

          ! Evaluate the radial derivatives we will need:

          ! B and Lambda are on the half mesh, so their radial derivatives are on the full mesh.
          ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.

          d_B_d_s_mns(2:ns-1) = (bmns(imn_nyq,3:ns) - bmns(imn_nyq,2:ns-1)) / ds
          ! Simplistic extrapolation at the endpoints:
          d_B_d_s_mns(1) = d_B_d_s_mns(2)
          d_B_d_s_mns(ns) = d_B_d_s_mns(ns-1)

          if (non_Nyquist_mode_available) then
             ! R is on the full mesh:
             d_R_d_s_mns(2:ns) = (rmns(imn,2:ns) - rmns(imn,1:ns-1)) / ds
             d_R_d_s_mns(1) = 0

             ! Z is on the full mesh:
             d_Z_d_s_mnc(2:ns) = (zmnc(imn,2:ns) - zmnc(imn,1:ns-1)) / ds
             d_Z_d_s_mnc(1) = 0

             ! Lambda is on the half mesh:
             d_Lambda_d_s_mnc(2:ns-1) = (lmnc(imn_nyq,3:ns) - lmnc(imn_nyq,2:ns-1)) / ds
             ! Simplistic extrapolation at the endpoints:
             d_Lambda_d_s_mnc(1) = d_Lambda_d_s_mnc(2)
             d_Lambda_d_s_mnc(ns) = d_Lambda_d_s_mnc(ns-1)
          else
             d_R_d_s_mns = 0
             d_Z_d_s_mnc = 0
             d_Lambda_d_s_mnc = 0
          end if

          ! End of evaluating radial derivatives.

          do ia = 1,Eq_Obj%nx2
             do iz = -Eq_Obj%nx3, Eq_Obj%nx3
                angle = m * theta_vmec(ia,iz) - n * nfp * zeta(iz)
                cos_angle = cos(angle)
                sin_angle = sin(angle)

                do isurf = 1,2

                   ! Handle |B|:
                   temp = bmns(imn_nyq,vmec_radial_index_half(1)) * vmec_radial_weight_half(1)
                   temp = temp * scale_factor
                   B(ia,iz) = B(ia,iz) + temp * sin_angle
                   d_B_d_theta_vmec(ia,iz) = d_B_d_theta_vmec(ia,iz) + m * temp * cos_angle
                   d_B_d_zeta(ia,iz) = d_B_d_zeta(ia,iz) - n * nfp * temp * cos_angle

                   ! Handle Jacobian:
                   temp = gmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   sqrt_g(ia,iz) = sqrt_g(ia,iz) + temp * sin_angle

                   ! Handle B sup theta:
                   temp = bsupumns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor	
                   B_sup_theta_vmec(ia,iz) = B_sup_theta_vmec(ia,iz) + temp * sin_angle

                   ! Handle B sup zeta:
                   temp = bsupvmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sup_zeta(ia,iz) = B_sup_zeta(ia,iz) + temp * sin_angle

                   ! Handle B sub theta:
                   temp = bsubumns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sub_theta_vmec(ia,iz) = B_sub_theta_vmec(ia,iz) + temp * sin_angle

                   ! Handle B sub zeta:
                   temp = bsubvmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                   temp = temp*scale_factor
                   B_sub_zeta(ia,iz) = B_sub_zeta(ia,iz) + temp * sin_angle

                   ! Handle B sub psi.
                   ! Unlike the other components of B, this one is on the full mesh.
                   temp = bsubsmnc(imn_nyq,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
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
                      temp = rmns(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                      temp = temp*scale_factor
                      R(ia,iz) = R(ia,iz) + temp * sin_angle
                      d_R_d_theta_vmec(ia,iz) = d_R_d_theta_vmec(ia,iz) + temp * m * cos_angle
                      d_R_d_zeta(ia,iz)  = d_R_d_zeta(ia,iz)  - temp * n * nfp * cos_angle

                      ! Handle Z, which is on the full mesh
                      temp = zmnc(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                      temp = temp*scale_factor
                      Z(ia,iz) = Z(ia,iz) + temp * cos_angle   ! We don't actually need Z itself, only derivatives of Z.
                      d_Z_d_theta_vmec(ia,iz) = d_Z_d_theta_vmec(ia,iz) - temp * m * sin_angle
                      d_Z_d_zeta(ia,iz)  = d_Z_d_zeta(ia,iz)  + temp * n * nfp * sin_angle

                      ! Handle Lambda, which is on the half mesh
                      temp = lmnc(imn,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                      temp = temp*scale_factor
                      ! We don't actually need Lambda itself, only derivatives of Lambda.
                      d_Lambda_d_theta_vmec(ia,iz) = d_Lambda_d_theta_vmec(ia,iz) - temp * m * sin_angle
                      d_Lambda_d_zeta(ia,iz)  = d_Lambda_d_zeta(ia,iz)  + temp * n * nfp * sin_angle

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
                end do
             end do
          end do
       end if
    end do

    !*********************************************************************
    ! Sanity check: If the conversion to theta_pest has been done 
    ! correctly, we should find that 
    ! (B dot grad theta_pest) / (B dot grad zeta) = iota.
    ! Let's verify this:
    !*********************************************************************
    
    if (test) then
      allocate(B_dot_grad_theta_pest_over_B_dot_grad_zeta(Eq_Obj%nx2, -Eq_Obj%nx3:Eq_Obj%nx3))
      ! Compute (B dot grad theta_pest) / (B dot grad zeta):
      B_dot_grad_theta_pest_over_B_dot_grad_zeta = &
        & (B_sup_theta_vmec * (1 + d_Lambda_d_theta_vmec) + B_sup_zeta * d_Lambda_d_zeta) / B_sup_zeta 
      temp2D = iota
      call test_arrays(B_dot_grad_theta_pest_over_B_dot_grad_zeta, temp2D, .false., 0.01, 'iota')
      deallocate(B_dot_grad_theta_pest_over_B_dot_grad_zeta)
    end if

    !*********************************************************************
    ! Using R(theta,zeta) and Z(theta,zeta), compute the Cartesian
    ! components of the gradient basis vectors using the dual relations:
    !*********************************************************************

    sqrt_s = sqrt(normalized_toroidal_flux_used)

    do iz = -Eq_Obj%nx3, Eq_Obj%nx3
       cos_angle = cos(zeta(iz))
       sin_angle = sin(zeta(iz))

       ! X = R * cos(zeta)
       d_X_d_theta_vmec(:,iz) = d_R_d_theta_vmec(:,iz) * cos_angle
       d_X_d_zeta(:,iz) = d_R_d_zeta(:,iz) * cos_angle - R(:,iz) * sin_angle
       d_X_d_s(:,iz) = d_R_d_s(:,iz) * cos_angle

       ! Y = R * sin(zeta)
       d_Y_d_theta_vmec(:,iz) = d_R_d_theta_vmec(:,iz) * sin_angle
       d_Y_d_zeta(:,iz) = d_R_d_zeta(:,iz) * sin_angle + R(:,iz) * cos_angle
       d_Y_d_s(:,iz) = d_R_d_s(:,iz) * sin_angle

    end do

    ! Use the dual relations to get the Cartesian components of grad s, grad theta_vmec, and grad zeta:
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
      ! Sanity check: grad_zeta_X should be -sin(zeta) / R:
      do iz = -Eq_Obj%nx3,Eq_Obj%nx3
         temp2D(:,iz) = -sin(zeta(iz)) / R(:,iz)
      end do
      call test_arrays(grad_zeta_X, temp2D, .false., 1.0e-2, 'grad_zeta_X')
      grad_zeta_X = temp2D ! We might as well use the exact value, which is in temp2D.

      ! Sanity check: grad_zeta_Y should be cos(zeta) / R:
      do iz = -Eq_Obj%nx3,Eq_Obj%nx3
         temp2D(:,iz) = cos(zeta(iz)) / R(:,iz)
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

    ! Form grad alpha = grad (theta_vmec + Lambda - iota * zeta)
    do iz = -Eq_Obj%nx3,Eq_Obj%nx3
       grad_alpha_X(:,iz) = (d_Lambda_d_s(:,iz) - zeta(iz) * d_iota_d_s) * grad_s_X(:,iz)
       grad_alpha_Y(:,iz) = (d_Lambda_d_s(:,iz) - zeta(iz) * d_iota_d_s) * grad_s_Y(:,iz)
       grad_alpha_Z(:,iz) = (d_Lambda_d_s(:,iz) - zeta(iz) * d_iota_d_s) * grad_s_Z(:,iz)
    end do
    grad_alpha_X = grad_alpha_X + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_X + (-iota + d_Lambda_d_zeta) * grad_zeta_X
    grad_alpha_Y = grad_alpha_Y + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Y + (-iota + d_Lambda_d_zeta) * grad_zeta_Y
    grad_alpha_Z = grad_alpha_Z + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Z + (-iota + d_Lambda_d_zeta) * grad_zeta_Z

    grad_B_X = d_B_d_s * grad_s_X + d_B_d_theta_vmec * grad_theta_vmec_X + d_B_d_zeta * grad_zeta_X
    grad_B_Y = d_B_d_s * grad_s_Y + d_B_d_theta_vmec * grad_theta_vmec_Y + d_B_d_zeta * grad_zeta_Y
    grad_B_Z = d_B_d_s * grad_s_Z + d_B_d_theta_vmec * grad_theta_vmec_Z + d_B_d_zeta * grad_zeta_Z

    B_X = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_X_d_zeta + (iota - d_Lambda_d_zeta) * d_X_d_theta_vmec) / sqrt_g
    B_Y = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_Y_d_zeta + (iota - d_Lambda_d_zeta) * d_Y_d_theta_vmec) / sqrt_g
    B_Z = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_Z_d_zeta + (iota - d_Lambda_d_zeta) * d_Z_d_theta_vmec) / sqrt_g

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

    jac_gist_inv = (L_reference**3)/(2*safety_factor_q)*(1 + d_Lambda_d_theta_vmec)/sqrt_g

    d_B_d_par = -safety_factor_q/L_reference *jac_gist_inv / B * d_B_d_zeta 

        
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
         - B_sub_theta_vmec * (d_Lambda_d_zeta - iota) ) / sqrt_g

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

    do iz = -Eq_Obj%nx3,Eq_Obj%nx3
       B_cross_grad_B_dot_grad_alpha(:,iz) = 0 &
            + (B_sub_s(:,iz) * d_B_d_theta_vmec(:,iz) * (d_Lambda_d_zeta(:,iz) - iota) &
            + B_sub_theta_vmec(:,iz) * d_B_d_zeta(:,iz) * (d_Lambda_d_s(:,iz) - zeta(iz) * d_iota_d_s) &
            + B_sub_zeta(:,iz) * d_B_d_s(:,iz) * (1 + d_Lambda_d_theta_vmec(:,iz)) &
            - B_sub_zeta(:,iz) * d_B_d_theta_vmec(:,iz) * (d_Lambda_d_s(:,iz) - zeta(iz) * d_iota_d_s) &
            - B_sub_theta_vmec(:,iz) * d_B_d_s(:,iz) * (d_Lambda_d_zeta(:,iz) - iota) &
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

    !*********************************************************************
    ! Finally, assemble the quantities needed for gs2.
    !*********************************************************************

    ! See the latex note gs2_full_surface_stellarator_geometry in the "doc" directory for a derivation of the formulae that follow.

    bmag = B / B_reference

    gradpar = L_reference * B_sup_zeta / B

    gds2 = (grad_alpha_X * grad_alpha_X + grad_alpha_Y * grad_alpha_Y + grad_alpha_Z * grad_alpha_Z) &
         * L_reference * L_reference * normalized_toroidal_flux_used

    gds21 = (grad_alpha_X * grad_psi_X + grad_alpha_Y * grad_psi_Y + grad_alpha_Z * grad_psi_Z) &
         * sign_toroidal_flux * shat / B_reference

    gds22 = (grad_psi_X * grad_psi_X + grad_psi_Y * grad_psi_Y + grad_psi_Z * grad_psi_Z) &
         * shat * shat / (L_reference * L_reference * B_reference * B_reference * normalized_toroidal_flux_used)

    gbdrift = sign_toroidal_flux * 2 * B_reference * L_reference * L_reference * sqrt_s * B_cross_grad_B_dot_grad_alpha &
         / (B * B * B)

    gbdrift0 = abs(edge_toroidal_flux_over_2pi) &
      & * (B_sub_theta_vmec * d_B_d_zeta - B_sub_zeta * d_B_d_theta_vmec) / sqrt_g &
      & * 2 * shat / (B * B * B * sqrt_s)
    ! In the above 2-line expression for gbdrift0, the first line is \vec{B} \times \nabla B \cdot \nabla \psi.

    cvdrift = gbdrift + sign_toroidal_flux * 2 * B_reference * L_reference * L_reference * sqrt_s * mu_0 * d_pressure_d_s &
         * B_cross_grad_s_dot_grad_alpha / (B * B * B * B)

    cvdrift0 = gbdrift0

    !*********************************************************************
    ! Copy the surface quantities
    !*********************************************************************

    Rsurf = R
    Zsurf = Z

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

    if (verbose) print *,"Leaving compute_sfl_quantities."


  contains

    subroutine test_arrays(array1, array2, should_be_0, tolerance, name)
      ! This subroutine is used for verifying the geometry arrays.
      ! When should_be_0 = .true., the subroutine verifies that |array1| = 0 to within 
      !     an absolute tolerance specified by 'tolerance'. array2 is ignored in this case.
      ! When should_be_0 = .false., the subroutine verifies that array1 = array2 
      !     to within a relative tolerance specified by 'tolerance'.

      implicit none

      real(dp), dimension(Eq_Obj%nx2,-Eq_Obj%nx3:Eq_Obj%nx3) :: array1, array2
      real(dp) :: tolerance
      character(len=*) :: name
      logical :: should_be_0
      real(dp) :: max_value, max_difference

      if (should_be_0) then
         max_value = maxval(abs(array1))
         if (verbose) print *,"  maxval(abs(",trim(name),")):",max_value,"(should be << 1.)"
         if (max_value > tolerance) then
            print *,"Error! ",trim(name)," should be 0, but instead it is:"
            do ia = 1,Eq_Obj%nx2
               print *,array1(ia,:)
            end do
            stop
         end if
      else
         max_difference = maxval(abs(array1 - array2)) / maxval(abs(array1) + abs(array2))
         if (verbose) print *,"  Relative difference between two methods for computing ",trim(name),":",max_difference,"(should be << 1.)"
         if (max_difference > tolerance) then
            print *,"Error! Two methods for computing ",trim(name)," disagree. Here comes method 1:"
            do ia = 1,Eq_Obj%nx2
               print *,array1(ia,:)
            end do
            print *,"Here comes method 2:"
            do ia = 1,Eq_Obj%nx2
               print *,array2(ia,:)
            end do
            print *,"Here comes the difference:"
            do ia = 1,Eq_Obj%nx2
               print *,array1(ia,:) - array2(ia,:)
            end do
            stop
         end if
      end if

    end subroutine test_arrays

  end subroutine compute_sfl_quantities

  ! --------------------------------------------------------------------------



  function fzero_residual(theta_vmec_try)

    use read_wout_mod, only: xm, xn, mnmax, lmns, lmnc, lasym
    ! Note that lmns and lmnc use the non-Nyquist xm, xn, and mnmax.
    ! Also note that lmns and lmnc are on the HALF grid.

    implicit none
    
    real(dp) :: theta_vmec_try, fzero_residual
    real(dp) :: angle, sinangle, cosangle
    integer :: imn, which_surface

    ! residual = (theta_pest based on theta_vmec_try) - theta_pest_target = theta_vmec_try + Lambda - theta_pest_target
    fzero_residual = theta_vmec_try - theta_pest_target

    do imn = 1, mnmax
       angle = xm(imn)*theta_vmec_try - xn(imn)*zeta0
       sinangle = sin(angle)
       cosangle = cos(angle)
       do which_surface = 1,2
          fzero_residual = fzero_residual + vmec_radial_weight_half(which_surface) * lmns(imn,vmec_radial_index_half(which_surface)) * sinangle
          if (lasym) then
             fzero_residual = fzero_residual + vmec_radial_weight_half(which_surface) * lmnc(imn,vmec_radial_index_half(which_surface)) * cosangle
          end if
       end do
    end do

  end function fzero_residual

end module vmec2sfl_mod
