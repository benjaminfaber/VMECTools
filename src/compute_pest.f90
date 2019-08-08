! This code is based off of the 
! full_surface_vmec_to_gs2.f90 code
! Written by Matt Landreman, University of Maryland
! Initial code written August 2017.


!******************************************************************************
! This module contains the actual code to compute the VMEC to PEST
! transformation 
! Modified by: B.J. Faber University of Wisconsin-Madison (bfaber@wisc.edu)
! Modification date: June 2019
!******************************************************************************
module compute_pest 

  use types, only: dp, pi, mu_0
  use vmec_object, only: VMEC_Obj
  use pest_object, only: PEST_Obj
  use io_core, only: verbose, test
  implicit none

  public :: compute_pest_geometry
  
  private


contains

  subroutine compute_pest_geometry(pest,x3_center,&
    & number_of_field_periods_to_include,vmec_surface_option)
    implicit none
    !*********************************************************************
    ! Input parameters
    !*********************************************************************

    ! If vmec_surface_option = 0, the magnetic surface specified by desired_normalized_toroidal_flux will be used,
    ! by intedpolating between the surfaces available in the vmec file.
    ! If vmec_surface_option = 1, the magnetic surface on vmec's HALF radial mesh will be used that is closest to desired_normalized_toroidal_flux.
    ! If vmec_surface_option = 2, the magnetic surface on vmec's FULL radial mesh will be used that is closest to desired_normalized_toroidal_flux.    
    ! Other values of vmec_surface_option will cause the program to abort with an error.
    integer, intent(in) :: vmec_surface_option

    ! The pest%x3 domain is centered at x3_center. Setting x3_center = 2*pi*N/nfp for any integer N should
    ! yield identical results to setting x3_center = 0, where nfp is the number of field periods (as in VMEC).
    real(dp), intent(in) :: x3_center

    ! If number_of_field_periods_to_include is > 0, then this parameter does what you think:
    ! the extent of the toroidal in pest%x3 will be 2*pi*number_of_field_periods_to_include/nfp.
    ! If number_of_field_periods_to_include is <= 0, the entire 2*pi toroidal domain will be included.
    real(dp), intent(in) :: number_of_field_periods_to_include

    type(PEST_Obj), intent(inout) :: pest



    !*********************************************************************
    ! Variables used internally by this subroutine
    !*********************************************************************

    integer :: j, idx, idx1, idx3, iztemp, idx2, which_surface, isurf, m, n, imn, imn_nyq, ns
    real(dp) :: angle, sin_angle, cos_angle, temp 
    integer :: fzero_flag, sign_toroidal_flux
    real(dp) :: number_of_field_periods_to_include_final
    real(dp) :: ds, d_pressure_d_s, d_iota_d_s, scale_factor, dphi, min_dr2
    real(dp) :: theta_vmec_min, theta_vmec_max, sqrt_s
    real(dp) :: root_solve_absolute_tolerance, root_solve_relative_tolerance
    real(dp) :: theta_pest_target, zeta0, edge_toroidal_flux_over_2pi
    logical :: non_Nyquist_mode_available, found_imn

    real(dp), dimension(2) :: vmec_radial_weight_full, vmec_radial_weight_half
    integer, dimension(2) :: vmec_radial_index_full, vmec_radial_index_half
    real(dp), dimension(:), allocatable :: dr2, normalized_toroidal_flux_full_grid, normalized_toroidal_flux_half_grid
    real(dp), dimension(:), allocatable :: d_pressure_d_s_on_half_grid, d_iota_d_s_on_half_grid


    real(dp), dimension(:,:), allocatable :: theta_vmec
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
    real(dp), dimension(:,:), allocatable :: grad_theta_X, grad_theta_Y, grad_theta_Z
    real(dp), dimension(:,:), allocatable :: grad_psi_X, grad_psi_Y, grad_psi_Z
    real(dp), dimension(:,:), allocatable :: grad_alpha_X, grad_alpha_Y, grad_alpha_Z
    real(dp), dimension(:,:), allocatable :: B_cross_grad_B_dot_grad_psi, B_cross_grad_B_dot_grad_alpha, B_cross_grad_B_dot_grad_alpha_alternate
    real(dp), dimension(:,:), allocatable :: B_cross_grad_psi_dot_grad_alpha, B_cross_grad_s_dot_grad_alpha, B_cross_grad_s_dot_grad_alpha_alternate
    real(dp), dimension(:,:), allocatable :: grad_B_X, grad_B_Y, grad_B_Z
    real(dp), dimension(:,:), allocatable :: B_X, B_Y, B_Z

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
    ns = pest%vmec%ns
    allocate(normalized_toroidal_flux_full_grid(ns))
    allocate(normalized_toroidal_flux_half_grid(ns-1))
    allocate(d_iota_d_s_on_half_grid(ns))
    allocate(d_pressure_d_s_on_half_grid(ns))

    allocate(theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(temp2D(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(sqrt_g(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(R(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_B_d_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_B_d_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_B_d_s(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_R_d_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_R_d_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_R_d_s(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Z_d_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Z_d_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Z_d_s(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Lambda_d_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Lambda_d_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Lambda_d_s(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_sub_s(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_sub_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_sub_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_sup_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_sup_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))

    allocate(d_B_d_s_mnc(ns))
    allocate(d_B_d_s_mns(ns))
    allocate(d_R_d_s_mnc(ns))
    allocate(d_R_d_s_mns(ns))
    allocate(d_Z_d_s_mnc(ns))
    allocate(d_Z_d_s_mns(ns))
    allocate(d_Lambda_d_s_mnc(ns))
    allocate(d_Lambda_d_s_mns(ns))

    allocate(d_X_d_s(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_X_d_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_X_d_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Y_d_s(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Y_d_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(d_Y_d_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))

    allocate(grad_s_X(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_s_Y(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_s_Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_theta_vmec_X(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_theta_vmec_Y(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_theta_vmec_Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_zeta_X(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_zeta_Y(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_zeta_Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_theta_X(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_theta_Y(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_theta_Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_psi_X(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_psi_Y(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_psi_Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_alpha_X(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_alpha_Y(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_alpha_Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    
    allocate(B_X(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_Y(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_B_X(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_B_Y(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(grad_B_Z(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_cross_grad_B_dot_grad_psi(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_cross_grad_B_dot_grad_alpha(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_cross_grad_B_dot_grad_alpha_alternate(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_cross_grad_psi_dot_grad_alpha(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_cross_grad_s_dot_grad_alpha(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
    allocate(B_cross_grad_s_dot_grad_alpha_alternate(pest%ix21:pest%ix22, pest%ix31:pest%ix32))


    if (verbose) print *,"Entering subroutine compute_pest_geometry."
    !*********************************************************************
    ! Do some validation.
    !*********************************************************************
    do idx1 = pest%ix11,pest%ix12 
      if (pest%x1(idx1) <= 0) then
         print *,"Error! desired_normalized_toroidal_flux must be >0. Instead it is", pest%x1(idx1)
         stop
      end if

      if (pest%x1(idx1) > 1) then
         print *,"Error! desired_normalized_toroidal_flux must be <= 1. Instead it is",pest%x1(idx1)
         stop
      end if
    end do

    if (pest%nx2<1) then
       print *,"Error! pest%nx2 must be >= 1. Instead it is",pest%nx2
       stop
    end if

    if (pest%nx3<1) then
       print *,"Error! pest%nx3 must be >= 1. Instead it is",pest%nx3
       stop
    end if

    !*********************************************************************
    ! Set up the coordinate grids.
    !*********************************************************************

    pest%x2 = [( (j*2*pi) / pest%nx2, j=pest%ix21, pest%ix22 )]


    number_of_field_periods_to_include_final = number_of_field_periods_to_include
    if (number_of_field_periods_to_include <= 0) then
       number_of_field_periods_to_include_final = pest%vmec%nfp
       if (verbose) print *,"  Since number_of_field_periods_to_include was <= 0, it is being reset to nfp =",pest%vmec%nfp
    end if





    edge_toroidal_flux_over_2pi = pest%vmec%phi(ns) / (2*pi) * pest%vmec%isigng ! isigns is called signgs in the wout*.nc file. Why is this signgs here?

    ! this gives the sign of the edge toroidal flux
    sign_toroidal_flux = int(sign(1.1,edge_toroidal_flux_over_2pi))
    !write (*,*) 'sign_toroidal_flux', sign_toroidal_flux


    normalized_toroidal_flux_full_grid = [( real(j-1)/(ns-1), j=1,ns )]

    ! Build an array of the half grid points:
    do j = 1,ns-1
       normalized_toroidal_flux_half_grid(j) = (normalized_toroidal_flux_full_grid(j) + normalized_toroidal_flux_full_grid(j+1))*(0.5d+0)
    end do

    d_iota_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    if (verbose) print *,"  ds =",ds
    d_iota_d_s_on_half_grid(2:ns) = (pest%vmec%iotaf(2:ns) - pest%vmec%iotaf(1:ns-1)) / ds

    d_pressure_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    d_pressure_d_s_on_half_grid(2:ns) = (pest%vmec%presf(2:ns) - pest%vmec%presf(1:ns-1)) / ds

    !*********************************************************************
    ! Determine which flux surface to use, based on 
    ! desired_normalized_toroidal_flux and vmec_surface_option.
    !*********************************************************************

    ! Possible values of vmec_surface_option:
    ! 0 = Use the exact radius requested.
    ! 1 = Use the nearest value of the VMEC half grid.
    ! 2 = Use the nearest value of the VMEC full grid.

    do idx1 = pest%ix11,pest%ix12
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


      select case (vmec_surface_option)
      case (0)
         ! Use exact radius requested.
         pest%x1(idx1) = pest%x1(idx1)

      case (1)
         ! Use nearest value of the VMEC half grid

         ! Compute differences
         allocate(dr2(ns-1))
         dr2 = (normalized_toroidal_flux_half_grid - pest%x1(idx1)) ** 2

         idx = 1
         min_dr2 = dr2(1)
         ! Find the index of minimum error:
         do j=2,ns-1
            if (dr2(j)<min_dr2) then
               idx = j
               min_dr2 = dr2(j)
            end if
         end do

         pest%x1(idx1) = normalized_toroidal_flux_half_grid(idx)
         deallocate(dr2)

      case (2)
         ! Use nearest value of the VMEC full grid

         ! Compute differences
         allocate(dr2(ns))
         dr2 = (normalized_toroidal_flux_full_grid - pest%x1(idx1)) ** 2

         idx = 1
         min_dr2 = dr2(1)
         ! Find the index of minimum error:
         do j=2,ns
            if (dr2(j)<min_dr2) then
               idx = j
               min_dr2 = dr2(j)
            end if
         end do

         pest%x1(idx1) = normalized_toroidal_flux_full_grid(idx)
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
      if (pest%x1(idx1)>1) then
         stop "Error! pest%x1(idx1) cannot be >1"
      elseif (pest%x1(idx1)<0) then
         stop "Error! pest%x1(idx1) cannot be <0"
      elseif (pest%x1(idx1)==1) then
         vmec_radial_index_full(1) = ns-1
         vmec_radial_index_full(2) = ns
         vmec_radial_weight_full(1) = 0.0d0
      else
         ! pest%x1(pest%ix11) is >= 0 and <1
         ! This is the most common case.
         vmec_radial_index_full(1) = floor(pest%x1(idx1)*(ns-1))+1
         vmec_radial_index_full(2) = vmec_radial_index_full(1) + 1
         vmec_radial_weight_full(1) = vmec_radial_index_full(1) - pest%x1(idx1)*(ns-1.0d0)
      end if
      vmec_radial_weight_full(2) = 1.0d0 - vmec_radial_weight_full(1)

      ! Handle quantities for the half grid
      if (pest%x1(idx1) < normalized_toroidal_flux_half_grid(1)) then
         print *,"Warning: extrapolating beyond the end of VMEC's half grid."
         print *,"(Extrapolating towards the magnetic axis.) Results are likely to be inaccurate."

         ! We start at element 2 since element 1 is always 0 for quantities on the half grid.
         vmec_radial_index_half(1) = 2
         vmec_radial_index_half(2) = 3
         vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(2) - pest%x1(idx1)) / (normalized_toroidal_flux_half_grid(2) - normalized_toroidal_flux_half_grid(1))

      elseif (pest%x1(idx1) > normalized_toroidal_flux_half_grid(ns-1)) then
         print *,"Warning: extrapolating beyond the end of VMEC's half grid."
         print *,"(Extrapolating towards the last closed flux surface.) Results may be inaccurate."
         vmec_radial_index_half(1) = ns-1
         vmec_radial_index_half(2) = ns
         vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(ns-1) - pest%x1(idx1)) &
              / (normalized_toroidal_flux_half_grid(ns-1) - normalized_toroidal_flux_half_grid(ns-2))

      elseif (pest%x1(idx1) == normalized_toroidal_flux_half_grid(ns-1)) then
         ! We are exactly at the last point of the half grid
         vmec_radial_index_half(1) = ns-1
         vmec_radial_index_half(2) = ns
         vmec_radial_weight_half(1) = 0.0d0
      else
         ! pest%x1(pest%ix11) is inside the half grid.
         ! This is the most common case.
         vmec_radial_index_half(1) = floor(pest%x1(idx1)*(ns-1) + 0.5d+0)+1
         if (vmec_radial_index_half(1) < 2) then
            ! This can occur sometimes due to roundoff error.
            vmec_radial_index_half(1) = 2
         end if
         vmec_radial_index_half(2) = vmec_radial_index_half(1) + 1
         vmec_radial_weight_half(1) = vmec_radial_index_half(1) - pest%x1(idx1)*(ns-1.0d0) - (0.5d+0)
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

      pest%iota(idx1) = pest%vmec%iotas(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
           + pest%vmec%iotas(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
      if (verbose) print *,"  iota =",pest%iota(idx1)
      pest%safety_factor_q(idx1) = 1/pest%iota(idx1)

      d_iota_d_s =  &
           d_iota_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
           + d_iota_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
      if (verbose) print *,"  d iota / d s =",d_iota_d_s
      ! shat = (r/q)(dq/dr) where r = a sqrt(s).
      !      = - (r/pest%iota(idx1)) (d pest%iota(idx1) / d r) = -2 (s/pest%iota(idx1)) (d pest%iota(idx1) / d s)
      pest%shat(idx1) = (-2 * pest%x1(idx1) / pest%iota(idx1)) * d_iota_d_s

      d_pressure_d_s =  &
           d_pressure_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
           + d_pressure_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
      if (verbose) print *,"  d pressure / d s =",d_pressure_d_s


      if (verbose) print *,"  Number of field periods (nfp):",pest%vmec%nfp

      !*********************************************************************
      ! Do some validation.
      !*********************************************************************

      !*********************************************************************
      ! We know theta_pest = alpha + pest%iota(idx1) * pest%x3, but we need to determine
      ! theta_vmec = theta_pest - Lambda.
      !*********************************************************************
      select case(trim(pest%x3_coord))
        case('zeta')
          pest%x3(:,idx1) = [( x3_center + 2.0*(pi*j*number_of_field_periods_to_include_final)/(pest%vmec%nfp*(pest%nx3-1)), j=pest%ix31,pest%ix32 )]
        case('theta')
          pest%x3(:,idx1) = [( pest%safety_factor_q(idx1)*(x3_center + 2.0*(pi*j*number_of_field_periods_to_include_final)/(pest%vmec%nfp*(pest%nx3-1))), j=pest%ix31,pest%ix32 )]
        case default
          pest%x3(:,idx1) = [( x3_center + 2.0*(pi*j*number_of_field_periods_to_include_final)/(pest%vmec%nfp*(pest%nx3-1)), j=pest%ix31,pest%ix32 )]
      end select


      if (verbose) print *,"  Beginning root solves to determine theta_vmec."
      root_solve_absolute_tolerance = 1.0d-10
      root_solve_relative_tolerance = 1.0d-10
      do idx3 = pest%ix31,pest%ix32
         zeta0 = pest%x3(idx3,idx1)
         do idx2 = pest%ix21,pest%ix22
            theta_pest_target = pest%x2(idx2) + pest%iota(idx1) * zeta0
            ! Guess that theta_vmec will be within 0.3 radians of theta_pest:
            theta_vmec_min = theta_pest_target - 0.5
            theta_vmec_max = theta_pest_target + 0.5

            ! In the 4th argument, we are telling the root-finder (fzero) to use theta_pest as the initial guess for theta_vmec.
            call fzero(fzero_residual, theta_vmec_min, theta_vmec_max, theta_pest_target, &
                 root_solve_relative_tolerance, root_solve_absolute_tolerance, fzero_flag)
            ! Note: fzero returns its answer in theta_vmec_min.
            theta_vmec(idx2,idx3) = theta_vmec_min
            if (fzero_flag == 4) then
               stop "ERROR: fzero returned error 4: no sign change in residual"
            else if (fzero_flag > 2) then
               print *,"WARNING: fzero returned an error code:",fzero_flag
            end if
         end do
      end do
      if (verbose) then
         print *,"  Done with root solves. Here comes theta_vmec:"
         do idx2 = pest%ix21,pest%ix22
            print *,theta_vmec(idx2,:)
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


      !*********************************************************************
      ! Now that we know the grid points in theta_vmec, we can evaluate
      ! all the geometric quantities on the grid points.
      !*********************************************************************
      
      do imn_nyq = 1, pest%vmec%mnmax_nyq ! All the quantities we need except R, Z, and Lambda use the _nyq mode numbers.
         m = pest%vmec%xm_nyq(imn_nyq)
         n = pest%vmec%xn_nyq(imn_nyq)/pest%vmec%nfp

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

         do idx3 = pest%ix31,pest%ix32
            do idx2 = pest%ix21,pest%ix22
               angle = m * theta_vmec(idx2,idx3) - n * pest%vmec%nfp * pest%x3(idx3,idx1)
               cos_angle = cos(angle)
               sin_angle = sin(angle)

               do isurf = 1,2

                  ! Handle |B|:
                  temp = pest%vmec%bmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                  temp = temp*scale_factor
                  B(idx2,idx3) = B(idx2,idx3) + temp * cos_angle
                  d_B_d_theta_vmec(idx2,idx3) = d_B_d_theta_vmec(idx2,idx3) - m * temp * sin_angle
                  d_B_d_zeta(idx2,idx3) = d_B_d_zeta(idx2,idx3) + n * pest%vmec%nfp * temp * sin_angle       

                  ! Handle Jacobian:
                  temp = pest%vmec%gmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                  temp = temp*scale_factor
                  sqrt_g(idx2,idx3) = sqrt_g(idx2,idx3) + temp * cos_angle

                  ! Handle B sup theta:
                  temp = pest%vmec%bsupumnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                  temp = temp*scale_factor
                  B_sup_theta_vmec(idx2,idx3) = B_sup_theta_vmec(idx2,idx3) + temp * cos_angle

                  ! Handle B sup pest%x3:
                  temp = pest%vmec%bsupvmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                  temp = temp*scale_factor
                  B_sup_zeta(idx2,idx3) = B_sup_zeta(idx2,idx3) + temp * cos_angle

                  ! Handle B sub theta:
                  temp = pest%vmec%bsubumnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                  temp = temp*scale_factor
                  B_sub_theta_vmec(idx2,idx3) = B_sub_theta_vmec(idx2,idx3) + temp * cos_angle

                  ! Handle B sub pest%x3:
                  temp = pest%vmec%bsubvmnc(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                  temp = temp*scale_factor
                  B_sub_zeta(idx2,idx3) = B_sub_zeta(idx2,idx3) + temp * cos_angle

                  ! Handle B sub psi.
                  ! Unlike the other components of B, this one is on the full mesh.
                  temp = pest%vmec%bsubsmns(imn_nyq,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                  temp = temp*scale_factor
                  B_sub_s(idx2,idx3) = B_sub_s(idx2,idx3) + temp * sin_angle

                  ! Handle d B / d s
                  ! Since bmnc is on the half mesh, its radial derivative is on the full mesh.
                  temp = d_B_d_s_mnc(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                  temp = temp*scale_factor
                  d_B_d_s(idx2,idx3) = d_B_d_s(idx2,idx3) + temp * cos_angle

                  ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                  if (non_Nyquist_mode_available) then

                     ! Handle R, which is on the full mesh
                     temp = pest%vmec%rmnc(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                     temp = temp*scale_factor
                     R(idx2,idx3) = R(idx2,idx3) + temp * cos_angle
                     d_R_d_theta_vmec(idx2,idx3) = d_R_d_theta_vmec(idx2,idx3) - temp * m * sin_angle
                     d_R_d_zeta(idx2,idx3)  = d_R_d_zeta(idx2,idx3)  + temp * n * pest%vmec%nfp * sin_angle

                     ! Handle Z, which is on the full mesh
                     temp = pest%vmec%zmns(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                     temp = temp*scale_factor
                     Z(idx2,idx3) = Z(idx2,idx3) + temp * sin_angle  ! We don't actually need Z itself, only derivatives of Z.
                     d_Z_d_theta_vmec(idx2,idx3) = d_Z_d_theta_vmec(idx2,idx3) + temp * m * cos_angle
                     d_Z_d_zeta(idx2,idx3)  = d_Z_d_zeta(idx2,idx3)  - temp * n * pest%vmec%nfp * cos_angle

                     ! Handle Lambda:
                     temp = pest%vmec%lmns(imn,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                     temp = temp*scale_factor
                     ! We don't need Lambda itself, just its derivatives.
                     d_Lambda_d_theta_vmec(idx2,idx3) = d_Lambda_d_theta_vmec(idx2,idx3) + m * temp * cos_angle
                     d_Lambda_d_zeta(idx2,idx3) = d_Lambda_d_zeta(idx2,idx3) - n * pest%vmec%nfp * temp * cos_angle       

                     ! Handle d R / d s
                     ! Since R is on the full mesh, its radial derivative is on the half mesh.
                     temp = d_R_d_s_mnc(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                     temp = temp*scale_factor
                     d_R_d_s(idx2,idx3) = d_R_d_s(idx2,idx3) + temp * cos_angle

                     ! Handle d Z / d s
                     ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                     temp = d_Z_d_s_mns(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                     temp = temp*scale_factor
                     d_Z_d_s(idx2,idx3) = d_Z_d_s(idx2,idx3) + temp * sin_angle

                     ! Handle d Lambda / d s
                     ! Since Lambda is on the half mesh, its radial derivative is on the full mesh.
                     temp = d_Lambda_d_s_mns(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                     temp = temp*scale_factor
                     d_Lambda_d_s(idx2,idx3) = d_Lambda_d_s(idx2,idx3) + temp * sin_angle

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

            do idx3 = pest%ix31,pest%ix32
               do idx2 = pest%ix21,pest%ix22
                  angle = m * theta_vmec(idx2,idx3) - n * pest%vmec%nfp * pest%x3(idx3,idx1)
                  cos_angle = cos(angle)
                  sin_angle = sin(angle)

                  do isurf = 1,2

                     ! Handle |B|:
                     temp = pest%vmec%bmns(imn_nyq,vmec_radial_index_half(1)) * vmec_radial_weight_half(1)
                     temp = temp * scale_factor
                     B(idx2,idx3) = B(idx2,idx3) + temp * sin_angle
                     d_B_d_theta_vmec(idx2,idx3) = d_B_d_theta_vmec(idx2,idx3) + m * temp * cos_angle
                     d_B_d_zeta(idx2,idx3) = d_B_d_zeta(idx2,idx3) - n * pest%vmec%nfp * temp * cos_angle

                     ! Handle Jacobian:
                     temp = pest%vmec%gmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                     temp = temp*scale_factor
                     sqrt_g(idx2,idx3) = sqrt_g(idx2,idx3) + temp * sin_angle

                     ! Handle B sup theta:
                     temp = pest%vmec%bsupumns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                     temp = temp*scale_factor
                     B_sup_theta_vmec(idx2,idx3) = B_sup_theta_vmec(idx2,idx3) + temp * sin_angle

                     ! Handle B sup pest%x3:
                     temp = pest%vmec%bsupvmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                     temp = temp*scale_factor
                     B_sup_zeta(idx2,idx3) = B_sup_zeta(idx2,idx3) + temp * sin_angle

                     ! Handle B sub theta:
                     temp = pest%vmec%bsubumns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                     temp = temp*scale_factor
                     B_sub_theta_vmec(idx2,idx3) = B_sub_theta_vmec(idx2,idx3) + temp * sin_angle

                     ! Handle B sub pest%x3:
                     temp = pest%vmec%bsubvmns(imn_nyq,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                     temp = temp*scale_factor
                     B_sub_zeta(idx2,idx3) = B_sub_zeta(idx2,idx3) + temp * sin_angle

                     ! Handle B sub psi.
                     ! Unlike the other components of B, this one is on the full mesh.
                     temp = pest%vmec%bsubsmnc(imn_nyq,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                     temp = temp*scale_factor
                     B_sub_s(idx2,idx3) = B_sub_s(idx2,idx3) + temp * cos_angle

                     ! Handle d B / d s.
                     ! Since bmns is on the half mesh, its radial derivative is on the full mesh.
                     temp = d_B_d_s_mns(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                     temp = temp*scale_factor
                     d_B_d_s(idx2,idx3) = d_B_d_s(idx2,idx3) + temp * sin_angle

                     ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                     if (non_Nyquist_mode_available) then

                        ! Handle R, which is on the full mesh
                        temp = pest%vmec%rmns(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                        temp = temp*scale_factor
                        R(idx2,idx3) = R(idx2,idx3) + temp * sin_angle
                        d_R_d_theta_vmec(idx2,idx3) = d_R_d_theta_vmec(idx2,idx3) + temp * m * cos_angle
                        d_R_d_zeta(idx2,idx3)  = d_R_d_zeta(idx2,idx3)  - temp * n * pest%vmec%nfp * cos_angle

                        ! Handle Z, which is on the full mesh
                        temp = pest%vmec%zmnc(imn,vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                        temp = temp*scale_factor
                        Z(idx2,idx3) = Z(idx2,idx3) + temp * cos_angle   ! We don't actually need Z itself, only derivatives of Z.
                        d_Z_d_theta_vmec(idx2,idx3) = d_Z_d_theta_vmec(idx2,idx3) - temp * m * sin_angle
                        d_Z_d_zeta(idx2,idx3)  = d_Z_d_zeta(idx2,idx3)  + temp * n * pest%vmec%nfp * sin_angle

                        ! Handle Lambda, which is on the half mesh
                        temp = pest%vmec%lmnc(imn,vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                        temp = temp*scale_factor
                        ! We don't actually need Lambda itself, only derivatives of Lambda.
                        d_Lambda_d_theta_vmec(idx2,idx3) = d_Lambda_d_theta_vmec(idx2,idx3) - temp * m * sin_angle
                        d_Lambda_d_zeta(idx2,idx3)  = d_Lambda_d_zeta(idx2,idx3)  + temp * n * pest%vmec%nfp * sin_angle

                        ! Handle d R / d s.
                        ! Since R is on the full mesh, its radial derivative is on the half mesh.
                        temp = d_R_d_s_mns(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                        temp = temp*scale_factor
                        d_R_d_s(idx2,idx3) = d_R_d_s(idx2,idx3) + temp * sin_angle

                        ! Handle d Z / d s.
                        ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                        temp = d_Z_d_s_mnc(vmec_radial_index_half(isurf)) * vmec_radial_weight_half(isurf)
                        temp = temp*scale_factor
                        d_Z_d_s(idx2,idx3) = d_Z_d_s(idx2,idx3) + temp * cos_angle

                        ! Handle d Lambda / d s.
                        ! Since Lambda is on the half mesh, its radial derivative is on the full mesh.
                        temp = d_Lambda_d_s_mnc(vmec_radial_index_full(isurf)) * vmec_radial_weight_full(isurf)
                        temp = temp*scale_factor
                        d_Lambda_d_s(idx2,idx3) = d_Lambda_d_s(idx2,idx3) + temp * cos_angle
                     end if
                  end do ! End surface loop
               end do ! End x3 loop
            end do ! End x2 loop
         end if ! End lasym
      end do ! End imn_nyq loop 

      !*********************************************************************
      ! Sanity check: If the conversion to theta_pest has been done 
      ! correctly, we should find that 
      ! (B dot grad theta_pest) / (B dot grad pest%x3) = pest%iota(idx1).
      ! Let's verify this:
      !*********************************************************************
      
      if (test) then
        allocate(B_dot_grad_theta_pest_over_B_dot_grad_zeta(pest%ix21:pest%ix22, pest%ix31:pest%ix32))
        ! Compute (B dot grad theta_pest) / (B dot grad pest%x3):
        B_dot_grad_theta_pest_over_B_dot_grad_zeta = &
          & (B_sup_theta_vmec * (1 + d_Lambda_d_theta_vmec) + B_sup_zeta * d_Lambda_d_zeta) / B_sup_zeta 
        temp2D = pest%iota(idx1)
        call test_arrays(B_dot_grad_theta_pest_over_B_dot_grad_zeta, temp2D, .false., 0.01, 'pest%iota(idx1)')
        deallocate(B_dot_grad_theta_pest_over_B_dot_grad_zeta)
      end if

      !*********************************************************************
      ! Using R(theta,zeta) and Z(theta,zeta), compute the Cartesian
      ! components of the gradient basis vectors using the dual relations:
      !*********************************************************************

      sqrt_s = sqrt(pest%x1(pest%ix11))

      do idx3 = pest%ix31,pest%ix32
         cos_angle = cos(pest%x3(idx3,idx1))
         sin_angle = sin(pest%x3(idx3,idx1))

         ! X = R * cos(pest%x3)
         d_X_d_theta_vmec(:,idx3) = d_R_d_theta_vmec(:,idx3) * cos_angle
         d_X_d_zeta(:,idx3) = d_R_d_zeta(:,idx3) * cos_angle - R(:,idx3) * sin_angle
         d_X_d_s(:,idx3) = d_R_d_s(:,idx3) * cos_angle

         ! Y = R * sin(pest%x3)
         d_Y_d_theta_vmec(:,idx3) = d_R_d_theta_vmec(:,idx3) * sin_angle
         d_Y_d_zeta(:,idx3) = d_R_d_zeta(:,idx3) * sin_angle + R(:,idx3) * cos_angle
         d_Y_d_s(:,idx3) = d_R_d_s(:,idx3) * sin_angle

      end do

      ! Use the dual relations to get the Cartesian components of grad s, grad theta_vmec, and grad pest%x3:
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
        ! Sanity check: grad_zeta_X should be -sin(pest%x3) / R:
        do idx3 = pest%ix31,pest%ix32
           temp2D(:,idx3) = -sin(pest%x3(idx3,idx1)) / R(:,idx3)
        end do
        call test_arrays(grad_zeta_X, temp2D, .false., 1.0e-2, 'grad_zeta_X')
        grad_zeta_X = temp2D ! We might as well use the exact value, which is in temp2D.

        ! Sanity check: grad_zeta_Y should be cos(pest%x3) / R:
        do idx3 = pest%ix31,pest%ix32
           temp2D(:,idx3) = cos(pest%x3(idx3,idx1)) / R(:,idx3)
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

      ! Form grad alpha = grad (theta_vmec + Lambda - pest%iota(idx1) * pest%x3)
      do idx3 = pest%ix31,pest%ix32
         grad_alpha_X(:,idx3) = (d_Lambda_d_s(:,idx3) - pest%x3(idx3,idx1) * d_iota_d_s) * grad_s_X(:,idx3)
         grad_alpha_Y(:,idx3) = (d_Lambda_d_s(:,idx3) - pest%x3(idx3,idx1) * d_iota_d_s) * grad_s_Y(:,idx3)
         grad_alpha_Z(:,idx3) = (d_Lambda_d_s(:,idx3) - pest%x3(idx3,idx1) * d_iota_d_s) * grad_s_Z(:,idx3)
      end do
      grad_alpha_X = grad_alpha_X + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_X + (-pest%iota(idx1) + d_Lambda_d_zeta) * grad_zeta_X
      grad_alpha_Y = grad_alpha_Y + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Y + (-pest%iota(idx1) + d_Lambda_d_zeta) * grad_zeta_Y
      grad_alpha_Z = grad_alpha_Z + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Z + (-pest%iota(idx1) + d_Lambda_d_zeta) * grad_zeta_Z

      grad_theta_X = d_Lambda_d_s * grad_s_X + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_X + d_Lambda_d_zeta * grad_zeta_X 
      grad_theta_Y = d_Lambda_d_s * grad_s_Y + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Y + d_Lambda_d_zeta * grad_zeta_Y 
      grad_theta_Z = d_Lambda_d_s * grad_s_Z + (1 + d_Lambda_d_theta_vmec) * grad_theta_vmec_Z + d_Lambda_d_zeta * grad_zeta_Z 

      grad_B_X = d_B_d_s * grad_s_X + d_B_d_theta_vmec * grad_theta_vmec_X + d_B_d_zeta * grad_zeta_X
      grad_B_Y = d_B_d_s * grad_s_Y + d_B_d_theta_vmec * grad_theta_vmec_Y + d_B_d_zeta * grad_zeta_Y
      grad_B_Z = d_B_d_s * grad_s_Z + d_B_d_theta_vmec * grad_theta_vmec_Z + d_B_d_zeta * grad_zeta_Z

      B_X = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_X_d_zeta + (pest%iota(idx1) - d_Lambda_d_zeta) * d_X_d_theta_vmec) / sqrt_g
      B_Y = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_Y_d_zeta + (pest%iota(idx1) - d_Lambda_d_zeta) * d_Y_d_theta_vmec) / sqrt_g
      B_Z = edge_toroidal_flux_over_2pi * ((1 + d_Lambda_d_theta_vmec) * d_Z_d_zeta + (pest%iota(idx1) - d_Lambda_d_zeta) * d_Z_d_theta_vmec) / sqrt_g

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

      B_cross_grad_psi_dot_grad_alpha = 0 + &
        & (B_Y * grad_psi_Z - B_Z * grad_psi_Y) * grad_alpha_X + &         
        & (B_Z * grad_psi_X - B_X * grad_psi_Z) * grad_alpha_Y + &
        & (B_X * grad_psi_Y - B_Y * grad_psi_X) * grad_alpha_Z

      B_cross_grad_B_dot_grad_psi = 0 + &
        & (B_Y * grad_B_Z - B_Z * grad_B_Y) * grad_psi_X + &         
        & (B_Z * grad_B_X - B_X * grad_B_Z) * grad_psi_Y + &
        & (B_X * grad_B_Y - B_Y * grad_B_X) * grad_psi_Z

      B_cross_grad_B_dot_grad_alpha = 0 + &
        & (B_Y * grad_B_Z - B_Z * grad_B_Y) * grad_alpha_X + &         
        & (B_Z * grad_B_X - B_X * grad_B_Z) * grad_alpha_Y + &
        & (B_X * grad_B_Y - B_Y * grad_B_X) * grad_alpha_Z


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
           - B_sub_theta_vmec * (d_Lambda_d_zeta - pest%iota(idx1)) ) / sqrt_g

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

      do idx3 = pest%ix31,pest%ix32
         B_cross_grad_B_dot_grad_alpha(:,idx3) = 0 &
              + (B_sub_s(:,idx3) * d_B_d_theta_vmec(:,idx3) * (d_Lambda_d_zeta(:,idx3) - pest%iota(idx1)) &
              + B_sub_theta_vmec(:,idx3) * d_B_d_zeta(:,idx3) * (d_Lambda_d_s(:,idx3) - pest%x3(idx3,idx1) * d_iota_d_s) &
              + B_sub_zeta(:,idx3) * d_B_d_s(:,idx3) * (1 + d_Lambda_d_theta_vmec(:,idx3)) &
              - B_sub_zeta(:,idx3) * d_B_d_theta_vmec(:,idx3) * (d_Lambda_d_s(:,idx3) - pest%x3(idx3,idx1) * d_iota_d_s) &
              - B_sub_theta_vmec(:,idx3) * d_B_d_s(:,idx3) * (d_Lambda_d_zeta(:,idx3) - pest%iota(idx1)) &
              - B_sub_s(:,idx3) * d_B_d_zeta(:,idx3) * (1 + d_Lambda_d_theta_vmec(:,idx3))) / sqrt_g(:,idx3)
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

      pest%bmag(:,:,idx1) = B 
      pest%jac(:,:,idx1) = sqrt_g
      pest%g11(:,:,idx1) = 4.0/((pest%B_ref**2)*(pest%L_ref**4))*(grad_psi_X * grad_psi_X + grad_psi_Y * grad_psi_Y + grad_psi_Z * grad_psi_Z)
      pest%g12(:,:,idx1) = 2.0*sign_toroidal_flux/(pest%B_ref*(pest%L_ref**2))*(grad_psi_X * grad_alpha_X + grad_psi_Y * grad_alpha_Y + grad_psi_Z * grad_alpha_Z)
      pest%g22(:,:,idx1) = grad_alpha_X * grad_alpha_X + grad_alpha_Y * grad_alpha_Y + grad_alpha_Z * grad_alpha_Z

      if (trim(pest%x3_coord) == 'zeta') then
        pest%g13(:,:,idx1) = 2.0*sign_toroidal_flux/(pest%B_ref*(pest%L_ref**2))*(grad_psi_X * grad_zeta_X + grad_psi_Y * grad_zeta_Y + grad_psi_Z * grad_zeta_Z)
        pest%g23(:,:,idx1) = grad_alpha_X * grad_zeta_X + grad_alpha_Y * grad_zeta_Y + grad_alpha_Z * grad_zeta_Z
        pest%g33(:,:,idx1) = grad_zeta_X * grad_zeta_X + grad_zeta_Y * grad_zeta_Y + grad_zeta_Z * grad_zeta_Z
        pest%d_B_d_x3(:,:,idx1) = d_B_d_zeta
      end if

      if (trim(pest%x3_coord) == 'theta') then
        pest%g13(:,:,idx1) = 2.0*sign_toroidal_flux/(pest%B_ref*(pest%L_ref**2))*(grad_psi_X * grad_theta_X + grad_psi_Y * grad_theta_Y + grad_psi_Z * grad_theta_Z)
        pest%g23(:,:,idx1) = grad_alpha_X * grad_theta_X + grad_alpha_Y * grad_theta_Y + grad_alpha_Z * grad_theta_Z
        pest%g33(:,:,idx1) = grad_theta_X * grad_theta_X + grad_theta_Y * grad_theta_Y + grad_theta_Z * grad_theta_Z
        pest%d_B_d_x3(:,:,idx1) = pest%safety_factor_q(idx1)*d_B_d_zeta
      end if

      pest%gradB_drift_x2(:,:,idx1) = 2.0*sign_toroidal_flux/(pest%B_ref*(pest%L_ref**2))*B_cross_grad_B_dot_grad_psi 
      pest%gradB_drift_x1(:,:,idx1) = B_cross_grad_B_dot_grad_alpha

      pest%curv_drift_x1(:,:,idx1) = sign_toroidal_flux*(B_cross_grad_B_dot_grad_alpha / B + 2.0/(pest%B_ref*(pest%L_ref**2))* mu_0 * d_pressure_d_s * B_cross_grad_psi_dot_grad_alpha / (B*B)) 
      pest%curv_drift_x2(:,:,idx1) = 2.0*sign_toroidal_flux/(pest%B_ref*(pest%L_ref**2))*B_cross_grad_B_dot_grad_psi / (B)

      pest%d_Lambda_d_theta_vmec(:,:,idx1) = d_Lambda_d_theta_vmec
      !*********************************************************************
      ! Finally, assemble the quantities needed for gs2.
      !*********************************************************************

      ! See the latex note gs2_full_surface_stellarator_geometry in the "doc" directory for a derivation of the formulae that follow.

  !    bmag = B / pest%B_ref
  !
  !    gradpar = pest%L_ref * B_sup_zeta / B
  !
  !    gds2 = (grad_alpha_X * grad_alpha_X + grad_alpha_Y * grad_alpha_Y + grad_alpha_Z * grad_alpha_Z) &
  !         * pest%L_ref * pest%L_ref * pest%x1(pest%ix11)
  !
  !    gds21 = (grad_alpha_X * grad_psi_X + grad_alpha_Y * grad_psi_Y + grad_alpha_Z * grad_psi_Z) &
  !         * sign_toroidal_flux * shat / pest%B_ref
  !
  !    gds22 = (grad_psi_X * grad_psi_X + grad_psi_Y * grad_psi_Y + grad_psi_Z * grad_psi_Z) &
  !         * shat * shat / (pest%L_ref * pest%L_ref * pest%B_ref * pest%B_ref * pest%x1(pest%ix11))
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

      pest%Rsurf(:,:,idx1) = R
      pest%Zsurf(:,:,idx1) = Z

      !*********************************************************************
      ! Free all arrays that were allocated.
      !*********************************************************************
    end do ! End global x1 loop

    deallocate(normalized_toroidal_flux_full_grid)
    deallocate(normalized_toroidal_flux_half_grid)
    deallocate(d_iota_d_s_on_half_grid)
    deallocate(d_pressure_d_s_on_half_grid)
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
    deallocate(grad_theta_X)
    deallocate(grad_theta_Y)
    deallocate(grad_theta_Z)
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
    deallocate(B_cross_grad_B_dot_grad_psi)
    deallocate(B_cross_grad_B_dot_grad_alpha)
    deallocate(B_cross_grad_B_dot_grad_alpha_alternate)
    deallocate(B_cross_grad_psi_dot_grad_alpha)
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

      real(dp), dimension(pest%ix21:pest%ix22, pest%ix31:pest%ix32) :: array1, array2
      real(dp) :: tolerance
      character(len=*) :: name
      logical :: should_be_0
      real(dp) :: max_value, max_difference

      if (should_be_0) then
         max_value = maxval(abs(array1))
         if (verbose) print *,"  maxval(abs(",trim(name),")):",max_value,"(should be << 1.)"
         if (max_value > tolerance) then
            print *,"Error! ",trim(name)," should be 0, but instead it is:"
            do idx2 = pest%ix21,pest%ix22
               print *,array1(idx2,:)
            end do
            stop
         end if
      else
         max_difference = maxval(abs(array1 - array2)) / maxval(abs(array1) + abs(array2))
         if (verbose) print *,"  Relative difference between two methods for computing ",trim(name),":",max_difference,"(should be << 1.)"
         if (max_difference > tolerance) then
            print *,"Error! Two methods for computing ",trim(name)," disagree. Here comes method 1:"
            do idx2 = pest%ix21,pest%ix22
               print *,array1(idx2,:)
            end do
            print *,"Here comes method 2:"
            do idx2 = pest%ix21,pest%ix22
               print *,array2(idx2,:)
            end do
            print *,"Here comes the difference:"
            do idx2 = pest%ix21,pest%ix22
               print *,array1(idx2,:) - array2(idx2,:)
            end do
            stop
         end if
      end if

    end subroutine test_arrays

    real(dp) function fzero_residual(theta_vmec_try) result(residual)
      ! Note that lmns and lmnc use the non-Nyquist xm, xn, and mnmax.
      ! Also note that lmns and lmnc are on the HALF grid.

      implicit none
      
      real(dp), intent(in) :: theta_vmec_try
      real(dp) :: angle, sinangle, cosangle
      integer :: imn, which_surface
      !type(VMEC_Obj), intent(in) :: vmec

      ! residual = (theta_pest based on theta_vmec_try) - theta_pest_target = theta_vmec_try + Lambda - theta_pest_target
      residual = theta_vmec_try - theta_pest_target

      do imn = 1, pest%vmec%mnmax
         angle = pest%vmec%xm(imn)*theta_vmec_try - pest%vmec%xn(imn)*zeta0
         sinangle = sin(angle)
         cosangle = cos(angle)
         do which_surface = 1,2
            residual = residual + vmec_radial_weight_half(which_surface) * pest%vmec%lmns(imn,vmec_radial_index_half(which_surface)) * sinangle
            if (pest%vmec%lasym) then
               residual = residual + vmec_radial_weight_half(which_surface) * pest%vmec%lmnc(imn,vmec_radial_index_half(which_surface)) * cosangle
            end if
         end do
      end do

    end function


  end subroutine

  ! --------------------------------------------------------------------------


end module
