!*******************************************************************************
! This module defines the public types used in the vmec2sfl calculation
! Author: B.J. Faber (bfaber@wisc.edu)
! Creation date: 30.08.2018
!******************************************************************************

module pest_object
  use types, only: dp, pi
  use vmec_object, only: VMEC_Obj, create_VMEC_Obj, destroy_VMEC_Obj
  implicit none

  public PEST_Obj, create_PEST_Obj, destroy_PEST_Obj

  private

  type :: PEST_Obj
    ! Pointer to the VMEC equilibrium
    type(VMEC_Obj) :: vmec

    ! List of all the calculated quantities
    ! On exit, s0 holds the flux surface that was actually used for the geometry,
    ! as measured by psi_toroidal / psi_{toroidal,edge}
    real(dp), allocatable :: s0

    ! The rotational transform iota
    real(dp) :: iota

    ! Safety factor q = 1/iota
    real(dp) :: safety_factor_q

    ! Magnetic shear shat = (x/q) * (d q / d x) where x = Aminor_p * sqrt(psi_toroidal / psi_{toroidal,edge})
    ! and Aminor_p is the minor radius calculated by VMEC.
    real(dp) :: shat

    ! L_ref is the reference length used for gs2's normalization, in meters.
    real(dp) :: L_ref

    ! The major radius in meters 
    real(dp) :: major_R

    ! The minor radius in meters
    real(dp) :: minor_r
  
    ! B_ref is the reference magnetic field strength used for gs2's normalization, in Tesla.
    real(dp) :: B_ref

    ! Integers representing the number of points in each SFL direction 
    integer :: nx1, nx2, nx3

    ! Integers representing the starting and ending array indices
    integer :: ix11, ix12, ix21, ix22, ix31, ix32

    ! On exit, s holds the grid points in the first coordinate (flux surface label)
    ! Typically, this is psi_t = normalized_toroidal_flux
    real(dp), dimension(:), allocatable :: x1(:)

    ! On exit, alpha holds the grid points in the second coordinate (field line label)
    ! Typically, this is alpha = theta_p - iota * zeta, where zeta is the PEST (geometric) toroidal angle
    ! and theta_p is the PEST poloidal angle
    real(dp), dimension(:), allocatable :: x2(:)

    ! On exit, zeta holds the grid points in the third coordinate (field line following coordinate)
    ! Typically, this is zeta = zeta, where zeta is the PEST (geometric) toroidal angle
    real(dp), dimension(:), allocatable :: x3(:)

    character(len=5) :: x3_coord

    ! Arrays that contain the geometric elements on the surface
    real(dp), dimension(:,:,:), allocatable :: bmag  ! Magnitude of B (|B|)
    real(dp), dimension(:,:,:), allocatable :: jac ! Jacobian (grad s . grad alpha x grad zeta)^(-1)
    real(dp), dimension(:,:,:), allocatable :: g11 ! Metric element gss
    real(dp), dimension(:,:,:), allocatable :: g12 ! Metric element gsa
    real(dp), dimension(:,:,:), allocatable :: g22 ! Metric element gaa
    real(dp), dimension(:,:,:), allocatable :: g13 ! Metric element gsz
    real(dp), dimension(:,:,:), allocatable :: g23 ! Metric element gaz
    real(dp), dimension(:,:,:), allocatable :: g33 ! Metric element gzz
    real(dp), dimension(:,:,:), allocatable :: d_B_d_x1 ! Derivative of |B| w.r.t. the 1 coordinate
    real(dp), dimension(:,:,:), allocatable :: d_B_d_x2 ! Derivative of |B| w.r.t. the 2 coordinate
    real(dp), dimension(:,:,:), allocatable :: d_B_d_x3 ! Derivative of |B| w.t.t. the 3 coordinate

    real(dp), dimension(:,:,:), allocatable :: Rsurf ! The R coordinate of the surfaces
    real(dp), dimension(:,:,:), allocatable :: Zsurf ! The Z coordinate of the surfaces
    real(dp), dimension(:,:,:), allocatable :: d_Lambda_d_theta_vmec ! The Z coordinate of the surfaces

  end type ! end type PEST_Obj


  interface create_PEST_Obj
    module procedure create_from_VMEC_Obj
    module procedure create_from_VMEC_file
  end interface
      
 
!    real(dp), pointer :: gds2(:,:) ! Metric element gss
!    real(dp), pointer :: gds21(:,:) ! Metric element gsa
!    real(dp), pointer :: gdsaa(:,:) ! Metric element gaa
!    real(dp), pointer :: gbdrift(:,:) ! grad B drift in s
!    real(dp), pointer :: gbdrift0(:,:) ! grad B drift in a
!    real(dp), pointer :: cvdrift(:,:) ! curv drift in s
!    real(dp), pointer :: cvdrift0(:,:) ! curv drift in a
!    real(dp), pointer :: jac_gist_inv(:,:) ! jacobian
!    real(dp), pointer :: d_B_d_par(:,:) ! parallel derivative of B
!
!
!
!  public :: dp, s0, safety_factor_q, shat, L_ref, B_ref, &
!    & alpha, zeta, bmag, gradpar, gds2, gds21, gdsaa, &
!    & gbdrift, gbdrift0, cvdrift, cvdrift0, jac_gist_inv, d_B_d_par, verbose
!
!  !*********************************************************************
!  ! Input quantities
!  !*********************************************************************
!  ! If verbose is .true., lots of diagnostic information is printed.
!  logical :: verbose, test
!
!  ! VMEC supplies the normalized toroidal flux, thus it is natural to use the
!  ! toroidal zeta coordinate as the angle parameterizing the field line.
!  ! However some codes use the poloidal angle theta, related to zeta by theta =
!  ! iota*zeta.  The maximum extent of the flux tube and the coordinate used in
!  ! output can be chosen to be either theta or zeta, but the calculation itself
!  ! is performed using the zeta coordinate
!  !character(len=5) :: zcoord 
!  !*********************************************************************
!  ! Output quantities
!  !*********************************************************************
!
!  ! Arrays that contain the geometric elements along field lines
!  real(dp), dimension(:,:), allocatable :: bmag  ! Magnitude of B
!  real(dp), dimension(:,:), allocatable :: gradpar ! Grad parallel
!  real(dp), dimension(:,:), allocatable :: gds2 ! Metric element gss
!  real(dp), dimension(:,:), allocatable :: gds21 ! Metric element gsa
!  real(dp), dimension(:,:), allocatable :: gdsaa ! Metric element gaa
!  real(dp), dimension(:,:), allocatable :: gbdrift ! grad B drift in s
!  real(dp), dimension(:,:), allocatable :: gbdrift0 ! grad B drift in a
!  real(dp), dimension(:,:), allocatable :: cvdrift ! curv drift in s
!  real(dp), dimension(:,:), allocatable :: cvdrift0 ! curv drift in a
!  real(dp), dimension(:,:), allocatable :: jac_gist_inv ! jacobian
!  real(dp), dimension(:,:), allocatable :: d_B_d_par ! parallel derivative of B
!
!  ! Surface quantities
!  real(dp), dimension(:,:), allocatable :: Rsurf ! R coordinate as a function of straight field line angles
!  real(dp), dimension(:,:), allocatable :: Zsurf ! Z coordinate as a function of straight field line angles
contains

  type(PEST_Obj) function create_from_VMEC_Obj(vmec,surfaces,n_field_lines,n_parallel_pts) result(pest)
    type(VMEC_Obj), intent(in) :: vmec
    integer, intent(in) :: n_field_lines, n_parallel_pts
    real(dp), dimension(:), intent(in) :: surfaces
    integer :: j, nsurf
    real(dp) :: edge_toroidal_flux_over_2pi
    character(len=7) :: norm_type
    character(len=5) :: x3_coord
    logical :: verbose
    norm_type = 'minor_r'
    x3_coord = 'theta'
    verbose = .true.
    nsurf = size(surfaces)

    pest%vmec = vmec
    pest%nx1 = nsurf
    pest%nx2 = n_field_lines
    pest%nx3 = n_parallel_pts+1
    pest%x3_coord = x3_coord

    pest%ix11 = 0
    pest%ix12 = nsurf-1
    pest%ix21 = 0
    pest%ix22 = n_field_lines-1
    pest%ix31 = -n_parallel_pts/2
    pest%ix32 = n_parallel_pts/2

    if(allocated(pest%x1)) deallocate(pest%x1)
    if(allocated(pest%x2)) deallocate(pest%x2)
    if(allocated(pest%x3)) deallocate(pest%x3)
    if(allocated(pest%bmag)) deallocate(pest%bmag)
    if(allocated(pest%jac)) deallocate(pest%jac)
    if(allocated(pest%g11)) deallocate(pest%g11)
    if(allocated(pest%g12)) deallocate(pest%g12)
    if(allocated(pest%g22)) deallocate(pest%g22)
    if(allocated(pest%g13)) deallocate(pest%g13)
    if(allocated(pest%g23)) deallocate(pest%g23)
    if(allocated(pest%g33)) deallocate(pest%g33)
    if(allocated(pest%d_B_d_x1)) deallocate(pest%d_B_d_x1)
    if(allocated(pest%d_B_d_x2)) deallocate(pest%d_B_d_x2)
    if(allocated(pest%d_B_d_x3)) deallocate(pest%d_B_d_x3)
    if(allocated(pest%Rsurf)) deallocate(pest%Rsurf)
    if(allocated(pest%Zsurf)) deallocate(pest%Zsurf)
    if(allocated(pest%d_Lambda_d_theta_vmec)) deallocate(pest%d_Lambda_d_theta_vmec)

    allocate(pest%x1(pest%ix11:pest%ix12))
    allocate(pest%x2(pest%ix21:pest%ix22))
    allocate(pest%x3(pest%ix31:pest%ix32))
    allocate(pest%bmag(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%jac(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%g11(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%g12(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%g22(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%g13(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%g23(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%g33(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%d_B_d_x1(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%d_B_d_x2(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%d_B_d_x3(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%Rsurf(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%Zsurf(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%d_Lambda_d_theta_vmec(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))

    pest%x1 = surfaces

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set reference values
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    edge_toroidal_flux_over_2pi = pest%vmec%phi(pest%vmec%ns) / (2*pi) * pest%vmec%isigng ! isigns is called signgs in the wout*.nc file. Why is this signgs here?

    ! Set reference length and magnetic field for GS2's normalization, 
    ! using the choices made by Pavlos Xanthopoulos in GIST:
    pest%major_R = pest%vmec%Rmajor
    pest%minor_r = pest%vmec%Aminor ! Note that 'Aminor' in read_wout_mod is called 'Aminor_p' in the wout*.nc file.

    if (norm_type == 'minor_r') then
      pest%L_ref = pest%minor_r
    else
      pest%L_ref = pest%major_R
    end if

    pest%B_ref = 2 * abs(edge_toroidal_flux_over_2pi) / (pest%L_ref * pest%L_ref)
    if (verbose) then
       print *,"  Reference length for normalization:",pest%L_ref," meters."
       print *,"  Reference magnetic field strength normalization:",pest%B_ref," Tesla."
    end if


  end function

  type(PEST_Obj) function create_from_VMEC_file(VMEC_file,surfaces,n_alpha,n_parallel) result(pest)
    character(len=2000), intent(in) :: VMEC_file
    integer, intent(in) :: n_alpha, n_parallel
    real(dp), dimension(:), intent(in) :: surfaces
    type(VMEC_Obj) :: vmec

    vmec = create_VMEC_Obj(VMEC_file)
    pest = create_from_VMEC_Obj(vmec,surfaces,n_alpha,n_parallel)
  end function
 
  subroutine destroy_PEST_Obj(pest)
    type(PEST_Obj), intent(inout) :: pest
    
    call destroy_VMEC_Obj(pest%vmec)
    deallocate(pest%x1)
    deallocate(pest%x2)
    deallocate(pest%x3)
    deallocate(pest%bmag)
    deallocate(pest%jac)
    deallocate(pest%g11)
    deallocate(pest%g12)
    deallocate(pest%g22)
    deallocate(pest%g13)
    deallocate(pest%g23)
    deallocate(pest%g33)
    deallocate(pest%d_B_d_x1)
    deallocate(pest%d_B_d_x2)
    deallocate(pest%d_B_d_x3)
    deallocate(pest%Rsurf)
    deallocate(pest%Zsurf)
    deallocate(pest%d_Lambda_d_theta_vmec)
  end subroutine

end module 
