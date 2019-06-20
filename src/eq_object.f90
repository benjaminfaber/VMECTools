!*******************************************************************************
! This module defines the public types used in the vmec2sfl calculation
! Author: B.J. Faber (bfaber@wisc.edu)
! Creation date: 30.08.2018
!******************************************************************************

module eq_object
  use types, only: dp
  use vmec_core, only: VMEC_Obj, create_VMEC_Obj, destroy_VMEC_Obj
  implicit none

  public Eq_Obj, create_Eq_Obj, destroy_Eq_Obj

  private

  type :: Eq_Obj
    ! Pointer to the VMEC equilibrium
    type(VMEC_Obj) :: vmec

    ! List of all the calculated quantities
    ! On exit, normalized_toroidal_flux_used holds the flux surface that was actually used for the geometry,
    ! as measured by psi_toroidal / psi_{toroidal,edge}
    real(dp) :: normalized_toroidal_flux_used

    ! The rotational transform iota
    real(dp) :: iota

    ! Safety factor q = 1/iota
    real(dp) :: safety_factor_q

    ! Magnetic shear shat = (x/q) * (d q / d x) where x = Aminor_p * sqrt(psi_toroidal / psi_{toroidal,edge})
    ! and Aminor_p is the minor radius calculated by VMEC.
    real(dp) :: shat

    ! L_reference is the reference length used for gs2's normalization, in meters.
    real(dp) :: L_reference

    ! The major radius in meters 
    real(dp) :: major_R

    ! The minor radius in meters
    real(dp) :: minor_r
  
    ! B_reference is the reference magnetic field strength used for gs2's normalization, in Tesla.
    real(dp) :: B_reference

    ! Integers representing the number of points in each SFL direction 
    integer :: nx1, nx2, nx3

    ! On exit, x1 holds the grid points in the first coordinate (flux surface label)
    ! Typically, this is psi_t = normalized_toroidal_flux
    real(dp), dimension(:), allocatable :: x1(:)

    ! On exit, x2 holds the grid points in the second coordinate (field line label)
    ! Typically, this is alpha = theta_p - iota * zeta, where zeta is the PEST (geometric) toroidal angle
    ! and theta_p is the PEST poloidal angle
    real(dp), dimension(:), allocatable :: x2(:)

    ! On exit, x3 holds the grid points in the third coordinate (field line following coordinate)
    ! Typically, this is x3 = zeta, where zeta is the PEST (geometric) toroidal angle
    real(dp), dimension(:), allocatable :: x3(:)

    ! Arrays that contain the geometric elements on the surface
    real(dp), dimension(:,:,:), allocatable :: bmag  ! Magnitude of B (|B|)
    real(dp), dimension(:,:,:), allocatable :: jac ! Jacobian (grad x1 . grad x2 x grad x3)^(-1)
    real(dp), dimension(:,:,:), allocatable :: g11 ! Metric element g11
    real(dp), dimension(:,:,:), allocatable :: g12 ! Metric element g12
    real(dp), dimension(:,:,:), allocatable :: g22 ! Metric element g22
    real(dp), dimension(:,:,:), allocatable :: g13 ! Metric element g13
    real(dp), dimension(:,:,:), allocatable :: g23 ! Metric element g23
    real(dp), dimension(:,:,:), allocatable :: g33 ! Metric element g33
    real(dp), dimension(:,:,:), allocatable :: d_B_d_x1 ! Derivative of |B| w.r.t. the 1 coordinate
    real(dp), dimension(:,:,:), allocatable :: d_B_d_x2 ! Derivative of |B| w.r.t. the 2 coordinate
    real(dp), dimension(:,:,:), allocatable :: d_B_d_x3 ! Derivative of |B| w.t.t. the 3 coordinate

  end type ! end type Eq_Obj


  interface create_Eq_Obj
    module procedure create_from_VMEC_Obj
    module procedure create_from_VMEC_file
  end interface
      
 
!    real(dp), pointer :: gds2(:,:) ! Metric element gss
!    real(dp), pointer :: gds21(:,:) ! Metric element gsa
!    real(dp), pointer :: gds22(:,:) ! Metric element gaa
!    real(dp), pointer :: gbdrift(:,:) ! grad B drift in s
!    real(dp), pointer :: gbdrift0(:,:) ! grad B drift in a
!    real(dp), pointer :: cvdrift(:,:) ! curv drift in s
!    real(dp), pointer :: cvdrift0(:,:) ! curv drift in a
!    real(dp), pointer :: jac_gist_inv(:,:) ! jacobian
!    real(dp), pointer :: d_B_d_par(:,:) ! parallel derivative of B
!
!
!
!  public :: dp, normalized_toroidal_flux_used, safety_factor_q, shat, L_reference, B_reference, &
!    & alpha, zeta, bmag, gradpar, gds2, gds21, gds22, &
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
!  real(dp), dimension(:,:), allocatable :: gds22 ! Metric element gaa
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

  type(Eq_Obj) function create_from_VMEC_Obj(vmec) result(eq)
    type(VMEC_Obj), intent(in) :: vmec

    eq%vmec = vmec
    if(allocated(eq%x1)) deallocate(eq%x1)
    if(allocated(eq%x2)) deallocate(eq%x2)
    if(allocated(eq%x3)) deallocate(eq%x3)
    if(allocated(eq%bmag)) deallocate(eq%bmag)
    if(allocated(eq%jac)) deallocate(eq%jac)
    if(allocated(eq%g11)) deallocate(eq%g11)
    if(allocated(eq%g12)) deallocate(eq%g12)
    if(allocated(eq%g22)) deallocate(eq%g22)
    if(allocated(eq%g13)) deallocate(eq%g13)
    if(allocated(eq%g23)) deallocate(eq%g23)
    if(allocated(eq%g33)) deallocate(eq%g33)
    if(allocated(eq%d_B_d_x1)) deallocate(eq%d_B_d_x1)
    if(allocated(eq%d_B_d_x2)) deallocate(eq%d_B_d_x2)
    if(allocated(eq%d_B_d_x3)) deallocate(eq%d_B_d_x3)

    allocate(eq%x1(eq%nx1))
    allocate(eq%x2(eq%nx2))
    allocate(eq%x3(eq%nx3))
    allocate(eq%bmag(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%jac(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g11(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g12(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g22(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g13(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g23(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g33(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%d_B_d_x1(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%d_B_d_x2(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%d_B_d_x3(eq%nx3,eq%nx2,eq%nx1))
  end function

  type(Eq_Obj) function create_from_VMEC_file(VMEC_file) result(eq)
    character(len=2000), intent(in) :: VMEC_file

    eq%vmec = create_VMEC_Obj(VMEC_file)
    if(allocated(eq%x1)) deallocate(eq%x1)
    if(allocated(eq%x2)) deallocate(eq%x2)
    if(allocated(eq%x3)) deallocate(eq%x3)
    if(allocated(eq%bmag)) deallocate(eq%bmag)
    if(allocated(eq%jac)) deallocate(eq%jac)
    if(allocated(eq%g11)) deallocate(eq%g11)
    if(allocated(eq%g12)) deallocate(eq%g12)
    if(allocated(eq%g22)) deallocate(eq%g22)
    if(allocated(eq%g13)) deallocate(eq%g13)
    if(allocated(eq%g23)) deallocate(eq%g23)
    if(allocated(eq%g33)) deallocate(eq%g33)
    if(allocated(eq%d_B_d_x1)) deallocate(eq%d_B_d_x1)
    if(allocated(eq%d_B_d_x2)) deallocate(eq%d_B_d_x2)
    if(allocated(eq%d_B_d_x3)) deallocate(eq%d_B_d_x3)

    allocate(eq%x1(eq%nx1))
    allocate(eq%x2(eq%nx2))
    allocate(eq%x3(eq%nx3))
    allocate(eq%bmag(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%jac(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g11(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g12(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g22(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g13(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g23(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%g33(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%d_B_d_x1(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%d_B_d_x2(eq%nx3,eq%nx2,eq%nx1))
    allocate(eq%d_B_d_x3(eq%nx3,eq%nx2,eq%nx1))
  end function
 
  subroutine destroy_Eq_Obj(eq)
    type(Eq_Obj), intent(inout) :: eq
    
    call destroy_VMEC_Obj(eq%vmec)
    deallocate(eq%x1)
    deallocate(eq%x2)
    deallocate(eq%x3)
    deallocate(eq%bmag)
    deallocate(eq%jac)
    deallocate(eq%g11)
    deallocate(eq%g12)
    deallocate(eq%g22)
    deallocate(eq%g13)
    deallocate(eq%g23)
    deallocate(eq%g33)
    deallocate(eq%d_B_d_x1)
    deallocate(eq%d_B_d_x2)
    deallocate(eq%d_B_d_x3)
  end subroutine

end module 
