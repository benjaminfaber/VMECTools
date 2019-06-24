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
    real(dp) :: s0

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
    integer :: ns, nalpha, nzeta

    ! Integers representing the starting and ending array indices
    integer :: is1, is2, ia1, ia2, iz1, iz2

    ! On exit, s holds the grid points in the first coordinate (flux surface label)
    ! Typically, this is psi_t = normalized_toroidal_flux
    real(dp), dimension(:), allocatable :: s(:)

    ! On exit, alpha holds the grid points in the second coordinate (field line label)
    ! Typically, this is alpha = theta_p - iota * zeta, where zeta is the PEST (geometric) toroidal angle
    ! and theta_p is the PEST poloidal angle
    real(dp), dimension(:), allocatable :: alpha(:)

    ! On exit, zeta holds the grid points in the third coordinate (field line following coordinate)
    ! Typically, this is zeta = zeta, where zeta is the PEST (geometric) toroidal angle
    real(dp), dimension(:), allocatable :: zeta(:)

    ! Arrays that contain the geometric elements on the surface
    real(dp), dimension(:,:,:), allocatable :: bmag  ! Magnitude of B (|B|)
    real(dp), dimension(:,:,:), allocatable :: jac ! Jacobian (grad s . grad alpha x grad zeta)^(-1)
    real(dp), dimension(:,:,:), allocatable :: gss ! Metric element gss
    real(dp), dimension(:,:,:), allocatable :: gsa ! Metric element gsa
    real(dp), dimension(:,:,:), allocatable :: gaa ! Metric element gaa
    real(dp), dimension(:,:,:), allocatable :: gsz ! Metric element gsz
    real(dp), dimension(:,:,:), allocatable :: gaz ! Metric element gaz
    real(dp), dimension(:,:,:), allocatable :: gzz ! Metric element gzz
    real(dp), dimension(:,:,:), allocatable :: d_B_d_s ! Derivative of |B| w.r.t. the 1 coordinate
    real(dp), dimension(:,:,:), allocatable :: d_B_d_alpha ! Derivative of |B| w.r.t. the 2 coordinate
    real(dp), dimension(:,:,:), allocatable :: d_B_d_zeta ! Derivative of |B| w.t.t. the 3 coordinate

    real(dp), dimension(:,:,:), allocatable :: Rsurf ! The R coordinate of the surfaces
    real(dp), dimension(:,:,:), allocatable :: Zsurf ! The Z coordinate of the surfaces
    real(dp), dimension(:,:,:), allocatable :: d_L_d_theta_v ! The Z coordinate of the surfaces

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

  type(PEST_Obj) function create_from_VMEC_Obj(vmec,ns,n_alpha,n_zeta) result(pest)
    type(VMEC_Obj), intent(in) :: vmec
    integer, intent(in) :: ns, n_alpha, n_zeta
    integer :: j
    real(dp) :: edge_toroidal_flux_over_2pi
    character(len=7) :: norm_type
    logical :: verbose
    norm_type = 'minor_r'
    verbose = .false.

    pest%vmec = vmec
    pest%ns = ns
    pest%nalpha = n_alpha
    pest%nzeta = n_zeta+1

    pest%is1 = 1
    pest%is2 = 1
    pest%ia1 = 0
    pest%ia2 = n_alpha-1
    pest%iz1 = -n_zeta/2
    pest%iz2 = n_zeta/2

    if(allocated(pest%s)) deallocate(pest%s)
    if(allocated(pest%alpha)) deallocate(pest%alpha)
    if(allocated(pest%zeta)) deallocate(pest%zeta)
    if(allocated(pest%bmag)) deallocate(pest%bmag)
    if(allocated(pest%jac)) deallocate(pest%jac)
    if(allocated(pest%gss)) deallocate(pest%gss)
    if(allocated(pest%gsa)) deallocate(pest%gsa)
    if(allocated(pest%gaa)) deallocate(pest%gaa)
    if(allocated(pest%gsz)) deallocate(pest%gsz)
    if(allocated(pest%gaz)) deallocate(pest%gaz)
    if(allocated(pest%gzz)) deallocate(pest%gzz)
    if(allocated(pest%d_B_d_s)) deallocate(pest%d_B_d_s)
    if(allocated(pest%d_B_d_alpha)) deallocate(pest%d_B_d_alpha)
    if(allocated(pest%d_B_d_zeta)) deallocate(pest%d_B_d_zeta)
    if(allocated(pest%Rsurf)) deallocate(pest%Rsurf)
    if(allocated(pest%Zsurf)) deallocate(pest%Zsurf)
    if(allocated(pest%d_L_d_theta_v)) deallocate(pest%d_L_d_theta_v)

    allocate(pest%s(pest%ns))
    allocate(pest%alpha(pest%ia1:pest%ia2))
    allocate(pest%zeta(pest%iz1:pest%iz2))
    allocate(pest%bmag(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%jac(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%gss(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%gsa(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%gaa(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%gsz(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%gaz(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%gzz(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%d_B_d_s(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%d_B_d_alpha(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%d_B_d_zeta(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%Rsurf(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%Zsurf(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))
    allocate(pest%d_L_d_theta_v(pest%ia1:pest%ia2,pest%iz1:pest%iz2,pest%ns))

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

  type(PEST_Obj) function create_from_VMEC_file(VMEC_file,ns,n_alpha,n_zeta) result(pest)
    character(len=2000), intent(in) :: VMEC_file
    integer, intent(in) :: ns, n_alpha, n_zeta
    type(VMEC_Obj) :: vmec

    vmec = create_VMEC_Obj(VMEC_file)
    pest = create_from_VMEC_Obj(vmec,ns,n_alpha,n_zeta)
  end function
 
  subroutine destroy_PEST_Obj(pest)
    type(PEST_Obj), intent(inout) :: pest
    
    call destroy_VMEC_Obj(pest%vmec)
    deallocate(pest%s)
    deallocate(pest%alpha)
    deallocate(pest%zeta)
    deallocate(pest%bmag)
    deallocate(pest%jac)
    deallocate(pest%gss)
    deallocate(pest%gsa)
    deallocate(pest%gaa)
    deallocate(pest%gsz)
    deallocate(pest%gaz)
    deallocate(pest%gzz)
    deallocate(pest%d_B_d_s)
    deallocate(pest%d_B_d_alpha)
    deallocate(pest%d_B_d_zeta)
    deallocate(pest%Rsurf)
    deallocate(pest%Zsurf)
    deallocate(pest%d_L_d_theta_v)
  end subroutine

end module 
