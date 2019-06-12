!*******************************************************************************
! This module defines the public variables used in the vmec2sfl calculation
! Author: B.J. Faber (bfaber@wisc.edu)
! Creation date: 30.08.2018
!******************************************************************************

module vmec2sfl_vars_mod

  implicit none

  public :: rp, normalized_toroidal_flux_used, safety_factor_q, shat, L_reference, B_reference, &
    & alpha, zeta, bmag, gradpar, gds2, gds21, gds22, &
    & gbdrift, gbdrift0, cvdrift, cvdrift0, jac_gist_inv, d_B_d_par, verbose

  integer, parameter :: rp = selected_real_kind(15,300)
  !*********************************************************************
  ! Input quantities
  !*********************************************************************
  ! If verbose is .true., lots of diagnostic information is printed.
  logical :: verbose, test

  ! VMEC supplies the normalized toroidal flux, thus it is natural to use the
  ! toroidal zeta coordinate as the angle parameterizing the field line.
  ! However some codes use the poloidal angle theta, related to zeta by theta =
  ! iota*zeta.  The maximum extent of the flux tube and the coordinate used in
  ! output can be chosen to be either theta or zeta, but the calculation itself
  ! is performed using the zeta coordinate
  !character(len=5) :: zcoord 
  !*********************************************************************
  ! Output quantities
  !*********************************************************************

  ! On exit, normalized_toroidal_flux_used holds the flux surface that was actually used for the geometry,
  ! as measured by psi_toroidal / psi_{toroidal,edge}
  real(rp) :: normalized_toroidal_flux_used

  ! Safety factor q = 1/iota
  real(rp) :: safety_factor_q

  ! Magnetic shear shat = (x/q) * (d q / d x) where x = Aminor_p * sqrt(psi_toroidal / psi_{toroidal,edge})
  ! and Aminor_p is the minor radius calculated by VMEC.
  real(rp) :: shat

  ! L_reference is the reference length used for gs2's normalization, in meters.
  real(rp) :: L_reference
  real(rp) :: R0

  ! B_reference is the reference magnetic field strength used for gs2's normalization, in Tesla.
  real(rp) :: B_reference

  ! On exit, alpha holds the grid points in alpha = theta_p - iota * zeta, where theta_p is the PEST toroidal angle
  real(rp), dimension(:), allocatable :: alpha

  ! On exit, zeta holds the grid points in the toroidal angle zeta
  real(rp), dimension(:), allocatable:: zeta

  ! Arrays that contain the geometric elements along field lines
  real(rp), dimension(:,:), allocatable :: bmag  ! Magnitude of B
  real(rp), dimension(:,:), allocatable :: gradpar ! Grad parallel
  real(rp), dimension(:,:), allocatable :: gds2 ! Metric element gss
  real(rp), dimension(:,:), allocatable :: gds21 ! Metric element gsa
  real(rp), dimension(:,:), allocatable :: gds22 ! Metric element gaa
  real(rp), dimension(:,:), allocatable :: gbdrift ! grad B drift in s
  real(rp), dimension(:,:), allocatable :: gbdrift0 ! grad B drift in a
  real(rp), dimension(:,:), allocatable :: cvdrift ! curv drift in s
  real(rp), dimension(:,:), allocatable :: cvdrift0 ! curv drift in a
  real(rp), dimension(:,:), allocatable :: jac_gist_inv ! jacobian
  real(rp), dimension(:,:), allocatable :: d_B_d_par ! parallel derivative of B

  ! Surface quantities
  real(rp), dimension(:,:), allocatable :: Rsurf ! R coordinate as a function of straight field line angles
  real(rp), dimension(:,:), allocatable :: Zsurf ! Z coordinate as a function of straight field line angles

end module 
