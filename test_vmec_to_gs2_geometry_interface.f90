program test_vmec_to_gs2_geometry_interface

  use vmec_to_gs2_geometry_interface_mod

  implicit none

  !*********************************************************************
  ! Input parameters
  !*********************************************************************

  character(len=2000) :: vmec_filename = 'equilibria/wout_w7x_standardConfig.nc'
  integer, parameter :: nalpha = 7
  integer, parameter :: nzgrid = 5
  real :: zeta_center = 0.0
  integer :: number_of_field_periods_to_include = 1
  real :: desired_normalized_toroidal_flux = 0.5d+0
  integer :: vmec_surface_option = 0
  logical :: verbose = .true.

  !*********************************************************************
  ! Output arrays
  !*********************************************************************
  
  real :: normalized_toroidal_flux_used, safety_factor_q, shat
  real, dimension(nalpha) :: alpha
  real, dimension(nzgrid*2+1) :: zeta
  real, dimension(nalpha, nzgrid*2+1) :: bmag, gradpar, gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, cvdrift0
  ! This code uses normalizations in which kxfac is always 1, so kxfac is not presently returned.

  !*********************************************************************
  ! Variables used internally by this program
  !*********************************************************************

  integer :: j

  !*********************************************************************
  ! Beginning of executable statements
  !*********************************************************************

  call vmec_to_gs2_geometry_interface(vmec_filename, nalpha, nzgrid, zeta_center, number_of_field_periods_to_include, &
       desired_normalized_toroidal_flux, vmec_surface_option, verbose, &
       normalized_toroidal_flux_used, safety_factor_q, shat, &
       alpha, zeta, bmag, gradpar, gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, cvdrift0)

  print *,"-------------- Input parameters ------------------"
  print *,"vmec_filename: ",trim(vmec_filename)
  print *,"nalpha:",nalpha
  print *,"nzgrid:",nzgrid
  print *,"zeta_center:",zeta_center
  print *,"number_of_field_periods_to_include:",number_of_field_periods_to_include
  print *,"desired_normalized_toroidal_flux:",desired_normalized_toroidal_flux
  print *,"vmec_surface_option:",vmec_surface_option

  print *,"-------------- Output parameters -----------------"
  print *,"normalized_toroidal_flux_used:",normalized_toroidal_flux_used
  print *,"safety_factor_q:",safety_factor_q
  print *,"shat:",shat
  print *,"alpha:"
  print *,alpha
  print *,"zeta:"
  print *,zeta

  print *,"bmag:"
  do j=1,nalpha
     print *,bmag(j,:)
  end do

  print *,"gradpar:"
  do j=1,nalpha
     print *,gradpar(j,:)
  end do

  print *,"gds2:"
  do j=1,nalpha
     print *,gds2(j,:)
  end do

  print *,"gds21:"
  do j=1,nalpha
     print *,gds21(j,:)
  end do

  print *,"gds22:"
  do j=1,nalpha
     print *,gds22(j,:)
  end do

  print *,"gbdrift:"
  do j=1,nalpha
     print *,gbdrift(j,:)
  end do

  print *,"gbdrift0:"
  do j=1,nalpha
     print *,gbdrift0(j,:)
  end do

  print *,"cvdrfit:"
  do j=1,nalpha
     print *,cvdrift(j,:)
  end do

  print *,"cvdrift0:"
  do j=1,nalpha
     print *,cvdrift0(j,:)
  end do

end program test_vmec_to_gs2_geometry_interface
