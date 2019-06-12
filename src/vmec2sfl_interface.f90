program vmec2sfl_interface

  use vmec2sfl_vars_mod
  use vmec2sfl_io_mod
  use vmec2sfl_mod
  use read_wout_mod, only: nfp

  implicit none
  real(rp), parameter :: pi = 4.0*atan(1.0)
  real(rp), parameter :: pi2 = 8.0*atan(1.0)
  integer :: nt_grid, nz_grid, j
  real(rp) :: nfpi, zeta_center
  character(len=256) :: infile

  real(rp) :: time1, time2
  real(rp) :: dtheta, max_theta, periods
  integer :: global_npol, max_point
  ! Output quantities
  
  call cpu_time(time1)
  ! Make sure all allocatable arrays are dellocated before proceeding
  if (allocated(alpha)) deallocate(alpha)
  if (allocated(zeta)) deallocate(zeta)
  if (allocated(bmag)) deallocate(bmag)
  if (allocated(gradpar)) deallocate(gradpar)
  if (allocated(gds2)) deallocate(gds2)
  if (allocated(gds21)) deallocate(gds21)
  if (allocated(gds22)) deallocate(gds22)
  if (allocated(gbdrift)) deallocate(gbdrift)
  if (allocated(gbdrift0)) deallocate(gbdrift0)
  if (allocated(cvdrift)) deallocate(cvdrift)
  if (allocated(cvdrift0)) deallocate(cvdrift0)
  if (allocated(jac_gist_inv)) deallocate(jac_gist_inv)
  if (allocated(d_B_d_par)) deallocate(d_B_d_par)
  
  infile = 'vmec2sfl.inp'
  call read_input(infile)

  call read_vmec_file(geom_file)
  call get_surface_quantities(s0,surf_opt)

!  dtheta = pi2/real(points_per_turn) 
  max_theta = kx_max/(abs(shat)*ky_min) + 24.0*pi2
  global_npol = nint(max_theta/pi)
  if (modulo(global_npol,2) == 1) global_npol = global_npol+1
  max_theta = real(global_npol)*pi
  max_point = global_npol*points_per_turn

  !global_npol = n_turns
  !periods = real(global_npol)*nfp
  !nz_grid = max_point
  !nz_grid = (points_per_turn - 1)/2 * global_npol
  !nz_grid = points_per_turn/2 * global_npol
  periods = 4.0
  nz_grid = points_per_turn/2
  nfpi = periods*safety_factor_q
  nfpi = periods

  allocate(alpha(n_alpha))
  allocate(zeta(-nz_grid:nz_grid))
  allocate(bmag(n_alpha,-nz_grid:nz_grid))
  allocate(gradpar(n_alpha,-nz_grid:nz_grid))
  allocate(gds2(n_alpha,-nz_grid:nz_grid))
  allocate(gds21(n_alpha,-nz_grid:nz_grid))
  allocate(gds22(n_alpha,-nz_grid:nz_grid))
  allocate(gbdrift(n_alpha,-nz_grid:nz_grid))
  allocate(gbdrift0(n_alpha,-nz_grid:nz_grid))
  allocate(cvdrift(n_alpha,-nz_grid:nz_grid))
  allocate(cvdrift0(n_alpha,-nz_grid:nz_grid))
  allocate(jac_gist_inv(n_alpha,-nz_grid:nz_grid))
  allocate(d_B_d_par(n_alpha,-nz_grid:nz_grid)) 
  allocate(Rsurf(n_alpha,-nz_grid:nz_grid))
  allocate(Zsurf(n_alpha,-nz_grid:nz_grid))
   
  call vmec2sfl(n_alpha,nz_grid,zeta_center,nfpi)

!print *, gds22(1,-nz_grid), gds21(1,-nz_grid),gds2(1,-nz_grid), Bmag(1,-nz_grid),jac_gist_inv(1,-nz_grid),cvdrift(1,-nz_grid),cvdrift0(1,-nz_grid),d_B_d_par(1,-nz_grid)
!  DO j = 0,nz_grid
!    WRITE(6, "(9(F15.10,2x))") gds22(1,j), gds21(1,j), gds2(1,j), Bmag(1,j), &
!      & jac_gist_inv(1,j), cvdrift(1,j), cvdrift0(1,j), d_B_d_par(1,j), zeta(j)/safety_factor_q
!  END DO
  !call write_sfl_file(nz_grid,global_npol)
  nt_grid = n_alpha
  call write_RZ_surface(nt_grid,nz_grid,periods,periods)

  deallocate(alpha)
  deallocate(zeta)
  deallocate(bmag)
  deallocate(gradpar)
  deallocate(gds2)
  deallocate(gds21)
  deallocate(gds22)
  deallocate(gbdrift)
  deallocate(gbdrift0)
  deallocate(cvdrift)
  deallocate(cvdrift0)
  deallocate(jac_gist_inv)
  deallocate(d_B_d_par)
  deallocate(Rsurf)
  deallocate(Zsurf)
  call cpu_time(time2)
write(6,"(A,F8.5)") "Time for vmec2sfl: ", time2-time1
end program vmec2sfl_interface 
