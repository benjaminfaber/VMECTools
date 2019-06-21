program vmec2sfl
  use types, only: dp, pi
  use pest_object, only: PEST_Obj, create_PEST_Obj, destroy_PEST_Obj
  use io_core, only: read_input, write_pest_file, geom_file, s0, surf_opt, z_center, n_alpha, n_zeta
  use compute_pest, only: compute_pest_surface, compute_pest_sfl
  use normalizations, only: gene_normalizations
  implicit none

  real(dp), parameter :: pi2 = pi*pi
  integer :: nt_grid, nz_grid, j
  real(dp) :: nfpi, zeta_center
  character(len=2000) :: infile

  real(dp) :: time1, time2
  real(dp) :: dtheta, max_theta, periods
  integer :: global_npol, max_point
 
  type(PEST_Obj) :: pest

  nfpi = 4.0

  call cpu_time(time1)
  infile = 'vmec2sfl.inp'
  call read_input(infile)
  pest = create_PEST_Obj(geom_file,1,n_alpha,n_zeta)
  call compute_pest_surface(s0,surf_opt,pest)
  call compute_pest_sfl(z_center,nfpi,pest)
  call gene_normalizations(pest) 
  call write_pest_file(pest)
  call destroy_PEST_Obj(pest)

end program
