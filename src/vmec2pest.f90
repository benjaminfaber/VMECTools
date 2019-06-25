program vmec2sfl
  use types, only: dp, pi
  use pest_object, only: PEST_Obj, create_PEST_Obj, destroy_PEST_Obj
  use io_core, only: read_input, write_pest_file, geom_file, surfaces, surf_opt, grid_center, n_field_lines, n_parallel_pts, n_field_periods
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

  call cpu_time(time1)
  infile = 'vmec2pest.inp'
  call read_input(infile)
  pest = create_PEST_Obj(geom_file,surfaces,n_field_lines,n_parallel_pts)
  call compute_pest_surface(pest%x1(pest%ix11),surf_opt,pest)
  call compute_pest_sfl(grid_center,n_field_periods,pest)
  call gene_normalizations(pest) 
  call write_pest_file(pest)
  call destroy_PEST_Obj(pest)

end program
