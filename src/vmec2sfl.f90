program vmec2sfl
  use types, only: dp, pi
  use Eq_Object, only: Eq_Obj, create_Eq_Obj, destroy_Eq_Obj
  use io_core, only: read_input, geom_file
  implicit none

  real(dp), parameter :: pi2 = pi*pi
  integer :: nt_grid, nz_grid, j
  real(dp) :: nfpi, zeta_center
  character(len=2000) :: infile

  real(dp) :: time1, time2
  real(dp) :: dtheta, max_theta, periods
  integer :: global_npol, max_point
 
  type(Eq_Obj) :: eq 


  call cpu_time(time1)
  infile = 'vmec2sfl.inp'
  call read_input(infile)
  eq = create_Eq_Obj(geom_file)
  print *, eq%vmec%mnmax
  print *, eq%vmec%rmnc(2,51)
  call destroy_Eq_Obj(eq)

end program
