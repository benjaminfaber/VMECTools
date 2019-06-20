module io_core
  use types, only: dp, pi
  implicit none

  public:: read_input
  public:: tag, geom_file, outdir, zcoord, n_alpha, points_per_turn, z_center, &
    & n_turns, max_angle, kx_max, ky_min, s0, surf_opt

  character(len=2000) :: tag, geom_file, outdir
  character(len=5) :: zcoord
  integer :: n_alpha, points_per_turn, n_turns, surf_opt
  real(dp) :: z_center, max_angle, kx_max, ky_min, s0
  logical :: verbose, test

contains
  
  subroutine read_input(filename)
    ! Reads the parameters from the input file
    implicit none
    character(len=256), intent(in) :: filename
    integer :: iunit

    namelist /parameters/ tag, geom_file, outdir, zcoord, n_alpha, points_per_turn, z_center, n_turns, &
      max_angle, kx_max, ky_min, s0, surf_opt, verbose, test
    ! Set default values of parameters
    tag = ''
    zcoord = 'zeta'
    outdir = ''
    
    n_alpha = 1
    points_per_turn = 128
    z_center = 0.0  

    n_turns = 1 ! Number physical turns around the device
    max_angle = pi
    kx_max = 1.0 ! In units of rho_s
    ky_min = 0.05 ! In units of rho_s 
   
    s0 = 0.5
    surf_opt = 0 
    verbose = .false.
    test = .false.

    iunit=273 
    ! Add line to ensure iunit is already not open
    open(iunit,file=trim(filename),status='old',action='read')
    read(iunit,NML=parameters)
    close(iunit)

  end subroutine read_input
end module
