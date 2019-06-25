module io_core
  use types, only: dp, pi
  use pest_object, only: PEST_Obj
  implicit none

  public:: read_input, write_pest_file
  public:: tag, geom_file, outdir, x3_coord, norm_type, &
    & n_field_lines, n_parallel_pts, x3_center, &
    & n_field_periods, surfaces, surf_opt

  character(len=2000) :: tag, geom_file, outdir
  character(len=5) :: x3_coord
  character(len=7) :: norm_type
  integer :: n_field_lines, n_parallel_pts, surf_opt
  real(dp) :: x3_center, n_field_periods
  real(dp), dimension(:), allocatable :: surfaces
  logical :: verbose, test

contains
  
  subroutine read_input(filename)
    ! Reads the parameters from the input file
    implicit none
    character(len=256), intent(in) :: filename
    integer :: iunit

    namelist /parameters/ tag, geom_file, outdir, x3_coord, norm_type, &
      & n_field_lines, n_parallel_pts, x3_center, n_field_periods, &
      & surfaces, surf_opt, verbose, test
    ! Set default values of parameters
    tag = ''
    x3_coord = 'zeta'
    norm_type = 'minor_r'
    outdir = ''
    
    n_field_lines = 1
    n_parallel_pts = 128
    x3_center = 0.0  

    n_field_periods = 1.0 ! Number of field periods to calculate
   
    surfaces = (/0.5/)
    surf_opt = 0 
    verbose = .false.
    test = .false.

    ! Add line to ensure iunit is already not open
    open(newunit=iunit,file=trim(filename),status='old',action='read')
    read(iunit,NML=parameters)
    close(iunit)

  end subroutine read_input

  subroutine write_pest_file(pest)
    implicit none
    type(PEST_Obj), intent(in) :: pest
    integer :: j, k, iunit
    character(len=2000) :: filename

    filename = "pest_"//trim(tag)//".dat"
    open(newunit=iunit,file=trim(filename))
    write (iunit,'(A)') '&parameters'
    write (iunit,'(A,F12.7)') '!s0 = ', pest%x1(pest%ix11)
    write (iunit,'(A,F12.7)') '!minor_r = ', pest%L_ref
    write (iunit,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit,'(A,F12.7)') 'q0 = ', pest%safety_factor_q
    write (iunit,'(A,F12.7)') 'shat = ', pest%shat
    write (iunit,'(A,I6)') 'gridpoints = ', pest%nx3-1
    write (iunit,'(A,I6)') 'n_pol = ', 1 
    write (iunit,'(A)') '/'
    write (iunit,'(9(A23,2x))') '#theta','g11','g12','g22','g13','g23','g33','|B|','sqrt(g)' 
    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        write(iunit,'(11(ES23.12E3,2x))') pest%x3(k)/pi, pest%g11(j,k,pest%ix11), pest%g12(j,k,pest%ix11), pest%g22(j,k,pest%ix11), &
          & pest%g13(j,k,pest%ix11), pest%g23(j,k,pest%ix11), pest%g33(j,k,pest%ix11), pest%bmag(j,k,pest%ix11), pest%jac(j,k,pest%ix11), pest%curv_drift_x1(j,k,pest%ix11), pest%curv_drift_x2(j,k,pest%ix11) 
      end do
    end do
    close(iunit)

  end subroutine


end module
