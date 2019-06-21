module io_core
  use types, only: dp, pi
  use pest_object, only: PEST_Obj
  implicit none

  public:: read_input, write_pest_file
  public:: tag, geom_file, outdir, zcoord, n_alpha, n_zeta, z_center, &
    & n_turns, max_angle, kx_max, ky_min, s0, surf_opt

  character(len=2000) :: tag, geom_file, outdir
  character(len=5) :: zcoord
  integer :: n_alpha, n_zeta, n_turns, surf_opt
  real(dp) :: z_center, max_angle, kx_max, ky_min, s0
  logical :: verbose, test

contains
  
  subroutine read_input(filename)
    ! Reads the parameters from the input file
    implicit none
    character(len=256), intent(in) :: filename
    integer :: iunit

    namelist /parameters/ tag, geom_file, outdir, zcoord, n_alpha, n_zeta, z_center, n_turns, &
      max_angle, kx_max, ky_min, s0, surf_opt, verbose, test
    ! Set default values of parameters
    tag = ''
    zcoord = 'zeta'
    outdir = ''
    
    n_alpha = 1
    n_zeta = 128
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

  subroutine write_pest_file(pest)
    implicit none
    type(PEST_Obj), intent(in) :: pest
    integer :: j, k, iunit
    character(len=2000) :: filename

    iunit = 500
    filename = "pest_"//trim(tag)//".dat"
    open(file=trim(filename),unit=iunit)
    write (iunit,'(A)') '&parameters'
    write (iunit,'(A,F12.7)') '!s0 = ', pest%s0 
    write (iunit,'(A,F12.7)') '!minor_r = ', pest%L_ref
    write (iunit,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit,'(A,F12.7)') 'q0 = ', pest%safety_factor_q
    write (iunit,'(A,F12.7)') 'shat = ', pest%shat
    write (iunit,'(A,I6)') 'gridpoints = ', pest%nzeta+1 
    write (iunit,'(A,I6)') 'n_pol = ', 1 
    write (iunit,'(A)') '/'
    write (iunit,'(9(A23,2x))') '#theta','gxx','gxy','gyy','gxz','gyz','gzz','|B|','sqrt(g)' 
    do k=pest%iz1,pest%iz2
      do j=pest%ia1,pest%ia2
        write(iunit,'(9(ES23.12E3,2x))') pest%iota*pest%zeta(k)/pi, pest%gss(j,k,1), pest%gsa(j,k,1), pest%gaa(j,k,1), &
          & pest%gsz(j,k,1), pest%gaz(j,k,1), pest%gzz(j,k,1), pest%bmag(j,k,1), pest%jac(j,k,1) 
      end do
    end do
    close(iunit)

  end subroutine


end module
