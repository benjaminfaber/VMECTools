module vmec2sfl_io_mod

  use vmec2sfl_vars_mod
  use read_wout_mod
  implicit none

  public:: read_input, write_sfl_file, write_RZ_surface
  public:: tag, geom_file, outdir, zcoord, n_alpha, points_per_turn, z_center, &
    & n_turns, max_angle, kx_max, ky_min, s0, surf_opt

  character(len=2000) :: tag, geom_file, outdir
  character(len=5) :: zcoord
  integer :: n_alpha, points_per_turn, n_turns, surf_opt
  real(rp) :: z_center, max_angle, kx_max, ky_min, s0

contains
  
  subroutine read_input(filename)
    ! Reads the parameters from the input file
    implicit none
    real(rp), parameter :: pi = 4.0*atan(1.0)
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

  subroutine read_vmec_file(vmec_filename) 
    use read_wout_mod
    implicit none
    
    ! vmec_filename is the vmec wout_* file that will be read.
    character(len=2000), intent(in) :: vmec_filename
    integer :: ierr, iopen

    if (verbose) print *,"  About to read VMEC wout file ",trim(vmec_filename)
    call read_wout_file(vmec_filename, ierr, iopen)
    if (iopen .ne. 0) stop 'error opening wout file'
    if (ierr .ne. 0) stop 'error reading wout file'
    if (verbose) print *,"  Successfully read VMEC data from ",trim(vmec_filename)

    if (verbose) print *,"  Number of field periods (nfp):",nfp
    if (verbose) print *,"  Stellarator-asymmetric? (lasym):",lasym
  end subroutine read_vmec_file

  subroutine write_sfl_file(nzgrid,nturns)
    implicit none
    real(rp), parameter :: pi = 4.0*atan(1.0)
    integer, intent(in) :: nzgrid, nturns
    real(rp) :: iota
    integer :: j, k, iunit
    character(len=2000) :: filename

    iota = 1.0/safety_factor_q
    iunit = 500
    filename = "gist_"//trim(tag)//".dat"
    open(file=trim(filename),unit=iunit)
    write (iunit,'(A)') '&parameters'
    write (iunit,'(A,F12.7)') '!s0 = ', normalized_toroidal_flux_used 
    write (iunit,'(A,F12.7)') '!minor_r = ', L_reference
    write (iunit,'(A,F12.7)') '!Bref = ', B_reference
    write (iunit,'(A,F12.7)') 'q0 = ', safety_factor_q
    write (iunit,'(A,F12.7)') 'shat = ', shat
    write (iunit,'(A,I6)') 'gridpoints = ', nzgrid*2
    write (iunit,'(A,I6)') 'n_pol = ', nturns 
    write (iunit,'(A)') '/'
    do j=1,n_alpha
      do k=-nzgrid, nzgrid-1
        write(iunit,'(9(ES22.12E3,2x))')  gds22(j,k)/(shat**2), gds21(j,k)/shat, gds2(j,k) ,bmag(j,k), &
        & abs(1.0/jac_gist_inv(j,k)), bmag(j,k)/2*cvdrift(j,k), -bmag(j,k)/(2*shat)*cvdrift0(j,k), d_B_d_par(j,k), iota*zeta(k)/pi
      end do
    end do

  end subroutine write_sfl_file

  subroutine write_RZ_surface(ntgrid,nzgrid,nfpi,nfp)
    implicit none
    integer, intent(in) :: ntgrid, nzgrid
    real(rp), intent(in) :: nfpi, nfp 
    real(rp) :: iota, theta, theta_j, theta_jp1, theta_interp, dt1, dt2, delta, Rsurf_interp, Zsurf_interp
    integer i, j, jp1, k, iunit_R, iunit_Z, iunit_cyc, theta_index
    character(len=2000) :: filename_R, filename_Z, filename_cyc
    real(rp), parameter :: pi = 4.0*atan(1.0)
    real(rp), parameter :: pi2 = 8.0*atan(1.0)
    real(rp), parameter :: eps = 1e-8

    delta = (real(ntgrid)*nfpi)/(nfp*real(nzgrid))

    iota = 1.0/safety_factor_q
    filename_R = trim(outdir)//"R_surface_"//trim(tag)//".dat"
    filename_Z = trim(outdir)//"Z_surface_"//trim(tag)//".dat"
    filename_cyc = trim(outdir)//"cyc_surface_"//trim(tag)//".dat"
    iunit_R = 500
    iunit_Z = 600
    iunit_cyc = 700
    open(file=trim(filename_R),unit=iunit_R)
    open(file=trim(filename_Z),unit=iunit_Z)
    open(file=trim(filename_cyc),unit=iunit_cyc)
    do k=-nzgrid,nzgrid
      do j=1,n_alpha
        if (j .lt. n_alpha) then 
          jp1 = j+1 
        else 
          jp1 = 1
        end if
        theta = 0.
        theta_j = alpha(j) + iota*zeta(k)
        theta_jp1 = alpha(jp1) + iota*zeta(k)
        if (theta_j .lt. (0.0 - eps)) then
          theta_j = pi2 + theta_j
        end if
        if (theta_jp1 .lt. (0.0 - eps)) then 
          theta_jp1 = pi2 + theta_jp1
        end if
        i = 0
        do while (theta .lt. (theta_j - eps))
          i = i + 1
          theta = real(i)*pi2/real(ntgrid)
        end do
        theta_j = mod(theta_j,pi2)
        theta_jp1 = mod(theta_jp1,pi2)

        theta_index = mod(i,ntgrid)
        theta_interp = real(theta_index)*pi2/real(ntgrid)
        dt1 = theta_interp - theta_j
        if (dt1 .lt. -eps) then
          dt1 = dt1 + pi2
        end if
        dt2 = theta_jp1 - theta_interp
        if (dt2 .lt. -eps) then
          dt2 = dt2 + pi2
        end if       
        Rsurf_interp = real(ntgrid)/pi2*(dt2*Rsurf(j,k) + dt1*Rsurf(jp1,k))
        Zsurf_interp = real(ntgrid)/pi2*(dt2*Zsurf(j,k) + dt1*Zsurf(jp1,k))
        write (iunit_R,'(2(I5,2x),F12.7)') theta_index+1, k+nzgrid+1, Rsurf_interp
        write (iunit_Z,'(2(I5,2x),F12.7)') theta_index+1, k+nzgrid+1, Zsurf_interp
        write (iunit_cyc,'(3(F12.7,2x))') Rsurf(j,k), Zsurf(j,k), -zeta(k)
      end do
      write (iunit_R,'(A)') " "
      write (iunit_Z,'(A)') " "
      write (iunit_cyc,'(A)') " "
    end do
    close(iunit_R)
    close(iunit_Z)
    close(iunit_cyc)
  end subroutine write_RZ_surface
end module vmec2sfl_io_mod
