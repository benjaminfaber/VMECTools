!*******************************************************************************
! This module defines the I/O routines used in the vmec2pest calculation
! Author: B.J. Faber (bfaber@wisc.edu)
! Creation date: 30.08.2018
!******************************************************************************

module io_core
  use types, only: dp, pi
  use pest_object, only: PEST_Obj
  implicit none

  public:: read_vmec2pest_input, write_pest_file, write_cylindrical_surface, write_RZ_theta_zeta_grid, &
    & write_gene_geometry_file, write_surface_quantity_cyl, write_surface_quantity_xyz, write_surface_quantity_theta_zeta
  public:: tag, geom_file, outdir, x3_coord, norm_type, &
    & n_field_lines, n_parallel_pts, x3_center, &
    & n_field_periods, surfaces, surf_opt, verbose, test, &
    & output_files, surface_quantities, n_surface_quantities, geom_id

  private 
    character(len=:), pointer :: geom_id
    character(len=2000), target :: geom_file
    character(len=2000) :: tag, outdir
    character(len=5) :: x3_coord
    character(len=7) :: norm_type
    integer :: n_field_lines, n_parallel_pts, surf_opt, n_surface_quantities
    real(dp) :: x3_center, n_field_periods
    real(dp), dimension(501) :: surfaces
    character(len=4), dimension(4) :: output_files
    character(len=32), dimension(24) :: surface_quantities
    logical :: verbose, test

contains
  
  subroutine read_vmec2pest_input(filename)
    ! Reads the parameters from the input file
    implicit none
    character(len=256), intent(in) :: filename
    character(len=:), allocatable, target :: temp
    integer :: iunit, j

    namelist /parameters/ tag, geom_file, outdir, x3_coord, norm_type, &
      & n_field_lines, n_parallel_pts, x3_center, n_field_periods, &
      & surfaces, surf_opt, verbose, test, output_files, surface_quantities
    ! Set default values of parameters
    tag = ''
    x3_coord = 'zeta'
    norm_type = 'minor_r'
    outdir = ''
    
    n_field_lines = 1
    n_parallel_pts = 128
    x3_center = 0.0  

    n_field_periods = 1.0 ! Number of field periods to calculate
   
    surfaces(1) = 0.5
    surf_opt = 0 
    verbose = .false.
    test = .false.
    output_files(1) = 'pest'
    surface_quantities(1) = "" 

    ! Add line to ensure iunit is already not open
    open(newunit=iunit,file=trim(filename),status='old',action='read')
    read(iunit,NML=parameters)
    close(iunit)

    j = 0
    do while(len(trim(surface_quantities(j+1))) .lt. 32) 
print *, len(trim(surface_quantities(j+1)))
      j = j + 1
    end do

    n_surface_quantities = j
    print *, n_surface_quantities
print *, len(trim(geom_file))

    geom_id => geom_file(1:len(trim(geom_file)))
print *, geom_id

  end subroutine

  subroutine write_pest_file(pest,idx1)
    implicit none
    integer, intent(in) :: idx1
    type(PEST_Obj), intent(in) :: pest
    integer :: j, k, iunit
    real(dp) :: prefac
    character(len=2000) :: filename
    character(len=2000) :: filenumber

    prefac = 1.0
    if (trim(pest%x3_coord) .eq. 'theta') then
      prefac = pest%iota(idx1)
    end if

    write(filenumber,"(I0.3)") idx1
    filename = trim(outdir)//"pest_"//trim(tag)//"_surf_"//trim(filenumber)//".dat"
    open(newunit=iunit,file=trim(filename))
    write (iunit,'(A)') '&parameters'
    write (iunit,'(A,F12.7)') '!s0 = ', pest%x1(idx1)
    write (iunit,'(A,F12.7)') '!minor_r = ', pest%L_ref
    write (iunit,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit,'(A,I6)') 'gridpoints = ', pest%nx3-1
    write (iunit,'(A,I6)') 'n_pol = ', ceiling(n_field_periods/pest%vmec%nfp)
    write (iunit,'(A)') '/'
    write (iunit,'(9(A23,2x))') '#theta','g11','g12','g22','g13','g23','g33','|B|','sqrt(g)', 'K2', 'K1', 'dBdz' 
    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        write(iunit,'(12(ES23.12E3,2x))') prefac*pest%x3(k,idx1)/pi, pest%g11(j,k,idx1), pest%g12(j,k,idx1), pest%g22(j,k,idx1), &
          & pest%g13(j,k,idx1), pest%g23(j,k,idx1), pest%g33(j,k,idx1), pest%bmag(j,k,idx1), pest%jac(j,k,idx1), pest%curv_drift_x1(j,k,idx1), pest%curv_drift_x2(j,k,idx1), pest%d_B_d_x3(j,k,idx1) 
      end do
    end do
    close(iunit)

  end subroutine

  subroutine write_gene_geometry_file(pest,idx1)
    implicit none
    integer, intent(in) :: idx1
    type(PEST_Obj), intent(in) :: pest
    integer :: j, k, iunit
    character(len=2000) :: filename
    character(len=2000) :: filenumber

    write(filenumber,"(I0.3)") idx1
    filename = trim(outdir)//"gene_"//trim(tag)//"_surf_"//trim(filenumber)//".dat"
    open(newunit=iunit,file=trim(filename))
    write (iunit,'(A)') '&parameters'
    write (iunit,'(A,F12.7)') '!s0 = ', pest%x1(idx1)
    write (iunit,'(A,F12.7)') '!minor_r = ', pest%L_ref
    write (iunit,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit,'(A,I6)') 'gridpoints = ', pest%nx3-1
    write (iunit,'(A,I6)') 'n_pol = ', ceiling(n_field_periods/pest%vmec%nfp)
    write (iunit,'(A)') '/'
    write (iunit,'(9(A23,2x))') '#theta','g11','g12','g22','g13','g23','g33','|B|','sqrt(g)', 'K2', 'K1', 'dBdz' 
    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        write(iunit,'(9(ES23.12E3,2x))') pest%g11(j,k,idx1), pest%g12(j,k,idx1), pest%g22(j,k,idx1), &
          & pest%bmag(j,k,idx1), pest%jac(j,k,idx1), pest%curv_drift_x1(j,k,idx1), pest%curv_drift_x2(j,k,idx1), pest%d_B_d_x3(j,k,idx1), 0.0 
      end do
    end do
    close(iunit)

  end subroutine

  subroutine write_cylindrical_surface(pest,idx1)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    integer :: j, k, iunit_cyl
    real(dp) :: prefac
    character(len=2000) :: filename_cyl, filenumber
    write(filenumber,"(I0.3)") idx1

    prefac = 1.0
    if (trim(pest%x3_coord) .eq. 'theta') then
      prefac = pest%safety_factor_q(idx1)
    end if

    filename_cyl = trim(outdir)//"cyl_surface_"//trim(tag)//"_surf_"//trim(filenumber)//".dat"
    open(file=trim(filename_cyl),newunit=iunit_cyl)
    write (iunit_cyl,'(A)') '&parameters'
    write (iunit_cyl,'(A,F12.7)') '!s0 = ', pest%x1(idx1) 
    write (iunit_cyl,'(A,F12.7)') '!minor_r = ', pest%minor_r
    write (iunit_cyl,'(A,F12.7)') '!major_R = ', pest%major_R
    write (iunit_cyl,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit_cyl,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit_cyl,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit_cyl,'(A)') '/'
    write (iunit_cyl,'(3(A12))') '# R', 'Z', 'Phi'

    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        write (iunit_cyl,'(3(F12.7,2x))') pest%Rsurf(j,k,idx1), pest%Zsurf(j,k,idx1), -prefac*pest%x3(k,idx1)
      end do
      write (iunit_cyl,'(A)') " "
    end do
    close(iunit_cyl)
  end subroutine

  subroutine write_surface_quantity_cyl(pest,idx1,data_name,surf_data)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    character(len=32), intent(in) :: data_name
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32), intent(in) :: surf_data
    integer :: j, k, iunit_cyl
    real(dp) :: prefac
    character(len=2000) :: filename_cyl, filenumber
    write(filenumber,"(I0.3)") idx1

    prefac = 1.0
    if (trim(pest%x3_coord) .eq. 'theta') then
      prefac = pest%safety_factor_q(idx1)
    end if

    filename_cyl = trim(outdir)//"cyl_surface_"//trim(tag)//"_"//trim(data_name)//"_surf_"//trim(filenumber)//".dat"
    open(file=trim(filename_cyl),newunit=iunit_cyl)
    write (iunit_cyl,'(A)') '&parameters'
    write (iunit_cyl,'(A,F12.7)') '!s0 = ', pest%x1(idx1) 
    write (iunit_cyl,'(A,F12.7)') '!minor_r = ', pest%minor_r
    write (iunit_cyl,'(A,F12.7)') '!major_R = ', pest%major_R
    write (iunit_cyl,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit_cyl,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit_cyl,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit_cyl,'(A)') '/'
    write (iunit_cyl,'(3(A12))') '# R', 'Z', 'Phi'

    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        write (iunit_cyl,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1), pest%Zsurf(j,k,idx1), -prefac*pest%x3(k,idx1), surf_data(j,k)
      end do
      write (iunit_cyl,'(A)') " "
      
      write (iunit_cyl,'(4(F12.7,2x))') pest%Rsurf(pest%ix21,k,idx1), pest%Zsurf(pest%ix21,k,idx1), -prefac*pest%x3(k,idx1), surf_data(pest%ix21,k)
    end do
    close(iunit_cyl)
  end subroutine

  subroutine write_surface_quantity_xyz(pest,idx1,data_name,surf_data)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    character(len=32), intent(in) :: data_name
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32), intent(in) :: surf_data
    integer :: j, k, iunit_xyz
    real(dp) :: prefac
    character(len=2000) :: filename_xyz, filenumber
    write(filenumber,"(I0.3)") idx1

    prefac = 1.0
    if (trim(pest%x3_coord) .eq. 'theta') then
      prefac = pest%safety_factor_q(idx1)
    end if

    filename_xyz = trim(outdir)//"xyz_surface_"//trim(tag)//"_"//trim(data_name)//"_surf_"//trim(filenumber)//".dat"
    open(file=trim(filename_xyz),newunit=iunit_xyz)
    write (iunit_xyz,'(A)') '&parameters'
    write (iunit_xyz,'(A,F12.7)') '!s0 = ', pest%x1(idx1) 
    write (iunit_xyz,'(A,F12.7)') '!minor_r = ', pest%minor_r
    write (iunit_xyz,'(A,F12.7)') '!major_R = ', pest%major_R
    write (iunit_xyz,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit_xyz,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit_xyz,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit_xyz,'(A)') '/'
    write (iunit_xyz,'(3(A12))') '# ', 'X', 'Y', 'Z' 

    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        write (iunit_xyz,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1)*cos(prefac*pest%x3(k,idx1)), pest%Rsurf(j,k,idx1)*sin(prefac*pest%x3(k,idx1)), pest%Zsurf(j,k,idx1), surf_data(j,k)
      end do
      j=pest%ix21
      write (iunit_xyz,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1)*cos(prefac*pest%x3(k,idx1)), pest%Rsurf(j,k,idx1)*sin(prefac*pest%x3(k,idx1)), pest%Zsurf(j,k,idx1), surf_data(j,k)
      write (iunit_xyz,'(A)') " "
    end do
    k=pest%ix31
    do j=pest%ix21,pest%ix22
      write (iunit_xyz,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1)*cos(prefac*pest%x3(k,idx1)), pest%Rsurf(j,k,idx1)*sin(prefac*pest%x3(k,idx1)), pest%Zsurf(j,k,idx1), surf_data(j,k)
    end do
    j=pest%ix21
    write (iunit_xyz,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1)*cos(prefac*pest%x3(k,idx1)), pest%Rsurf(j,k,idx1)*sin(prefac*pest%x3(k,idx1)), pest%Zsurf(j,k,idx1), surf_data(j,k)
    write (iunit_xyz,'(A)') " "

    close(iunit_xyz)
  end subroutine

  subroutine write_surface_quantity_theta_zeta(pest,idx1,data_name,surf_data)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    character(len=32), intent(in) :: data_name
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32), intent(in) :: surf_data
    integer :: i, j, jp1, k, iunit_tz, theta_index
    real(dp) :: prefac, surf_interp, theta, theta_j, theta_jp1, theta_interp, dt1, dt2, delta
    character(len=2000) :: filename_tz, filenumber
    real(dp), parameter :: pi2 = pi+pi
    real(dp), parameter :: eps = 1e-8

    write(filenumber,"(I0.3)") idx1

    prefac = 1.0
    if (trim(pest%x3_coord) .eq. 'theta') then
      prefac = pest%safety_factor_q(idx1)
    end if

    filename_tz = trim(outdir)//"theta_zeta_surface_"//trim(tag)//"_"//trim(data_name)//"_surf_"//trim(filenumber)//".dat"
    open(file=trim(filename_tz),newunit=iunit_tz)
    write (iunit_tz,'(A)') '&parameters'
    write (iunit_tz,'(A,F12.7)') '!s0 = ', pest%x1(idx1) 
    write (iunit_tz,'(A,F12.7)') '!minor_r = ', pest%minor_r
    write (iunit_tz,'(A,F12.7)') '!major_R = ', pest%major_R
    write (iunit_tz,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit_tz,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit_tz,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit_tz,'(A)') '/'
    write (iunit_tz,'(5(A12))') '#theta_index', 'theta', 'zeta_index', 'zeta', trim(data_name)

    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        if (j .lt. pest%ix22) then 
          jp1 = j+1 
        else 
          jp1 = pest%ix21
        end if
        theta = 0.
        theta_j = pest%x2(j) + pest%iota(idx1)*prefac*pest%x3(k,idx1)
        theta_jp1 = pest%x2(jp1) + pest%iota(idx1)*prefac*pest%x3(k,idx1)
        if (theta_j .lt. (0.0 - eps)) then
          theta_j = pi2 + theta_j
        end if
        if (theta_jp1 .lt. (0.0 - eps)) then 
          theta_jp1 = pi2 + theta_jp1
        end if
        i = 0
        do while (theta .lt. (theta_j - eps))
          i = i + 1
          theta = real(i)*pi2/real(pest%nx2)
        end do
        theta_j = mod(theta_j,pi2)
        theta_jp1 = mod(theta_jp1,pi2)

        theta_index = mod(i,pest%nx2)
        theta_interp = real(theta_index)*pi2/real(pest%nx2)
        dt1 = theta_interp - theta_j
        if (dt1 .lt. -eps) then
          dt1 = dt1 + pi2
        end if
        dt2 = theta_jp1 - theta_interp
        if (dt2 .lt. -eps) then
          dt2 = dt2 + pi2
        end if       
        surf_interp = real(pest%nx2)/pi2*(dt2*surf_data(j,k) + dt1*surf_data(jp1,k))
        write (iunit_tz,'(2(I5,2x,F12.7,2x),F12.7)') theta_index+1, theta_interp, k+pest%nx3/2+1, prefac*pest%x3(k,idx1), surf_interp
      end do
      write (iunit_tz,'(A)') " "
    end do

    close(iunit_tz)
  end subroutine

  subroutine write_RZ_theta_zeta_grid(pest,idx1,nfpi)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1, nfpi
    real(dp) :: prefac, theta, theta_j, theta_jp1, theta_interp, dt1, dt2, delta, Rsurf_interp, Zsurf_interp
    integer i, j, jp1, k, iunit_R, iunit_Z, theta_index
    character(len=2000) :: filename_R, filename_Z, filenumber
    real(dp), parameter :: pi2 = pi*pi
    real(dp), parameter :: eps = 1e-8

    if (pest%x3_coord .eq. 'theta') then
      prefac = pest%safety_factor_q(idx1)
    else
      prefac = 1
    end if
    write(filenumber,"(I0.3)") idx1
    delta = (real(pest%nx2)*nfpi)/(pest%vmec%nfp*real(pest%nx3/2))

    filename_R = trim(outdir)//"R_surface_"//trim(tag)//"_surf_"//trim(filenumber)//".dat"
    filename_Z = trim(outdir)//"Z_surface_"//trim(tag)//"_surf_"//trim(filenumber)//".dat"
    open(file=trim(filename_R),newunit=iunit_R)
    open(file=trim(filename_Z),newunit=iunit_Z)
    write (iunit_R,'(A)') '&parameters'
    write (iunit_R,'(A,F12.7)') '!s0 = ', pest%x1(idx1) 
    write (iunit_R,'(A,F12.7)') '!minor_r = ', pest%minor_r
    write (iunit_R,'(A,F12.7)') '!major_R = ', pest%major_R
    write (iunit_R,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit_R,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit_R,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit_R,'(A)') '/'
    write (iunit_Z,'(A)') '&parameters'
    write (iunit_Z,'(A,F12.7)') '!s0 = ', pest%x1(idx1) 
    write (iunit_Z,'(A,F12.7)') '!minor_r = ', pest%minor_r
    write (iunit_Z,'(A,F12.7)') '!major_R = ', pest%major_R
    write (iunit_Z,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit_Z,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit_Z,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit_Z,'(A)') '/'


    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        if (j .lt. pest%ix22) then 
          jp1 = j+1 
        else 
          jp1 = pest%ix21
        end if
        theta = 0.
        theta_j = pest%x2(j) + pest%iota(idx1)*prefac*pest%x3(k,idx1)
        theta_jp1 = pest%x2(jp1) + pest%iota(idx1)*prefac*pest%x3(k,idx1)
        if (theta_j .lt. (0.0 - eps)) then
          theta_j = pi2 + theta_j
        end if
        if (theta_jp1 .lt. (0.0 - eps)) then 
          theta_jp1 = pi2 + theta_jp1
        end if
        i = 0
        do while (theta .lt. (theta_j - eps))
          i = i + 1
          theta = real(i)*pi2/real(pest%nx2)
        end do
        theta_j = mod(theta_j,pi2)
        theta_jp1 = mod(theta_jp1,pi2)

        theta_index = mod(i,pest%nx2)
        theta_interp = real(theta_index)*pi2/real(pest%nx2)
        dt1 = theta_interp - theta_j
        if (dt1 .lt. -eps) then
          dt1 = dt1 + pi2
        end if
        dt2 = theta_jp1 - theta_interp
        if (dt2 .lt. -eps) then
          dt2 = dt2 + pi2
        end if       
        Rsurf_interp = real(pest%nx2)/pi2*(dt2*pest%Rsurf(j,k,idx1) + dt1*pest%Rsurf(jp1,k,idx1))
        Zsurf_interp = real(pest%nx2)/pi2*(dt2*pest%Zsurf(j,k,idx1) + dt1*pest%Zsurf(jp1,k,idx1))
        write (iunit_R,'(2(I5,2x),F12.7)') theta_index+1, k+pest%nx3/2+1, Rsurf_interp
        write (iunit_Z,'(2(I5,2x),F12.7)') theta_index+1, k+pest%nx3/2+1, Zsurf_interp
      end do
      write (iunit_R,'(A)') " "
      write (iunit_Z,'(A)') " "
    end do
    close(iunit_R)
    close(iunit_Z)
  end subroutine
end module
