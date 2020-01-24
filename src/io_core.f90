!*******************************************************************************
! This module defines the I/O routines used in the vmec2pest calculation
! Author: B.J. Faber (bfaber@wisc.edu)
! Creation date: 30.08.2018
!******************************************************************************

module io_core
  use types, only: dp, pi
  use pest_object, only: PEST_Obj
  use hdf5
  implicit none

  public:: read_vmec2pest_input, write_pest_file, write_cylindrical_surface, write_RZ_theta_zeta_grid, &
    & write_gene_geometry_file, write_surface_quantity_cyl, write_surface_quantity_xyz, write_surface_quantity_theta_zeta, &
    & create_mapping_file_hdf5, close_mapping_file_hdf5, write_mapping_file_info_hdf5, write_mapping_file_surf_data_hdf5
  public:: tag, geom_file, outdir, x3_coord, norm_type, &
    & n_surf, n_field_lines, n_parallel_pts, n_pol, x2_center, x3_center, &
    & n_field_periods, surfaces, surf_opt, verbose, test, &
    & output_files, surface_quantities, n_surface_quantities, geom_id

  private 
    character(len=:), pointer :: geom_id
    character(len=2000), target :: geom_file
    character(len=2000) :: tag, outdir
    character(len=5) :: x3_coord
    character(len=7) :: norm_type
    integer :: n_surf, n_field_lines, n_parallel_pts, surf_opt, n_surface_quantities, n_pol
    real(dp) :: x2_center, x3_center, n_field_periods
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
    character(len=32) :: temp_str
    integer :: iunit, j
    logical :: end_found

    namelist /parameters/ tag, geom_file, outdir, x3_coord, norm_type, &
      & n_field_lines, n_parallel_pts, x2_center, x3_center, n_field_periods, n_pol, &
      & surfaces, surf_opt, verbose, test, output_files, surface_quantities
    ! Set default values of parameters
    tag = ''
    x3_coord = 'zeta'
    norm_type = 'minor_r'
    outdir = ''
    
    n_field_lines = 1
    n_parallel_pts = 128
    x2_center = 0.0
    x3_center = 0.0  

    n_field_periods = -1.0 ! Number of field periods to calculate
    n_pol = -1
   
    surfaces(1) = 0.5
    surf_opt = 0 
    verbose = .false.
    test = .false.
    output_files(1) = 'pest'
    surface_quantities = "stop" 

    ! Add line to ensure iunit is already not open
    open(newunit=iunit,file=trim(filename),status='old',action='read')
    read(iunit,NML=parameters)
    close(iunit)

    do j = 1,size(surfaces)
      if (isnan(surfaces(j))) then
        exit
      else
        if (surfaces(j) .lt. 1e-8) then 
          exit
        else if (surfaces(j) .gt. 1.0) then
          exit
        end if
      end if
    end do
    n_surf = j-1

    j = 0

    block
      character(:), allocatable :: temp_str
      end_found = .false.
      do while(end_found .eqv. .false.)
        j = j + 1
        temp_str = trim(surface_quantities(j))
        if (temp_str == 'stop') then
          end_found = .true.
        end if
      end do
    end block

    n_surface_quantities = j-1
    write(*,*) n_surface_quantities

    geom_id => geom_file(1:len(trim(geom_file)))

  end subroutine

  subroutine write_pest_file(pest,idx1,n_field_periods_final)
    implicit none
    integer, intent(in) :: idx1
    type(PEST_Obj), intent(in) :: pest
    real(dp), intent(in) :: n_field_periods_final
    integer :: j, k, iunit
    real(dp) :: prefac
    character(len=2000) :: filename
    character(len=2000) :: filenumber
    character(len=6) :: x3_string

    prefac = 1.0
    if (trim(pest%x3_coord) .eq. 'theta') then
      prefac = pest%iota(idx1)
    end if
    x3_string = '#'//trim(pest%x3_coord)

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
    write (iunit,'(A,F8.3)') 'n_pol = ', pest%iota(idx1)*n_field_periods_final/pest%vmec%nfp
    write (iunit,'(A,F8.3)') 'n_tor = ', n_field_periods_final/pest%vmec%nfp
    write (iunit,'(A)') '/'
    write (iunit,'(9(A23,2x))') x3_string,'g11','g12','g22','g13','g23','g33','|B|','sqrt(g)', 'K2', 'K1', 'dBdz' 
    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        write(iunit,'(12(ES23.12E3,2x))') prefac*pest%x3(k,idx1)/pi, pest%g11(j,k,idx1), pest%g12(j,k,idx1), pest%g22(j,k,idx1), &
          & pest%g13(j,k,idx1), pest%g23(j,k,idx1), pest%g33(j,k,idx1), pest%bmag(j,k,idx1), pest%jac(j,k,idx1), pest%curv_drift_x1(j,k,idx1), pest%curv_drift_x2(j,k,idx1), pest%d_B_d_x3(j,k,idx1) 
      end do
    end do
    close(iunit)

  end subroutine

  subroutine write_gene_geometry_file(pest,idx1,n_field_periods_final)
    implicit none
    integer, intent(in) :: idx1
    type(PEST_Obj), intent(in) :: pest
    real(dp), intent(in) :: n_field_periods_final
    integer :: j, k, iunit
    character(len=2000) :: filename
    character(len=2000) :: filenumber

    write(filenumber,"(I0.3)") idx1
    filename = trim(outdir)//"gene_"//trim(tag)//"_surf_"//trim(filenumber)//".dat"
    open(newunit=iunit,file=trim(filename))
    write (iunit,'(A)') '&parameters'
    write (iunit,'(A,F12.7)') '!s0 = ', pest%x1(idx1)
    write (iunit,'(A,F12.7)') '!alpha0 = ', pest%x2(pest%ix21)
    write (iunit,'(A,F12.7)') '!minor_r = ', pest%L_ref
    write (iunit,'(A,F12.7)') '!Bref = ', pest%B_ref
    write (iunit,'(A,F12.7)') 'q0 = ', pest%safety_factor_q(idx1)
    write (iunit,'(A,F12.7)') 'shat = ', pest%shat(idx1)
    write (iunit,'(A,I6)') 'gridpoints = ', pest%nx3-1
    write (iunit,'(A,I6)') 'n_pol = ', ceiling(pest%iota(idx1)*n_field_periods_final/pest%vmec%nfp)
    write (iunit,'(A)') '/'
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
    character(len=2000) :: filename_cyl, filenumber
    write(filenumber,"(I0.3)") idx1

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
        write (iunit_cyl,'(3(F12.7,2x))') pest%Rsurf(j,k,idx1), pest%Zsurf(j,k,idx1), -pest%x3(k,idx1)
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
    character(len=2000) :: filename_cyl, filenumber
    write(filenumber,"(I0.3)") idx1


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
        write (iunit_cyl,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1), pest%Zsurf(j,k,idx1), pest%x3(k,idx1), surf_data(j,k)
      end do
      write (iunit_cyl,'(A)') " "
      
      write (iunit_cyl,'(4(F12.7,2x))') pest%Rsurf(pest%ix21,k,idx1), pest%Zsurf(pest%ix21,k,idx1), pest%x3(k,idx1), surf_data(pest%ix21,k)
    end do
    close(iunit_cyl)
  end subroutine

  subroutine write_surface_quantity_xyz(pest,idx1,data_name,surf_data)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    character(len=32), intent(in) :: data_name
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32), intent(in) :: surf_data
    integer :: j, k, iunit_xyz
    character(len=2000) :: filename_xyz, filenumber
    write(filenumber,"(I0.3)") idx1

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
        write (iunit_xyz,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1)*cos(pest%x3(k,idx1)), pest%Rsurf(j,k,idx1)*sin(pest%x3(k,idx1)), pest%Zsurf(j,k,idx1), surf_data(j,k)
      end do
      j=pest%ix21
      write (iunit_xyz,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1)*cos(pest%x3(k,idx1)), pest%Rsurf(j,k,idx1)*sin(pest%x3(k,idx1)), pest%Zsurf(j,k,idx1), surf_data(j,k)
      write (iunit_xyz,'(A)') " "
    end do
    k=pest%ix31
    do j=pest%ix21,pest%ix22
      write (iunit_xyz,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1)*cos(pest%x3(k,idx1)), pest%Rsurf(j,k,idx1)*sin(pest%x3(k,idx1)), pest%Zsurf(j,k,idx1), surf_data(j,k)
    end do
    j=pest%ix21
    write (iunit_xyz,'(4(F12.7,2x))') pest%Rsurf(j,k,idx1)*cos(pest%x3(k,idx1)), pest%Rsurf(j,k,idx1)*sin(pest%x3(k,idx1)), pest%Zsurf(j,k,idx1), surf_data(j,k)
    write (iunit_xyz,'(A)') " "

    close(iunit_xyz)
  end subroutine

  subroutine write_surface_quantity_theta_zeta(pest,idx1,data_name,surf_data)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    character(len=32), intent(in) :: data_name
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32), intent(in) :: surf_data
    integer :: i, j, jp1, k, iunit_tz, theta_index
    real(dp) :: surf_interp, theta, theta_j, theta_jp1, theta_interp, dt1, dt2, delta
    character(len=2000) :: filename_tz, filenumber
    real(dp), parameter :: pi2 = pi+pi
    real(dp), parameter :: eps = 1e-8

    write(filenumber,"(I0.3)") idx1


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
      if (pest%nx2 .gt. 1) then
        do j=pest%ix21,pest%ix22
          if (j .lt. pest%ix22) then 
            jp1 = j+1 
          else 
            jp1 = pest%ix21
          end if
          theta = 0.
          theta_j = pest%x2(j) + pest%iota(idx1)*pest%x3(k,idx1)
          theta_jp1 = pest%x2(jp1) + pest%iota(idx1)*pest%x3(k,idx1)
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
          write (iunit_tz,'(2(I5,2x,F12.7,2x),F12.7)') theta_index+1, theta_interp, k+pest%nx3/2+1, pest%x3(k,idx1), surf_interp
        end do
      else
        j = 1
        theta_index = 0
        theta_interp = pest%x2(j) + pest%iota(idx1)*pest%x3(k,idx1)

        write (iunit_tz,'(2(I5,2x,F12.7,2x),F12.7)') theta_index+1, theta_interp, k+pest%nx3/2+1, pest%x3(k,idx1),surf_data(j,k)
      end if
      write (iunit_tz,'(A)') " "
    end do

    close(iunit_tz)
  end subroutine

  subroutine write_RZ_theta_zeta_grid(pest,idx1,nfpi)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1, nfpi
    real(dp) :: theta, theta_j, theta_jp1, theta_interp, dt1, dt2, delta, Rsurf_interp, Zsurf_interp
    integer i, j, jp1, k, iunit_R, iunit_Z, theta_index
    character(len=2000) :: filename_R, filename_Z, filenumber
    real(dp), parameter :: pi2 = pi*pi
    real(dp), parameter :: eps = 1e-8

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
        theta_j = pest%x2(j) + pest%iota(idx1)*pest%x3(k,idx1)
        theta_jp1 = pest%x2(jp1) + pest%iota(idx1)*pest%x3(k,idx1)
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

  subroutine create_mapping_file_hdf5(filename,pest,file_id)
    character(len=*), intent(in) :: filename
    type(PEST_Obj), intent(in) :: pest
    integer(hid_t), intent(out) :: file_id

    integer(hid_t) :: space_id, data_id
    integer :: ierr
    integer, parameter :: one = 1
    integer(hsize_t), dimension(1) :: one_size = (/ 1 /)

    call h5open_f(ierr)

    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, ierr)

    call h5screate_simple_f(one, one_size, space_id, ierr)
    call h5dcreate_f(file_id, "B_ref", H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, pest%B_ref, one_size, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "Major R", H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, pest%major_R, one_size, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "Minor r", H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, pest%minor_r, one_size, ierr)
    call h5dclose_f(data_id, ierr)
    call h5sclose_f(space_id, ierr)

  end subroutine

  subroutine close_mapping_file_hdf5(file_id)
    integer(hid_t), intent(in) :: file_id
    integer :: ierr

    call h5fclose_f(file_id,ierr)
  end subroutine

  subroutine write_mapping_file_info_hdf5(file_id,pest,idx1)
    integer(hid_t), intent(in) :: file_id
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32) :: theta_coords
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32) :: zeta_coords
    integer, dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32) :: theta_indices
    integer, dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32) :: zeta_indices
    integer :: i, j, jp1, k, iunit_tz, theta_index, ierr
    real(dp) :: theta, theta_j, theta_jp1, theta_interp
    real(dp), parameter :: pi2 = pi+pi
    real(dp), parameter :: eps = 1e-8

    integer(hid_t) :: data_id, space_id
    integer, parameter :: rank = 2
    integer(hsize_t), dimension(2) :: dims
    integer, parameter :: one = 1
    integer(hsize_t), dimension(1) :: one_size = (/ 1 /)

    call h5screate_simple_f(one, one_size, space_id, ierr)
    call h5dcreate_f(file_id, "s0", H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, pest%x1(idx1), one_size, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "iota", H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, pest%iota(idx1), one_size, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "shat", H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, pest%shat(idx1), one_size, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "M_theta", H5T_NATIVE_INTEGER, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_INTEGER, pest%nx2, one_size, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "N_zeta", H5T_NATIVE_INTEGER, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_INTEGER, pest%nx3-1, one_size, ierr)
    call h5dclose_f(data_id, ierr)
    call h5sclose_f(space_id, ierr)

    dims(1) = pest%nx2
    dims(2) = pest%nx3-1

    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        if (j .lt. pest%ix22) then 
          jp1 = j+1 
        else 
          jp1 = pest%ix21
        end if
        theta = 0.
        theta_j = pest%x2(j) + pest%iota(idx1)*pest%x3(k,idx1)
        theta_jp1 = pest%x2(jp1) + pest%iota(idx1)*pest%x3(k,idx1)
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
        theta_coords(j,k) = theta_interp
        zeta_coords(j,k) = pest%x3(k,idx1)
        theta_indices(j,k) = theta_index+1
        zeta_indices(j,k) = k+pest%nx3/2+1
      end do
    end do

    call h5screate_simple_f(rank,dims,space_id,ierr) 
   
    call h5dcreate_f(file_id, "Theta", H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, theta_coords, dims, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "Zeta", H5T_NATIVE_DOUBLE, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, zeta_coords, dims, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "Theta Indices", H5T_NATIVE_INTEGER, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_INTEGER, theta_indices, dims, ierr)
    call h5dclose_f(data_id, ierr)

    call h5dcreate_f(file_id, "Zeta Indices", H5T_NATIVE_INTEGER, space_id, data_id, ierr)
    call h5dwrite_f(data_id, H5T_NATIVE_INTEGER, zeta_indices, dims, ierr)
    call h5dclose_f(data_id, ierr)

    call h5sclose_f(space_id,ierr)

  end subroutine


  subroutine write_mapping_file_surf_data_hdf5(file_id,pest,idx1,data_name,surf_data)
    integer(hid_t), intent(in) :: file_id
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    character(len=32), intent(in) :: data_name
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32), intent(in) :: surf_data
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32) :: surf_interp
    integer :: i, j, jp1, k, iunit_tz, theta_index, ierr
    real(dp) :: theta, theta_j, theta_jp1, theta_interp, dt1, dt2, delta
    real(dp), parameter :: pi2 = pi+pi
    real(dp), parameter :: eps = 1e-8

    integer(hid_t) :: data_id, space_id
    integer, parameter :: rank = 2
    integer(hsize_t), dimension(2) :: dims

    dims(1) = pest%nx2
    dims(2) = pest%nx3-1
    do k=pest%ix31,pest%ix32-1
      do j=pest%ix21,pest%ix22
        if (j .lt. pest%ix22) then 
          jp1 = j+1 
        else 
          jp1 = pest%ix21
        end if
        theta = 0.
        theta_j = pest%x2(j) + pest%iota(idx1)*pest%x3(k,idx1)
        theta_jp1 = pest%x2(jp1) + pest%iota(idx1)*pest%x3(k,idx1)
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
        surf_interp(j,k) = real(pest%nx2)/pi2*(dt2*surf_data(j,k) + dt1*surf_data(jp1,k))

      end do
    end do
    call h5screate_simple_f(rank,dims,space_id,ierr) 
   
    call h5dcreate_f(file_id, trim(data_name), H5T_NATIVE_DOUBLE, space_id, data_id, ierr)

    call h5dwrite_f(data_id, H5T_NATIVE_DOUBLE, surf_interp, dims, ierr)

    call h5dclose_f(data_id, ierr)

    call h5sclose_f(space_id,ierr)
  end subroutine

end module
