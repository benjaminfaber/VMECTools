!******************************************************************************
! This is the top-level program to convert a VMEC equilibrium to a PEST
! coordinate system.
! Written by B.J. Faber, University of Wisconsin-Madison
! Created June 2019
!******************************************************************************

program vmec2pest
  use hdf5
  use types, only: dp, pi
  use pest_object, only: PEST_Obj, create_PEST_Obj, destroy_PEST_Obj, set_PEST_reference_values, &
    & get_PEST_data
  use io_core, only: read_vmec2pest_input, write_pest_file, geom_file, surfaces, n_surf, surf_opt, &
    & x2_center, x3_center, n_field_lines, n_parallel_pts, n_field_periods, n_pol, x3_coord, norm_type, output_files, &
    & write_gene_geometry_file, write_cylindrical_surface, surface_quantities, n_surface_quantities, &
    & write_surface_quantity_cyl, write_surface_quantity_xyz, write_surface_quantity_theta_zeta, geom_id, &
    & create_mapping_file_hdf5, close_mapping_file_hdf5, write_mapping_file_info_hdf5, write_mapping_file_surf_data_hdf5
  use compute_pest, only: compute_pest_geometry
  use normalizations, only: set_normalizations
  implicit none

  real(dp), parameter :: pi2 = pi*pi
  integer :: i, j, arg_count
  character(len=2000), target :: infile
  character(len=2000) :: ext_file
  character(len=8) :: grid_type

  real(dp) :: time1, time2, n_field_periods_final
  !real(dp), dimension(:), allocatable :: line_data
  real(dp), dimension(:,:), allocatable :: surf_data
  !real(dp), dimension(:,:,:), allocatable :: vol_data
 
  type(PEST_Obj) :: pest

  integer(hid_t) :: fid

  call cpu_time(time1)

  arg_count = command_argument_count()
  if (arg_count .eq. 1) then
    call get_command_argument(1,ext_file)
    infile = trim(ext_file)
  else
    infile = 'vmec2pest.inp'
  end if
  call read_vmec2pest_input(infile)

  pest = create_PEST_Obj(geom_id,surfaces,n_surf,n_field_lines,n_parallel_pts)
  call set_PEST_reference_values(pest,norm_type)
  pest%x3_coord = x3_coord
  if (trim(pest%x3_coord) .eq. 'theta') then
    pest%x3_max_interval = real(n_pol)
  else
    pest%x3_max_interval = n_field_periods
  end if
  call compute_pest_geometry(pest,x2_center,x3_center,n_field_periods_final,surf_opt)
  grid_type = 'gene'
  call set_normalizations(pest,grid_type) 
  do i=1,4 
    select case (trim(output_files(i)))
      case('pest')
        do j=pest%ix11,pest%ix12
          call write_pest_file(pest,j,n_field_periods_final)
        end do
      case('gene')
        do j=pest%ix11,pest%ix12
          call write_gene_geometry_file(pest,j,n_field_periods_final)
        end do
      case('surf')
        do j=pest%ix11,pest%ix12
          call write_cylindrical_surface(pest,j)
        end do
      case('gs2')
        print *, "Error! gs2 geometry options not yet implemented"
        stop
      case default
    end select 
  end do
  if (n_surface_quantities .gt. 0) then
    call create_mapping_file_hdf5("test.h5",pest,fid)
    allocate(surf_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32))
    do j=pest%ix11,pest%ix12
      call write_mapping_file_info_hdf5(fid,pest,j)
    end do
    do i=1,n_surface_quantities
      if (trim(surface_quantities(i)) .ne. "") then
        do j=pest%ix11,pest%ix12
          call get_PEST_data(pest,j,surface_quantities(i),surf_data)
          call write_surface_quantity_cyl(pest,j,surface_quantities(i),surf_data)
          call write_surface_quantity_xyz(pest,j,surface_quantities(i),surf_data)
          call write_surface_quantity_theta_zeta(pest,j,surface_quantities(i),surf_data)
          call write_mapping_file_surf_data_hdf5(fid,pest,j,surface_quantities(i),surf_data)
        end do
      end if
    end do
    deallocate(surf_data)
    call close_mapping_file_hdf5(fid)
  end if 
!  if (n_volume_quantities .gt. 0) then
!    allocate(vol_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
!    do i=1,n_surface_quantities
!      do j=pest%ix11,pest%ix12
!        call get_PEST_data(pest,j,surface_quantities(i),vol_data)
!        call write_surface_quantity(pest,j,surface_quantities(i),vol_data)
!      end do
!    end do
!    deallocate(vol_data)
!  end if

  call destroy_PEST_Obj(pest)
  call cpu_time(time2)
  write(6,"(A,F8.4,A)") "vmec2sfl completed in ",time2-time1," seconds"

end program

