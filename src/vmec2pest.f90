!******************************************************************************
! This is the top-level program to convert a VMEC equilibrium to a PEST
! coordinate system.
! Written by B.J. Faber, University of Wisconsin-Madison
! Created June 2019
!******************************************************************************

program vmec2sfl
  use types, only: dp, pi
  use pest_object, only: PEST_Obj, create_PEST_Obj, destroy_PEST_Obj, set_PEST_reference_values, &
    & get_PEST_data
  use io_core, only: read_vmec2pest_input, write_pest_file, geom_file, surfaces, surf_opt, &
    & x3_center, n_field_lines, n_parallel_pts, n_field_periods, x3_coord, norm_type, output_files, &
    & write_gene_geometry_file, write_cylindrical_surface, surface_quantities, n_surface_quantities, &
    & write_surface_quantity_cyl, write_surface_quantity_xyz, geom_id
  use compute_pest, only: compute_pest_geometry
  use normalizations, only: set_normalizations
  implicit none

  real(dp), parameter :: pi2 = pi*pi
  integer :: i, j
  character(len=2000), target :: infile
  character(len=8) :: grid_type

  real(dp) :: time1, time2
  !real(dp), dimension(:), allocatable :: line_data
  real(dp), dimension(:,:), allocatable :: surf_data
  !real(dp), dimension(:,:,:), allocatable :: vol_data
 
  type(PEST_Obj) :: pest

  call cpu_time(time1)
  infile = 'vmec2pest.inp'
  call read_vmec2pest_input(infile)
print *,geom_id(1:50)
  pest = create_PEST_Obj(geom_id,surfaces,n_field_lines,n_parallel_pts)
  call set_PEST_reference_values(pest,norm_type)
  pest%x3_coord = x3_coord
  call compute_pest_geometry(pest,x3_center,n_field_periods,surf_opt)
  grid_type = 'gene'
  call set_normalizations(pest,grid_type) 
  do i=1,4 
    select case (trim(output_files(i)))
      case('pest')
        do j=pest%ix11,pest%ix12
          call write_pest_file(pest,j)
        end do
      case('gene')
        do j=pest%ix11,pest%ix12
          call write_gene_geometry_file(pest,j)
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
    allocate(surf_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32))
    do i=1,n_surface_quantities
      do j=pest%ix11,pest%ix12
!        call get_PEST_data(pest,j,surface_quantities(i),surf_data)
!        call write_surface_quantity_cyl(pest,j,surface_quantities(i),surf_data)
!        call write_surface_quantity_xyz(pest,j,surface_quantities(i),surf_data)
      end do
    end do
    deallocate(surf_data)
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
