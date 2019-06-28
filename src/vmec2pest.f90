!******************************************************************************
! This is the top-level program to convert a VMEC equilibrium to a PEST
! coordinate system.
! Written by B.J. Faber, University of Wisconsin-Madison
! Created June 2019
!******************************************************************************

program vmec2sfl
  use types, only: dp, pi
  use pest_object, only: PEST_Obj, create_PEST_Obj, destroy_PEST_Obj, set_PEST_reference_values, &
    & get_PEST_surface_data
  use io_core, only: read_input, write_pest_file, geom_file, surfaces, surf_opt, &
    & x3_center, n_field_lines, n_parallel_pts, n_field_periods, x3_coord, norm_type, output_files, &
    & write_gene_geometry_file, write_cylindrical_surface, surface_quantities, n_surface_quantities, &
    & write_surface_quantity
  use compute_pest, only: compute_pest_geometry
  use normalizations, only: gene_normalizations
  implicit none

  real(dp), parameter :: pi2 = pi*pi
  integer :: i, j
  character(len=2000) :: infile

  real(dp) :: time1, time2
  real(dp), dimension(:,:), allocatable :: surf_data
 
  type(PEST_Obj) :: pest

  call cpu_time(time1)
  infile = 'vmec2pest.inp'
  call read_input(infile)
  pest = create_PEST_Obj(geom_file,surfaces,n_field_lines,n_parallel_pts)
print *, surface_quantities
  call set_PEST_reference_values(pest,norm_type)
  pest%x3_coord = x3_coord
  allocate(surf_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32))
  call compute_pest_geometry(pest,x3_center,n_field_periods,surf_opt)
  call gene_normalizations(pest) 
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
    do i=1,n_surface_quantities
      do j=pest%ix11,pest%ix12
        call get_PEST_surface_data(pest,j,surface_quantities(i),surf_data)
        call write_surface_quantity(pest,j,pest%nx2,pest%nx3,surface_quantities(i),surf_data)
      end do
    end do
  end if 
  call destroy_PEST_Obj(pest)
  call cpu_time(time2)
  write(6,"(A,F8.4,A)") "vmec2sfl completed in ",time2-time1," seconds"

end program
