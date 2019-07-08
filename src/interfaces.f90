!******************************************************************************
! This module contains routines that allow other codes to interface with 
! routines in VMECTools.
! Author: B.J. Faber, University of Wisconsin-Madsion (bfaber@wisc.edu)
! Creation date: June 2019
!******************************************************************************

module interfaces
  use, intrinsic :: iso_c_binding
  use types, only: dp
  use pest_object, only: PEST_Obj, create_PEST_Obj, set_PEST_reference_values, get_PEST_data
  use compute_pest, only: compute_pest_geometry
  use normalizations, only: set_normalizations 
  implicit none

  public :: pest2vmec_c_interface, pest2vmec_stellopt_interface, pest2vmec_setup_interface

  private
    type(PEST_Obj) :: pest ! Make this pest_object shared within the module so data can be read using get_PEST_data

    type, bind(C) :: vmec2pest_options
      ! Currently this must be a NetCDF file
      type(c_ptr) :: geom_file
      type(c_ptr) :: grid_type
      type(c_ptr) :: x3_coord
      type(c_ptr) :: norm_type
      integer(c_int) :: nx1, nx2, nx3
      type(c_ptr) :: x1
      real(c_double) :: x3_center, nfpi
    end type

  interface c2f
    module procedure c2f_string_array_1d
    module procedure c2f_real_array_1d
  end interface 

  interface pest2vmec_setup_interface
    module procedure pest2vmec_setup_c_interface
  end interface

contains
  !****************************************************************************
  ! Interface to set up pest2vmec from C/C++
  !****************************************************************************
  subroutine pest2vmec_setup_c_interface(options) bind(C,name='pest2vmec_setup_c_interface')
    use, intrinsic :: iso_c_binding
    type(vmec2pest_options), target :: options
    character(len=:), allocatable, target :: geom_file
    character(len=:), pointer :: geom_id
    real(c_double), dimension(:), allocatable :: surfaces
    integer :: nx2, nx3
    character(len=2000), target :: file_name
    
    nx2 = options%nx2
    nx3 = options%nx3
    call c2f(options%geom_file,geom_file)
print *, len(geom_file)
    geom_id => geom_file(1:len(geom_file))
print *, geom_id
!file_name(1:size(geom_file)) = geom_file(1:size(geom_file))
    call c2f(options%x1,surfaces,options%nx1)
    pest = create_PEST_Obj(geom_id,surfaces,nx2,nx3)
print *, options%nx2, pest%nx2

    deallocate(geom_file,surfaces)
   
  end subroutine

  !****************************************************************************
  ! Interface to call pest2vmec from C/C++
  !****************************************************************************
  subroutine pest2vmec_c_interface(nx1,nx2,nx3,nfpi,x1,vmec_file,grid_type) bind(C,name='pest2vmec_c_interface')
    use, intrinsic :: iso_c_binding
    integer(c_int), intent(in) :: nx1, nx2, nx3
    real(c_double), intent(in) :: nfpi
    real(c_double), dimension(nx1) :: x1 
    type(c_ptr), target :: vmec_file
    type(c_ptr), target :: grid_type
    character(len=2000), pointer :: vmec_ptr
    character(len=2000), pointer :: grid_ptr
    character(len=2000), target :: vmec_string
    character(len=2000), target :: grid_string
    character(len=7) :: norm_type
    character(len=5) :: x3_coord
    integer :: i, surf_opt, str_len
    surf_opt = 0

print *, nx1, nx2, nx3, nfpi
print *, x1
    call c_f_pointer(c_loc(vmec_file),vmec_ptr)
    str_len = index(vmec_ptr,c_null_char)
    vmec_string = trim(vmec_ptr(1:str_len))
    print *, trim(vmec_string)

    call c_f_pointer(c_loc(grid_type),grid_ptr)
    str_len = index(grid_ptr,c_null_char)
    grid_string = trim(grid_ptr(1:str_len))
    print *, trim(grid_string)

!    pest = create_PEST_Obj(trim(vmec_string),x1,nx2,nx3)
    if (trim(grid_string) .eq. "tok") then 
      norm_type = 'major_R'
    else
      norm_type = 'minor_r'
      x3_coord = 'theta'
    end if 

    call set_PEST_reference_values(pest,norm_type)
    pest%x3_coord = x3_coord
    call compute_pest_geometry(pest,0.0,nfpi,surf_opt)
    call set_normalizations(pest,trim(grid_string))
print *, pest%g11(pest%ix21,:,pest%ix11)
    
  end subroutine

  subroutine get_pest_data_c_interface(x1,x2,data_name,n_dims,data_size,c_data) bind(C,name="get_pest_data_c_interface")
    use, intrinsic :: iso_c_binding
    integer(c_int), intent(in) :: n_dims, data_size, x1, x2
    !character, dimension(*) :: data_name
    type(c_ptr), target :: data_name
    real(c_double), dimension(data_size), intent(out) :: c_data
    ! Depending on the size of n_dims, one of these will be allocated
    real(dp), dimension(:), allocatable :: line_data
    real(dp), dimension(:,:), allocatable :: surf_data
    real(dp), dimension(:,:,:), allocatable :: vol_data 
    character(len=1000), pointer :: data_string
    character(len=1000), target :: temp
    integer idx, idx1, idx2, idx3, str_len;

    temp = 'a'
    data_string => temp 
    call c_f_pointer(c_loc(data_name),data_string)
    str_len = index(data_string,c_null_char)-1
    temp = data_string(1:str_len)
print *, str_len, data_string(1:str_len)
print *, data_string
print *, trim(temp)
!    data_string = ""
!    do idx1=1,size(data_name)
!      data_string = trim(data_string)//trim(data_name(idx1))
!    end do


    select case(n_dims)
      case(1)
print *, trim(temp)
        if (data_size .ne. pest%nx3) then
          print *, "Error! C array does not have the same size as ",pest%nx3,"!"
          stop
        end if
        allocate(line_data(pest%ix31:pest%ix32))
        call get_PEST_data(pest,x1,x2,trim(temp),line_data)
        do idx3 = pest%ix31,pest%ix32
          idx = idx3 - pest%ix31 + 1
print *, idx, pest%bmag(x2,idx3,x1), pest%g11(x2,idx3,x1), line_data(idx3), pest%L_ref
          c_data(idx) = line_data(idx3)
        end do
        
        deallocate(line_data)
      case(2)
        if (data_size .ne. pest%nx2*pest%nx3) then
          print *, "Error! C array does not have the same size as ",pest%nx2*pest%nx3,"!"
          stop
        end if
        allocate(surf_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32))
        call get_PEST_data(pest,x1,trim(temp),surf_data)
        do idx3 = pest%ix31,pest%ix32
          do idx2 = pest%ix21,pest%ix22
            idx = pest%nx2*(idx3 - pest%ix31) + idx2 - pest%ix21 + 1 
            c_data(idx) = surf_data(idx2,idx3)
          end do
        end do
        deallocate(surf_data)
      case(3)
        if (data_size .ne. pest%nx1*pest%nx2*pest%nx3) then
          print *, "Error! C array does not have the same size as ",pest%nx1*pest%nx2*pest%nx3,"!"
          stop
        end if
        allocate(vol_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
        call get_PEST_data(pest,trim(data_string),vol_data)
        do idx1 = pest%ix11,pest%ix12
          do idx3 = pest%ix31,pest%ix32
            do idx2 = pest%ix21,pest%ix22
              idx = pest%nx2*pest%nx3*(idx1 - pest%ix11) + pest%nx2*(idx3 - pest%ix31) + idx2 - pest%ix21 + 1 
              c_data(idx) = vol_data(idx2,idx3,idx1)
            end do
          end do
        end do
        deallocate(vol_data)
      case default
        print *, "Error! n_dims must be specified to 1, 2, or 3!"
        stop
    end select

  end subroutine

  subroutine c2f_string_array_1d(c_pointer,f_string)
    use, intrinsic :: iso_c_binding
    type(c_ptr), intent(in) :: c_pointer
    character(len=:), allocatable, intent(out) :: f_string
    character, dimension(:), pointer :: f_pointer
    character(len=:), pointer :: f_char_pointer

    integer :: i
    logical :: null_char_found
    i = 0
    null_char_found = .false.
    do while (null_char_found .eqv. .false.)
      i = i + 1
      call c_f_pointer(c_pointer,f_pointer,[i])
      if (f_pointer(i) == c_null_char) then
        null_char_found = .true.
      end if
    end do
    i = i - 1
    
    allocate(character(len=i)::f_string)
    call c_f_pointer(c_loc(f_pointer),f_char_pointer)

    f_string = f_char_pointer(1:i)
print *, len(f_string)

  end subroutine

  subroutine c2f_real_array_1d(c_pointer,f_real,data_size)
    use, intrinsic :: iso_c_binding
    type(c_ptr), intent(in) :: c_pointer
    real(dp), dimension(:), allocatable, intent(out) :: f_real
    integer, intent(in) :: data_size
    real(dp), dimension(:), pointer :: f_pointer

    call c_f_pointer(c_pointer,f_pointer,[data_size])

    allocate(f_real(data_size))
    f_real = f_pointer(1:data_size)
  end subroutine
 
  subroutine pest2vmec_stellopt_interface
  end subroutine

end module
