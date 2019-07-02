!******************************************************************************
! This module contains routines that allow other codes to interface with 
! routines in VMECTools.
! Author: B.J. Faber, University of Wisconsin-Madsion (bfaber@wisc.edu)
! Creation date: June 2019
!******************************************************************************

module interfaces
  use types, only: dp
  use pest_object, only: PEST_Obj, create_PEST_Obj, set_PEST_reference_values, get_PEST_data
  use compute_pest, only: compute_pest_geometry
  use normalizations, only: set_normalizations 
  implicit none

  public :: pest2vmec_c_interface, pest2vmec_stellopt_interface

  private
    type(PEST_Obj) :: pest ! Make this pest_object shared within the module so data can be read using get_PEST_data
contains

  !****************************************************************************
  ! Interface to call pest2vmec from C/C++
  !****************************************************************************
  subroutine pest2vmec_c_interface(nx1,nx2,nx3,nfpi,x1,vmec_file,ncvf,grid_type,ncgt) bind(C,name='pest2vmec_c_interface')
    use, intrinsic :: iso_c_binding
    integer(c_int), intent(in) :: nx1, nx2, nx3, ncvf, ncgt
    real(c_double), intent(in) :: nfpi
    real(c_double), dimension(nx1) :: x1 
    character(c_char), dimension(ncvf), intent(in) :: vmec_file
    character(c_char), dimension(ncgt), intent(in) :: grid_type
    character(len=2000) :: vmec_string
    character(len=8) :: grid_string
    character(len=7) :: norm_type
    character(len=5) :: x3_coord
    integer :: i, surf_opt
    surf_opt = 0

print *, nx1, nx2, nx3, nfpi
print *, x1
print *, grid_type(1:4)
    vmec_string = ""
    do i = 1,ncvf
      vmec_string = trim(vmec_string)//trim(vmec_file(i))
    end do
print *, trim(vmec_string)

    grid_string = ""
    do i = 1,ncgt
      grid_string = trim(grid_string)//trim(grid_type(i))
    end do

    pest = create_PEST_Obj(trim(vmec_string),x1,nx2,nx3)
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

  subroutine get_pest_data_c_interface(x1,x2,data_name,ncdn,n_dims,data_size,c_data) bind(C)
    use, intrinsic :: iso_c_binding
    integer(c_int), intent(in) :: n_dims, data_size, x1, x2, ncdn
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
print *, data_string
        if (data_size .ne. pest%nx3) then
          print *, "Error! C array does not have the same size as ",pest%nx3,"!"
          stop
        end if
        allocate(line_data(pest%ix31:pest%ix32))
        call get_PEST_data(pest,x1,x2,trim(data_string),line_data)
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

  subroutine pest2vmec_stellopt_interface
  end subroutine

end module
