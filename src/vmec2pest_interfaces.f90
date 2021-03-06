!******************************************************************************
! This module contains routines that allow other codes to interface with 
! routines in VMECTools.
! Author: B.J. Faber, University of Wisconsin-Madsion (bfaber@wisc.edu)
! Creation date: June 2019
!******************************************************************************

module vmec2pest_interfaces
  use, intrinsic :: iso_c_binding
  use types, only: dp
  use pest_object, only: PEST_Obj, create_PEST_Obj, set_PEST_reference_values, get_PEST_data
  use compute_pest, only: compute_pest_geometry
  use normalizations, only: set_normalizations 
  implicit none

  public :: vmec2pest_interface, vmec2pest_c_interface, vmec2pest_stellopt_interface, get_pest_data_interface

  private
    type(PEST_Obj) :: pest ! Make this pest_object shared within the module so data can be read using get_PEST_data

    type, bind(C) :: vmec2pest_c_options
      ! Currently this must be a NetCDF file
      type(c_ptr) :: geom_file
      type(c_ptr) :: grid_type
      type(c_ptr) :: x3_coord
      type(c_ptr) :: norm_type
      integer(c_int) :: nx1, nx2, nx3
      type(c_ptr) :: x1
      real(c_double) :: x2_center, x3_center, nfpi
      integer(c_int) :: npol
    end type

  interface c2f
    module procedure c2f_string_array_1d
    module procedure c2f_real_array_1d
  end interface 

contains
  !****************************************************************************
  ! Interface to call vmec2pest from C/C++
  !****************************************************************************
  subroutine vmec2pest_c_interface(options) bind(C,name='vmec2pest_c_interface')
    use, intrinsic :: iso_c_binding
    type(vmec2pest_c_options), intent(in) :: options
    character(len=:), allocatable, target :: geom_file
    character(len=:), pointer :: geom_id
    real(c_double), dimension(:), allocatable :: surfaces
    integer :: surf_opt, i
    character(len=:), allocatable :: grid_type, norm_type, x3_coord
    real(dp) :: nfpif
    surf_opt = 0

    call c2f(options%geom_file,geom_file)
    geom_id => geom_file(1:len(geom_file))
    call c2f(options%x1,surfaces,options%nx1)
    pest = create_PEST_Obj(geom_id,surfaces,options%nx1,options%nx2,options%nx3)
 
    call c2f(options%grid_type,grid_type)
    call c2f(options%x3_coord,x3_coord)
    call c2f(options%norm_type,norm_type)

    call set_PEST_reference_values(pest,norm_type)
    pest%x3_coord = x3_coord
    if (trim(pest%x3_coord) .eq. 'theta') then
      pest%x3_max_interval = options%npol
    else
      pest%x3_max_interval = options%nfpi
    endif
    call compute_pest_geometry(pest,options%x2_center,options%x3_center,nfpif,surf_opt)
    call set_normalizations(pest,trim(grid_type))

    deallocate(geom_file,grid_type,norm_type,x3_coord,surfaces)
    
  end subroutine

  subroutine get_pest_data_c_interface(x1,x2,data_name,n_dims,data_size,c_data) bind(C,name="get_pest_data_c_interface")
    use, intrinsic :: iso_c_binding
    integer(c_int), intent(in) :: n_dims, data_size, x1, x2
    type(c_ptr), intent(in) :: data_name
    real(c_double), dimension(data_size), intent(out) :: c_data
    real(dp) :: vmec_data, rad_data
    real(dp), dimension(:), allocatable :: line_data
    real(dp), dimension(:,:), allocatable :: surf_data
    real(dp), dimension(:,:,:), allocatable :: vol_data 
    character(len=:), allocatable :: data_string
    integer idx, idx1, idx2, idx3, str_len;

    call c2f(data_name,data_string)

    ! TODO: implement bounds checking on x1 and x2 to ensure they aren't outside the dimension of
    ! the PEST object
    select case(n_dims)
      case(-1)
        if (data_size .ne. 1) then
          print *, "VMEC2PEST Error! C array does not have the same size as 1!"
          stop
        endif
        call get_PEST_data(pest,trim(data_string),vmec_data)
        c_data(1) = vmec_data
      case(0)
        if (data_size .ne. 1) then
          print *, "VMEC2PEST Error! C array does not have the same size as 1!"
          stop
        endif
        call get_PEST_data(pest,x1,trim(data_string),rad_data)
        c_data(1) = rad_data
      case(1)
        if (data_size .ne. pest%nx3) then
          print *, "VMEC2PEST Error! C array does not have the same size as ",pest%nx3,"!"
          stop
        end if
        allocate(line_data(pest%ix31:pest%ix32))
        call get_PEST_data(pest,x1,x2,trim(data_string),line_data)
        do idx3 = pest%ix31,pest%ix32
          idx = idx3 - pest%ix31 + 1
          c_data(idx) = line_data(idx3)
        end do
        deallocate(line_data)
      case(2)
        if (data_size .ne. pest%nx2*pest%nx3) then
          print *, "VMEC2PEST Error! C array does not have the same size as ",pest%nx2*pest%nx3,"!"
          stop
        end if
        allocate(surf_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32))
        call get_PEST_data(pest,x1,trim(data_string),surf_data)
        do idx3 = pest%ix31,pest%ix32
          do idx2 = pest%ix21,pest%ix22
            idx = pest%nx2*(idx3 - pest%ix31) + idx2 - pest%ix21 + 1 
            c_data(idx) = surf_data(idx2,idx3)
          end do
        end do
        deallocate(surf_data)
      case(3)
        if (data_size .ne. pest%nx1*pest%nx2*pest%nx3) then
          print *, "VMEC2PEST Error! C array does not have the same size as ",pest%nx1*pest%nx2*pest%nx3,"!"
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
        print *, "VMEC2PEST Error! n_dims must be specified to 1, 2, or 3!"
        stop
    end select

    deallocate(data_string)

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
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Interface for calling vmec2pest from STELLOPT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine vmec2pest_stellopt_interface(surfaces,nx2,nx3,x2_center,x3_center,x3_coord,nfpi,npol,norm_type,grid_type)
    character(*), intent(in) :: x3_coord, norm_type, grid_type
    real(dp), dimension(:), intent(in) :: surfaces
    real(dp), intent(in) :: x2_center, x3_center, nfpi
    integer, intent(in) :: nx2, nx3, npol
    character(len=:), pointer :: geom_id
    real(dp) :: nfpif
    
    geom_id => NULL()
    pest = create_PEST_Obj(geom_id,surfaces,size(surfaces),nx2,nx3)
    call set_PEST_reference_values(pest,norm_type)
    pest%x3_coord = x3_coord
    pest%x3_coord = x3_coord
    if (trim(pest%x3_coord) .eq. 'theta') then
      pest%x3_max_interval = real(npol)
    else
      pest%x3_max_interval = nfpi
    end if
    call compute_pest_geometry(pest,x2_center,x3_center,nfpif,0)
    call set_normalizations(pest,grid_type) 
    
  end subroutine

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Interface for calling vmec2pest from Fortran
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine vmec2pest_interface(geom_file,surfaces,nx2,nx3,x2_center,x3_center,x3_coord,nfpi,npol,norm_type,grid_type)
    character(*), intent(in), target :: geom_file
    character(*), intent(in) :: x3_coord, norm_type, grid_type
    real(dp), dimension(:), intent(in) :: surfaces
    real(dp), intent(in) :: x2_center, x3_center, nfpi
    integer, intent(in) :: nx2, nx3, npol
    character(len=:), pointer :: geom_id
    real(dp) :: nfpif
    
    geom_id => geom_file
    pest = create_PEST_Obj(geom_id,surfaces,size(surfaces),nx2,nx3)
    call set_PEST_reference_values(pest,norm_type)
    pest%x3_coord = x3_coord
    if (trim(pest%x3_coord) .eq. 'theta') then
      pest%x3_max_interval = real(npol)
    else
      pest%x3_max_interval = nfpi
    end if
    call compute_pest_geometry(pest,x2_center,x3_center,nfpif,0)
    call set_normalizations(pest,grid_type) 
    

  end subroutine

  subroutine get_pest_data_interface(x1,x2,data_string,n_dims,data_size,data_arr)
    integer, intent(in) :: n_dims, data_size, x1, x2
    character(*), intent(in) :: data_string
    real(dp), dimension(data_size), intent(out) :: data_arr
    real(dp) :: rad_data, vmec_data
    real(dp), dimension(:), allocatable :: line_data
    real(dp), dimension(:,:), allocatable :: surf_data
    real(dp), dimension(:,:,:), allocatable :: vol_data 
    integer idx, idx1, idx2, idx3


    ! TODO: implement bounds checking on x1 and x2 to ensure they aren't outside the dimension of
    ! the PEST object
    select case(n_dims)
      case(-1)
        if (data_size .ne. 1) then
          print *, "VMEC2PEST Error! Array does not have the same size as 1!"
          stop
        endif
        call get_PEST_data(pest,trim(data_string),vmec_data)
        data_arr(1) = vmec_data
      case(0)
        if (data_size .ne. 1) then
          print *, "VMEC2PEST Error! Array does not have the same size as 1!"
          stop
        endif
        call get_PEST_data(pest,x1,trim(data_string),rad_data)
        data_arr(1) = rad_data
      case(1)
        if (data_size .ne. pest%nx3) then
          print *, "VMEC2PEST Error! Array does not have the same size as ",pest%nx3,"!"
          stop
        end if
        allocate(line_data(pest%ix31:pest%ix32))
        call get_PEST_data(pest,x1,x2,trim(data_string),line_data)
        do idx3 = pest%ix31,pest%ix32
          idx = idx3 - pest%ix31 + 1
          data_arr(idx) = line_data(idx3)
        end do
        deallocate(line_data)
      case(2)
        if (data_size .ne. pest%nx2*pest%nx3) then
          print *, "VMEC2PEST Error! Array does not have the same size as ",pest%nx2*pest%nx3,"!"
          stop
        end if
        allocate(surf_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32))
        call get_PEST_data(pest,x1,trim(data_string),surf_data)
        do idx3 = pest%ix31,pest%ix32
          do idx2 = pest%ix21,pest%ix22
            idx = pest%nx2*(idx3 - pest%ix31) + idx2 - pest%ix21 + 1 
            data_arr(idx) = surf_data(idx2,idx3)
          end do
        end do
        deallocate(surf_data)
      case(3)
        if (data_size .ne. pest%nx1*pest%nx2*pest%nx3) then
          print *, "VMEC2PEST Error! Array does not have the same size as ",pest%nx1*pest%nx2*pest%nx3,"!"
          stop
        end if
        allocate(vol_data(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12))
        call get_PEST_data(pest,trim(data_string),vol_data)
        do idx1 = pest%ix11,pest%ix12
          do idx3 = pest%ix31,pest%ix32
            do idx2 = pest%ix21,pest%ix22
              idx = pest%nx2*pest%nx3*(idx1 - pest%ix11) + pest%nx2*(idx3 - pest%ix31) + idx2 - pest%ix21 + 1 
              data_arr(idx) = vol_data(idx2,idx3,idx1)
            end do
          end do
        end do
        deallocate(vol_data)
      case default
        print *, "VMEC2PEST Error! n_dims must one of the following:"
        print *, "n_dims = -1 - VMEC data"
        print *, "n_dims = 0 - radial data"
        print *, "n_dims = 1 - field line data"
        print *, "n_dims = 2 - flux surface data"
        print *, "n_dims = 3 - volume data"
        stop
    end select

  end subroutine

end module
