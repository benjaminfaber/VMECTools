!******************************************************************************
! This module defines the PEST object used in the vmec2pest calculation
! Author: B.J. Faber (bfaber@wisc.edu)
! Creation date: 30.08.2018
!******************************************************************************

module pest_object
  use types, only: dp, pi
  use vmec_object, only: VMEC_Obj, create_VMEC_Obj, destroy_VMEC_Obj
  implicit none

  public PEST_Obj, create_PEST_Obj, destroy_PEST_Obj, set_PEST_reference_values, get_PEST_data

  private

  type :: PEST_Obj
    ! Pointer to the VMEC equilibrium
    type(VMEC_Obj) :: vmec

    ! The rotational transform iota
    real(dp), dimension(:), allocatable :: iota

    ! Safety factor q = 1/iota
    real(dp), dimension(:), allocatable :: safety_factor_q

    ! Magnetic shear shat = (x/q) * (d q / d x) where x = Aminor_p * sqrt(psi_toroidal / psi_{toroidal,edge})
    ! and Aminor_p is the minor radius calculated by VMEC.
    real(dp), dimension(:), allocatable :: shat

    ! L_ref is the reference length used for gs2's normalization, in meters.
    real(dp) :: L_ref

    ! The major radius in meters 
    real(dp) :: major_R

    ! The minor radius in meters
    real(dp) :: minor_r
  
    ! B_ref is the reference magnetic field strength used for gs2's normalization, in Tesla.
    real(dp) :: B_ref

    ! Integers representing the number of points in each SFL direction 
    integer :: nx1, nx2, nx3

    ! Integers representing the starting and ending array indices
    integer :: ix11, ix12, ix21, ix22, ix31, ix32

    ! On exit, s holds the grid points in the first coordinate (flux surface label)
    ! in units of normalized flux
    ! Typically, this is x1 = s = psi_toroidal/psi_edge 
    real(dp), dimension(:), allocatable :: x1

    ! On exit, x2 holds the grid points in the second coordinate (field line label)
    ! Typically, this is x2 = alpha = theta_p - iota * zeta, where zeta is the PEST (geometric) toroidal angle
    ! and theta_p is the PEST poloidal angle
    real(dp), dimension(:), allocatable :: x2

    ! On exit, x3 holds the grid points in the third coordinate (field line following coordinate)
    ! Typically, this is x3 = zeta, where zeta is the PEST (geometric) toroidal angle
    real(dp), dimension(:,:), allocatable :: x3

    ! VMEC supplies the normalized toroidal flux, thus it is natural to use the
    ! toroidal zeta coordinate as the angle parameterizing the field line.
    ! However some codes use the poloidal angle theta, related to zeta by theta =
    ! iota*zeta.  The maximum extent of the flux tube and the coordinate used in
    ! output can be chosen to be either theta or zeta, but the calculation itself
    ! is performed using the zeta coordinate

    character(len=5) :: x3_coord
    real(dp) :: x3_max_interval

    ! Arrays that contain the geometric elements on the surface
    real(dp), dimension(:,:,:), allocatable :: bmag  ! Magnitude of B (|B|)
    real(dp), dimension(:,:,:), allocatable :: jac ! Jacobian (grad s . grad alpha x grad zeta)^(-1)
    real(dp), dimension(:,:,:), allocatable :: g11 ! Metric element gss
    real(dp), dimension(:,:,:), allocatable :: g12 ! Metric element gsa
    real(dp), dimension(:,:,:), allocatable :: g22 ! Metric element gaa
    real(dp), dimension(:,:,:), allocatable :: g13 ! Metric element gsz
    real(dp), dimension(:,:,:), allocatable :: g23 ! Metric element gaz
    real(dp), dimension(:,:,:), allocatable :: g33 ! Metric element gzz
    real(dp), dimension(:,:,:), allocatable :: gradB_drift_x1 ! Derivative of |B| w.r.t. the 1 coordinate
    real(dp), dimension(:,:,:), allocatable :: gradB_drift_x2 ! Derivative of |B| w.r.t. the 2 coordinate
    real(dp), dimension(:,:,:), allocatable :: curv_drift_x1 ! Derivative of |B| w.r.t. the 1 coordinate
    real(dp), dimension(:,:,:), allocatable :: curv_drift_x2 ! Derivative of |B| w.r.t. the 2 coordinate
    real(dp), dimension(:,:,:), allocatable :: d_B_d_x3 ! Derivative of |B| w.t.t. the 3 coordinate

    real(dp), dimension(:,:,:), allocatable :: Rsurf ! The R coordinate of the surfaces
    real(dp), dimension(:,:,:), allocatable :: Zsurf ! The Z coordinate of the surfaces
    real(dp), dimension(:,:,:), allocatable :: d_Lambda_d_theta_vmec ! The Z coordinate of the surfaces

  end type ! end type PEST_Obj


  interface create_PEST_Obj
    module procedure create_from_VMEC
    module procedure create_from_field_line_file
  end interface

  interface get_PEST_data
    module procedure get_PEST_VMEC_data
    module procedure get_PEST_radial_data
    module procedure get_PEST_field_line_data
    module procedure get_PEST_surface_data
    module procedure get_PEST_volume_data
  end interface

contains

  type(PEST_Obj) function create_from_VMEC(VMEC_id,surfaces,n_surf,n_field_lines,n_parallel_pts) result(pest)
    character(len=:), intent(in), pointer :: VMEC_id
    real(dp), dimension(:), intent(in) :: surfaces
    integer, intent(in) :: n_surf, n_field_lines, n_parallel_pts
    logical :: verbose

    type(VMEC_Obj) :: vmec
    vmec = create_VMEC_Obj(VMEC_id)

!    j = 0
!    do while (surfaces(j+1) .gt. 1e-8)
!      j = j+1
!      print *,surfaces(j), j
!    end do
!
    !if (j .eq. 0) then
    !  print *, "VMEC2PEST error! Cannot compute a surface at s = 0"
    !  stop
    !end if 
    !nsurf = j

    pest%vmec = vmec
    pest%nx1 = n_surf
    pest%nx2 = n_field_lines
    pest%nx3 = n_parallel_pts+1

    pest%ix11 = 0
    pest%ix12 = n_surf-1
    pest%ix21 = 0
    pest%ix22 = n_field_lines-1
    pest%ix31 = -n_parallel_pts/2
    pest%ix32 = n_parallel_pts/2

    if(allocated(pest%iota)) deallocate(pest%iota)
    if(allocated(pest%shat)) deallocate(pest%shat)
    if(allocated(pest%safety_factor_q)) deallocate(pest%safety_factor_q)

    if(allocated(pest%x1)) deallocate(pest%x1)
    if(allocated(pest%x2)) deallocate(pest%x2)
    if(allocated(pest%x3)) deallocate(pest%x3)
    if(allocated(pest%bmag)) deallocate(pest%bmag)
    if(allocated(pest%jac)) deallocate(pest%jac)
    if(allocated(pest%g11)) deallocate(pest%g11)
    if(allocated(pest%g12)) deallocate(pest%g12)
    if(allocated(pest%g22)) deallocate(pest%g22)
    if(allocated(pest%g13)) deallocate(pest%g13)
    if(allocated(pest%g23)) deallocate(pest%g23)
    if(allocated(pest%g33)) deallocate(pest%g33)
    if(allocated(pest%gradB_drift_x1)) deallocate(pest%curv_drift_x1)
    if(allocated(pest%gradB_drift_x2)) deallocate(pest%curv_drift_x2)
    if(allocated(pest%curv_drift_x1)) deallocate(pest%curv_drift_x1)
    if(allocated(pest%curv_drift_x2)) deallocate(pest%curv_drift_x2)
    if(allocated(pest%d_B_d_x3)) deallocate(pest%d_B_d_x3)
    if(allocated(pest%Rsurf)) deallocate(pest%Rsurf)
    if(allocated(pest%Zsurf)) deallocate(pest%Zsurf)
    if(allocated(pest%d_Lambda_d_theta_vmec)) deallocate(pest%d_Lambda_d_theta_vmec)

    allocate(pest%iota(pest%ix11:pest%ix12))
    allocate(pest%shat(pest%ix11:pest%ix12))
    allocate(pest%safety_factor_q(pest%ix11:pest%ix12))
    allocate(pest%x1(pest%ix11:pest%ix12))
    allocate(pest%x2(pest%ix21:pest%ix22))
    allocate(pest%x3(pest%ix31:pest%ix32,pest%ix11:pest%ix12))
    allocate(pest%bmag(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%jac(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%g11(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%g12(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%g22(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%g13(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%g23(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%g33(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%curv_drift_x1(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%curv_drift_x2(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%gradB_drift_x1(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%gradB_drift_x2(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%d_B_d_x3(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%Rsurf(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%Zsurf(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))
    allocate(pest%d_Lambda_d_theta_vmec(pest%ix21:pest%ix22, pest%ix31:pest%ix32, pest%ix11:pest%ix12))

    pest%x1(pest%ix11:pest%ix12) = surfaces(1:n_surf)
  end function

  type(PEST_Obj) function create_from_field_line_file(fl_file) result(pest)
    character(len=2000), intent(in) :: fl_file

  end function

  subroutine set_PEST_reference_values(pest,norm_type)
    type(PEST_Obj), intent(inout) :: pest
    character(len=7), intent(in) :: norm_type
    real(dp) :: edge_toroidal_flux_over_2pi
    logical :: verbose
    verbose = .false.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Set reference values
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    edge_toroidal_flux_over_2pi = pest%vmec%phi(pest%vmec%ns) / (2*pi) * pest%vmec%isigng ! isigns is called signgs in the wout*.nc file. Why is this signgs here?

    ! Set reference length and magnetic field for GS2's normalization, 
    ! using the choices made by Pavlos Xanthopoulos in GIST:
    pest%major_R = pest%vmec%Rmajor
    pest%minor_r = pest%vmec%Aminor ! Note that 'Aminor' in read_wout_mod is called 'Aminor_p' in the wout*.nc file.

    if (norm_type == 'minor_r') then
      pest%L_ref = pest%minor_r
    else
      pest%L_ref = pest%major_R
    end if

    pest%B_ref = 2 * abs(edge_toroidal_flux_over_2pi) / (pest%L_ref * pest%L_ref)
    if (verbose) then
       print *,"  Reference length for normalization:",pest%L_ref," meters."
       print *,"  Reference magnetic field strength normalization:",pest%B_ref," Tesla."
    end if
  end subroutine

  subroutine destroy_PEST_Obj(pest)
    type(PEST_Obj), intent(inout) :: pest
    
    call destroy_VMEC_Obj(pest%vmec)
    deallocate(pest%iota)
    deallocate(pest%shat)
    deallocate(pest%safety_factor_q)
    deallocate(pest%x1)
    deallocate(pest%x2)
    deallocate(pest%x3)
    deallocate(pest%bmag)
    deallocate(pest%jac)
    deallocate(pest%g11)
    deallocate(pest%g12)
    deallocate(pest%g22)
    deallocate(pest%g13)
    deallocate(pest%g23)
    deallocate(pest%g33)
    deallocate(pest%gradB_drift_x1)
    deallocate(pest%gradB_drift_x2)
    deallocate(pest%curv_drift_x1)
    deallocate(pest%curv_drift_x2)
    deallocate(pest%d_B_d_x3)
    deallocate(pest%Rsurf)
    deallocate(pest%Zsurf)
    deallocate(pest%d_Lambda_d_theta_vmec)
  end subroutine

  subroutine get_PEST_surface_data(pest,idx1,data_name,surf_data)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    character(*), intent(in) :: data_name
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32), intent(out) :: surf_data
    integer :: str_len

    select case(trim(data_name))
      case('g11')
        surf_data = pest%g11(:,:,idx1)
      case('g12')
        surf_data = pest%g12(:,:,idx1)
      case('g22')
        surf_data = pest%g22(:,:,idx1)
      case('g13')
        surf_data = pest%g13(:,:,idx1)
      case('g23')
        surf_data = pest%g23(:,:,idx1)
      case('g33')
        surf_data = pest%g33(:,:,idx1)
      case('bmag')
        surf_data = pest%bmag(:,:,idx1)
      case('jac')
        surf_data = pest%jac(:,:,idx1)
      case('curv_drift_x1')
        surf_data = pest%curv_drift_x1(:,:,idx1)
      case('curv_drift_x2')
        surf_data = pest%curv_drift_x2(:,:,idx1)
      case('d_B_d_x3')
        surf_data = pest%d_B_d_x3(:,:,idx1)
      case('R')
        surf_data = pest%Rsurf(:,:,idx1)
      case('Z')
        surf_data = pest%Zsurf(:,:,idx1)
      case default
        write(6,"(A)") "VMEC2PEST error! Data field with name ",trim(data_name)," does not exist!"
        stop
    end select
  end subroutine 

  subroutine get_PEST_volume_data(pest,data_name,vol_data)
    type(PEST_Obj), intent(in) :: pest
    character(*), intent(in) :: data_name
    real(dp), dimension(pest%ix21:pest%ix22,pest%ix31:pest%ix32,pest%ix11:pest%ix12), intent(out) :: vol_data 

    select case(trim(data_name))
      case('g11')
        vol_data = pest%g11
      case('g12')
        vol_data = pest%g12
      case('g22')
        vol_data = pest%g22
      case('g13')
        vol_data = pest%g13
      case('g23')
        vol_data = pest%g23
      case('g33')
        vol_data = pest%g33
      case('bmag')
        vol_data = pest%bmag
      case('jac')
        vol_data = pest%jac
      case('curv_drift_x1')
        vol_data = pest%curv_drift_x1
      case('curv_drift_x2')
        vol_data = pest%curv_drift_x2
      case('d_B_d_x3')
        vol_data = pest%d_B_d_x3
      case default
        write(6,"(A)") "VMEC2PEST error! Data field with name ",trim(data_name)," does not exist!"
        stop
    end select
  end subroutine

  subroutine get_PEST_field_line_data(pest,idx1,idx2,data_name,line_data)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1, idx2
    character(*), intent(in) :: data_name
    real(dp), dimension(pest%ix31:pest%ix32), intent(out) :: line_data

    select case(trim(data_name))
      case('g11')
        line_data = pest%g11(idx2,:,idx1)
      case('g12')
        line_data = pest%g12(idx2,:,idx1)
      case('g22')
        line_data = pest%g22(idx2,:,idx1)
      case('g13')
        line_data = pest%g13(idx2,:,idx1)
      case('g23')
        line_data = pest%g23(idx2,:,idx1)
      case('g33')
        line_data = pest%g33(idx2,:,idx1)
      case('bmag')
        line_data = pest%bmag(idx2,:,idx1)
      case('jac')
        line_data = pest%jac(idx2,:,idx1)
      case('curv_drift_x1')
        line_data = pest%curv_drift_x1(idx2,:,idx1)
      case('curv_drift_x2')
        line_data = pest%curv_drift_x2(idx2,:,idx1)
      case('d_B_d_x3')
        line_data = pest%d_B_d_x3(idx2,:,idx1)
      case('x3')
        line_data = pest%x3(:,idx1)
      case default
        write(6,"(A)") "VMEC2PEST error! Data field with name ",trim(data_name)," does not exist!"
        stop
    end select
  end subroutine

  subroutine get_PEST_radial_data(pest,idx1,data_name,rad_data)
    type(PEST_Obj), intent(in) :: pest
    integer, intent(in) :: idx1
    character(*), intent(in) :: data_name
    real(dp), intent(out) :: rad_data

    select case(trim(data_name))
      case('iota')
        rad_data = pest%iota(idx1)
      case('safety_factor_q')
        rad_data = pest%safety_factor_q(idx1)
      case('shat')
        rad_data = pest%shat(idx1)
      case('x1')
        rad_data = pest%x1(idx1)
      case default
        write(6,"(A)") "VMEC2PEST error! Data field with name ",trim(data_name)," does not exist!"
        stop
    end select
  end subroutine

  subroutine get_PEST_VMEC_data(pest,data_name,vmec_data)
    type(PEST_Obj), intent(in) :: pest
    character(*), intent(in) :: data_name
    real(dp), intent(out) :: vmec_data

    select case(trim(data_name))
      case('nfp')
        vmec_data = pest%vmec%nfp
      case default
        write(6,"(A)") "VMEC2PEST error! Data field with name ",trim(data_name)," does not exist!"
        stop
    end select
  end subroutine
end module 
