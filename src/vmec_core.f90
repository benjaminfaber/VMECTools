module vmec_core
  use types: dp
  
  public VMEC_Object, create_VMEC_Object, destroy_VMEC_Object

  type :: VMEC_Object
    integer :: ns, nfp

    real(dp), dimension(:), allocatable :: phi ! Toroidal flux 
  end type

contains

  type(VMEC_Object) function create_VMEC_Object(VMEC_file) result(VMEC_Obj)
    character(len=2000), intent(in) :: VMEC_file
    integer :: ierr, iopen

    if (verbose) print *,"  About to read VMEC wout file ",trim(vmec_filename)
    call read_wout_file(VMEC_file, ierr, iopen)
    if (iopen .ne. 0) stop 'error opening wout file'
    if (ierr .ne. 0) stop 'error reading wout file'
    if (verbose) print *,"  Successfully read VMEC data from ",trim(vmec_filename)

    if (verbose) print *,"  Number of field periods (nfp):",nfp
    if (verbose) print *,"  Stellarator-asymmetric? (lasym):",lasym
    
    ! Assign the VMEC object entries
    VMEC_Obj%nfp = nfp
    VMEC_Obj%lasym = lasym
    VMEC_Obj%ns = ns
    VMEC_Obj%isigng = isigng
    VMEC_Obj%Aminor = Aminor
    VMEC_Obj%Rmajor = Rmajor

    if(allocated(VMEC_Obj%phi)) deallocate(VMEC_Obj%phi)

    allocate(VMEC_Obj%phi(ns))
    VMEC_Obj%phi = phi
  end function

  subroutine read_VMEC_file(VMEC_file,VMEC_Obj)
    use read_wout_mod
  end subroutine

  subroutine destroy_VMEC_Object
    deallocate(phi)
  end subroutine
end module
