module vmec_core
  use types, only: dp
  implicit none
  
  public VMEC_Obj, create_VMEC_Obj, destroy_VMEC_Obj

  type :: VMEC_Obj
    integer :: ns, nfp, isigng, mpol, ntor, mnmax, mnmax_nyq

    real(dp) :: Aminor, Rmajor

    logical :: lasym

    real(dp), dimension(:), allocatable :: phi ! Toroidal flux 
    real(dp), dimension(:), allocatable :: xm
    real(dp), dimension(:), allocatable :: xn
    real(dp), dimension(:), allocatable :: xm_nyq
    real(dp), dimension(:), allocatable :: xn_nyq
    real(dp), dimension(:), allocatable :: iotas
    real(dp), dimension(:), allocatable :: iotaf

    real(dp), dimension(:,:), allocatable :: rmnc
    real(dp), dimension(:,:), allocatable :: zmnc
    real(dp), dimension(:,:), allocatable :: lmnc
    real(dp), dimension(:,:), allocatable :: bmnc
    real(dp), dimension(:,:), allocatable :: gmnc
    real(dp), dimension(:,:), allocatable :: rmns
    real(dp), dimension(:,:), allocatable :: zmns
    real(dp), dimension(:,:), allocatable :: lmns
    real(dp), dimension(:,:), allocatable :: bmns
    real(dp), dimension(:,:), allocatable :: gmns
    real(dp), dimension(:,:), allocatable :: bsubumnc
    real(dp), dimension(:,:), allocatable :: bsubvmnc
    real(dp), dimension(:,:), allocatable :: bsubsmnc
    real(dp), dimension(:,:), allocatable :: bsupumnc
    real(dp), dimension(:,:), allocatable :: bsupvmnc
    real(dp), dimension(:,:), allocatable :: bsubumns
    real(dp), dimension(:,:), allocatable :: bsubvmns
    real(dp), dimension(:,:), allocatable :: bsubsmns
    real(dp), dimension(:,:), allocatable :: bsupumns
    real(dp), dimension(:,:), allocatable :: bsupvmns
    real(dp), dimension(:,:), allocatable :: currvmnc
    real(dp), dimension(:,:), allocatable :: currumnc
    real(dp), dimension(:,:), allocatable :: currvmns
    real(dp), dimension(:,:), allocatable :: currumns

  end type ! end type VMEC_Obj


contains

  type(VMEC_Obj) function create_VMEC_Obj(VMEC_file) result(vmec)
    use read_wout_mod
    character(len=2000), intent(in) :: VMEC_file
    integer :: ierr, iopen
    logical :: verbose

    verbose = .true.

    if (verbose) print *,"  About to read VMEC wout file ",trim(VMEC_file)
    call read_wout_file(VMEC_file, ierr, iopen)
    if (iopen .ne. 0) stop 'error opening wout file'
    if (ierr .ne. 0) stop 'error reading wout file'
    if (verbose) print *,"  Successfully read VMEC data from ",trim(VMEC_file)

    if (verbose) print *,"  Number of field periods (nfp):",nfp
    if (verbose) print *,"  Stellarator-asymmetric? (lasym):",lasym

    
    ! Assign the VMEC object entries
    vmec%nfp = nfp
    vmec%lasym = lasym
    vmec%ns = ns
    vmec%mnmax = mnmax
    vmec%mnmax_nyq = mnmax_nyq
    vmec%isigng = isigng
    vmec%Aminor = Aminor
    vmec%Rmajor = Rmajor
        

    if(allocated(vmec%phi)) deallocate(vmec%phi)
    if(allocated(vmec%xm)) deallocate(vmec%xm)
    if(allocated(vmec%xn)) deallocate(vmec%xn)
    if(allocated(vmec%xm_nyq)) deallocate(vmec%xm_nyq)
    if(allocated(vmec%xn_nyq)) deallocate(vmec%xn_nyq)
    if(allocated(vmec%iotas)) deallocate(vmec%iotas)
    if(allocated(vmec%iotaf)) deallocate(vmec%iotaf)

    if(allocated(vmec%rmnc)) deallocate(vmec%rmnc)
    if(allocated(vmec%bmnc)) deallocate(vmec%bmnc)
    if(allocated(vmec%gmnc)) deallocate(vmec%gmnc)
    if(allocated(vmec%zmns)) deallocate(vmec%zmns)
    if(allocated(vmec%lmns)) deallocate(vmec%lmns)
    if(allocated(vmec%bsubumnc)) deallocate(vmec%bsubumnc)
    if(allocated(vmec%bsubumnc)) deallocate(vmec%bsubvmnc)
    if(allocated(vmec%bsupumnc)) deallocate(vmec%bsupumnc)
    if(allocated(vmec%bsupvmnc)) deallocate(vmec%bsupvmnc)
    if(allocated(vmec%bsubsmns)) deallocate(vmec%bsubsmns)
    if(allocated(vmec%currumnc)) deallocate(vmec%currumnc)
    if(allocated(vmec%currvmnc)) deallocate(vmec%currvmnc)

    allocate(vmec%phi(ns))
    vmec%phi = phi

    allocate(vmec%xm(mnmax))
    vmec%xm = xm

    allocate(vmec%xn(mnmax))
    vmec%xn = xn 

    allocate(vmec%xm_nyq(mnmax_nyq))
    vmec%xm_nyq = xm_nyq

    allocate(vmec%xn_nyq(mnmax_nyq))
    vmec%xn_nyq = xn_nyq

    allocate(vmec%iotas(ns))
    vmec%iotas = iotas

    allocate(vmec%iotaf(ns))
    vmec%iotaf = iotaf

    allocate(vmec%rmnc(mnmax,ns))
    vmec%rmnc = rmnc
 
    allocate(vmec%bmnc(mnmax,ns))
    vmec%bmnc = bmnc

    allocate(vmec%gmnc(mnmax,ns))
    vmec%gmnc = gmnc

    allocate(vmec%zmns(mnmax,ns))
    vmec%zmns = zmns

    allocate(vmec%lmns(mnmax,ns))
    vmec%lmns = lmns

    allocate(vmec%bsubumnc(mnmax_nyq,ns))
    vmec%bsubumnc = bsubumnc

    allocate(vmec%bsubvmnc(mnmax_nyq,ns))
    vmec%bsubvmnc = bsubvmnc

    allocate(vmec%bsupumnc(mnmax_nyq,ns))
    vmec%bsupumnc = bsupumnc

    allocate(vmec%bsupvmnc(mnmax_nyq,ns))
    vmec%bsupvmnc = bsupvmnc

    allocate(vmec%bsubsmns(mnmax_nyq,ns))
    vmec%bsubsmns = bsubsmns

    allocate(vmec%currumnc(mnmax_nyq,ns))
    vmec%currumnc = currumnc

    allocate(vmec%currvmnc(mnmax_nyq,ns))
    vmec%currvmnc = currvmnc

    if (lasym) then
      if(allocated(vmec%zmnc)) deallocate(vmec%zmnc)
      if(allocated(vmec%lmnc)) deallocate(vmec%lmnc)
      if(allocated(vmec%rmns)) deallocate(vmec%rmns)
      if(allocated(vmec%bmns)) deallocate(vmec%bmns)
      if(allocated(vmec%gmns)) deallocate(vmec%gmns)
      if(allocated(vmec%bsubsmnc)) deallocate(vmec%bsubsmnc)
      if(allocated(vmec%bsubumns)) deallocate(vmec%bsubumns)
      if(allocated(vmec%bsubvmns)) deallocate(vmec%bsubvmns)
      if(allocated(vmec%bsupumns)) deallocate(vmec%bsupumns)
      if(allocated(vmec%bsupvmns)) deallocate(vmec%bsupvmns)
      if(allocated(vmec%currumns)) deallocate(vmec%currumns)
      if(allocated(vmec%currvmns)) deallocate(vmec%currvmns)
      allocate(vmec%zmnc(mnmax,ns))
      vmec%zmnc = zmnc
    
      allocate(vmec%lmnc(mnmax,ns))
      vmec%lmnc = lmnc

      allocate(vmec%rmns(mnmax,ns))
      vmec%rmns = rmns

      allocate(vmec%bmns(mnmax,ns))
      vmec%bmns = bmns

      allocate(vmec%gmns(mnmax,ns))
      vmec%gmns = gmns

      allocate(vmec%bsubumns(mnmax_nyq,ns))
      vmec%bsubumns = bsubumns

      allocate(vmec%bsubvmns(mnmax_nyq,ns))
      vmec%bsubvmns = bsubvmns

      allocate(vmec%bsupumns(mnmax_nyq,ns))
      vmec%bsupumns = bsupumns

      allocate(vmec%bsupvmns(mnmax_nyq,ns))
      vmec%bsupvmns = bsupvmns

      allocate(vmec%bsubsmnc(mnmax_nyq,ns))
      vmec%bsubsmnc = bsubsmnc

      allocate(vmec%currumns(mnmax_nyq,ns))
      vmec%currumns = currumns

      allocate(vmec%currvmns(mnmax_nyq,ns))
      vmec%currvmns = currvmns   

    end if

  end function

!  subroutine read_VMEC_file(VMEC_file,vmec)
!    use read_wout_mod
!  end subroutine

  subroutine destroy_VMEC_Obj(vmec)
    type(VMEC_Obj), intent(inout) :: vmec
    deallocate(vmec%phi)
    deallocate(vmec%xm)
    deallocate(vmec%xn)
    deallocate(vmec%xm_nyq)
    deallocate(vmec%xn_nyq)
    deallocate(vmec%iotas)
    deallocate(vmec%iotaf)

    deallocate(vmec%rmnc)
    deallocate(vmec%bmnc)
    deallocate(vmec%gmnc)
    deallocate(vmec%zmns)
    deallocate(vmec%lmns)
    deallocate(vmec%bsubumnc)
    deallocate(vmec%bsubvmnc)
    deallocate(vmec%bsupumnc)
    deallocate(vmec%bsupvmnc)
    deallocate(vmec%bsubsmns)
    deallocate(vmec%currumnc)
    deallocate(vmec%currvmnc)

    if (vmec%lasym) then
      deallocate(vmec%zmnc)
      deallocate(vmec%lmnc)
      deallocate(vmec%rmns)
      deallocate(vmec%bmns)
      deallocate(vmec%gmns)
      deallocate(vmec%bsubsmnc)
      deallocate(vmec%bsubumns)
      deallocate(vmec%bsubvmns)
      deallocate(vmec%bsupumns)
      deallocate(vmec%bsupvmns)
      deallocate(vmec%currumns)
      deallocate(vmec%currvmns)
    end if
  end subroutine

end module
