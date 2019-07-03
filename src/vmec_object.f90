!*******************************************************************************
! This module defines the VMEC object used in the vmec2pest calculation
! Author: B.J. Faber - University of Wisconsin-Madison (bfaber@wisc.edu)
! Creation date: June 2019
!******************************************************************************

module vmec_object
  use types, only: dp, pi
  implicit none
  
  public VMEC_Obj, create_VMEC_Obj, destroy_VMEC_Obj

  type :: VMEC_Obj
    integer :: ns, nfp, isigng, mpol, ntor, mnmax, mnmax_nyq

    real(dp) :: Aminor, Rmajor

    logical :: lasym

    real(dp), dimension(:), allocatable :: phi ! Toroidal flux 
    real(dp), dimension(:), allocatable :: phip
    real(dp), dimension(:), allocatable :: xm
    real(dp), dimension(:), allocatable :: xn
    real(dp), dimension(:), allocatable :: xm_nyq
    real(dp), dimension(:), allocatable :: xn_nyq
    real(dp), dimension(:), allocatable :: iotas
    real(dp), dimension(:), allocatable :: iotaf
    real(dp), dimension(:), allocatable :: presf

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

  type(VMEC_Obj) function create_VMEC_Obj(VMEC_id) result(vmec)
    use read_wout_mod
    character(len=2000), intent(in) :: VMEC_id
    integer :: ierr, iopen, j
    real(dp) :: dphi
    logical :: verbose

    verbose = .false.

    ! If we are reading a file, then we are running outside the scope of STELLOPT,
    ! which supplies the read_wout_mod module
    if (trim(VMEC_id) .ne. "") then
      if (verbose) print *,"  About to read VMEC wout file ",trim(VMEC_id)
      call read_wout_file(trim(VMEC_id), ierr, iopen)
      if (iopen .ne. 0) stop 'error opening wout file'
      if (ierr .ne. 0) stop 'error reading wout file'
      if (verbose) print *,"  Successfully read VMEC data from ",trim(VMEC_id)
    end if

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
    vmec%mpol = mpol
    vmec%ntor = ntor
        

    if(allocated(vmec%phi)) deallocate(vmec%phi)
    if(allocated(vmec%phip)) deallocate(vmec%phip)
    if(allocated(vmec%xm)) deallocate(vmec%xm)
    if(allocated(vmec%xn)) deallocate(vmec%xn)
    if(allocated(vmec%xm_nyq)) deallocate(vmec%xm_nyq)
    if(allocated(vmec%xn_nyq)) deallocate(vmec%xn_nyq)
    if(allocated(vmec%iotas)) deallocate(vmec%iotas)
    if(allocated(vmec%iotaf)) deallocate(vmec%iotaf)
    if(allocated(vmec%presf)) deallocate(vmec%presf)

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

    !*********************************************************************
    ! Read in everything from the vmec wout file using libstell.
    !*********************************************************************
 
    allocate(vmec%phi(ns))
    vmec%phi = phi

    allocate(vmec%phip(ns))
    vmec%phip = phip

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

    allocate(vmec%presf(ns))
    vmec%presf = presf

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

   
    ! There is a bug in libstell read_wout_file for ASCII-format wout files, in which the xm_nyq and xn_nyq arrays are sometimes
    ! not populated. The next few lines here provide a workaround:
    if (maxval(abs(vmec%xm_nyq)) < 1 .and. maxval(abs(vmec%xn_nyq)) < 1) then
       if (vmec%mnmax_nyq == vmec%mnmax) then
          if (verbose) print *,"xm_nyq and xn_nyq arrays are not populated in the wout file. Using xm and xn instead."
          vmec%xm_nyq = vmec%xm
          vmec%xn_nyq = vmec%xn
       else
          print *,"Error! xm_nyq and xn_nyq arrays are not populated in the wout file, and mnmax_nyq != mnmax."
          stop
       end if
    end if

    ! --------------------------------------------------------------------------------
    ! Do some sanity checking to ensure the VMEC arrays have some expected properties.
    ! --------------------------------------------------------------------------------

    ! 'phi' is vmec's array of the toroidal flux (not divided by 2pi!) on vmec's radial grid.
    if (abs(vmec%phi(1)) > 1d-14) then
       print *,"Error! VMEC phi array does not begin with 0."
       print *,"phi:",vmec%phi
       stop
    end if

    dphi = vmec%phi(2) - vmec%phi(1)
    do j=3,vmec%ns
       if (abs(vmec%phi(j)-vmec%phi(j-1)-dphi) > 1d-11) then
          print *,"Error! VMEC phi array is not uniformly spaced."
          print *,"phi:",vmec%phi
          stop
       end if
    end do

    ! The variable called 'phips' in the wout file is called just 'phip' in read_wout_mod.F.
    ! phips is on the half-mesh, so skip first point.
    do j=2,vmec%ns
       if (abs(vmec%phip(j)+vmec%phi(ns)/(2*pi)) > 1d-11) then
          print *,"Error! VMEC phips array is not constant and equal to -phi(ns)/(2*pi)."
          print *,"phip(s):",vmec%phip
          stop
       end if
    end do

    ! The first mode in the m and n arrays should be m=n=0:
    if (vmec%xm(1) .ne. 0) stop "First element of xm in the wout file should be 0."
    if (vmec%xn(1) .ne. 0) stop "First element of xn in the wout file should be 0."
    if (vmec%xm_nyq(1) .ne. 0) stop "First element of xm_nyq in the wout file should be 0."
    if (vmec%xn_nyq(1) .ne. 0) stop "First element of xn_nyq in the wout file should be 0."

    ! Lambda should be on the half mesh, so its value at radial index 1 should be 0 for all (m,n)
    if (maxval(abs(vmec%lmns(:,1))) > 0) then
       print *,"Error! Expected lmns to be on the half mesh, but its value at radial index 1 is nonzero."
       print *,"Here comes lmns(:,1):", vmec%lmns(:,1)
       stop
    end if
    if (vmec%lasym) then
       if (maxval(abs(vmec%lmnc(:,1))) > 0) then
          print *,"Error! Expected lmnc to be on the half mesh, but its value at radial index 1 is nonzero."
          print *,"Here comes lmnc(:,1):", vmec%lmnc(:,1)
          stop
       end if
    end if

    ! --------------------------------------------------------------------------------
    ! End of sanity checks.
    ! --------------------------------------------------------------------------------

  end function

!  subroutine read_VMEC_file(VMEC_file,vmec)
!    use read_wout_mod
!  end subroutine

  subroutine destroy_VMEC_Obj(vmec)
    type(VMEC_Obj), intent(inout) :: vmec
    deallocate(vmec%phi)
    deallocate(vmec%phip)
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
