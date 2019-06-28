!******************************************************************************
! This module defines normalizations for different simuation grid types
! that are computed through the  vmec2pest calculation
! Author: B.J. Faber (bfaber@wisc.edu)
! Creation date: June 2019
!******************************************************************************
module normalizations
  use types, only: dp, pi
  use pest_object, only: PEST_Obj
  implicit none

  public gene_normalizations, gs2_normalizations

contains
  
  subroutine gene_normalizations(pest)
    type(PEST_Obj), intent(inout) :: pest
    integer :: idx1, idx2, idx3
    
    ! Adjust the jacobian
    do idx1 = pest%ix11,pest%ix12
      do idx2 = pest%ix31,pest%ix32
        do idx3 = pest%ix21, pest%ix22
          pest%jac(idx2,idx3,idx1) = abs(pest%jac(idx2,idx3,idx1)/(1.0+pest%d_Lambda_d_theta_vmec(idx2,idx3,idx1)) )
        end do
      end do
      ! GENE uses theta as the parallel coordinate 
      pest%g11(:,:,idx1) = pest%L_ref * pest%L_ref/(4.0*pest%x1(idx1))*pest%g11(:,:,idx1) 
      pest%g12(:,:,idx1) = 0.5*pest%L_ref*pest%L_ref*pest%g12(:,:,idx1)
      pest%g22(:,:,idx1) = pest%x1(pest%ix11) * pest%L_ref * pest%L_ref * pest%g22(:,:,idx1)
      pest%g13(:,:,idx1) = 0.5*pest%L_ref*pest%L_ref/sqrt(pest%x1(idx1))*pest%g13(:,:,idx1)
      pest%g23(:,:,idx1) = pest%L_ref*pest%L_ref*sqrt(pest%x1(idx1))*pest%g23(:,:,idx1)
      pest%g33(:,:,idx1) = pest%L_ref*pest%L_ref*pest%g33(:,:,idx1)
      pest%curv_drift_x1(:,:,idx1) = pest%L_ref*pest%L_ref*sqrt(pest%x1(idx1))*pest%curv_drift_x1(:,:,idx1)/pest%bmag(:,:,idx1)
      pest%curv_drift_x2(:,:,idx1) = pest%L_ref*pest%L_ref/(2.0*sqrt(pest%x1(idx1)))*pest%curv_drift_x2(:,:,idx1)/pest%bmag(:,:,idx1)
      pest%d_B_d_x3(:,:,idx1) = pest%L_ref*pest%d_B_d_x3(:,:,idx1)/(pest%vmec%nfp*pest%jac(:,:,idx1)*pest%bmag(:,:,idx1))
      pest%bmag(:,:,idx1) = 1.0/pest%B_ref*pest%bmag(:,:,idx1)
      pest%jac(:,:,idx1) = 2.0*pest%safety_factor_q(idx1)/(pest%L_ref*pest%L_ref*pest%L_ref)*pest%jac(:,:,idx1)
    end do
  end subroutine

  subroutine gs2_normalizations(pest)
    type(PEST_Obj), intent(inout) :: pest

    pest%bmag = 1.0/pest%B_ref*pest%bmag

  end subroutine

end module
