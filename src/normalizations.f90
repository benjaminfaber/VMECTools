module normalizations
  use types, only: dp, pi
  use pest_object, only: PEST_Obj
  implicit none

  public gene_normalizations, gs2_normalizations

contains
  
  subroutine gene_normalizations(pest)
    type(PEST_Obj), intent(inout) :: pest
    integer :: j, k, s
    
    ! Adjust the jacobian
    do s = pest%ix11,pest%ix12
      do j = pest%ix31,pest%ix32
        do k = pest%ix21, pest%ix22
          pest%jac(k,j,s) = abs(pest%jac(k,j,s)/(1.0+pest%d_Lambda_d_theta_vmec(k,j,s)) )
        end do
      end do
    end do
    ! GENE uses theta as the parallel coordinate 
    pest%x3 = pest%iota*pest%x3
    pest%g11 = pest%L_ref * pest%L_ref/(4.0*pest%x1(pest%ix11))*pest%g11 
    pest%g12 = 0.5*pest%L_ref*pest%L_ref*pest%g12
    pest%g22 = pest%x1(pest%ix11) * pest%L_ref * pest%L_ref * pest%g22
    pest%g13 = 0.5*pest%L_ref*pest%L_ref/sqrt(pest%x1(pest%ix11))*pest%g13
    pest%g23 = pest%L_ref*pest%L_ref*sqrt(pest%x1(pest%ix11))*pest%g23
    pest%g33 = pest%L_ref*pest%L_ref*pest%g33
    pest%curv_drift_x1 = pest%L_ref*pest%L_ref*sqrt(pest%x1(pest%ix11))*pest%curv_drift_x1/pest%bmag 
    pest%curv_drift_x2 = pest%L_ref*pest%L_ref/(2.0*sqrt(pest%x1(pest%ix11)))*pest%curv_drift_x2/pest%bmag
    pest%bmag = 1.0/pest%B_ref*pest%bmag 
    pest%jac = 2.0*pest%safety_factor_q/(pest%L_ref*pest%L_ref*pest%L_ref)*pest%jac

  end subroutine

  subroutine gs2_normalizations(pest)
    type(PEST_Obj), intent(inout) :: pest

    pest%bmag = 1.0/pest%B_ref*pest%bmag

  end subroutine

end module
