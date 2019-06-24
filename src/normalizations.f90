module normalizations
  use types, only: dp, pi
  use pest_object, only: PEST_Obj
  implicit none

  public gene_normalizations, gs2_normalizations

contains
  
  subroutine gene_normalizations(pest)
    type(PEST_Obj), intent(inout) :: pest
    real(dp) :: d_iota_d_s
    integer :: j, k
    
    ! GENE uses theta as the parallel coordinate 
    d_iota_d_s = - 0.5*pest%iota*pest%shat/pest%s0
    print *, d_iota_d_s, pest%L_ref
    do j=pest%iz1,pest%iz2
      do k=pest%ia1,pest%ia2
        pest%gsz(k,j,1) = pest%gsa(k,j,1) - pest%zeta(j)*d_iota_d_s*pest%gss(k,j,1) + pest%iota*pest%gsz(k,j,1)
        pest%gaz(k,j,1) = pest%gaa(k,j,1) - pest%zeta(j)*d_iota_d_s*pest%gsa(k,j,1) + pest%iota*pest%gaz(k,j,1)
        pest%gzz(k,j,1) = pest%gaa(k,j,1) + pest%zeta(j)*pest%zeta(j)*d_iota_d_s*d_iota_d_s*pest%gss(k,j,1) + &
          & pest%iota*pest%iota*pest%gzz(k,j,1) - 2.0*(pest%iota*pest%gaz(j,k,1) + &
          & pest%zeta(j)*d_iota_d_s*pest%gsa(k,j,1) + pest%zeta(j)*d_iota_d_s*pest%gsz(j,k,1))
        pest%jac(k,j,1) = abs(pest%jac(k,j,1)/(1.0+pest%d_L_d_theta_v(k,j,1)))
      end do
    end do
    pest%jac = 2.0*pest%safety_factor_q/(pest%L_ref*pest%L_ref*pest%L_ref)*pest%jac 
    pest%gss = pest%L_ref * pest%L_ref/(4.0*pest%s0)*pest%gss 
    pest%gsa = 0.5*pest%iota*pest%L_ref*pest%L_ref*pest%gsa
    pest%gaa = pest%s0 * pest%L_ref * pest%L_ref * pest%iota * pest%iota * pest%gaa
    pest%gsz = 0.5*pest%L_ref*pest%L_ref/sqrt(pest%s0)*pest%gsz
    pest%gaz = pest%L_ref*pest%L_ref*sqrt(pest%s0)*pest%iota*pest%gaz
    pest%gzz = pest%L_ref*pest%L_ref*pest%gzz
    pest%bmag = 1.0/pest%B_ref*pest%bmag
    print  *, pest%jac(1,1,1)

  end subroutine

  subroutine gs2_normalizations(pest)
    type(PEST_Obj), intent(inout) :: pest
    integer :: j, k

!    do j=pest%iz1,pest%iz2
!      do k=pest%ia1,pest%ia2
!          
!      end do
!    end do
  end subroutine

end module
