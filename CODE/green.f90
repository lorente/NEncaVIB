module green
Use declarations
Use convolution
CONTAINS

!
!D
!
   subroutine D_green_R (omega, SelfD, Dr)
   implicit none
   real (q) :: omega (:)
   complex (qc) :: SelfD(:), Dr (:)

       Dr (:) = 1._q/(omega(:)-SelfD(:))
 
    return
   end subroutine D_green_R
!
   subroutine D_lesser (Dr, SelfDl, Dl)
   implicit none
   complex (qc) :: SelfDl(:), Dl (:), Dr (:)

      Dl = Dr * SelfDl * conjg(Dr)
 
    return
   end subroutine D_lesser
!
   subroutine D_bigger (Dr, SelfDb, Db)
   implicit none
   complex (qc) :: SelfDb(:), Db (:), Dr (:)

      Db = Dr * SelfDb * conjg(Dr)
 
    return
   end subroutine D_bigger
!
!G
!
   subroutine G_green_R (N_spin, omega, Hamiltonian, SelfG, Gr)
   implicit none
   integer :: i_spin, N_spin
   real (q) :: omega (:)
   real (q) :: Hamiltonian (:)
   complex (qc) :: SelfG(:), Gr (:,:)

      do i_spin = 1, N_spin
        Gr (i_spin, :) = 1._q / (omega (:) - Hamiltonian (i_spin) - SelfG (:))
      enddo

    return
   end subroutine G_green_R
!
   subroutine G_bigger (N_spin, Gr, SelfGb, Gb)
   implicit none
   integer :: i_spin, N_spin
   complex (qc) :: SelfGb(:,:), Gr (:,:), Gb (:,:)

      do i_spin = 1, N_spin
        Gb (i_spin, :) = Gr (i_spin, :)*SelfGb (i_spin, :)*conjg(Gr (i_spin, :))
      enddo

    return
   end subroutine G_bigger
!
   subroutine G_lesser (N_spin, Gr, SelfGl, Gl)
   implicit none
   integer :: i_spin, N_spin
   complex (qc) :: SelfGl(:,:), Gl (:,:), Gr (:,:)

      do i_spin = 1, N_spin
        Gl (i_spin, :) = Gr (i_spin, :)*SelfGl (i_spin, :)*conjg(Gr (i_spin, :))
      enddo

    return
   end subroutine G_lesser
!
! normalization
!
    subroutine normalization (N_omega, Dl, Gl, constant, omega)
    implicit none
    integer :: i, N_omega
    real (q) :: step_omega
    real (q) :: omega (:)
    complex (qc) :: constant
    complex (qc) :: Dl (:), Gl (:,:)

      step_omega = omega (2)-omega (1)
      constant = zero
       do i = 1, N_omega
         constant = constant + Dl (i)*step_omega &
      &           - sum (Gl (:,i))*step_omega
       enddo
      constant = 1._q/constant
    return
    end subroutine normalization

!
! Physical Greens function
!
    subroutine GreenR (N_omega, N_spin, constant, Dl, Db, Gl, Gb, GGR, omega)
    implicit none
    integer :: i_spin, N_spin, N_omega
    real (q) :: step_omega
    real (q) :: omega (:)
    complex (qc) :: constant
    complex (qc) :: Dl (:), Db (:), Gl (:,:), Gb (:,:), GGR (:,:)
    complex (qc) :: convol (N_omega)
    complex (qc) :: h (N_omega)
    complex (qc) :: g (N_omega)


         GGR = zero
         step_omega = omega (2) - omega (1)

      do i_spin = 1, N_spin
                g (:) = Gl (i_spin, :)

                h = Db

                call correlation_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

         GGR (i_spin,:) = GGR (i_spin,:) + convol (:)

                g (:) = Gb (i_spin, :)

                h = Dl

                call correlation_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

         GGR (i_spin,:) = GGR (i_spin,:) - convol (:)


      enddo

         GGR = - constant * GGR / 2._q

      return
    

    end subroutine GreenR

end module green
