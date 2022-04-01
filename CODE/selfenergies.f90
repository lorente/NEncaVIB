module selfenergies
Use declarations
Use convolution
Contains

! Pi retarded
      subroutine Self_Pi (N_omega, N_spin, Gr, Gamma_t, Gamma_s, Fermi_t, Fermi_s, omega, SelfD)
      implicit none
      integer :: N_omega, i_spin, N_spin
      real (q) :: step_omega
      real (q), intent (in) :: omega (:)
      real (q), dimension (:), intent (in) :: Gamma_t, Gamma_S, Fermi_t, Fermi_s
      complex (qc), intent (in) :: Gr (:, :)
      complex (qc), intent (out) :: SelfD (:)
      complex (qc) :: convol (N_omega) 
      complex (qc) :: h (N_omega) 
      complex (qc) :: g (N_omega) 

         SelfD = zero
         step_omega = omega (2) - omega (1)

      do i_spin = 1, N_spin 

                g (:) = Gr (i_spin, :)

                h = Fermi_t * Gamma_t + Fermi_s * Gamma_s

                call correlation_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

         SelfD = SelfD + convol

      enddo

      return
      end subroutine Self_Pi
! Pi lesser than
      subroutine Pi_lesser (N_omega, N_spin, omega, Gamma_t, Gamma_s, &
      &      fermi_t, fermi_s, Gl, SelfDl)
      implicit none
      integer :: N_omega, N_spin, i_spin
      real (q) :: step_omega
      real (q), dimension (:), intent (in) :: Gamma_t, Gamma_S, fermi_t, fermi_s
      real (q), intent (in) :: omega (:)
      complex (qc), intent (in) :: Gl (:, :)
      complex (qc), intent (out) :: SelfDl (:)
      complex (qc) :: convol (N_omega)
      complex (qc) :: h (N_omega)
      complex (qc) :: g (N_omega)

         SelfDl = zero
         step_omega = omega (2) - omega (1)

      do i_spin = 1, N_spin
                g (:) = Gl (i_spin, :)

                h = (Fermi_t-1) * Gamma_t + (Fermi_s-1) * Gamma_s

                call correlation_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

         SelfDl = SelfDl + convol

      enddo

      return
      end subroutine Pi_lesser
! Pi bigger than
      subroutine Pi_bigger (N_omega, N_spin, omega, Gamma_t, Gamma_s, &
      &      fermi_t, fermi_s, Gb, SelfDb)
      implicit none
      integer :: N_omega, N_spin, i_spin
      real (q) :: step_omega
      real (q), dimension (:), intent (in) :: Gamma_t, Gamma_S, fermi_t, fermi_s
      real (q), intent (in) :: omega (:)
      complex (qc), intent (in) :: Gb (:, :)
      complex (qc), intent (out) :: SelfDb (:)
      complex (qc) :: convol (N_omega)
      complex (qc) :: h (N_omega)
      complex (qc) :: g (N_omega)

         SelfDb = zero
         step_omega = omega (2) - omega (1)

      do i_spin = 1, N_spin
                g (:) = Gb (i_spin, :)

                h = Fermi_t * Gamma_t + Fermi_s * Gamma_s

                call correlation_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

         SelfDb = SelfDb + convol
      
      enddo

      return
      end subroutine Pi_bigger
! Sigma retarded
      subroutine Self_Sigma (N_omega, N_spin, Dr, Gamma_t, Gamma_s, fermi_t, fermi_s, omega, SelfG)
      implicit none
      integer :: N_omega, i_spin, N_spin
      real (q) :: step_omega
      real (q), intent (in) :: omega (:)
      real (q), dimension (:), intent (in) :: Gamma_t, Gamma_S, fermi_t, fermi_s
      complex (qc), intent (in) :: Dr ( :)
      complex (qc), intent (out) :: SelfG (:)
      complex (qc) :: convol (N_omega)
      complex (qc) :: h (N_omega)     
      complex (qc) :: g (N_omega)

         step_omega = omega (2) - omega (1)

                g (:) = Dr (:)

                h = (1-Fermi_t) * Gamma_t + (1-Fermi_s) * Gamma_s

                call convolution_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

         SelfG = convol


      return
      end subroutine Self_Sigma
! Sigma lesser
      subroutine Sigma_lesser (N_omega, N_spin, Dl, Gamma_t, Gamma_s, &
      &     fermi_t, fermi_s, omega, N_vib, indexvib, W, Bose, Gl,  SelfGl)
      implicit none
      integer :: i_omega, N_omega, i_spin, N_spin
      integer :: i_vib, N_vib
      integer, dimension (:):: indexvib
      real (q), dimension (:) :: W, Bose
      real (q) :: step_omega
      real (q), intent (in) :: omega (:)
      real (q), dimension (:), intent (in) :: Gamma_t, Gamma_S, fermi_t, fermi_s
      complex (qc), intent (in) :: Dl ( :)
      complex (qc), intent (in) :: Gl (:,:)
      complex (qc), intent (out) :: SelfGl (:,:)
      complex (qc) :: convol (N_omega)
      complex (qc) :: h (N_omega)
      complex (qc) :: g (N_omega)

          step_omega = omega (2) - omega (1)
          
 
                 g (:) = Dl (:)
 
                 h = (-Fermi_t) * Gamma_t + (-Fermi_s) * Gamma_s
 
                 call convolution_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

        do i_spin = 1, N_spin

         SelfGl (i_spin, :) = convol (:)

         do i_vib = 1, N_vib

           do i_omega = indexvib(i_vib)+1, N_omega-indexvib(i_vib)

         SelfGl (i_spin, i_omega) = SelfGl (i_spin, i_omega) + W(i_vib)**2*( Bose(i_vib) *  &
      &    Gl (i_spin, i_omega-indexvib(i_vib)) + (Bose(i_vib) + 1) *  &
      &    Gl (i_spin, i_omega+indexvib(i_vib))) 

           enddo

         enddo

        enddo


      return
      end subroutine Sigma_lesser
! Sigma lesser
      subroutine Sigma_lesser0 (N_omega, N_spin, Dl, Gamma_t, Gamma_s, &
      &     fermi_t, fermi_s, omega, SelfGl)
      implicit none
      integer :: N_omega, i_spin, N_spin
      real (q) :: step_omega
      real (q), intent (in) :: omega (:)
      real (q), dimension (:), intent (in) :: Gamma_t, Gamma_S, fermi_t, fermi_s
      complex (qc), intent (in) :: Dl ( :)
      complex (qc), intent (out) :: SelfGl (:,:)
      complex (qc) :: convol (N_omega)
      complex (qc) :: h (N_omega)
      complex (qc) :: g (N_omega)

         step_omega = omega (2) - omega (1)

                g (:) = Dl (:)

                h = (-Fermi_t) * Gamma_t + (-Fermi_s) * Gamma_s

                call convolution_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

         do i_spin = 1, N_spin
         SelfGl (i_spin,:) = convol (:)
         enddo


      return
      end subroutine Sigma_lesser0

! Sigma BIGGER
      subroutine Sigma_bigger (N_omega, N_spin, Db, Gamma_t, Gamma_s, &
      &     fermi_t, fermi_s, omega, N_vib, indexvib, W, Bose, Gb,  SelfGb)
      implicit none
      integer :: i_omega, N_omega, i_spin, N_spin
      integer :: i_vib, N_vib
      integer, dimension (:) :: indexvib
      real (q), dimension (:) ::  W, Bose
      real (q) :: step_omega
      real (q), intent (in) :: omega (:)
      real (q), dimension (:), intent (in) :: Gamma_t, Gamma_S, fermi_t, fermi_s
      complex (qc), intent (in) :: Db ( :), Gb (:,:)
      complex (qc), intent (out) :: SelfGb (:,:)
      complex (qc) :: convol (N_omega)
      complex (qc) :: h (N_omega)
      complex (qc) :: g (N_omega)

          step_omega = omega (2) - omega (1)
 
                 g (:) = Db (:)
 
                 h = (1-Fermi_t) * Gamma_t + (1-Fermi_s) * Gamma_s
 
                 call convolution_fft (N_omega, 2*N_omega, step_omega, g, h, convol)
 
         do i_spin = 1, N_spin
 
         SelfGb (i_spin, :) = convol (:)

         do i_vib = 1, N_vib

           do i_omega = indexvib(i_vib)+1, N_omega-indexvib(i_vib)

         SelfGb (i_spin, i_omega) = SelfGb (i_spin, i_omega) + W(i_vib)**2*(Bose(i_vib) *  &
      &    Gb (i_spin, i_omega+indexvib(i_vib)) + (Bose(i_vib) + 1) *  &
      &    Gb (i_spin, i_omega-indexvib(i_vib))) 

           enddo


         enddo

        enddo

      return
      end subroutine Sigma_bigger

      subroutine Sigma_bigger0 (N_omega, N_spin, Db, Gamma_t, Gamma_s, &
      &    fermi_t, fermi_s, omega, SelfGb)
      implicit none
      integer :: N_omega, i_spin, N_spin
      real (q) :: step_omega
      real (q), intent (in) :: omega (:)
      real (q), dimension (:), intent (in) :: Gamma_t, Gamma_S, fermi_t, fermi_s
      complex (qc), intent (in) :: Db ( :)
      complex (qc), intent (out) :: SelfGb (:)
      complex (qc) :: convol (N_omega)
      complex (qc) :: h (N_omega)
      complex (qc) :: g (N_omega)

         step_omega = omega (2) - omega (1)

                g (:) = Db (:)

                h = (1-Fermi_t) * Gamma_t + (1-Fermi_s) * Gamma_s

                call convolution_fft (N_omega, 2*N_omega, step_omega, g, h, convol)

         SelfGb = convol

      return
      end subroutine Sigma_bigger0


end module selfenergies 
