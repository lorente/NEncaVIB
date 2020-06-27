module electronic
Use declarations
Use used_files
CONTAINS
!
! check for accuracy
!
      subroutine checkSumRules (omega, N_spin, D, G)
      implicit none
      integer :: i, N_spin
      real (q) :: omega (:)
      real (q) :: somme, step_omega
      complex (qc), intent (in):: D (:), G (:,:)

      step_omega=omega (2)-omega (1)
      somme = -sum (aimag(D))/pi_d

      print *, " "
      print *, "____________________________________________"
      print *, "Règles de somme pour les retardées : "
      print *, " (doit être 1) "
      print *, "              D :", somme*step_omega
      print *, " "

      do i= 1, N_spin
      somme = -sum (aimag(G (i,:)))/pi_d
      print *, i, "      G :", somme*step_omega
      print *, "____________________________________________"
      enddo

      return
      end subroutine checkSumRules
!
! check for convergence
!
      subroutine convergence (N_omega, N_spin, D, D_store,  &
                               G, G_store, tol, converged)

      implicit none

      real (q),intent (in):: tol
      real (q):: somme
      integer :: N_omega, i, j, N_spin, counter
      logical, intent(out) :: converged
      complex (qc), intent (in) :: D (:), D_store (:), G (:,:), G_store (:,:)

        somme = 0._q
        counter = 0
        converged = .false.

      do i= 1, N_omega
        somme = somme + abs (D(i)-D_store(i))
        counter = counter + 1
      enddo

      do i = 1, N_spin
      do j = 1,N_omega
        somme = somme + abs (G(i,j)-G_store(i,j))
        counter = counter + 1
      enddo
      enddo

      somme = somme / counter
   
      print *,"Deviation in the sum of all elements:",somme

          if (somme < tol ) then
              converged = .true.
          endif

      return
      end subroutine convergence
!
! Electronic current
!
      subroutine current (Gamma_t, Gamma_s, GGR, fermi_t, fermi_s, &
     &       curr, omega, N_omega, N_spin)
      implicit none
      integer :: N_omega, i, j, N_spin
      real (q) :: curr, step_omega
      real (q) :: omega (:)
      real (q) :: Gamma_t (:), Gamma_s (:), fermi_t (:), fermi_s (:)
      complex (qc) :: GGR (:,:)

       step_omega = omega (2) - omega (1)
       curr = 0._q
        
      do j=1, N_spin
      do i = 1, N_omega
        curr = curr + (fermi_t (i)-fermi_s(i))*Gamma_t(i)*Gamma_s(i)  &
      &   *aimag (GGR(j,i))/ (Gamma_t(i)+Gamma_s(i))
      enddo 
      enddo 

      if (N_spin == 1) then
        curr = -2*curr*step_omega/pi_d
      else
        curr = -curr*step_omega/pi_d
      endif

      return
      end subroutine current
end module electronic
