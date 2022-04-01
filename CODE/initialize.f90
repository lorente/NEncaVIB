module initialize
Use declarations
Use used_files
CONTAINS
   subroutine initialize_quantities (N_omega, N_spin, omega, gamma_file,  &
      & Gamma_t, Gamma_s, bias_fraction, &
      & Fermi_t, Fermi_s, bias, Temperature, impurity_hamiltonian_file,  &
      & Hamiltonian, eta, Gr, Dr)  
        implicit none
    integer :: i_omega, N_omega, N_spin, i_spin
    real (q) :: epsilon_0
    real (q) :: Gamma_0, omega_0, Delta_0
    real (q), intent (out) :: Gamma_t (:), Gamma_s (:)
    real (q), intent (out) :: Fermi_t (:), Fermi_s (:)
    real (q), intent (out) :: Hamiltonian (:)
    real (q), intent (in) :: bias, Temperature, bias_fraction
    real (q), intent (in) :: omega (:)

    complex (qc), intent (in) :: eta
    complex (qc), intent (out) :: Gr (:,:), Dr (:)

    character (len=1) :: type_fit
    character ( len = 100 ) :: gamma_file
    character ( len = 100 ) :: impurity_hamiltonian_file



! Set up Gamma's 

  open (unit_gamma, file=gamma_file, status='old')
      print *, '___________________________________________'
                  print *, '  '
         print *, "REMEMBER: Gamma is defined without 2 pi!!"
      print *, '___________________________________________'
  read (unit_gamma, *) type_fit
         if (type_fit == 'L' ) then
                  print *, '  '
                  print *, ' Gamma fitted to a Lorentzian '
                  print *, ' '
      print *, 'Parameters for tip Gamma are:'
      read (unit_gamma, *) Gamma_0, omega_0, Delta_0
      print *,'Gamma_0, omega_0, Delta_0', Gamma_0, omega_0, Delta_0
      print *, '___________________________________________'

      Gamma_0 = Gamma_0 / hartree
      omega_0 = omega_0 / hartree
      Delta_0 = Delta_0 / hartree

      Gamma_t (:) = Gamma_0*Delta_0**2 &
     &                       /((omega (:)-omega_0)**2 + Delta_0**2 )



      print *, 'Parameters for substrate Gamma are:'
      read (unit_gamma, *) Gamma_0, omega_0, Delta_0
      print *,'Gamma_0, omega_0, Delta_0', Gamma_0, omega_0, Delta_0
                  print *, ' '
      Gamma_0 = Gamma_0 / hartree
      omega_0 = omega_0 / hartree
      Delta_0 = Delta_0 / hartree

      Gamma_s (:) = Gamma_0*Delta_0**2 &
     &                       /((omega (:)-omega_0)**2 + Delta_0**2 )
      else

                  print *, '  '
                  print *, ' Gamma fitted to a Gaussian '
                  print *, ' '
      print *, 'Parameters for tip Gamma are:'
      read (unit_gamma, *) Gamma_0, omega_0, Delta_0
      print *,'Gamma_0, omega_0, Delta_0', Gamma_0, omega_0, Delta_0
      print *, '___________________________________________'

      Gamma_0 = Gamma_0 / hartree
      omega_0 = omega_0 / hartree
      Delta_0 = Delta_0 / hartree

      Gamma_t (:) = Gamma_0 &
     &                   *exp(-((omega (:)-omega_0 )/Delta_0)**2)



      print *, 'Parameters for substrate Gamma are:'
      read (unit_gamma, *) Gamma_0, omega_0, Delta_0
      print *,'Gamma_0, omega_0, Delta_0', Gamma_0, omega_0, Delta_0
                  print *, ' '

      Gamma_0 = Gamma_0 / hartree
      omega_0 = omega_0 / hartree
      Delta_0 = Delta_0 / hartree

      Gamma_s (:) = Gamma_0 &
     &                   *exp(-((omega (:)-omega_0 )/Delta_0)**2)
   endif
   close (unit_gamma)

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                               !
!  Definition of array in omega containing      !
!  the Fermi occupation factor at Temperature   !
!  including biases of tip and grounded sample  !
!  although the sign criteria correspond to the !
!  experimentally grounding the tip             !
!  (positive bias-> empty sample states)        !
!                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i_omega = 1, N_omega
          Fermi_t (i_omega) = Fermi (omega (i_omega)-bias*bias_fraction, Temperature)
          Fermi_s (i_omega) = Fermi (omega (i_omega)+bias*(1-bias_fraction), Temperature)
      enddo

    
! Hamiltonians on both charge states:
       open (unit_impurity, file=impurity_hamiltonian_file, status='old')

       read (unit_impurity,*) epsilon_0; epsilon_0=epsilon_0/hartree
       close (unit_impurity)

! add Magnetic field, only spin 1/2... boring...
       do i_spin = 1, N_spin
       Hamiltonian (i_spin) = epsilon_0 + B_au*(2*i_spin-1)*0.5 
       enddo


! Initialize Gr to start self-consistence loop

       do i_spin = 1, N_spin
       Gr (i_spin, :) = 1./(omega(:) - Hamiltonian (i_spin) + eta)
       enddo


! Initialize Dr although not needed (just for the
! convergence test)

      Dr (:) = 1./(omega(:) +  eta)

   return
   end subroutine initialize_quantities
!
! Fermi occupation function
!
    function Fermi (e, T)
      real (q) :: Fermi, e, T
      if ( e > 0._q ) then
          Fermi = exp(-e/T)/(exp(-e/T)+1._q)
       else
         Fermi = 1._q/(exp(e/T)+1._q)
       endif
      return
    end function Fermi
end module initialize
