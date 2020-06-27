program WNncaVIB
! This is a Keldysh implementation
! of nca for a single electronic level
! only spin degeneracion
!
! we add vibrational modes
!
! Notation:
!        D< is Dl; D> is Db; Dr is Dr
!
! Following Wingreen and Meir PRB 49, 11040 (1994)
!
! Date
!      13 May 2020

Use declarations
Use used_files
Use present_version
Use algebra
Use electronic
Use selfenergies
Use green
Use initialize
Use io

                  implicit none

      print *, " "
      print *, "RUNNING non-equilibrium NCA wth vibrations"
      print *, " "
      print *, "Version::", version, "Date::", DateVersion
      print *, " "
! read input file (WNnca.input):
      call  input (omega, step_omega, Temperature, eta, N_spin, &
      &        tol, Max_loop_index, B_au, gamma_file, output_file, &
      &        impurity_hamiltonian_file, &
      &        Bias_ini, Bias_fin, N_bias, &
      &        bias_fraction, N_vib, Freq, W, vib_loop )

! Compute Bose-Einstein distribution
          allocate (Bose (N_vib))
        do i = 1, N_vib
          Bose (i) = 1._q/(exp(Freq(i)/Temperature)-1._q)
        enddo
! change Frequencies into indices of the Omega array
         allocate (indexvib (N_vib))
        do i = 1, N_vib
          indexvib (i) = int ((Freq(i)-omega(1))/step_omega)+1 &
      &                - int (-omega(1)/step_omega)-1  !from zero energy
        enddo 


    open (unit_pdos, file = 'PDOS_bias.dat')
    open (unit_ddos, file = 'dPDOS_bias.dat')
    open (unit_current, file = output_file)
    open (unit_conductance, file = 'conductance.dat')

! Allocation of main arrays that have not been allocated on input
      N_omega = size (omega, 1)
     call memory (N_omega, N_spin, Gamma_t, Gamma_s, Fermi_t, Fermi_s,   &
      &   Dr, Gr, SelfD, selfG, Dl, Gl, Db, Gb, Hamiltonian, &
      &   SelfDl, SelfGl, &
      &   SelfDb, SelfGb, GGR, Gr_store, Dr_store, Dl_store, Gl_store) 

BIASLOOP: do i_V=1, N_bias

        if (N_bias == 1) then
        bias = 0
        else
        bias = Bias_ini + (i_V-1._q)*(Bias_fin-Bias_ini) / (N_bias-1._q)
        endif
  
        if (i_V == 1) then 
        curr_store = 0._q
        else
        curr_store = curr
        endif

        print *, ' '
        print *, 'Bias iteration::', i_V, 'BIAS (V)::', bias*hartree
        print *, ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For each Bias                                       !
! initialize calculation of retarded greens functions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     call initialize_quantities (N_omega, N_spin, omega, gamma_file,  &
      & Gamma_t, Gamma_s, bias_fraction, &
      & Fermi_t, Fermi_s, bias, Temperature, impurity_hamiltonian_file,  &
      & Hamiltonian, eta, Gr, Dr) 
 
! Use previous bias converged GF for better convergence
!removed because this gives convergence issues with bias!!
!     if (i_V /= 1) then

!        Dr=Dr_store
!        Gr=Gr_store

!     endif

      print *, ' '
      print *, '1. Starting convergence of retarded GFs'
      print *, ' '
Self_consistency_R: do self_index = 1, Max_loop_index

      call Self_Pi (N_omega, N_spin, Gr, Gamma_t, Gamma_s,  &
        &                   Fermi_t, Fermi_s, omega, SelfD)

         Dr_store = Dr

      call D_green_R (omega, SelfD, Dr)

      call Self_Sigma (N_omega, N_spin, Dr, Gamma_t, Gamma_s, &
        &                      Fermi_t, Fermi_s, omega, SelfG)

         Gr_store = Gr

      call G_green_R (N_spin, omega, Hamiltonian, SelfG, Gr)
!
     if (N_bias == 1) then
     open (unit_screen, file="PDOS_fermion.dat")
     do i_omega = 1, N_omega
     write (unit_screen, *) omega (i_omega)*hartree, -aimag(sum (Gr (:,i_omega)))/(hartree*pi_d)
     enddo
     close (unit_screen)
     endif

! Are we converged?

      call convergence (N_omega, N_spin, Dr, Dr_store,  &
                               Gr, Gr_store, tol, converged)


               if (converged) exit



      enddo Self_consistency_R

              if (.not.converged) then
                  print *, "ATTENTION : le calcul n'est pas convergé !"
              endif

      call checkSumRules (omega, N_spin, Dr, Gr)

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  calculation of lesser and bigger greens functions  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      Db =2*ui*aimag (Dr)
      Gb =2*ui*aimag (Gr)

! initialize Dl , Gl  
! to something like the thermal one
! or use previous bias converged GF for better convergence

!Changed 6 Juin pour éviter "accumulation d'erreurs" entre V
! Change back removing !#
 !#    if (i_V == 1) then

           Dl = -2*ui*aimag(Dr)
         do i_spin =1, N_spin
           Gl (i_spin,:) =-2*ui*aimag(Gr(i_spin,:)) 
         enddo

 !#    else

!#         Dl=Dl_store
!#         Gl=Gl_store

 !#    endif

      print *, ' '
      print *, '2. Starting preconvergence of lesser GFs'
      print *, ' '

! PRE-converging


      do self_index = 1, Max_loop_index

       call Sigma_lesser0 (N_omega, N_spin, Dl, Gamma_t, Gamma_s, &
      &       fermi_t, fermi_s, omega, SelfGl)


           Gl_store = Gl
       call G_lesser (N_spin, Gr, SelfGl, Gl)

       call Pi_lesser (N_omega, N_spin, omega, Gamma_t, Gamma_s, &
      &      fermi_t, fermi_s, Gl, SelfDl)

           Dl_store = Dl
       call D_lesser (Dr, SelfDl, Dl)

      call convergence (N_omega, N_spin, Dl, Dl_store,  &
                               Gl, Gl_store, tol, converged)

               if (converged) exit
      enddo 


      print *, ' '
      print *, '3. Starting convergence of lesser GFs'
      print *, ' '

     if (vib_loop) then
Self_consistency_lesser: do self_index = 1, Max_loop_index

       call Sigma_lesser (N_omega, N_spin, Dl, Gamma_t, Gamma_s, &
      &     fermi_t, fermi_s, omega, N_vib, indexvib, W, Bose, Gl,  SelfGl)



           Gl_store = Gl
       call G_lesser (N_spin, Gr, SelfGl, Gl)

       call Sigma_bigger (N_omega, N_spin, Db, Gamma_t, Gamma_s, &
      &     fermi_t, fermi_s, omega, N_vib, indexvib, W, Bose, Gb, SelfGb)

       call G_bigger (N_spin, Gr, SelfGb, Gb)

! Changed 6 Juin
!      call Pi_lesser (N_omega, N_spin, omega, Gamma_t, Gamma_s, &
!     &      fermi_t, fermi_s, Gl, SelfDl)

!          Dl_store = Dl
!      call D_lesser (Dr, SelfDl, Dl)

!      call Pi_bigger (N_omega, N_spin, omega, Gamma_t, Gamma_s, &
!     &      fermi_t, fermi_s, Gb, SelfDb)

!      call D_bigger (Dr, SelfDb, Db)

      call convergence (N_omega, N_spin, Dl, Dl_store,  &
                               Gl, Gl_store, tol, converged)

               if (converged) exit

      enddo Self_consistency_lesser

              if (.not.converged) then
                  print *, "ATTENTION : le calcul n'est pas convergé"
              endif
      else !no consistence loops

       call Sigma_lesser (N_omega, N_spin, Dl, Gamma_t, Gamma_s, &
      &     fermi_t, fermi_s, omega, N_vib, indexvib, W, Bose, Gl,  SelfGl)

       call G_lesser (N_spin, Gr, SelfGl, Gl)

       call Sigma_bigger (N_omega, N_spin, Db, Gamma_t, Gamma_s, &
      &     fermi_t, fermi_s, omega, N_vib, indexvib, W, Bose, Gb, SelfGb)

       call G_bigger (N_spin, Gr, SelfGb, Gb)

            print *, ' '
            print *, '     No convegence loops for phonons, only 1st iteration!'
            print *, ' '
       endif
     if (N_bias == 1) then
     open (unit_screen, file="PDOS_Gl.dat")
     do i_omega = 1, N_omega
     write (unit_screen, *) omega (i_omega)*hartree, aimag(sum (Gl (:,i_omega)))/(2*pi_d*hartree)
     enddo
     close (unit_screen)
     endif

       call normalization (N_omega, Dl, Gl, constant, omega)

       call GreenR (N_omega, N_spin, constant, Dl, Db, Gl, Gb, GGR, omega)

! we write the PDOS for each bias
     if (N_spin == 1) then
       do i =1, N_omega
       write (unit_pdos, '(2g14.4)') omega (i)*hartree, &
     &         (-aimag (GGR(i_spin,i))/(pi_d*hartree), i_spin=1, N_spin)
       enddo
       do i =5, N_omega,5
       write (unit_ddos, '(2g14.4)') omega (i)*hartree, &
     &    -aimag (GGR(1,i)-GGR(1,i-5))/(pi_d*hartree*5*step_omega*hartree)
       enddo
     else
       do i =1, N_omega
       write (unit_pdos, '(3g14.4)') omega (i)*hartree, &
     &         (-aimag (GGR(i_spin,i))/(pi_d*hartree), i_spin=1, N_spin)
       enddo
       endif
       write (unit_pdos, *) 

       call current (Gamma_t, Gamma_s, GGR, fermi_t, fermi_s, &
     &       curr, omega, N_omega, N_spin)

       if (i_V == 1) curr_store=curr

! Final result

       write (unit_current, '(2g14.4)') bias*Hartree, curr*6.6236E6
       write (unit_conductance, '(2g14.4)') bias*Hartree,(curr-curr_store)/(Bias_fin-Bias_ini)*(N_bias-1._q)/Hartree*6.6236E6
       print *, "Bias (V) et Courant (nA)::", bias*Hartree, curr*6.6236E6
       print *, "__________________________________________________"

       enddo BIASLOOP

       stop

end program WNncaVIB
