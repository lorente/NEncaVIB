# NEncaVIB
Non-Equilibrium Non-Crossing Approximation with VIBrations

This code uses the equations by Ned S. Wingreen and Yigal Meir
Phys. Rev. B 49, 11040 – Published 15 April 1994 to account for the non-equilibrium Kondo effect between two electrodes. The impurity is in the present version treated as a single orbital and its total spin si fixed to 1/2. The code is very efficient and converges with great accuracy by using convolutions coded in FFT (needs fftw3).
The code uses the self-consisten Born Approximation to deal with a local vibration of the impurity. As is done, we can add as many vibrations as needed. This follwos the paper by P. Roura-Bas, L. Tosi, and A. A. Aligia
Phys. Rev. B 93, 115139 – Published 24 March 2016.
We have run the code using the Kondo efect with various e-vib couplind and it converges excedeengly well even in the strong coupling regime. Moreover the code is totally able to recover the vibronic regime in the presence of a bias drop, which is very interesting.

The input is in the file WMnca.input. And example is:

-5.0  !Omega_ini (eV) (the bandwidth extends -2*Omega_ini)

0.00001 !step_omega (eV)

0.0001 !Broadening eta (eV)

0. !B_field (Teslas)

2.0 !gyromagnetic factor

5.0 !Temperature (K)

0.01 !Convergence tolerance

20 !Maximum number of loops

-0.5 !Bias ini (Volts)

0.5 !Bias fin (Volts)

1 !Number of bias

1.0 !Fraction of bias drop at the tip (1-Tip drops at the sample)

Hamiltonian.dat !name of the file containing the level energy (rpt Fermi)

Gamma.dat !file with Go, Omega_0, Delta_0 for tip and substrate

Current.dat !file with V and Current OUTPUT

9 !N_vib number of modes

8.3 ! Freq meV

8.8 ! Freq meV

20.26 ! Freq meV

20.48 ! Freq meV

26.72 ! Freq meV

31.47 ! Freq meV

35.82 ! Freq meV

42.1 ! Freq meV

43.66 ! Freq meV

1.0 ! e-ph coupling meV

1.0 ! e-ph coupling meV

2.0 ! e-ph coupling meV

2.0 ! e-ph coupling meV

1.0 ! e-ph coupling meV

1.0 ! e-ph coupling meV

5.0 ! e-ph coupling meV

10.0 ! e-ph coupling meV

10.0 ! e-ph coupling meV

.true. ! convergence loops on phonon?

the code uses two more files:

The Gamma:dat file that reads as:

G  ! two options G (Gaussian= or L (Lorentzian) for the Hybridisation function
0.0 0 1.0  ! G0, E0, D0 for left electrode (eV)
0.1 0 2.0   ! same for right electrode 8eV)

where if it is a gaussina the hybridization function looks like G0*exp (-(omega-E0)^2/D0^2), a Lorentzian like G0/(omega-E0)^2+D0^2)
so that D0 is like an effective bandwidth and G0 is the value of the hybridisation (times 2 pi would be the width of the level) at
the energy E0.

The other file is the Hamiltonian that here is just a level referred to the Fermi energy.
Hamiltonian.dat ::

-0.65

in eV.


ENJOY
