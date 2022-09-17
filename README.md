# Huckel-method
This program allows you to study a generic π-conjugated molecule with extended Huckel model. In output (file output.dat) the program gives you the following information:
1. Huckel Hamiltonian;
2. Eigenvalues;
3. Eigenvectors;
4. Ground state energy;
5. Risonance energy for ground state;
6. First excited state energy;
7. Exitation energy from gs to first excited state;
8. Ground state atomic charges;
9. First excited state atomic charges;
10. π Bond orders;
11. Permanent dipole moment;
12. Transition dipole moment;
13. Single excitation transition with their energy.
14. Absorption spectra.
All this information are given to you as energetic units in fucntion of alpha and beta. You can obtain all the results in real energy value, by changing the parameters "ac" (carbon value of alpha) and "bc" (carbon value of beta). Usually values in always negative and the program already consider this feature, to put in the spectroscopic or thermodinic β only absolute value as to be considered.

Hamiltonian matrix parameter values are calculated with the following formulae [Quantum chemistry - Levine - 5th edition]:
α_X = α_c + h_x * β_cc
β_xy = k_xy * β_cc

To use this program you have to write an input file (name 'input.dat', i. e. uploaded file) whith the following information:
- In the firts row you have to give: numerber of atom in π system, number of p-electrons, numerber of bonds in the molecule. Between each information you have to leave a space and not to put a comma.
- In different lines you have to write these information in the following order: atom number, h_x value, number of p-electrons, x coordinate, y coordinate, z coordinate.
- In different lines you have to write these information in the following order: number of first atom in the bond, umber of first atom in the bond, k_xy value.

Some k_xy and h_x values can be found in Quantum chemistry - Levine - 5th edition.
