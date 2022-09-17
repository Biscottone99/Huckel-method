# Huckel-method
This program allows you to study a generic π-conjugated molecule with extended Huckel model. In output (file output.dat) the program gives you the following information:
1. Huckel Hamiltonian;
2. Eigenvalues;
3. Eigenvectors;
4. Ground state energy;
5. Ground state risonance energy;
6. First excited state energy;
7. Exitation energy from gs to first excited state;
8. Ground state atomic charges;
9. First excited state atomic charges;
10. π Bond orders;
11. Permanent dipole moment;
12. Transition dipole moment;
13. Single excitation transition with their energy.
14. Absorption spectra.

All this information are given to you as energetic units in function of alpha and beta. You can obtain all the results in real energy value by changing the parameters "ac" (carbon value of α) and "bc" (carbon value of β). β values are always negative and the program already consider it: to put in the spectroscopic or thermodinic β, only absolute value has to be considered.

Hamiltonian matrix parameter values are calculated with the following formulae [Quantum chemistry - Levine - 5th edition]:

α_X = α_c + h_x * β_cc
β_xy = k_xy * β_cc

To use this program you have to write an input file (name 'input.dat', i. e. uploaded file) whith the following information:
- In the firt row you have to put in: numerber of atoms in π system, number of p-electrons, number of bonds in the molecule. Between each information you have to leave a space and not to put a [,].
- For every atom of your molecule you have to put in these information in the given order: atom number, h_x value, number of p-electrons, x coordinate, y coordinate, z coordinate. One row contains the information of a single atom.
- For every bond you have to put in these information in the following order: number of first atom in the bond, number of first atom in the bond, k_xy value. One row contains the information of a single bond.

Some k_xy and h_x values can be found in Quantum chemistry - Levine - 5th edition.
