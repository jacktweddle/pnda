# Pulsed Neutron Die Away Experiment in the NILE bunker

## This is the primary project I am working on during my placement year in the Neutronics group at the ISIS Neutron and Muon Source 2021/22. 

In the field of neutron moderators, obtaining a strong understanding of the scattering mechanisms is essential for efficient and useful moderator design. A key aspect of this is the thermal neutron scattering kernel, *S(&alpha;, &beta;)*. This is a fundamental aspect of thermal neutron scattering laws and is produced by theory, but must be verified experimentally in order to obtain confidence in its validity. One such method of validation is via pulsed neutron die away (PNDA) experiments.

We utilise a neutron generator employing either the DD (2.5 MeV) or DT (14.1 MeV) reaction to produce a neutron pulse a few &mu;s or ms in length and strike a cylinder of moderating material. The neutrons establish spatial modes of flux, each with a characteristic decay time. After a certain waiting time determined, the fundamental mode is dominant and the pulse decays via a single exponential channel. By using least squares regression to fit to the decat curve, we can determine the decay constant, &alpha;, given by:

*&alpha;* = *v&Sigma;<sub>a</sub> + DB<sub>0</sub><sup>2</sup> + CB<sub>0</sub><sup>4</sup>*

Varying the geometric buckling, B<sub>0</sub><sup>2</sup>, by altering the geometry of the moderator allows for the determination of the parameters *D* and *C* and thus the valiation of thermal scattering laws. This experiment utilises a stackable polyethylene cylinder to easily vary the buckling and is being conducted to determine the feasibility of conducting such PNDA experiments in the NILE bunker. If successful, the facility could be used for future experiments to characterise, design and test moderators at ISIS.

This repository contains simulations of the PNDA experiment using both OpenMC and FETCH, in the hopes of also using the experiment to benchmark these new neutron transport codes. 

### OpenMC

Initial OpenMC simulations include recreating MCNP simulations of a PNDA testbed conducted by Daniel Siefman et al. from LLNL in order to determine whether OpenMC is capable of performing PNDA simulations. Subsequent simulations then include a polyethylene cylinder in vacuum as a reference test case and full scale simulations of the NILE bunker. These simulations provide predictions for what to expect during the experiment and can be used to investigate background contributions, room return and the potential for additional shielding. They also act as a test of the newly implemented time dependent functionality.

### FETCH

FETCH simulations have thus far been used to perform full scale simulations of the NILE bunker and used to obtain flux and dose profiles to compare with MCNP simulations. Current work is being undertaken to investigate time dependent functionality as well as implementing a point source to determine if it is suitable to perform PNDA simulations in FETCH

## In this repository

## References




