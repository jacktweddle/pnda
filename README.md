# Pulsed Neutron Die Away Experiment in the NILE bunker

## This is the primary project I am working on during my placement year in the Neutronics group at the ISIS Neutron and Muon Source 2021/22. 

In the field of neutron moderators, obtaining a strong understanding of the scattering mechanisms is essential for efficient and useful moderator design. A key aspect of this is the thermal neutron scattering kernel, *S(&alpha;, &beta;)*. This is a fundamental aspect of thermal neutron scattering laws and is produced by theory, but must be verified experimentally in order to obtain confidence in its validity. One such method of validation is via pulsed neutron die away (PNDA) experiments.

We utilise a neutron generator employing either the DD (2.5 MeV) or DT (14.1 MeV) reaction to produce a neutron pulse a few &mu;s or ms in length and strike a cylinder of moderating material. The neutrons establish spatial modes of flux, each with a characteristic decay time. After a certain waiting time determined, the fundamental mode is dominant and the pulse decays via a single exponential channel. By using least squares regression to fit to the decat curve, we can determine the decay constant, &alpha;, given by:

*&alpha;* = *v&Sigma;<sub>a</sub> + DB<sub>0</sub><sup>2</sup> + CB<sub>0</sub><sup>4</sup>*

By varying the geometric buckling, B<sub>0</sub><sup>2</sup>
