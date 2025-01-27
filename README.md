# **Microscopic many body calculations of spherical nuclei with modern NN-NNN forces - Spherical HF-(E)RPA solver package**

Repository with all codes & scripts related to my master thesis taken at CUNI Faculty of Mathematics & Physics in Prague at Insitute of Particle & Nuclear Physics (IPNP). Work is supervised by doc. Mgr. František Knapp PhD. from IPNP.

At present time all the scripts use Julia language. _Code should be work fine on all distributions v >= 1.4_

<br/>

Most of the code is tested & produces reasonable results. There may be some minor tweaks though.

<br/>

## **Running code requires these Julia packages:**

* _DelimitedFiles_

* _CGCoefficient_

* _LinearAlgebra_

## **Repository includes:**

* Spherical Hartree-Fock solver for Even-Even nuclei.

* HF based Leading Order Many-Body Perturbation Theory calculation.

* Spherical Hartree-Fock (Extended) Random-Phase Approximation & Tamm-Dancoff Approximation solver for Even-Even nuclei.

<br/>

Feel free to use & improve the code. Just please credit us whenever you explicitly use parts of the code.

<br/>

Later, when most of code is finished & I have more time, I will add more detailed description & a pdf manual on how to run & use the code, furthermore I will also include the full text of the thesis & likely include all data & scripts used for plots.

<br/>

Note there are small configuration space interactions files binaries included (DN2LOGO394 chiral potential, N_max = 3, N_3max = 9, hw = 16 MeV). In principle, you can run & play with the code in this small test space.
