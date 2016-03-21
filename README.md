## Resonances in quantum elastic scattering: a time-dependent numerical investigation
In my bachelor thesis, I investigate quantum scattering processes by numerically solving the time-dependent Schr ̈odinger equation, in the specific case of the hydrogen-krypton scattering. I discuss how scattering resonances are related to an increase of the probability, observed at some specific energies, of finding the scattered particle close to the interaction center.

##
*Thesis.pdf* contains the thesis text.

##
**Source code**
- most of the code is written in C (requires GSL lib)
- plotting is done with matplotlib

##
**Repo structure**
- *source/wkb* contains the algorithms described in Chapter 2 of the thesis and the WKB approximation
- *source/extract_solution* contains the algorithms described in Chapter 3 the time-dependent solution of the Schroedinger equation
- *source/basic* contains basic numerical routines in C (e.g. integral evaluation, roots finding, differential equations solving)
