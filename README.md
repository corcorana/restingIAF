# restingIAF
Source code for `restingIAF`, an automated resting-state individual alpha frequency (IAF) estimation routine implemented in EEGLAB. 
Repo also contains a manuscript (submitted) outlining rationale for programme development and its performance across simulated and non-simulated EEG datasets.
Materials required to replicate the simulation analyses reported in the manuscipt can be found in the [simulations](https://github.com/corcorana/restingIAF/tree/master/code/simulations) subdirectory.
A [tutorial](https://github.com/corcorana/restingIAF/tree/master/code/tutorial) analysis routine (including sample EEG data) is provided to familiarise you with various features of the package.
If you experience any difficulties implementing `restingIAF`, please attempt to replicate the analysis detailed in the [tutorial readme](https://github.com/corcorana/restingIAF/tree/master/code/tutorial/tute_README.md) before contacting us for assistance.

NOTE: All functions were developed and tested in MATLAB 2014/2015, and may not be backward compatible with older versions of MATLAB. The algorithm is also dependent on the Signal Processing Toolbox.

A preliminary implementation in Python for use with [MNE-Python](https://martinos.org/mne/) is available as part of the [`philistine` package](https://gitlab.com/palday/philistine).

Software developed in collaboration with Dr. Phillip Alday, Prof. Matthias Schlesewsky, & Prof. Ina Bornkessel-Schlesewsky, University of South Australia Cognitive Neuroscience Laboratory.

Copyright (c) 2016-2017 Andrew W. Corcoran.
(Any coding errors or inelegancies are solely his fault.)

This software was developed in the hope that it would be of some use to the EEG research community, and is freely available for redistribution and/or modification under the terms of the GNU General Public Licence. 
It is distributed WITHOUT WARRANTY; without even the implied warranty of merchantability or fitness for a particular purpose. 
See the GNU General Public License (LICENCE.md) for more details.

Constructive criticism, corrections, and potential improvements are most welcome.
E-mail: andrew.corcoran1\@monash.edu

## Citation
If this software was useful for your work, please cite the following paper: 

> Corcoran, A.W., Alday, P.M., Schlesewsky, M., & Bornkessel-Schlesewsky, I. (2017). Towards a reliable, automated method of individual alpha frequency (IAF) quantification. bioRxiv. doi: https://doi.org/10.1101/176792

## Releases
`v1.0` was released in conjunction with publication of the bioRxiv preprint.

Release `v0.1.0` formed the basis of Cross et al's (manuscript in preparation) resting-state IAF analysis. 
Additional analysis scripts customised for the purposes of this study are included in release `v0.1.1`.
