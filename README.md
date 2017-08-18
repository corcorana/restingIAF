# restingIAF
Source code for `restingIAF`, an automated resting-state individual alpha frequency (IAF) estimation routine implemented in EEGLAB. 
Repo also contains a manuscript (submitted) outlining rationale for programme development and its performance on simulated and non-simulated datasets.
An example analysis script (including real EEG data) is also [provided](https://github.com/corcorana/restingIAF/tree/master/code/tutorial).
We recommend ensuring you can replicate the analysis demonstrated in the tutorial before implementing the package in your own analyses.

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
Release `v0.1.0` formed the basis of Cross et al's (manuscript in preparation) resting-state IAF analysis. 
Additional analysis scripts customised for the purposes of this study are included in release `v0.1.1`.
