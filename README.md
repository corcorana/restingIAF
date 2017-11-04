# restingIAF
## General information
Source code for `restingIAF`, an automated resting-state individual alpha frequency (IAF) estimation routine implemented in EEGLAB. 
Repo also contains a manuscript (under revision) outlining rationale for programme development and its performance across simulated and non-simulated EEG datasets.
Materials required to replicate the simulation analyses reported in the manuscipt can be found in the [simulations](https://github.com/corcorana/restingIAF/tree/master/code/simulations) subdirectory.
A [tutorial](https://github.com/corcorana/restingIAF/tree/master/code/tutorial) analysis routine (including sample EEG data) is provided to familiarise you with various features of the package.
If you experience any difficulties implementing `restingIAF`, please attempt to replicate the analysis detailed in the [tutorial readme](https://github.com/corcorana/restingIAF/tree/master/code/tutorial/tute_README.md) before contacting us for assistance.

NOTE: All functions were developed and tested in MATLAB 2014b/2015a, and may not be backward compatible with older versions of MATLAB. The algorithm is also dependent on the Signal Processing Toolbox.

Software developed in collaboration with Dr. Phillip Alday, Prof. Matthias Schlesewsky, & Prof. Ina Bornkessel-Schlesewsky at the Centre for Cognitive and Systems Neuroscience, University of South Australia.

Constructive criticism, corrections, and potential improvements are most welcome.
E-mail: andrew.corcoran1\@monash.edu

## Python implementation
A preliminary implementation in Python for use with [MNE-Python](https://martinos.org/mne/) is available as part of the [`philistine` package](https://gitlab.com/palday/philistine).

## Citation
If this software was useful for your work, please consider citing our methods paper: 

> Corcoran, A.W., Alday, P.M., Schlesewsky, M., & Bornkessel-Schlesewsky, I. (2017). Towards a reliable, automated method of individual alpha frequency (IAF) quantification. bioRxiv. doi: https://doi.org/10.1101/176792

The current release version can be cited as:

> Corcoran, A.W., Alday, P.M., Schlesewsky, M., & Bornkessel-Schlesewsky, I. (2017). restingIAF (version number) [Software]. Retrieved from https://github.com/corcorana/restingIAF. [doi]

Each version release is issued a unique doi from [Zenodo](https://zenodo.org/record/846797#.WbEh4a2B2fQ).


## Releases
`v1.0.1` [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.888071.svg)](https://doi.org/10.5281/zenodo.888071)

Includes IAW estimates in summary output.
Also includes additional tutorial dataset showing output when 1 recording fails to satisfy $cMin$ for PAF estimation, tutorial script contingency for handling epoched data. 

`v1.0` [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.846797.svg)](https://doi.org/10.5281/zenodo.846797)

Released in conjunction with publication of Corcoran et al. bioRxiv preprint.

`v0.1.0` [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.268602.svg)](https://doi.org/10.5281/zenodo.268602)

Formed the basis of Cross et al's (manuscript in preparation) resting-state IAF analysis. 
Additional analysis scripts customised for the purposes of this study are included in release `v0.1.1`.
UPDATE (07/09/2017): The analysis for this paper was repeated using `v1.0` prior to submission.

## Licence & Copyright
This software was developed in the hope that it would be of some use to the EEG research community, and is freely available for redistribution and/or modification under the terms of the GNU General Public Licence. 
It is distributed WITHOUT WARRANTY; without even the implied warranty of merchantability or fitness for a particular purpose. 
See the GNU General Public License (LICENCE.md) for more details.

Copyright (c) 2016-2017 Andrew W. Corcoran.
(Any coding errors or inelegancies are solely his fault.)
