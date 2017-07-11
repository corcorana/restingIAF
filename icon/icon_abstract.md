---
title: 'restingIAF: A reliable, automated, open source method for quantifying individual alpha frequency'

author:
  - Andrew W. Corcoran
  - Phillip M. Alday
  - Matthias Schlesewsky
  - Ina Bornkessel-Schlesewsky

output: 
  pdf_document

bibliography: 
  libraryAC.bib
---

## Conference abstract
### Background:
Individual alpha frequency (IAF) is a promising electrophysiological marker of interindividual differences in cognitive function [@grandy2013; @grandy2013a].
In particular, IAF has been shown to predict performance across a variety of psychophysical and cognitive tasks [@bazanova2014; @cecere2015], and may underpin trait-like differences in information processing [@klimesch1996b] and general intelligence [@grandy2013].
IAF has also been cited as a useful anchor point for determining individually-tailored frequency band analyses [@klimesch2012].
Despite this large body of literature, however, there seems to be no clear consensus on the optimal method for estimating IAF.
We therefore sought to develop a reliable, automated method for IAF estimation that could be easily integrated within existing analysis pipelines.

### Method:
We implemented a method of calculating two common IAF estimators (peak and gravity frequency) in MATLAB and Python (both available on GitHub; Python implentation also in the `philistine` package on PyPi).
This `restingIAF` routine locates the bounds of the dominant alpha component according to the first and second derivatives of its Savitzky-Golay smoothed spectral density.
We evaluated its performance characteristics in both empirical and simulated EEG datasets.

### Results:
`restingIAF` generated 61 IAF estimates from 63 healthy adults, thus yielding a higher proportion of estimates than both standard analysis [e.g., @bornkessel-schlesewsky2015] and Gaussian curve-fitting techniques [e.g, @haegens2014].
The distribution of IAF estimates was consistent with that derived using more complex algorithms within large-scale datasets [e.g., @chiang2011].
Preliminary analysis of simulated data revealed that `restingIAF` accurately extracts the peak frequency of underlying alpha-band components, even when signal-to-noise ratio is highly degraded.

### Conclusion:
`restingIAF` is a novel method of extracting IAF estimates that can easily be applied to large datasets.
It is fast, reliable, open source, and available in two popular programming languages, and thus easily integrated in the most popular M/EEG toolsets (EEGLAB, FieldTrip and MNE-Python).
Widespread adoption might significantly improve the consistency and replicability of future IAF-related research.

## Theme: *Electrophysiology methods*

## References
