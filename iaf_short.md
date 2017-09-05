---
title: "Towards a reliable, automated method of individual alpha frequency (IAF) quantification"
header-includes:
- \usepackage{lineno}
- \linenumbers
- \usepackage{setspace}
- \doublespacing
geometry: margin = 2.5cm
output:
  pdf_document:
    toc: no
    number_sections: yes
bibliography: libraryAC.bib

---
\begin{center}
\textbf{Andrew W. Corcoran\textsuperscript{a,b}, Phillip M. Alday\textsuperscript{c,b},}
\linebreak
\textbf{Matthias Schlesewsky\textsuperscript{b}, Ina Bornkessel-Schlesewsky\textsuperscript{b}}
\linebreak
\end{center}

\begin{flushleft}
\textsuperscript{a}Cognition and Philosophy Laboratory, Monash University, Melbourne, Australia.
\linebreak 
\textsuperscript{b}Centre for Cognitive and Systems Neuroscience, University of South Australia, Adelaide, Australia.
\linebreak 
\textsuperscript{c}Max Planck Institute for Psycholinguistics, Nijmegen, 6500 AH, The Netherlands.
\linebreak
\linebreak
\textbf{Correspondence} Andrew Corcoran, Cognition and Philosophy Laboratory, Room E672, Menzies Building, 20 Chancellors Walk, Monash University, Clayton, VIC 3800, Australia.
\linebreak
Email: andrew.corcoran1@monash.edu
\linebreak
Telephone: +61 (3) 9905 9166 
\linebreak
\linebreak
\textbf{Funding} This work was partially supported by funding from the University of South Australia Ehrenberg-Bass Institute for Marketing Science. 
This funding supported AC while he collected and analysed the empirical EEG dataset reported in this manuscript. 
The Institute had no influence on the design, analysis, or interpretation of the reported study.
AC is supported by an Australian Government Research Training Program (RTP) scholarship.
IBS is supported by an Australian Research Council Future Fellowship (FT160100437).
\end{flushleft}

\pagebreak

\begin{center}
\textbf{Abstract}
\end{center}
Individual alpha frequency (IAF) is a promising electrophysiological marker of interindividual differences in cognitive function. 
IAF has been linked with trait-like differences in information processing and general intelligence, and provides an empirical basis for the definition of individualised frequency bands. 
Despite its widespread application, however, there is little consensus on the optimal method for estimating IAF, and many existing approaches are prone to various sources of bias and/or inconsistency. 
We describe an automated method for deriving the two most common estimators of IAF: peak alpha frequency (PAF) and centre of gravity (CoG).
These indices are calculated from resting-state power spectra that have been smoothed by a Savitzky-Golay filter (SGF). 
We evaluated the performance characteristics of this SGF analysis routine in both empirical and simulated EEG datasets. 
Application of the SGF technique to resting-state data from $n$=63 healthy adults resulted in 61 PAF, and 62 CoG estimates.
The statistical properties of these estimates were consistent with previous studies. 
Analysis of simulated electrophysiological signals revealed that the automated SGF routine reliably extracts target alpha components, even under relatively noisy spectral conditions. 
The routine consistently outperformed a simpler method of automated peak localisation that did not involve spectral smoothing. 
The SGF technique is fast, open-source, and available in two popular programming languages (MATLAB and Python), and thus can easily be integrated within the most popular M/EEG toolsets (EEGLAB, FieldTrip and MNE-Python). 
As such, it affords a convenient opportunity for improving the reliability and replicability of future IAF-related research. 
\linebreak
\begin{flushleft}
\textbf{Keywords:} Alpha Rhythm, EEG, Oscillation/Time Frequency Analyses, Savitzky-Golay Filter, Individual Alpha Frequency
\end{flushleft}
\pagebreak

# Introduction
Alpha is the dominant rhythm in the human EEG, and its importance for cognitive processing has been recognised since Hans Berger's seminal work in the early 20th century [@berger1929; cf. @adrian1934].
Interindividual differences in the predominant frequency of alpha band oscillations (i.e. individual alpha frequency; IAF) has been linked with variability in cognitive performance since the 1930s [see @vogel1964; for a more recent review, see @klimesch1999].
IAF is a trait-like characteristic of the human EEG [@grandy2013a], which shows high heritability [@lykken1974; @malone2014; @smit2006] and test-retest reliability [@gasser1985; @kondacs1999; @naepflin2007].
There exists a considerable body of evidence that IAF predicts performance across a variety of perceptual [e.g., @cecere2015; @samaha2015] and cognitive [e.g., @bornkessel2004; @klimesch2006] tasks.
Individuals with a low IAF process information more slowly [@klimesch1996b; @surwillo1961; @surwillo1963], and show reduced performance on memory tasks [@klimesch1999] and general intelligence measures [*g*; @grandy2013], in comparison to their high-IAF counterparts.
IAF decreases with age from young adulthood onwards [@chiang2011; @kopruner1984], hence lifelong changes in IAF accompany the decline of many cognitive abilities in older adulthood [e.g. @hedden2004; @salthouse2011].
Taken together, this evidence suggests that IAF constitutes a promising neurophysiological marker of certain fundamental properties of central nervous system information processing [@grandy2013; @grandy2013a].

Beyond the quantification of interindividual variation in the dominant alpha rhythm per se, IAF can also be used to derive individualised estimates of the canonical frequency bands [@klimesch2012].
Such empirically-driven approach to the definition of frequency band regions may help to sharpen the precision of frequency-domain analyses more broadly [@klimesch2012].
Indeed, taking the IAF as a marker for distinguishing subregions of the alpha band has revealed functional dissociations between distinct alpha-rhythms [@klimesch1997].
Despite the potential advantages of using the IAF as an anchor-point for individualised spectral analysis, however, no consensus currently exists as to the optimal method of IAF quantification.
This paper thus sets out to develop a rigorous, automated strategy for estimating two of the most widely reported indices of IAF: peak alpha frequency (PAF) and alpha frequency centre of gravity (CoG).
We begin by briefly describing common strategies for extracting these estimators, and some of their attendant problems.

## Peak alpha frequency
IAF estimation typically depends on the delineation of a singular, prominent spectral peak within the alpha bandwidth [standardly defined as 8-13 Hz; @noachtar2004].
In many cases, PAF can be easily identified on visual inspection of the power spectral density (PSD) derived from eyes-closed resting-state EEG [particularly from channels over the centro-posterior region of the scalp; @sadaghiani2016].
However, this strategy is complicated when participants manifest two (or more) peaks within the alpha band [so-called 'split-peaks'; @chiang2011], or fail to demonstrate any obvious deviation from the characteristic $1/f$-like power spectrum scaling of resting-state M/EEG [the 'inverse power-law'; @pritchard1992].
Under such circumstances, subjective PAF estimation may be prone to biased and inconsistent analysis [@chiang2008], thus posing a significant challenge to replicability.
Conservative approaches to subjective PAF estimation in the context of ambiguous spectral conditions may also result in high rates of attrition if peaks cannot be confidently delineated [see for e.g., @bornkessel-schlesewsky2015].

One approach for improving the rigour, objectivity, and (for larger datasets) practicality of PAF estimation is to implement a peak-detection algorithm.
While this strategy does not solve the basic problem of deciding the criteria by which valid PAF estimates are discriminated from split-peaks or spurious fluctuations, it at least applies such criteria consistently across all subjects.
Simple algorithms may however introduce new sources of bias.
For instance, a basic routine that searches for local maxima within the alpha band may arbitrarily assign the PAF to the lower bound of the search window in the absence of any distinct deviation from the inverse-power law.
A more sophisticated implementation [such as searching the first derivative of the PSD for downward-going zero-crossings; e.g. @grandy2013a] obviates this error, but is incapable of distinguishing a substantive peak from split-peak or arbitrarily small deviation from background spectral activity.
Such routines may therefore be too liberal in the estimates that class as spectral peaks.

## Alpha frequency centre of gravity
Klimesch, Schimke, and Pfurtscheller [-@klimesch1993; see also @klimesch1997] proposed the CoG estimator as an alternative to PAF that circumvents some of the difficulties posed by the absence of a dominant alpha peak.
CoG takes into account the distribution of PSD estimates within the defined alpha interval, and may therefore be subject to bias if the bounds of alpha-band are specified inaccurately.
Since individuals show variation in the span and location of alpha-band activity, Klimesch and colleagues [@klimesch1990] recommended computing CoG on the basis of bespoke frequency windows designed to capture this interval.
However, the definition of such individualised alpha-band windows (IAWs) poses a nontrivial challenge, and may depend on subjective assessments [e.g., @klimesch1990] or arbitrary criteria [e.g., @klimesch1999].
One principled solution is to this problem is to derive the IAW from reactivity-based contrasts between two conditions [e.g., eyes-closed vs. eyes-open resting-state EEG, @klimesch1999; pre- vs. peri-stimulus presentation, @goljahani2012].
This approach is also vulnerable to bias, however, since alpha rhythms are not always substantially attentuated by opening the eyes [@gaal2010; @kreitman1965].
The picture is even more complicated when relying on stimulus-induced reactivity, since alpha rhythms might only be selectively attenuated [e.g., @klimesch2006] -- or indeed, enhanced [e.g., @rihs2007] -- during experimental tasks.

## Curve-fitting approaches to alpha-rhythm quantification
One promising approach to spectral peak quantification exploits iterative curve-fitting techniques to parameterise the statistical properties of the PSD [@chiang2008; @lodder2011].
The practical advantages of such methods is clearly evident from their application to large $n$ datasets [e.g., @chiang2011; @van_albada2013], while comparison of Lodder and Putten's [-@lodder2011] algorithm with subjectively-assessed PAF estimates indicated a high degree of agreement.
It is puzzling then why such methods have not been taken up more widely in the IAF literature [cf. @haegens2014, for a notable exception].
One possibility is that investigators are generally unaware of such techniques, given that they have mostly been applied in the domain of spectral modeling [indeed, neither @goljahani2012, nor @bazanova2014, mention the existence of these methods in their reviews of IAF estimators].
Alternatively, investigators may be put off by the practical burden involved in accessing these programmes (which we have not been able to locate publically), and integrating them within existing analysis pipelines (which may be incompatible with the algorithm).
We suggest that one of the critical steps towards achieving a more widespread adoption of automated IAF estimation routines is to make these tools as openly available as possible, in formats that are easily assimilated within popular methods of EEG data analysis.

## Aims of the present study
In sum, common methodological approaches to IAF estimation are either time-consuming and vulnerable to inconsistencies arising from subjective interpretation, or at risk of producing spurious or biased estimates under certain plausible spectral conditions.
More recent innovations that address these problems via the application of sophisticated curve-fitting algorithms have so far found limited uptake within the broader IAF literature, perhaps on account of practical barriers pertaining to software access and implementation. 
Consequently, we sought to develop an automated method of alpha-band quantification that provides fast, reliable, and easily replicated estimates of resting-state IAF in two major programming languages: MATLAB^&reg;^ (The MathWorks, Inc., Natick, MA, USA) and Python&trade;. 
This goal is consistent with recent proposals to make the analysis of cognitive electrophysiological data as open, transparent, and amenable to replication as possible [@cohen2017].

# Method
Our approach aims to emulate Klimesch and colleagues' [-@klimesch1990] original attempt to characterise individual profiles of resting-state alpha-band activity by means of a relatively simple, non-parametric curve-fitting technique.
Our basic strategy runs as follows:
First, we extract PSD estimates from preprocessed and fast Fourier-transformed EEG signals.
Second, we apply a curve-fitting procedure (Savitzky-Golay filtering) to smooth the PSD function and estimate its first- and second-order derivatives.
Third, these derivatives are analysed for evidence of a distinct spectral peak (i.e. PAF) within the alpha band.
Finally, the first derivative of the PSD is reanalysed to delineate upper and lower bounds of the IAW, from which the CoG is calculated.

## Savitzky-Golay smoothing and differentiation



\begin{center}
\textbf{Acknowledgements}
\end{center}
\begin{flushleft}
We thank Jessica Gysin-Webster and Daniel A. Rogers for their assistance with data collection.
\end{flushleft}
\pagebreak
\begin{center}
\textbf{References}
\end{center}
