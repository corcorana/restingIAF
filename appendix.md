---
title: Appendix
output:
  pdf_document:
    number_sections: yes
bibliography: libraryAC.bib
---
This appendix gives a more formal account of the estimators implemented in the algorithm, and technical details about the operations we use to derive them. 
We also provide some additional advice regarding SGF parameter selection.

# IAF indices
Peak alpha frequency (PAF) typically involves the visual identification of a singular, prominant peak within a predefined alpha-band range.
PAF can be formalised in terms of the local (i.e. relative) maximum within the alpha band

$$ \text{PAF} =  
\begin{cases}
\text{arg}\max_{f \in  \text{alpha band}} \text{PSD}(f), & \text{PSD}(f) - \text{PSD}(f') \geq \varphi \, \forall f' \neq f \\ 
\text{undefined}, & \text{otherwise}
\end{cases}
$$

where $\text{arg} \max$ returns the frequency bin (or subset of bins) $f$ containing the maximal power value $\max \text{PSD}(f)$ registered within the set of frequency bins constituting the alpha band.
Note that, for the output of $\text{arg max}$ to qualify as an estimate of PAF, it must return a single frequency bin $f$  with a corresponding power spectral density $\geq \varphi$, where $\varphi$ defines the minimum threshold value differentiating a substantive spectral peak from background noise.


Klimesch and colleagues [@klimesch1990; @klimesch1993] suggest using the PSD-weighted mean alpha frequency (i.e. alpha centre of gravity; CoG) to estimate IAF.
Mathematically, we can express the COG as

$$ \text{CoG} = \frac{\int_{f_1}^{f_2} \text{PSD}(f) \cdot f \; df}{\int_{f_1}^{f_2} \text{PSD}(f) \; df} , $$

where $f_1$ and $f_2$ index the frequency bins bounding the alpha-band interval.

# Savitzky-Golay curve-fitting procedure
SGFs work by centring a sampling window frame of length $F_w$ on a portion of the input signal and computing the least-squares fit of a specified polynomial to each $i$^th^ data point spanned by $F_w$.
The window is then shifted one point along the input signal, and the polynomial fit recalculated accordingly.
The centre value of the polynomial fit is taken as the filter output at each iteration of the sliding window calculation, and these output values are concatenated to render the smoothed estimate of the input function.
For a more detailed treatment of the SGF and its technical performance properties, see Schafer [-@schafer2011].

Both polynomial degree $k$ and filter window frame width $F_w$ are required to define the least-squares minimisation operation.
$k$ must be $< F_w$, and $F_w$ must be odd to ensure an equal number of sample points either side of the centre coefficient.
Also note that no smoothing will occur if $k = F_w - 1$.
A convenient heuristic is to set the length of $F_w$ approximately 1 to 2 times the anticipated full width at half maximum (FWHM) of the PAF component [@enke1976; @press1992].
Relatively higher $F_w$ lengths are expected to result in more aggressive smoothing of the input function [@bromba1981], and hence may be desirable in cases of noisy spectral densities.
Note however that excessively flat peaks following application of the smoothing procedure are indicative of a suboptimally large $F_w$.
We favour higher-order polynomials (e.g., $k = 5$) due to their peak-height preserving properties, but acknowledge that they might render suboptimal fits (and less smoothing) in the context of relatively broad component structures [@press1992].
This may be of particular concern when dealing (for instance) with elderly populations, given that spectral power (and relative peak dominance) is typically diminished in older adults [e.g., @chiang2011; @dustman1999].

# First- and second-derivative tests
Derivatives describe the relative rate of change in the dependent variable or function $g(x)$ given some value of the independent variable $x$.
The first derivative of a vector of PSD estimates thus provides point estimates of the (instantaneous) rate of change in the amount of spectral power estimated for each frequency bin resolved in the analysis.
This relationship can be formalised as

$$ g'(x) = \lim_{ \Delta{x} \rightarrow{} 0} \frac{\Delta g(x)} {\Delta x} , $$

where $g'(x)$ is the first derivative of the relative change in the power estimate $g(x)$ at frequency bin $x$.

Another way to conceptualise this relationship is to construe the derivative as describing the slope of the tangent line to the PSD function $g(x)$ at any given frequency bin $x$.
From this perspective, it becomes clear that the first derivative will be zero (i.e. the slope of the tangent will be horizontal) at any point in the function corresponding to a peak or trough.
In the case of the former, the derivative will change from a positive value (as the function ascends towards its peak) to a negative value (once the function begins to descend) as the tangent traverses the local maximum.
As such, positive to negative sign changes (i.e. downward going zero crossings) within the first derivative offer a convenient index of local maxima.
Conversely, sign changes in the opposite direction (i.e. upward going zero crossings) can likewise be used to identify local minima.
^[A lack of sign change -- e.g., a positive derivative going to zero and then becoming strictly positive again -- corresponds to a plateau.]

Differentiating the second-derivative of the PSD is achieved via differentiation of its first derivative:

$$ g''(x) = \lim_{ \Delta{x} \rightarrow{} 0}  \frac {\Delta g'(x)} {\Delta x} , $$

where $g''(x)$ is the derivative of the first derivative $g'(x)$ at frequency bin $x$.
In other words, the second derivative is simply the rate of change of the first derivative of some function $g(x)$.
Second derivatives are useful for evaluating whether the curvature of a function is concave up (i.e. convex) or concave down at any given value of $x$.
The transition of a curve's direction between concave up and concave down is characterised by an inflection point, which registers a second derivative value of zero.
Consequently, inflection points can be identified by applying the same zero-crossing procedure described for locating optima within the first derivative.

# Quantifying peak quality
Having defined both the height and width of the putative alpha peak by means of the first- and second-derivative test, relative peak quality is quantified as

$$ Q = \frac{\int_{i_1}^{i_2} \text{PSD}(f)\; df } { i_2 - i_1 } , $$

where $Q$ is the scaled average power within the peak interval $[i_1,i_2]$.
^[Notice that the interval bounded by $[i_1,i_2]$ is distinct from that bounded by $[f_1,f_2]$, the estimated span of the individualised alpha-band window.
The former yields a narrower frequency range than the latter, and does not take into account secondary peaks within the alpha band.]
(In a very strict sense, $Q$ is the mean value of the power spectral density function on the peak interval as given by the Mean Value Theorem.)
Note that the inclusion of the denominator ensures that spectral width is taken into account when calculating $Q$.  
Given equal values of $\int_{i_1}^{i_2}\text{PSD}(f)$, the denominator adjusts the integrand such that narrower, sharper peaks are assigned a larger $Q$ value than their broader, flatter counterparts.
This formulation thus penalises 'less peaky' components by assigning a heavier weighting to those estimates containing evidence of a relatively more dominant spectral peak.

# Analysis procedure
For each channel subjected to analysis, the PSD is estimated, and extraneous frequency bins beyond the span of the filter passband excluded.
Spectral data are then normalised by dividing each power estimate by the mean power of the truncated spectrum.
A log-transformed version of the PSD is derived in order to fit the regression model that will be used to estimate an upper bound on background spectral activity (i.e. $minP$).
The SGF is applied to the (nontransformed) PSD to estimate its zeroth (i.e. smoothed function), first, and second derivatives.

Initially, the first derivative is searched for downward going zero crossings within the frequency domain defined by $W_\alpha$.
If no zero crossings are identified, or if candidate zero crossings fail to exceed the corresponding  power predicted by the regression fit by more than 1 standard deviation of the estimated prediction error, the channel is excluded from further PAF-based analysis (see supplementary materials for examples of this $minP$ criterion).
If more than one peak satisfies $minP$ within $W_\alpha$, these candidates are rank ordered according to their normalised power estimates, and the magnitude of the difference between the two largest peaks compared.
The primary peak must exceed the height of its closest competitor by more than the proportion defined by $pDiff$ in order to qualify as the channel PAF estimate.
If the primary peak satisfies this condition, the second derivative is examined to determine the location of its associated inflection points, and the $Q$ value subsequently computed.
If not, the channel is excluded from PAF analysis.

CoG calculation follows the standard procedure described by Klimesch and colleagues [-@klimesch1990; -@klimesch1993], with the exception that the bounds of each channel's alpha interval are detected automatically.
The analysis routine derives these bounds by taking the left- and right-most peaks within $W_\alpha$ (i.e. those peaks in the lowest and highest frequency bins, respectively; these may coincide with the PAF), and searching the first derivative for evidence of the nearest local minimum (1) prior to the left-most peak ($f_1$), and (2) following the right-most peak ($f_2$).
^[This contingency allows for the individualised alpha-band window, and thus the CoG, to be estimated even in cases where there is no clear PAF; e.g., in the presence of split-peaks.]
Since some spectra show a relatively shallow roll-off as the edges of the alpha peak diminish, and thus do not culminate in a local minimum several frequency bins away from the main body of the component structure, we relaxed the requirement for an upward going zero crossing (i.e. evidence of a local minimum) such that the transition into a prolonged shallow gradient is taken as sufficient evidence of the individual alpha bounds $f_1$ or $f_2$.
This criterion was formalised as

$$f_1 = \text{arg} \max_{f < \text{PAF}} |PSD'(f)| < 1 , $$

$$f_2 = \text{arg} \min_{f > \text{PAF}} |PSD'(f)| < 1 . $$

$f_1$ and $f_2$ estimates from each eligible channel are averaged to yield the individualised alpha-band window.
This window defines the index of summation (i.e. frequency band coverage) used to calculate the CoG across all available channels.

If a sufficient number of channels (as stipulated by $cMin$) furnish PAF and individualised alpha window estimates, channel PAF and CoG estimates are averaged to generate IAF summary statistics.
Mean PAF (PAF$_M$) is a weighted average that takes into account the $Q$ values associated with each peak estimate:

$$ \text{PAF}_M = \frac{\sum\limits_{c=1}^C \text{PAF}_c \cdot \lambda_c} {\sum\limits_{c=1}^C \lambda_c}, $$

where $c$ identifies the channel drawn from the set of all available channels $C$, and $\lambda_c$ is the channel weighting derived by dividing $Q_c$ by the maximum $Q_c$ in $C$ (such that $\sum \lambda_c = 1$).
In contrast to PAF$_M$, all CoG channel estimates contribute equally to the calculation of mean CoG (CoG$_M$).
If there are an insufficient number of channel estimates to satisfy $cMin$, no PAF$_M$ or CoG$_M$ estimates are returned (in some cases, $cMin$ will be satisfied for CoG$_M$, but not PAF$_M$, on account of the latter's more stringent criteria).

Since pre/post-experiment recordings may not be equivalent in terms of quantity of information rendered [e.g., @samaha2015], grand averaged IAF estimates (IAF$_{GA}$) are calculated using a weighted mean which takes into account the proportion of channels that contributed to each constituent summary statistic:

$$ \text{IAF}_{GA} = \frac{ \text{IAF}_1 \beta_1 + \text{IAF}_2 \beta_2 } {\beta_1 + \beta_2} , $$

where either PAF or CoG are substituted in place of IAF, $\beta$ constitutes the weighting afforded to the channel-wise mean estimates derived from each recording, and subscript indices indicate the identity of the EEG recording.
For PAF estimates, $\beta$ is the number of channels used to estimate PAF$_M$ divided by total number of channels included in the analysis.
For CoG estimates, $\beta$ is the number of channels used to estimate the mean individualised alpha-band window divided by total number of channels included in the analysis.

\textbf{References}
