---
title: Towards a reliable, automated method of individual alpha frequency (IAF) quantification

author:
  - Andrew W. Corcoran$^{a,b}$
  - Phillip M. Alday$^b$
  - Matthias Schlesewsky$^b$
  - name: Ina Bornkessel-Schlesewsky$^b$
    affiliation: $^a$Cognition and Philosophy Lab, School of Philosophical, Historical and International Studies, Monash University, Clayton, Victoria 3800, Australia $^b$Cognitive Neuroscience Lab, School of Psychology, Social Work and Social Policy, University of South Australia, Magill 5072, Australia.

abstract: 
  Great new IAF technique described here.

keywords:   
  Individual alpha frequency, peak frequency, centre of gravity, alpha rhythm, posterior dominant rhythm, Savitzky-Golay filter

output: 
  html_document:
    toc: true
    toc_float: true

bibliography: 
  libraryAC.bib

csl: 
  apa-old-doi-prefix.csl
---
## 1 Introduction
Oscillatory activity is an inherent property of neurons and neuronal assemblies, and the timing of oscillatory dynamics is thought to encode information [e.g. @buzsaki2004;@fries2005;@vanrullen2016].
Neuronal oscillations reflect fluctuations between states of high and low receptivity, such that communication between individual neurons and broader neuronal populations is optimised via the establishment of oscillatory coherence [@fries2005;@fries2015].
Complex cognitive tasks typically require coordination between distant brain regions and systems, thus requiring effective connectivity to be established within task-relevant neural networks at relatively short timescales [@fries2005;@palva2011].
Task-irrelevant and potentially interfering connections must concomitantly be inhibited, i.e. task-relevant neural networks are gated by inhibition [@jensen2010].
The alpha rhythm of the human EEG is thought be the primary carrier of this inhibitory function [@klimesch2007;@klimesch2012;@jensen2010;@jensen2012;@sadaghiani2016], with alpha synchronisation in task-irrelevant regions reflecting inhibition, and alpha desynchronisation in task-relevant regions reflecting release from inhibition [@pfurtscheller2003].
This account is gaining increasing acceptance over alternative accounts of the alpha rhythm such as the proposal that it reflects cognitive idling [@adrian1934;@pfurtscheller1996].

While the importance of the alpha rhythm for cognitive processing has been recognised since Hans Berger’s seminal work on the human EEG in the early 20th century [@berger1929; cf. @adrian1934], a more recent line of research has focused on the importance of inter-individual variability in resting alpha activity for cognitive processing [cf. @klimesch1999, for a review].
According to this body of literature, the frequency at which alpha-generating neural circuits predominantly oscillate (i.e. the individual alpha frequency; IAF) while one relaxes in a state of alert wakefulness predicts performance across a variety of perceptual [e.g., @cecere2015;@samaha2015] and cognitive [e.g., @bornkessel2004;@clark2004] tasks.
The IAF, which varies between approximately 9.5 and 11.5 Hz in healthy young adults [@klimesch1999], is a trait-like characteristic of the human EEG [@grandy2013a], which shows high heritability [@lykken1974;@malone2014;@posthuma2001;@smit2006] and test-retest reliability [@gasser1985;@kondacs1999;@naepflin2007], while remaining stable across cognitive training interventions [@grandy2013a].
Individuals with a low IAF process information more slowly [@klimesch1996b;@surwillo1961;@surwillo1963], possibly due to decreased efficiency of thalamo-cortical feedback loops [@klimesch1997;@steriade1990].
They also show a reduced performance on memory tasks [@klimesch1999] and general intelligence measures [*g*; @grandy2013] in comparison to their high-IAF counterparts. 
IAF decreases with age from young adulthood onwards [@chiang2011;@klimesch1999;@kopruner1984;@obrist1979], and the age-related slowing of the alpha rhythm thus accompanies the well-known decline of many cognitive abilities in older adulthood [e.g. @hedden2004;@salthouse2011].
Taken together, this evidence suggests that IAF constitutes a promising neurophysiological marker of certain fundamental properties of central nervous system functioning [@grandy2013;@grandy2013a].

<!-- IBS: not sure we need following paragraph -->IAF modulates early brain responses that accompany basic aspects of visual perception [@klimesch2004;@koch2008], while individual alpha power and alpha timing (entrainment) affect well-known attentional phenomena such as the attentional blink [@maclean2012;@zauner2012]. 
Moreover, IAF-informed brain stimulation has been shown to improve cognitive performance [@klimesch2006] and to shift the boundaries of integration windows for cross-modal perception [@cecere2015].
Accordingly, individual differences in alpha activity appear to influence the efficiency of attentional filtering [@klimesch2012].
This mechanism is especially important for increasing the signal-to-noise ratio when attentional demands increase, for example in cognitive processing scenarios requiring top-down control or prediction [@klimesch2011].
Adjusting frequency band analysis in relation to IAF has also been argued to enhance sensitivity to band-specific oscillatory dynamics [@citation], and has proved useful for dissociating distinct functional properties of alpha-band sub-regions [@citation].
Indeed, recent evidence [@van_albada2013] supporting the idea that the canonical frequency bands display an approximately harmonic relationship to the dominant alpha rhythm [@klimesch2012] suggests that IAF may offer a valuable tool for sharpening the precision of frequency domain analysis across the entire spectrum.

In spite of its promise as a marker of apparently enduring, trait-like individual differences in cognitive functioning, and its potential utility for tuning frequency band analysis to the individual properties of electrophysiological data, no consensus currently exists as to the optimal method for quantifying IAF.
This paper thus sets out to develop a rigorous, automated strategy for estimating two of the most widely reported indices in the IAF literature; namely, peak alpha frequency and alpha frequency centre of gravity.
We begin by surveying various ways in which these measures have been operationalised and implemented in previous research, and highlight some of the problematic aspects of these methods.

### 1.1 Peak alpha frequency
The classical method of estimating IAF relies on delineating the peak alpha frequency (PAF); a singular, prominent peak within the alpha-band frequency range [8-13 Hz; @noachtar2004] of the power spectral density (PSD) plot ([Fig_pafs](#pafs)).
This expression can be formalised in terms of the local (i.e. relative) maximum within the alpha band:

$$ PAF = \text{arg max}_{\substack{x\in A \subseteq PSD }} f(x) \iff \text{arg max}_{\substack{x\in A \subseteq PSD }} f(x) = \text{singleton} \land \text{max}_{\substack{x\in A \subseteq PSD }} f(x) \geq \phi , $$

$$ \text{arg max}_{\substack{x\in A \subseteq PSD }} f(x) := \{ x \mid x  \in A \land \forall y \in A : f(y) \leq f(x) \} , $$

where $\text{arg max}$ returns the frequency bin (or subset of bins) $x$ containing the maximal power value $\text{max } f(x)$ registered within that subset of frequency bins $A$ constituting the alpha band, and where $PSD$ denotes the complete set of frequency bins resolved by the spectral analysis. 
Note that, for the output of $\text{arg max}$ to qualify as an estimate of PAF, it must return a single frequency bin $x_k$ containing a power value $\geq \phi$, where $\phi$ defines the minimum threshold value differentiating a substantive spectral peak from background noise. 
The definition of both $A$ and $\phi$ pose non-trivial problems (more on which shortly).

![*Fig_pafs.* Power spectral density (PSD) plots displaying frequency component distribution of averaged signal variance across a 2 min eyes-closed resting-state EEG recording (POz). Light grey column indicates the standard alpha band interval, which constitutes the search window for the peak frequency. *Left panel*: Linear-scaled PSD ranging from 1 to 25 Hz. Strong alpha band activity is clearly evidenced by the sharp component spanning ~7.5 to 12.5 Hz, and peaking at ~9.75 Hz. *Central panel*: Alternative depiction of the PSD displayed in the left panel. Here, ordinate data have been log-transformed into decibels. Decibel-scaling accentuates the relatively minor peak detected in the beta range of the spectrum (this activity approximates the first harmonic of the dominant alpha rhythm). *Right panel*: Log-log plot of spectral density estimates across all frequency bins resolved within the range of 1 to 100 Hz (frequencies and power estimates both log~10~-transformed). The alpha peak represents a marked deviation from the $1/f$ inverse power-law (indicated by the broken line) characteristically approximated by log-transformed EEG power spectra.](figs/pafs.png?raw=true){#pafs}

PAF estimates are typically extracted from parieto-occipital EEG channels while the participant relaxes with their eyes closed.
This strategy exploits the classic observation that alpha oscillations dominate EEG recorded over centro-posterior scalp regions when visual sensory input is suppressed [@barry2007;@sadaghiani2016].
Although PAF can often be rapidly ascertained upon visual inspection of the PSD function ([Fig_pafs](#pafs)), this approach to IAF extraction is inefficient and potentially impractical when dealing with large datasets [@chiang2008;@goljahani2012]. 
Moreover, it is well documented that a sizeable proportion of individuals fail to manifest an unambiguous PAF, on account of their being either two (or more) peaks within the alpha band [so-called ‘split-peaks’; @chiang2011], or a general absence of prominant alpha-band activity [e.g., @anokhin1996] (see [Fig_bad_pafs](#bad_pafs)).
Under the former circumstances, the adjudicator<!-- AC: this term has proved controversial. i deliberately wanted to underscore the tacit subjectivity involved in judging whether a given PSD manifests a (singular) peak. i'm open to alteratives but would like to retain this connotation, if possible. --> must decide whether a single, primary peak can be justifiably discerned amidst competing peak candidates; under the latter, they must determine whether the signal is too noisy to derive reliable inferences pertaining to the IAF. 
Relating back to the PAF formulation outlined above, these scenarios can be construed as the problem of evaluating whether the $\phi$ criterion has been satisfied (or alternatively, of determining the correct threshold level where this criterion ought to be located).
Such cases may be prone to biased or inconsistent assessment [@chiang2008], pose significant challenges to replicability, and can result in the exclusion of a substantial subset of participants from IAF-related analyses [see for e.g., @bornkessel-schlesewsky2015].

![*Fig_bad_pafs.* Examples of problematic cases for peak alpha frequency identification. *Left panel*: Multiple peaks resolved within the alpha band search window (grey column). *Central panel*: Split- (or bimodal) alpha peak. *Right panel*: No discernable alpha peak. Note differences in ordinate scaling.](figs/bad_pafs.png){#bad_pafs}

Automated analysis routines offer one way of ensuring that PAF criteria are explicitly stipulated and consistently applied within and across datasets.
However, while simple automated peak detection methods are relatively straightforward to implement in many popular EEG processing software packages, such strategies are prone to various sources of error. 
For instance, searching for the maximal power estimate $\text{max } f(x)$ within the alpha-band search window $A$ can result in the arbitrary assignment of PAF at $f_1$ (i.e. the lower bound of $A$) in the absence of a clear peak.
That is, in the case where the PSD function declines approximately monotonically [e.g., conforms to the $1/f$-like power distribution without showing the characteristic deviation about 10 Hz; cf. @pritchard1992], the highest power value within $A$ will be the supremum encountered at the first frequency bin encompassed by the alpha-band interval.
One solution to this problem is to stipulate that $\text {arg max}_{\substack{x\in A \subseteq PSD }} f(x)$ qualifies as a viable PAF estimate $iff$ the power estimate of frequency bin $x_k$ exceeds that of its neighbouring frequency bins $x_{k-1}$ and $x_{k+1}$. 
While this approach ensures that the selected frequency component constitutes a local maximum rather than a supremum, it is still vulnerable to two of the problems identified above. 
First, it fails to distinguish spectra featuring singular, dominant peaks from those possessing prominent secondary peaks (i.e. where visual inspection of the plot would suggest two or more alpha peaks that differ in height by some – potentially trivial – magnitude; e.g., right and central panels, [Fig_bad_pafs](#bad_pafs)). 
Second, it fails to differentiate maximum power values at the apex of a genuine spectral peak from those at the apex of noisy fluctuations within the PSD function (i.e. where visual inspection of the plot would suggest the absence of any substantive alpha peak; e.g., right panel, [Fig_bad_pafs](#bad_pafs)). 
Automated routines of this sort might therefore render rapid and consistent estimates of PAF, but are likely to be too liberal in what they accept as an alpha peak.

### 1.2 Alpha centre of gravity and individualised frequency bands
Klimesch and colleagues [@klimesch1993; @klimesch1997] proposed using alpha mean or centre of gravity (CoG) frequency, an IAF estimator originally formulated by Klimesch, Schimke, Ladurner, and Pfurtscheller [-@klimesch1990], in order to circumvent some of the difficulties posed by the absence of a clear spectral peak. 
The CoG furnishes a weighted sum of spectral estimates divided by total power within the selected alpha frequency window $A$, as defined by the bounds $f_1$ and $f_2$:

$$ CoG = \frac { \sum\limits_{f_1}^{f_2} f(x) \ x } { \sum\limits_{f_1}^{f_2} f(x) }, $$

where $f(x)$ is the PSD estimate (ordinate) for frequency bin $x$.

Since the CoG is sensitive to the shape of the power distribution within the selected alpha band window, and the precise bandwidth of alpha-rhythm activity varies across individuals, Klimesch and colleagues [-@klimesch1990, see also @klimesch1997] discouraged calculating the CoG according to a fixed index of summation corresponding to some standard, a priori-defined alpha bandwidth (e.g., $f_1$ = 8 Hz, $f_2$ = 13 Hz). 
Rather, they recommended computing CoG on the basis of bespoke frequency windows that capture the entire range of the individual’s alpha-band activity. 
To this end, Klimesch and colleagues [-@klimesch1990] proposed the following procedure for estimating the IAF bandwidth: 
First, PSD plots are extracted from all EEG channels for each participant and examined for evidence of a clear alpha peak. 
Second, $f_1$ and $f_2$ are assigned to those frequency bins where the ascending and descending edges of the peak are deemed to start and end, respectively [[Fig_f1f2](#f1f2); cf. Table 1, @klimesch1997]. 
Finally, these channel-wise $f_1$ and $f_2$ values are averaged to render the bounds of the frequency interval that will be used to calculate the participant’s CoG. 
Notice that, even though EEG channels that fail to manifest distinctive peaks do not contribute to the definition of the IAF bandwidth, the CoG is computed on the basis of spectral data compiled from all frequency bins bounded by $f_1$ and $f_2$ across all available channels.

![*Fig_f1f2.* Fixed vs. individually adapted alpha frequncy band windows. *Left panel*: Individual bandwidth (indicated by shaded area under PSD) closely approximates the standard alpha band range indicated by the grey column. *Central panel*: Fixed bandwidths may fail to fully and selectively capture empirical alpha peak distributions, thus giving rise to biased estimates of alpha centre of gravity. *Right panel*: Spectral density from central panel reproduced with superposed eyes-open resting-state PSD (red function). Here, the transition frequency (TF) at which the ascending edge of the eyes-closed alpha peak intersects with (and surpasses) the corresponding eyes-open power estimates (indicated by the broken black line) closely approximates the lower bound of the alpha interval $f_1$ located via qualitative analysis of the plot. In this case, the empirical TF of ~6.25 Hz deviates from standard definitions of the theta/alpha boundary, thus locating the beginning of the alpha component within what would canonically be considered the theta band (indicated by the pink column). Power estimates log-scaled to aid visual identification of $f_1$ and $f_2$ [@klimesch1997].](figs/f1f2.png){#f1f2}

Klimesch [-@klimesch1999] later proposed an alternative method for defining individualised alpha-band windows that relies on a somewhat less subjective interpretation of the PSD. 
This technique exploits the typically observed anticorrelation between theta band [4–8 Hz; @noachtar2004] and alpha oscillatory dynamics in response to task demands or cognitive load [@doppelmayr1998a; @klimesch1996a; cf. @rugg1982]. 
Klimesch and colleagues [-@klimesch1996a] first used this approach to delineate the ‘transition frequency’ (TF) where the alpha-band activity that dominates the relaxed individual’s EEG gives way to theta oscillations induced by stimulus processing (see [Fig_f1f2](#f1f2), right panel). 
Adopting an event-related desynchronisation [ERD; @pfurtscheller1977; @pfurtscheller1999] paradigm, spectra from a pre-stimulus ‘reference’ interval and a peri-stimulus ‘test’ interval were averaged across trials and superimposed, and the TF designated as the frequency at which the reference and test interval PSD functions intersected (or where the difference between them was minimised, in cases where they failed to intersect). 
Klimesch [-@klimesch1999] proposed generalising this procedure to resting-state EEG recordings, where an analogous shift from prominent alpha- to theta-band activity is classically evoked by the visual stimulation experienced upon opening the eyes [this phenomenon, known variously as alpha blocking, desynchronisation, suppression, or attenuation, was first documented by @berger1929]. 
This method thus renders a systematic means of estimating the lower alpha frequency bound $f_1$. 

In the absence of any analogous means of inferring the upper bound of the alpha-band window, Klimesch [-@klimesch1999] recommended determining $f_2$ on the basis of $f_1$ and the IAF.
One suggestion was to set $f_2$ rather pragmatically at 1 or 2 Hz above the IAF ([Fig_findf2](#findf2), left column).
Alternatively, $f_2$ could be derived by subtracting the difference between IAF and $f_1$ (i.e. the lower alpha band) from some presupposed alpha bandwidth, and treating the remainder as the span of the upper alpha band.
For instance, if the alpha interval were assumed to span 5 Hz, and the PSD manifested a TF of 7 Hz and an IAF of 10.5 Hz, then $f_2 = 5 - (10.5 - 7) + 10.5 = 12$ Hz ([Fig_findf2](#findf2), central column). 
Although these approaches do indeed adjust the value of $f_2$ in relation to the IAF, they are insensitive to interindividual differences in alpha bandwidth [@doppelmayr1998; @goljahani2012]. 
A more promising solution that attempts to capture such variance is given by calculating $f_2$ as a proportion of the difference between IAF and $f_1$ [@doppelmayr1998;@klimesch1999]; e.g., $f_2 = ( (IAF - f_1) / 2 ) + IAF$ ([Fig_findf2](#findf2), right column). 
While this heuristic offers a more nuanced approach towards individually adapted estimates of $f_2$, it clearly depends on the assumption that the proportional difference between the IAF and the alpha bandwidth is consistent across individuals (such that the lower alpha band spans approximately double the frequency range of the upper alpha band).
Furthermore, all of these methods presuppose that IAF has already been estimated via the PAF (since CoG cannot be calculated prior to the definition of $f_1$ and $f_2$). 
This is obviously problematic given that one of the chief advantages of the CoG is its supposed capacity to deliver IAF estimates when the PAF is unavailable.

![*Fig_findf2.* Application of three procedures for deducing $f_2$ (right-most broken line) to the two channel spectra presented in [Fig_f1f2](#f1f2). In each case, $f_1$ (left-most broken line) is defined according to the transition frequency (TF), and IAF (solid line) is defined according to the peak frequency. *Top row*: Defining $f_2$ as $IAF + 2$ Hz (left panel) and $IAF + ((IAF - TF)/2)$ (right panel) render similar results, closely approximating the boundary located by visual inspection. However, defining $f_2$ as the residual of a fixed bandwidth (in this case, 6 Hz) following subtraction of the lower alpha band (blue shading) fails, collapsing $f_2$ into the IAF on account of the broadness of the lower alpha band. *Bottom row*: All three methods return similar estimates of $f_2$ when the lower alpha region conforms to the 3.5-4 Hz bandwidth assumed by Klimesch [-@klimesch1999]. Each of these attempts to calculate the span of the upper alpha band region (pink shading) appear to be suboptimally narrow.](figs/findf2.png){#findf2}

### 1.3 Peak attenuation and channel reactivity based (CRB) methods
We turn now to two interesting extensions of the TF approach that depend on the alpha blocking or desynchronisation phenomenon mentioned above. 
The first method, which we refer to as peak attenuation, was described by Posthuma, Neale, Boomsma, and de Geus [-@posthuma2001].

Similar in principle to Klimesch's [-@klimesch1999] delineation of the TF, this technique simply subtracts eyes-oped PSD estimates from the corresponding ordinates of the eyes-closed PSD within some a priori-defined target region (typically, this interval is extended beyond the canonical alpha bandwidth; e.g., 7–14 Hz).
The resulting peak, which constitutes the maximal power difference between the two spectra, is taken as the PAF ([Fig_blocking](#blocking), top row). 
To avoid trivial peaks arising from small fluctuations in the PSD, Posthuma and colleagues [-@posthuma2001] excluded spectra that featured consistent low power (< 1.5 $\mu$V/Hz) within the alpha-band.
It isn’t clear, however, whether this criterion (which is analagous to the $\phi$ parameter in the PAF equation defined above) was also applied in cases where the maximal difference between spectra was similarly small (i.e. where both spectra present a substantive, approximately overlapping alpha-band component).
Posthuma and colleagues [-@posthuma2001] also remarked that peak attenuation estimates deviated from eyes-closed resting-state PAFs in 21% of cases. 
This lack of convergence, coupled with the observation that peak attenuation may lead to distorted representations of the PAF when the assumption of significant eyes-open alpha desynchronisation fails to obtain [see [Fig_blocking](#blocking), central and bottom row; cf. @gaal2010; @kreitman1965], suggests peak attenuation may constitute a suboptimal IAF estimator (at least insofar as IAF is conceived in terms of the PAF).

![*Fig_blocking.* Illustration of the peak attenuation method, in which PSD estimates from a period of eyes-open EEG (EO) are subtracted from a corresponding eyes-closed recording (EC). The peak frequency of the resulting difference function (within the extended alpha band 7-14 Hz, grey column in difference spectra) is taken as the IAF estimate. *Top row*: Alpha-band activity registered during the EO recording (magenta shading) is relatively low compared to that of the EC spectrum (blue shading). In this case, the peak attenuation method returns the PAF of the EC spectrum as the IAF estimate. *Central row*: When alpha-band power is similar across EO and EC recordings, peak attenuation may be vulnerable to bias. In this case, the primary peak at ~9.75 Hz is present in both spectra, thus cancelling one another out in the difference spectrum. Peak attenuation returns 8.5 Hz (the frequency of the secondary peak) as the IAF. *Bottom row*: Data from another participant in which partial overlap of EC and EO spectral components shifts the estimated alpha peak down by ~1 Hz.](figs/blocking.png){#blocking}

Although Posthuma and colleagues [-@posthuma2001] did not attempt to locate the bounds of the individual alpha bandwidth on the basis of peak attenuation, the logic motivating this technique could be applied to infer $f_1$ and $f_2$ in much the same way as the TF (i.e. by taking those frequency bins either side of the difference peak where the difference between corresponding spectral estimates is *minimised*).
This is precisely the approach that Goljahani and colleagues formalised in their channel reactivity based (CRB) method [@goljahani2012; @goljahani2014]. 
This technique, which is conceptually reminiscent of Klimesch and colleagues’ [-@klimesch1996a] attempt to characterise phasic shifts in band power, quantifies the difference between reference and test PSDs in terms of the alpha responsiveness (or desynchronisation) region; i.e. the area between the PSD functions spanning frequency bins $f_1...f_2$ (where the spectra intersect, or the residual difference between corresponding spectral estimates is minimised). 
IAF is estimated by computing the CoG for the reference interval, taking the frequency bounds delimiting the responsiveness region as the index of summation.

The CRB method offers an elegant solution to the problem of finding individualised frequency bands for CoG estimation. 
Not only does it improve on previous attempts to characterise the bounds of the alpha range by means of a rapid, automated, and empirically driven appraisal of reference vs. test interval PSD functions, it does so without any assumptions about (or dependence on) the peak-like quality of the spectra. 
As such, CRB promises to maximise the potential utility of the CoG as an estimator of IAF that can be reliably computed irrespective of the presence or absence of a single, distinct peak component. 
We note however that, much like the peak attenuation technique, CRB estimates may be prone to bias in cases where the ERD task elicits asymmetrical patterns of alpha-band desynchronisation [@klimesch1996b; see also @klimesch1997, for evidence of the dissociation between lower and upper alpha-band reactivity]. 
This is because the partial overlap of reference and test spectra that will ensue should the test stimulus fail to evoke comprehensive alpha desynchronisation will restrict the responsiveness region to a limited segment of the alpha band, thus precipitating a narrower (and potentially shifted) index of summation. 
Indeed, this phenomenon might go some way to explaining the substantial discrepancies observed between some of the CRB CoGs reported by Goljahani and colleagues (some of which were extreme by conventional IAF standards; e.g., 14.9 Hz) and the corresponding PAFs derived from the same set of channel data [see Figure 6(d), @goljahani2012].

A related concern deriving from the CRB method’s reliance on phasic changes in rhythmic activity is the possibility that within-subject estimates of IAF might vary depending on the specific processing mode evoked by the event (e.g., target discrimination vs. memory retrieval; visual vs. auditory modality), and the relative timing of the reference/test intervals subjected to spectral analysis. <!-- PA: How do these concerns interact with the idea that IAF is trait-like? AC response: Several authors posit that there are a variety of alpha generators, the cumulative activity of which accounts for the characteristic resting peak (or peaks), and thus shifts in peak frequency mark shifts in the relative activity of these populations. The argument here then is that, while resting eyes closed IAF is stable/enduring and indicative of at least some cognitive capacities (presumably by tapping into specific underlying neuronal processes), measures of alpha activity during cognitive tasks (or the difference in alpha activity between different states of rest/anticipation/action) may not offer an analagous index of IAF as construed by research involving EC resting state data. I see no necessary contradiction between the notion that resting state IAF is stable/'trait-like' etc, and the observation that dominant alpha band activity dynamically shifts during 'non-restful' states (presumably as a consequency of the shifting division of labour across alpha generators). indeed, since you asked me this question, i've re-read a few articles and noticed a number of references to this idea of IAF as possessing both trait- and state-like qualities-->
ERD studies have revealed that both the qualitative profile and temporal course of alpha- and theta-band desynchronisation are contingent upon the particular nature of the task used to induce oscillatory power shifts [@klimesch2006; @klimesch2007]. 
If different paradigms do precipitate distinct patterns of ERD during the selected test interval [or indeed, *enhance* rather than attenuate the alpha rhythm; e.g., @kreitman1965; @rihs2007], then the ensuing responsiveness regions used to define the coverage of the CoG estimate will span non-identical frequency bands [cf. @haegens2014, for evidence of analogous intraindividual shifts in PAF as a function of varying task conditions]. 
While this property of the CRB method need not be a problem for ERD-type applications (indeed, sensitivity to such selective changes in band power might prove theoretically interesting and productive in this context), it renders the approach less suited to the task of estimating the IAF as a marker of stable, trait-like differences in information processing capacities.

### 1.4 Automated curve-fitting approaches to alpha rhythm quantification
Finally, we turn briefly to a promising line of research that attempts to quantify the spectral features of EEG data, and in particular spectral peaks, by means of statistical curve-fitting techniques. 
Chiang and colleagues [-@chiang2008] developed an algorithm (with corresponding implementation in C) that parameterises alpha-band peaks via a two-step procedure: Peaks are first identified and parameterised via the fitting of a Gaussian function, before being fine-tuned in relation to the spread of fitted estimates across multiple electrode sites. 
Chiang and colleagues [-@chiang2011] and van Albada and Robinson [-@van_albada2013] later demonstrated the potential utility of such automated routines by applying this general technique to datasets comprising 1498 and 1424 individuals, respectively. 
Another automated spectral analysis technique that similarly instantiates an iterative curve-fitting and clustering procedure has been proposed by Lodder and van Putten [-@lodder2011;-@lodder2013].
This technique has likewise been applied to a large dataset comprising resting-state EEG recordings from 1215 individuals [@lodder2011], where it produced estimates that were highly correlated with previously documented visual estimates of PAF (on average, visual and automated estimates differed by 0.52 Hz).
Lodder and van Putten reported that higher degrees of estimate accuracy (i.e. agreement with visually scored PAFs) could be achieved by eliminating low-certainty estimates, although at its extreme this strategy resulted in a substantial rate of attrition (41% of estimates discarded).

Given the obvious advantages of automated analysis tools, it is somewhat puzzling that these curve-fitting techniques do not yet appear to have gained widespread currency in the contemporary IAF literature [although cf. Haegens and colleagues, -@haegens2014, for a notable counterexample].
For instance, neither Goljahani and colleagues [-@goljahani2012] nor Bazanova and Vernon [-@bazanova2014] mention the development of such algorithms in their reviews of IAF methods. 
One possibility is that many researchers are simply unaware of the existence of these methods, since their application to date has been predominantly focused on spectral modeling, rather than the quantification of IAF per se. 
An alternative explanation is that researchers deem these methods too complex to be a worthwhile investment of their time (especially if quantifying IAF is only an intermediary step within a broader analysis framework, rather than the main topic of inquiry). 
This attitude might be reinforced by the additional burden involved in obtaining and implementing an algorithm that may have been written in an unfamiliar programming language, and which poses nontrivial challenges with respect to integration within existing analysis pipelines. 
We suggest then that one of the critical steps towards achieving a more widespread adoption of automated IAF estimation routines is to make these tools as openly available as possible, in formats that are easy to assimilate within popular methods of EEG data analysis.

### 1.5 Aims of the present study
In sum, common methodological approaches to IAF estimation are either (1) time-consuming and vulnerable to inconsistencies arising from qualitative interpretation, (2) at risk of producing spurious or biased estimates under certain plausible spectral conditions, (3) conflate trait-like alpha properties with variable phasic effects, or (4) show some combination of the above. 
More recent innovations designed to address these problems via the application of sophisticated curve-fitting algorithms have so far found limited uptake within the broader IAF literature, perhaps on account of practical barriers pertaining to software access and implementation.
Consequently, we sought to articulate an automated method of alpha-band quantification that provides fast, reliable, and easily replicated estimates of resting-state IAF in two major programming languages: MATLAB^&reg;^ (The MathWorks, Inc., Natick, MA, USA) and Python&trade;.
This goal is consistent with recent proposals to make the analysis of cognitive electrophysiological data as open, transparent, and amenable to replication as possible [@cohen2017].

Our approach aims to emulate Klimesch and colleagues’ [-@klimesch1990] original attempt to characterise individual profiles of resting-state oscillatory activity across the entirety of the alpha band by means of a relatively simple, non-parametric curve-fitting algorithm. 
Our strategy for accomplishing this task runs as follows: <!-- AC: am I correct in thinking it's non-parametric because it makes no assumptions about the underlying shape of the function (only that it can be fitted by polynomials) ? -->
First, we extract PSD estimates from preprocessed and fast Fourier transformed EEG signals. 
Second, we apply a least-squares curve-fitting procedure (i.e. Savitzky-Golay filtering) to accomplish the dual task of smoothing the PSD function and estimating its first- and second-order derivatives.
Third, these derivative functions are analysed to evaluate the quality of evidence that a distinct spectral peak (i.e. PAF) exists within the alpha frequency region.
Finally, the first derivative of the PSD is reanalysed to delineate upper and lower bounds of the individualised alpha interval, which are taken as the index of summation required to calculate the CoG. 
The effectiveness of this routine will first be demonstrated using empirical (i.e. non-simulated) EEG data. 
We will then turn to simulated data in order to assess how well our proposed technique performs under conditions of varied spectral composition and signal-to-noise ratio (SNR).

## 2 Method

### 2.1 Overview of methodological approach: Differentiation as a means of spectral peak quantification
In the following section, we show how differential calculus can be exploited for the purposes of alpha peak parameterisation. 
First, we outline how spectral peaks (and troughs) can be localised via differentiation of the first-order derivative. 
We then address the problem of multiple (potentially trivial) zero crossings, and propose Savitzky-Golay filtering as an elegant solution to this concern. 
Finally, we turn to the second-order derivative in order to arrive at a means of evaluating the relative quality of individual channel peak estimates.

#### 2.1.1 Local extrema and first derivative zero crossings 
<!-- AC: PA suggested deleting / relegating this section to the supp mats, as will be rather elementary for a proportion of readers. At the moment I'm inclined to keep it unless the reviewers kick up a fuss. When PA asked me to consider my target audience, I think that would be me ~6 months ago (and I would have needed this section to make sense of the method). I'm open to persuasion on this point though. -->As pointed out by Grandy and colleagues [-@grandy2013; -@grandy2013a], one solution to the problem of automated peak detection is to search for downward going zero crossings in the first derivative of the PSD. 
Derivatives describe the relative rate of change in the dependent variable or function $f(x)$ given some value of the independent variable $x$. 
The first derivative of a vector of PSD estimates thus provides point estimates of the (instantaneous) rate of change in the amount of spectral power estimated for each frequency bin resolved in the analysis. 
This relationship can be formalised as follows:

$$ f’(x) = \frac{\Delta f(x)} {\Delta x} , $$

where $f’(x)$ is the first derivative of the relative change in the power estimate $f(x)$ at frequency bin $x$.

Another way to conceptualise this relationship is to construe the derivative as describing the slope of the tangent line to the PSD function $f(x)$ at any given frequency bin $x$. 
From this perspective, it becomes clear that the first derivative will be zero (i.e. the slope of the tangent will be horizontal) at any point in the function corresponding to a peak or trough. 
In the case of the former, the derivative will change from a positive value (as the function ascends towards its peak) to a negative value (once the function begins to descend) as the tangent traverses the local maximum. 
As such, positive to negative sign changes (i.e. downward going zero crossings) within the first derivative offer a convenient index of local maxima. 
Conversely, sign changes in the opposite direction (i.e. upward going zero crossings) can likewise be used to identify local minima.

#### 2.1.2 Savitzky-Golay smoothing and differentiation
Although Grandy and colleagues [-@grandy2013; -@grandy2013a] correctly observe that searching for downward going zero crossings avoids the problem of arbitrary boundary effects in the absence of any clear alpha peak, they fail to articulate a systematic method for differentiating substantive peaks from trivial fluctuations in the PSD.
We suggest that the situation in which spectral analysis is degraded by signal noise can be substantially improved via the application of a smoothing procedure. 
The idea here is to attenuate noisy fluctuations about the true alpha peak such that the vast majority of zero crossings deriving from trivial variations are eliminated from the signal. 
However, since standard filtering techniques (such as the moving average) can result in marked distortions of the underlying peak structure [e.g., @press1992; @ziegler1981], the challenge is to find a smoothing operation that preserves the spectral characteristics of critical import to IAF analysis.

With this concern in mind, we turn to the Savitzky-Golay filter (SGF), a least-squares curve-fitting procedure specifically designed to aid in the detection of spectral peaks amidst noisy conditions [@savitzky1964]. 
The SGF has a number of properties that make it well suited to the task of smoothing PSD functions, not least of which being its capacity to render smoothed curves that conserve the height, width, position, area, and centre of gravity of the underlying spectral structure [see @ziegler1981].
SGFs work by centring a sampling window of length $F_w$ on a portion of the input signal and computing the least-squares fit of a specified polynomial to each $i$^th^ data point spanned by $F_w$.
The window is then shifted one point along the input signal, and the polynomial fit recalculated accordingly.
The centre value of the polynomial fit is taken as the filter output at each iteration of the sliding window calculation, and these output values are concatenated to render the smoothed estimate of the input function.
For a more detailed treatment of the SGF and its technical performance properties, the interested reader is referred to Schafer [-@schafer2011].

In addition to its smoothing capability, SGFs can also be applied to calculate the $n$^th^ order derivative of the input signal.
Indeed, commensurate with their desirable spectral smoothing characteristics, SGFs are optimal (or near optimal) digital differentiators [@luo2005]. 
The performance features thus qualify the SGF as a valuable tool for both (1) refining the precision of standard methods used to characterise the spectral profile of alpha-band rhythms, and (2) improving the reliability of first derivative zero crossing approaches to spectral peak (and trough) localisation. 
Before describing how the dual function of the SGF can be implemented for the purpose of IAF analysis, however, we turn to one final innovation involving the second derivative.

#### 2.1.3 Assessment of peak quality
Typically, resting-state EEG recordings afford data from multiple electrode channels, a selection of which may contribute to the final estimate of the IAF. 
Channels that are in close proximity to one another are expected to produce highly correlated data; hence, a set of channels concentrated on the centro-posterior region of the scalp should ideally render highly convergent estimates of spectral power. 
However, since channels may be differentially affected by various sources of signal noise (e.g., high or fluctuating levels of impedance between scalp and electrode), SNR might be degraded in analyses that treat all data sources uniformly.
We therefore propose an automated method of peak analysis that seeks to evaluate which channels provide the strongest evidence of a prominent alpha peak, so that these channels can be assigned heavier weights for cross-channel averaging. 
This procedure is intended to sharpen the precision of PAF estimation amidst cross-channel PSD variability without resorting to the exclusion of channels that are less conformant to the conceptual ideal of a high powered, narrowly defined spectral peak.

This approach relies on differentiation of the second derivative of the PSD function:

$$ f’’(x) = \frac {\Delta f'(x)} {\Delta x} , $$

where $f’’(x)$ is the derivative of the first derivative $f’(x)$ at frequency bin $x$.
In other words, the second derivative is simply the rate of change of the first derivative of some function $f(x)$. 
Second derivatives are useful for evaluating whether the curvature of a function is concave up (i.e. convex) or concave down at any given value of $x$. 
The transition of a curve’s direction between concave up and concave down is characterised by an inflection point, which registers a second derivative value of zero.

We suggest that the inflection points $i_1$ and $i_2$ on either side of $\text{max } f(x)$ offer a convenient objective standard for evaluating the relative quality of channel peak estimates. 
The basic idea here is to quantify the area under the peak in such a way that distinguishes stronger (i.e. containing a greater proportion of spectral power) and less variable (i.e. spanning fewer frequency bins) peaks from shallower, broader, or otherwise noisier components. 
Given a set of spectral data from a variety of electrode channels, those PSDs which give rise to higher quality peaks (as operationalised above) are weighted more heavily than their less prominent counterparts, and thus contribute relatively more information when calculating the cross-channel PAF estimate. 
Note that this procedure has no bearing on CoG estimation (since the CoG may be derived from spectra in which no clear evidence of a single alpha peak was detected, and thus for which no inflection point data will be available).

Having defined both the height and width of the putative alpha peak by means of the first and second derivative zero crossings, peak quality is quantified via the following formula:<!-- AC: I've used integration here, or more properly an approximation to integration rendered via MATLABs `trapz` function, but I wonder if it's better to use the sum of frequency bins approach as per the CoG formula in order to draw out the similarities with the CoG calculation ? -->

$$ Q = \frac{\int_{i_1}^{i_2} f(x) } { i_2 – i_1 } , $$
<!-- PA: From the Mean Value Theorem, this is exactly the “mean” power on that interval, which tells you a bit about how strong the peak is – peaks have higher mean values. AC response: All I can tell from wikipedia etc is that the MVT proves there is at least one point on a differentiable function that has a tangent that is parallel to the secant through the function's endpoints - I don't see how this relates to Q calculation -->
where $Q$ is the scaled quantity of normalised alpha peak power, and $f(x)$ is the estimated power at each frequency bin $x$ that falls within the index of integration $i_1... i_2$. 
Note that the inclusion of the denominator ensures that spectral width is taken into account when calculating $Q$.  
Given equal values of $\int_{i_1}^{i_2}f(x)$, the denominator adjusts the integrand such that narrower, sharper peaks are assigned a larger $Q$ value than their broader, flatter counterparts. 
This formulation thus penalises ‘less peaky’ components by assigning a heavier weighting to those estimates containing evidence of a relatively more dominant spectral peak. 
However, it is perhaps worth emphasising that this calculation only influences PAF estimation in cases where channel data produce divergent PAF estimates (i.e. channel estimate $Q$ weights have no impact on the mean PAF calculated from channels that furnish identical estimates of the peak frequency).

![*Fig_q_wts.* Power spectra from four individuals in which two channels from the same eyes-closed resting-state recording are superposed. Each spectrum is shaded within the region bounded by the inflection points either side of the  peak alpha frequency estimate. Respective $Q$ values are also presented. *Left column*: Blue channel peaks dominate their red counterparts, however cross-channel averages are unaffected due to the identity of peak estimates. *Right column*: Where peak estimates diverge across spectra, channels which manifest a greater concentration of power (holding peak width constant) will be assigned higher weightings compared to those with relatively less power. Note that secondary components are ignored on account of the inflection point delimitation of the main peak. Power values normalised within each channel according to mean spectral power.](figs/q_wts.png){#q_wts}

### 2.2 Implementation of proposed analysis techniques

#### 2.2.1 Software requirements
The afore-described approach to IAF estimation has been implemented via a set of customised functions programmed in MATLAB and Python.
The following report focusses on the MATLAB implementation of the programme, which is dependent upon the Signal Processing Toolbox&trade; and the EEGLAB toolbox [@delorme2004].
EEGLAB is necessary for data importation, since our analysis programme assumes that EEG data are structured according to EEGLAB conventions. 
Signal Processing Toolbox is required for the `pwelch` and `sgolay` functions, which are responsible for executing the PSD estimation and SGF design components of the programme, respectively. 
`sgolay` was preferred over `sgolayfilt` on the basis that it outputs the coefficients necessary for calculating higher-order derivative functions, as well as those required for the zero-order (i.e. smoothed) PSD function. 
All functions that were developed in order to conduct the following analyses are open source and can be accessed (along with sample datasets) via [GitHub](https://github.com/corcorana/restingIAF).

We chose to base our IAF estimation strategy on a Hamming-windowed implementation of Welch’s modified periodogram method of PSD estimation [@welch1967] on account of the prevalent application of this technique in the IAF literature. 
It should be noted however that our approach could quite readily be modified to enable the `pwelch` routine to be substituted with some alternative method of PSD estimation. 
We also stress that our selection of the Hamming window to taper time series data prior to FFT does not reflect any strong commitments concerning the optimal method of preparing EEG data for spectral analysis.
It stands rather as one common preference amongst several possible alternatives, each of which could reasonably have been applied. 
Indeed, while a rigorous comparison of how our IAF estimation technique performs in conjunction with various PSD estimation methods would be desirable, it is beyond the scope of this paper. 
Our focus here, rather, is to provide a convincing proof of concept for the general tenets of our smoothing-differentiation approach.

#### 2.2.2 Parameters for estimating PAF and CoG
A number of parameters must be specified in order to execute the IAF estimation programme. 
In relation to the SGF, both polynomial degree $k$ and a filter window frame width $F_w$ are required to define the least-squares minimisation operation.
$k$ must be $< F_w$, and $F_w$ must be odd to ensure an equal number of sample points either side of the centre coefficient. 
Also note that no smoothing will occur if $k = F_w - 1$. 
Although some degree of exploratory analysis may be desirable to find the optimal combination of $k$ and $F_w$ parameters for the dataset under analysis, a convenient heuristic is to set the length of $F_w$ approximately 1 to 2 times the anticipated FWHM of the PAF component [@enke1976; @press1992]. <!-- IBS: Doesn’t this detract somewhat from the ‘objectivity’ and replicability of the approach? Can we provide a more precise default choice (or an algorithm for finding one) to ensure that anyone using the procedure with this standard parameter choice will obtain the same results? The same applies to all other parameters needing to be set… but perhaps section 2.2.4 will clarify this. AC response: I completely agree this seems like a concession that undercuts some of the previous arguments. The 'cost' of using the SGF is that it comes with this flexibility in choice of Fw and k, the setting of which (general heuristics aside - which are only of limited use assuming peak width varies across individuals) can look a bit arbitrary. I also think some degree of exploratory fitting is justified if trying to extend the method out to populations that may not be expected to show the sort of classical peak components we associate with healthy young adults (older adult/clinical/children etc). Of course, so long as param settings are known, there is no problem on the replication side of things. I guess the only worry is people might fish around to get the peaks they want / exclude the ones they don't -->
Relatively higher $F_w$ lengths are expected to result in more aggressive smoothing of the input function [@bromba1981], and hence may be desirable in cases of noisy spectral densities.
We note however that excessively flat peaks following application of the smoothing procedure are indicative of a suboptimally large $F_w$.
We favour higher-order polynomials (e.g., $k = 5$) due to their peak-height preserving properties, but acknowledge that they might render suboptimal fits (and less smoothing) in the context of relatively broad component structures [@press1992]. 
This may be of particular concern for instance when dealing with elderly populations, given that spectral power (and relative peak dominance) is typically diminished in older adults [e.g., @chiang2011; @dustman1999].

Additional parameters include:

- $W_\alpha$, the domain of the PSD searched for evidence of peak activity (corresponds to the putative alpha bandwidth);
- $minP$, the minimum quantity of normalised power required to qualify as a potential PAF candidate (determined in relation to the statistical properties of the individual spectrum);
- $pDiff$, the minimum proportion of peak height by which the highest peak component must exceed all other competitors in order to qualify as the PAF;
- $cMin$, the minimum number of channel estimates necessary for computing cross-channel averages. 

Examples of what we consider to be reasonable parameter values are outlined in section 2.3.4.

#### 2.2.3 Overview of analysis procedure
The analysis pipeline can be summarised as follows<!-- AC: ?include a diagram of the procedural flow-->: Once parameters have been defined and preprocessed resting-state data have been imported into the workspace, the PSD is estimated for each channel to be included in the analysis. 
Each PSD is then subjected to the following procedure:
First, the PSD is truncated to exclude extraneous frequency bins located beyond the span of the filter passband.
The data are then normalised by dividing each ordinate value by the mean spectral power calculated across all frequency bins included in the spectrum. 
The SGF is then applied in order to estimate the zeroth (i.e. smoothed), first, and second derivatives of the PSD function.
Next, the first derivative is searched for downward going zero crossings within $W_\alpha$. 
If no zero crossings are identified, the channel is excluded from further PAF-based analysis (it may still be included in later CoG estimates). 

Identified zero crossings are first assessed as to whether they satisfy $minP$, which is calculated through fitting a least-squares linear regression to the log~10~-transformed power spectrum. 
In order to register as a candidate peak, the PSD estimate at the zero crossing must exceed the corresponding value predicted by the regression model by more than the standard deviation of the estimated prediction error. 
This parameter thus provides a convenient threshold for the exclusion of zero crossings emanating from trivial fluctuations in the power spectrum (see [Fig_minPow](#minPow)).
If more than one peak within the spectral domain defined by $W_\alpha$ is found to exceed $minP$, these candidates are rank ordered according to their normalised power estimates, and the magnitude difference between the two largest peaks is compared. 
If the primary peak exceeds the height of its closest competitor by more than the threshold defined by $pDiff$, it is assigned as the channel PAF. 
The second derivative is then examined to determine the location of the associated peak inflection points, and the $Q$ value computed.
Should the primary peak fail to satisfy the $pDiff$ criterion, this channel would fail to register an estimate of the PAF. 
It would however still qualify for inclusion in the alpha window localisation procedure required for estimating the CoG.

![*Fig_minPow.* Visualisation of smoothed power spectral density (PSD) plots with corresponding $minP$ thresholds superposed (red curve, inverse log~10~ transformation of the regression line). *Left panel*: PSD estimates for all frequency bins beyond the delta band fail to exceed $minP$; no peak registered for this channel. *Central panel*: Data from the same participant as in the left panel. In this channel, the spectral peak at ~10 Hz is sufficient to exceed the $minP$ threshold. *Right panel*: Data from another participant showing a marked alpha peak that comfortably exceeds $minP$. Note differences in ordinate scaling.](figs/minPow.png){#minPow}

CoG calculation follows the standard procedure described by Klimesch and colleagues [-@klimesch1990], with the exception that the bounds of each channel's alpha interval were detected automatically. 
The programme derives these bounds by taking the left- and right-most peaks within $W_\alpha$ (i.e. those peaks in the lowest and highest frequency bins, respectively; these may coincide in some cases with the PAF), and searching the first derivative for evidence of the nearest local minimum prior to the left-most peak ($f_1$) / following the right-most peak ($f_2$).
Since some spectra show a relatively shallow roll-off as the edges of the alpha peak diminish, and thus do not culminate in a local minimum for several frequency bins, we relaxed the requirement for an upward going zero crossing (i.e. evidence of a local minimum) such that the transition into a prolonged shallow function is taken as sufficient evidence of the individual alpha bounds $f_1$ or $f_2$.
This criterion was formalised as follows: 

$$f_1 = f’(x) < |1| \text{ for } f(x_{j}...x_{k}),$$ 
$$f_2 = f’(x) < |1| \text{ for } f(x_{k}...x_{l}),$$

where $f(x_k)$ is the first encountered frequency bin where $f’(x) < -1$, $f(x_{j})$ is the frequency bin whose centre frequency is ~1 Hz below that of $f(x_{k})$, and $f(x_{l})$ is the frequency bin whose centre frequency is ~1 Hz above that of $f(x_{k})$ (exact scaling will depend on the frequency resolution of the analysis). <!-- AC: ? better way to notate this -->
$f_1$ and $f_2$ estimates from each eligible channel are averaged to yield the individualised alpha window.
This window defines the index of summation (i.e. frequency band coverage) used to calculate the CoG across all available channels.

If a sufficient number of channels (as stipulated by $cMin$) furnish PAF and individual alpha window estimates, channel PAF and CoG estimates are averaged to generate IAF summary statistics.
Mean PAF ($PAF_M$) is a weighted average that takes into account the $Q$ values associated with each peak estimate:

$$ PAF_M = \frac{\sum\limits_{c=1}^C PAF_c \times \lambda_c} {\sum\limits_{c=1}^C \lambda_c}, $$

where $c$ identifies the channel drawn from the set of all available channels $C$, and $\lambda_c$ is the channel weighting derived by dividing $Q_c$ by the maximum $Q_c$ in $C$ (such that $\sum \lambda_c = 1$).
In contrast to $PAF_M$, all CoG channel estimates contribute equally to the calculation of mean CoG ($CoG_M$).
If there are an insufficient number of channel estimates to satisfy $cMin$, no $PAF_M$ or $CoG_M$ estimates are returned (in some cases, $cMin$ will be satisfied for $CoG_M$, but not $PAF_M$, on account of the latter's more stringent criteria).

Given that resting-state EEG is often recorded both before and after an experimental session, we also included the facility to compute repeated-measures comparisons and grand averages across IAF summary statistics. 
Indeed, since both the quality and quantity of intraindividual alpha-band activity may vary across recordings [see for e.g., @samaha2015], we recommend such cross-recording comparisons in order to maximise the likelihood of deriving reliable IAF estimates. 
Since separate EEG recordings may not be equivalent in terms of quantity of information rendered, grand averaged IAF estimates ($IAF_{GA}$) are calculated using a weighted average that takes into account the proportion of channels that contributed to each constituent summary statistic:

$$ IAF_{GA} = \frac{ IAF_1 \beta_1 + IAF_2 \beta_2 } {\beta_1 + \beta_2} , $$

where either $PAF$ or $CoG$ are substituted in place of $IAF$, $\beta$ constitutes the weighting afforded to the channel-wise mean estimates derived from each recording, and subscript indices indicate the identity of the EEG recording.
For PAF estimates, $\beta$ is the number of channels used to estimate $PAF_M$ divided by total number of channels included in the analysis. 
For CoG estimates, $\beta$ is the number of channels used to estimate the mean individual alpha bandwidth divided by total number of channels included in the analysis.<!-- AC: Let me know if there are convenient ways of capturing these details in equation format, to reduce need for in-text explication -->

### 2.3 Empirical EEG data

#### 2.3.1 Participants
Sixty-three right-handed [Edinburgh Handedness Inventory; @oldfield1971], native English-speaking adults (42 female, mean age = 35 years, age range = 18-74 years) with normal (or corrected-to-normal) vision and audition, and no history of psychiatric, neurological, or cognitive disorder, participated in the study. 
All participants provided written, informed consent, and received financial remuneration for their time. 
This study, which formed part of a larger research project investigating EEG responses to complex, naturalistic stimuli [@gysin-websterInprep], was approved by the University of South Australia Human Research Ethics Committee (Application ID: 0000035576).

#### 2.3.2 Procedure
Following screening and consent procedures, participants were seated in a dimly-lit, sound-attenuated room for the duration of the session. 
Two sets of resting-state EEG recordings were acquired approximately 90 min apart at the beginning and end of an experimental procedure. 
This procedure involved watching approximately 70 min of pre-recorded television programming, followed by an old/new cued recall task. 
As per our standard laboratory protocol, both sets of resting-state recordings comprised of approximately 2 min of eyes-open EEG followed by 2 min of eyes-closed EEG. 
Participants were instructed to sit still, relax, and avoid excessive eye movements during this time. 
In total, the entire session lasted between 2.5 and 3 hr. 
Note, only data from the eyes-closed component of the resting-state recordings will be reported in the analyses that follow.
Although our approach could easily be extended to the analysis of eyes-open EEG, we favour eyes-closed data on the basis that it demonstrates (1) greater interindividual variability in alpha power [@chen2008], and (2) improved within-session reliability and test-retest stability of IAF estimates [@grandy2013a].
Eyes-closed recordings are also advantageous in reducing the incidence of ocular artifact.

#### 2.3.3 EEG acquisition and preprocessing
EEG was recorded continuously from 64 cap-mounted Ag/AgCl electrodes via Scan 4.5 software for the SynAmpsRT amplifier (Compumedics^&reg;^ Neuroscan&trade;, Charlotte, NC, USA). 
The online recording was digitised at a rate of 1000 Hz, bandpass filtered (passband: 0.05-200 Hz), and referenced to the vertex electrode (AFz served as the ground electrode). 
Eye movements were also recorded from bipolar channels positioned above and below the left eye, and on the outer canthi of both eyes. 
Electrode impedances were maintained below 12.5 k$\Omega$.

EEG data acquired during eyes-closed resting-state recordings were preprocessed in MATLAB 2015a (version 8.5.0.197613). 
First, all EEG data channels were imported into the MATLAB workspace via EEGLAB version 13.6.5b and re-referenced to linked mastoids. 
Each dataset was then trimmed to retain only the EOG and the nine centro-posterior electrodes that constituted the region of interest for resting-state IAF analysis: Pz, P1/2, POz, PO3/4, Oz, O1/2. 
Trimmed datasets were subjected to zero-phase, finite impulse response (FIR) highpass (passband: 1 Hz, -6 dB cutoff: 0.5 Hz) and lowpass (passband: 40 Hz, -6 dB cutoff: 45 Hz), Hamming-windowed sinc filters.
Automated artifact detection routines where applied to identify regions of channel data (2 s segments) containing frequency components between 15 and 40 Hz that exceeded the upper spectral threshold by 10 dB.
Channels that exhibited improbable data distribution (kurtosis z-score > 5) were excluded from analysis.
EOG channels were then removed, and remaining channels downsampled to 250 Hz in preparation for spectral analysis.
If datasets exceeded 120 s following artifact rejection, they were trimmed down to this duration in order to reduce variability in the quantity of data contributed by each participant.

#### 2.3.4 IAF analysis parameters
Initial parameters for the IAF analysis were determined on the basis of preliminary testing with an independent set of resting-state data. (These data were collected as part of a separate EEG protocol).

The length of the Hamming window implemented in `pwelch` was set at 2048 samples (4 times the sampling rate raised to the next power of 2), which yielded a frequency resolution of ~0.24 Hz. This sliding window was applied with 50% overlap. The bounds of the alpha-band window $W_\alpha$ were set at 7 and 13 Hz. Width of the SGF sliding window $F_w$ was set at 11 samples (i.e. frequency bins), corresponding to a frequency span of ~2.44 Hz. 
A fifth-degree polynomial was selected as the curve-fitting parameter $k$.

The minimum power difference parameter was set at $pDiff = .20$. This meant that the largest peak detected within $W_\alpha$ had to register a PSD value at least 20% greater than that of the next largest peak to qualify as the PAF estimate for that channel. 
Although this criterion constitutes an arbitrarily defined threshold, we know of no objective standard (or theoretical rationale) that can be used to guide the determination of this necessary boundary condition. 
However, we consider it a strength of the current approach that this limit must be explicitly defined prior to data analysis, and uniformly applied across all encountered cases.

In addition to investigating the extent to which IAF estimates (and related alpha-band parameters) varied across the sample, we also performed within-subject comparisons of pre- and post-experiment spectral data. 
For both of these intra- and interindividual analyses, the minimum number of channel-wise PSD estimates required for PAF and CoG calculation was set at $cMin = 3$. 
In the event that only one of the paired recordings satisfied $cMin$, IAF estimates were derived solely from this data. 
From a methodological perspective, we were interested to observe how many cases failed to satisfy criteria necessary for PAF and/or CoG estimation; and in particular, the extent to which CoG provided additional information in cases where the PAF could not be reliably ascertained.

### 2.4 Simulated EEG data

#### 2.4.1 Overview of generic simulation procedure
Synthetic resting-state EEG data were generated by combining a sine wave time series oscillating at some alpha band frequency with a time series whose frequency distribution broadly conformed to the $1/f$ inverse power-law scaling characteristic of lower frequency resting-state M/EEG activity [@novikov1997;@pritchard1992].
This latter background or pink noise series was produced using the `pinknoise` MATLAB function (Zhivomirov, 2013).
This programme works by generating a vector of Gaussian-distributed random values, transforming this series into the frequency domain in order to rescale these samples according to the inverse power-law, and then converting back into the time domain after mean-centring and unity-based normalisation.
This pink noise time series was then multiplied with an alpha signal generated from sine waves oscillating at a specified centre frequency $f_c$.
The resulting composite signal was then subjected to PSD estimation and smoothing (see [Fig_sim_peaks](#sim_peaks) for an illustration of the signal generation procedure).
All synthetic signals were designed to replicate the 250 Hz sampling rate and 2 min duration of the empirical data reported above.

![*Fig_sim_peaks.* Illustration of the simuation general scheme for constructing synthetic resting-state EEG data. *Top and central rows*: 2 s portion of randomly synthesised alpha and pink noise signals, together with their respective power spectral densities (right panels; log-scaled). *Bottom row*: Time series produced by combining the above two signals together (point-by-point multiplication in the time domain). *a.u.*: arbitrary unit.](figs/sim_peaks.png){#sim_peaks}

#### 2.4.2 Preliminary analysis: Capability to recover PAF amidst varying levels of background noise
As an initial proof of concept, we investigated the capacity of the SGF technique to extract accurate PAF estimates from synthetic spectra.
These spectra were comprised of a single alpha-band component embedded within a time series whose frequency properties were scaled to approximate the $1/f$ inverse power-law.
We systematically varied the SNR of the target alpha component relative to the pink noise signal in order to assess how well the PAF estimation routine performed as a function of relative alpha peak power.
SNR was operationalised as the proportion of the composite (alpha plus pink noise) signal containing information from the alpha-band time series.
Hence, for an SNR of 1, each sample point of the composite signal was the product of the corresponding sample points in the alpha and pink noise time series; whereas only half of the alpha time series contributed to the composite signal when SNR was set to 0.5.

We examined PAF estimation at SNR = 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, and 0.50.
For each set of simulations, 1000 PSD estimates were generated using the afore-described method of embedding a singular alpha component within a pink noise signal.
The frequency of each alpha signal was randomly assigned by sampling (with replacement) a vector of frequency values ranging from 7.5 to 12.5 Hz in 0.1 Hz interations.
Given that the target alpha peak consisted of a singular frequency component, we expected the peak estimation routine to perform with a high degree of accuracy (within the limits of spectral resolution) when a peak could be discriminated from background spectral noise.
We were however particularly interested in establishing the proportion of spectral peaks that would be detected across various levels of SNR.

#### 2.4.3 Simulation of heterogeneous alpha-band activity within multi-channel datasets
To ascertain how the SGF method performs under more realistic conditions, we simulated datasets that comprised a randomly selected alpha peak situated at the centre of an approximately Gaussian-distributed set of neighbouring components.
This configuration was designed to replicate empirical observations of broader alpha peaks, which may arise on account of distributed alpha-generators oscillating at different frequencies to one another.
To achieve this effect, we generated a matrix of alpha-band oscillations from which the composite signal was drawn.
A Gaussian window spanning 5 Hz was centred on the target PAF, and those frequencies spanned by the window were sampled in proportion to the height of the window function (such that the centre frequency contributed most sample points, followed by the two adjacent signals, and so forth).
The area of the Gaussian window was varied by manipulating its $\alpha$ value (the reciprocal of the standard deviation), where higher $\alpha$ corresponds to a narrower distribution.
We stratified our simulations into three $\alpha$ levels: 4.0, 2.5, and 1.0. 

Once the composite alpha signal had been constructed, it was then combined with a randomly sampled pink noise signal at a given SNR.
This final step was repeated nine times for each simulated dataset, such that the alpha signal was replicated and combined with an independently sampled pink noise signal.
Consequently, each dataset contained nine synthetic 'channels' that differed from one another as a consequence of stochastic variance in the background characteristics of the spectrum (while comprising the same alpha component signal).
This enabled us to examine consistency of IAF estimation in the context of random noise fluctuations, as well as providing information about rates of channel rejection as a function of alpha-band distribution at differing levels of SNR.

## 3 Results

### 3.1 Empirical EEG data

#### 3.1.1 Global performance of the IAF estimation routine
Post-experiment resting-state EEG recordings were missing for three participants.
These individuals thus contributed only one set of channel estimates each to the following analysis.
A further 11 channels were excluded from analysis on the basis of excessive kurtosis (note, no more than 1 channel was rejected for any given recording).
This left a grand total of 1096 PSDs to estimate across the sample (561 pre- and 535 post-experiment spectra).
From these channel data, a total 948 PAF estimates (pre = 472, post = 476) and 1009 CoG estimates (pre = 506, post = 503) were extracted across all recordings.

[Fig_stacked_chans](#stacked_chans) shows the relative distribution of channel estimates derived from all participants.
This figure indicates that the estimation routine derived PAF and CoG estimates from most (if not all) available channels in the majority of cases.
Of those participants who registered a limited number of channel estimates, two (#29 and #50) failed surpass the minimal threshold (i.e. $cMin = 3$) required for $PAF_M$ estimation. 
This result indicates that these participants should be excluded from IAF-related analyses on account of the general lack of alpha-band activity in the majority of their channel spectra.
Visual inspection of these individuals' PSD plots confirmed that neither participant's data demonstrated compelling evidence of a distinct alpha peak component ([Fig_no_pafs](#no_pafs)).

![*Fig_stacked_chans.* Stacked bar chart displaying the number of channels from which PAF (lower half) and CoG (upper half) estimates were derived across participants. PAF and CoG estimates are further divided according to order of EEG recording. Totals normalised to take into account excluded channels. Note, post-experiment data where unavailable for participants #26, 56, and #57.](figs/stacked_chans.png){#stacked_chans}

![*Fig_no_pafs.* Representative channel PSDs (three superposed channel spectra) from two sets of EEG recordings in which the IAF estimation routine failed to detect evidence of substantive alpha rhythm activity (top row: participant #29; bottom row: participant #50). Shaded areas indicate range of the alpha peak search window $W_\alpha$; red line indicates representative $minP$ threshold.](figs/no_pafs.png){#no_pafs}

#### 3.1.2 Estimator distributions & correlation coefficients
Consistent with previous reports [@kopruner1984;@klimesch1990], mean estimates of IAF were centred about 10 Hz, with the majority of estimates falling in the range of 9 to 11 Hz.
$PAF_M$ and $CoG_M$ were similarly distributed across both sets of recordings (see histograms, [Fig_pre_post](#pre_post)).
Intraclass correlation coefficients ($ICC_{3,k}$: $PAF_M = .97; CoG_M = .94$) indicated that variance in $PAF_M$ and $CoG_M$ estimates can be predominantly attributed to interindividual differences across the sample, rather than intraindividual differences between recordings (see scatterplots, [Fig_pre_post](#pre_post)).
These data are therefore in accord with previous studies reporting evidence of the IAF's high temporal stability (at least within the same recording session) and interindividual variability (at least in the context of eyes-closed resting-state recordings).

![*Fig_pre_post.* *Top and central rows*: Histograms displaying distribution of $PAF_M$ (left panel) and $CoG_M$ (right panel) estimates across recordings. *Bottom row*: Scatterplots showing high degree of agreement between corresponding estimators of IAF calculated from pre- and post-experiment resting-state recordings. Broken line indicates perfect positive correlation between pre- and post-experiment IAF estimates.](figs/pre_post.png){#pre_post}

Final analyses were conducted on grand-averaged alpha peak and gravity estimates ($PAF_{GA}$ and $CoG_{GA}$, respectively).
Note that, in cases where only one channel mean estimate was registered, this value was included as the participant's grand average.
Kernal density estimation of the probability density function underlying $PAF_{GM}$ data suggested that sample estimates are well-characterised by a Gaussian distribution ([Fig_GA_distribs](#GA_distribs)).
This is in keeping with previous efforts to quantify the distribution of PAF [@klimesch1996].
$CoG_{GA}$ however showed some evidence of deviation from the Gaussian distribution on account of its somewhat heavier tails.
To the best of our knowledge, this is the first study to examine the distributional properties of the $CoG$.
It is therefore unclear how our estimation procedure compares to previous methods of calculating alpha-band $CoG$ in this regard.
In spite of this apparent difference in distributional characteristics, $PAF_{GA}$ and $CoG_{GA}$ produced remarkably consistent results ($ICC_{3,k} =.97$; $R^2 = .88$).
This finding, which extends that reported in a smaller sample by Jann, Koenig, Dierks, Boesch, and Federspiel [-@jann2010], lends support to the claim that these two estimators tap into the same fundamental oscillatory process.

![*Fig_GA_distribs. *Top row*: Probability density of $PAF_{GA}$ (left panel) and $CoG_{GA}$ (right panel) following application of a Gaussian kernal density estimator. *Bottom row*: Normal probability plots indicating how well $PAF_{GA}$ and $CoG_{GA}$ estimates approximate a normal distribution (indicated by the broken red line). While $PAF_{GA}$ estimates vary in a fashion broadly consistent with an underlying normal distribution, the $CoG_{GA}$ estimator shows evidence of a somewhat heavier-tailed distribution.*](figs/GA_distribs.png){#GA_distribs}


![*Fig_intercorrel.* Scatterplot displaying association between grand averaged PAF and CoG estimates across all participants for whom both estimators were calculated. Broken line indicates perfect positive correlation between IAF estimators.](figs/intercorrel.png){#intercorrel}

### 3.2 Simulated EEG data

#### 3.2.1 PAF estimation performance as a function of SNR
Preliminary analysis of synthetic channel data focused on the number of PAF estimates extracted for each level of SNR, and the extent to which these estimates approximated the ground truth as stipulated by the alpha frequency sine wave component of the simulated signal [see Table for summary statistics across SNR conditions].
While the SGF technique failed to extract PAF estimates from almost one-third of simulations when SNR = 0.05, the proportion of estimated alpha peaks rapidly approached ceiling as SNR increased beyond 0.10.
RMSE was generally low for all levels of SNR, suggesting that alpha peaks were on average estimated with a high degree of accuracy when detected by the analysis routine.
The SGF technique approached near-optimal performance at SNR = 0.20, in which 1% of simulated spectral peaks were either undetected or located within an adjacent frequency bin.
The routine performed optimally at SNR ≥ 0.30 (with the exception of one undetected peak at SNR = 0.30).
At this level of SNR, all alpha peak components were correctly estimated within the limits of the stipulated frequency resolution (~0.24 Hz bin intervals). 

Table: Summary statistics characterising peak alpha frequency (PAF) estimation as a function of signal-to-noise ratio (SNR). 
*n PAF*: total number of PAF estimates extracted from 1000 simulated time series; *RMSE*: root mean squared error of PAF estimates; *maxDiff*: maximum absolute difference between PAF estimates and corresponding frequency (Hz) of underlying alpha-band component; *n shift*: number of PAF estimates that diverged from the target alpha frequency by more than 0.24 Hz, the frequency resolution of the analysis.

| SNR | 0.05 | 0.10 | 0.15 | 0.20 | 0.25 | 0.30 | 0.40 | 0.50 |
|---|---|---|---|---|---|---|---|---|
| *n PAF* | 672 | 914 | 988 | 993 | 998 | 999 | 1000 | 1000 |
| *RMSE* | 0.12 | 0.10 | 0.09 | 0.07 | 0.08 | 0.07 | 0.07 | 0.07 |
| *maxDiff* | 1.13 | 0.78 | 0.68 | 0.40 | 0.40 | 0.15 | 0.16 | 0.15 |
| *n shift* | 25 | 19 | 10 | 3 | 5 | 0 | 0 | 0 |

Although we were interested in observing how the SGF technique performed in comparison to the `pwelch` function on which it depends, this was not of primary concern to the present analysis.
This was because one should expect the `pwelch` routine to perform well in cases where the alpha-band consists of only one dominiant frequency component, especially if the SNR is relatively high (recall that the primary motivation for the SGF approach is to resolve IAF estimates when more ambiguous - and thus problematic - spectral conditions obtain; e.g., split-peaks).
To give a flavour of how smoothing influenced `pwelch` PSD estimates, a random selection of the PSD estimates produced by both methods is illustrated in [Fig_snr](#snr).
We note that the SGF resulted in a distorted estimate of the true PAF in the SNR = 0.05 simulation depicted in the top left panel of [Fig_snr](#snr) (indeed, this is one case where the standard `pwelch` estimation routine would in fact be preferable in virtue of its superior PAF estimate).
The SNR = 0.10 simulation however shows the potential advantage of the SGF technique.
In this case, the SGF renders a veridical estimate of the underlying alpha component, which otherwise would have been falsely resolved (and rejected) as a split-peak using `pwelch`.
At higher SNRs, the two techniques produce identical PAF results (although the SGF attenuates peak height, as would be expected of a smoothing operation).

![*Fig_snr.* Synthetic channel spectra randomly selected from each simulated SNR condition. Blue functions represent PSD estimates generated by the standard MATLAB implementation of Welch's periodogram method (`pwelch`). Red functions signify outcome of subjecting `pwelch` estimates to the Savitzky-Golay filter (SGF) applied to the empirical data. $F\alpha$: The true frequency of the target alpha component. $PAF_{PW}$ and $PAF_{SG}$: Estimates of $F\alpha$ rendered by the `pwelch` and SGF methods, respectively.](figs/snr.png){#snr}

Ten percent of the simulation data that failed to produce PAF estimates in the SNR = 0.05 condition were randomly selected for visual inspection to evaluate the source of these omissions.
This analysis confirmed that the vast majority of sampled spectra lacked a distinct alpha peak (irrespective of the adopted PSD estimation method).
Five of the spectra within this subset were rejected due to the appearance of more than one spectral peak within the alpha window (none of which were clearly dominant).
Overall, it seems unlikely that any of these spectra would be considered to satisfy the criteria for PAF assignment as stipulated in section 1.1.<!-- AC: ? make plots available in supp mats -->
Visual inspection of the 25 PAF estimates that deviated most from their target frequencies within the 0.05 SNR condition confirmed that SGF estimates generally tracked their corresponding `pwelch` counterparts closely.
Where smoothed estimates were relatively distorted (i.e. displaced into a neighbouring frequency bin due to the influence of spurious local fluctuations), the more accurate `pwelch` estimate was by the same token disadvantaged by having to compete with these nearby spectral subpeaks.
For the most part, then, the least accurate estimates rendered via the SGF technique did not appear to be substantially worse than what would otherwise have been estimated via the `pwelch` technique. 

In sum, this preliminary analysis provides strong initial evidence that the SGF method generally furnishes highly accurate estimates of the PAF when a singular alpha component is present within the PSD.
This degree of accuracy is maintained even at relatively low levels of SNR, although reliable resolution of low powered spectral peaks amidst background noise becomes more challenging when SNR drops below 0.15 (at least when the SGF technique is implemented with the same parameters as those used in our analysis of the empirical EEG data).
Diminished sensitivity at low levels of SNR seems to derive from a generic reduction in the capacity of the underlying `pwelch` routine to accurately extract low powered alpha components under such conditions, rather than any specific deficiency of the SGF procedure per se.

#### 3.2.2 Multi-channel dataset simulations
Given that PAF estimation approached ceiling performance at moderate levels of SNR in the previous analysis, we limited our analysis of multi-channel simulaton data to a low (0.15) and a moderate (0.40) SNR condition. 100 datasets, each comprising 9 synthetic EEG channels, were simulated for each of the three levels of alpha component dispersal in both SNR conditions (resulting in a total 5400 PSD estimates).
The results of this analysis are summarised in Table.

Table: PAF and CoG estimation performance as a function of SNR (0.40 vs. 0.15) and alpha component distribution. Synthetic alpha peak sampled using a Gaussian windowing function whose inverse standard deviation was parametrically varied (where $\alpha$ = 4.0 corresponds to a narrow peak, $\alpha$ = 1.0 a broad peak). *$PAF_{PW}$*: PAF estimated via averaging of channel-wise PSD estimates furnished by `pwelch`; *$PAF_{SG}$*: PAF estimated via averaging of Savitzky-Golay smoothed, Q-weighted channel-wise PSD estimates; *$CoG$*: CoG estimated by SGF routine; *RMSE*: root mean squared error of PAF/CoG estimates; *maxDiff*: maximum absolute difference between PAF/CoG estimates and corresponding peak frequency (Hz) of underlying alpha-band component; *% Dev*: percentage of PAF/CoG estimates that diverged from the target alpha frequency by more than 0.5 Hz; *n chans*: median number of channel estimates per simulated dataset retained for PAF/CoG estimation (and associated standard deviation).

SNR |  | 0.40 |  |  | 0.15 |  |
|---|---|---|---|---|---|---|
$\alpha$ | 4.0 | 2.5 | 1.0 | 4.0 | 2.5 | 1.0 |
**RMSE** |  |  |  |  |  |  |
$PAF_{PW}$| 0.23 | 0.44 | 0.81 | 0.32 | 0.38 | 0.85 |
$PAF_{SG}$| 0.12 | 0.22 | 0.65 | 0.23 | 0.27 | 0.70 |
$CoG$     | 0.11 | 0.16 | 0.61 | 0.28 | 0.50 | 0.70 |
**maxDiff** | | | | | | |
$PAF_{PW}$| 0.53 | 1.14 | 2.02 | 0.93 | 0.82 | 2.33 |
$PAF_{SG}$| 0.35 | 0.68 | 1.75 | 0.61 | 0.57 | 2.34 |
$CoG$     | 0.35 | 0.46 | 1.20 | 0.81 | 1.50 | 1.56 |
**% Dev** | | | | | | |
$PAF_{PW}$| 4 | 32 | 44 | 12 | 22 | 54 |
$PAF_{SG}$| 0 | 1 | 43 | 4 | 1 | 49 |
$CoG$     | 0 | 0 | 39 | 7 | 32 | 58 |
**n chans (s.d.)** | | | | | | |
$PAF_{SG}$| 8 (0.86) | 7 (1.37) | 6 (1.58) | 7 (1.27) | 6 (1.67) | 6 (1.64) |
$CoG$| 9 (0.14) | 9 (0.14) | 9 (0.14) | 9 (0.43) | 9 (0.52) | 9 (0.90) |

For the moderate SNR condition, PAF and CoG estimates became less accurate on average (as indexed by RMSE) as the dispersal of the target component increased (i.e. as synthetic alpha peaks became broader and less prominent).
This overarching pattern was replicated in the low SNR condition, which also tended to produce estimates with higher average errors than those generated at corresponding levels of alpha dispersal in the moderate SNR condition.
The magnitude of maximum absolute estimation errors likewise tended to increase as a function of synthetic alpha peak dispersal (although this pattern was less consistently realised in the low SNR simulations).
In both SNR conditions, the least accurate estimates of broad peak components deviated from their target frequencies by more than 1 Hz. 
The proportion of PAF and CoG estimates that substantially erred from their target frequency (i.e. deviated by more than 0.5 Hz) likewise tended to increase as a function of alpha component dispersal.

Of the 600 sets of channel data simulated across all alpha distribution $\times$ SNR conditions, PAF estimates could not be derived from 6 datasets (1 from SNR = 0.30, $\alpha$ = 1.0; 4 from SNR = 0.15, $\alpha$ = 2.5; 1 from SNR = 0.15, $\alpha$ = 1.0).
The average number of simulated channels that contributed to the calculation of cross-channel PAF estimates generally diminished as alpha dispersal increased.
The number of channels in which a viable alpha peak was detected also demonstrated greater variability as function of both alpha distribution and SNR, with broader component structures and low SNR both associated with increased volatility of channel selection.
Together with the consistent trends towards increased estimation error described above, these observations support the notion that broader, flatter alpha peaks pose a greater challenge to automated PAF detection (even under moderate levels of SNR).
By contrast, the SGF routine was able to extract CoG estimates for all  simulated dataset, and retained a higher proportion of channel data across all conditions.
Indeed, channel inclusion was close to ceiling for the moderate SNR conditions, and remained high despite increasing variability as a funtion of increased alpha component dispersal in the low SNR conditions.

Amongst the two methods for estimating PAF, the Savitzky-Golay smoothed, Q-weighted average estimates were consistently less error-prone (i.e. registered lower RMSEs across conditions) than those derived by averaging across the 'raw' (i.e. unsmoothed) channel-wise PSDs generated by `pwelch`.
The latter approach also produced more extreme deviations from the target PAF in all but the least favourable simulation condition (i.e. SNR = 0.15, $\alpha$ = 1.0), in which the maximum error magnitude was roughly equal to that of the SGF technique.
Significantly, the SGF method produced far fewer estimates that deviated substantially (> 0.5 Hz) from the target peak frequency in all simulations except those featuring the most broadly distributed components (i.e. $\alpha$ = 1.0; where performance remained superior by a less substantial margin).
This pattern was observed irrespective of SNR, and is strongly suggestive of the SGF technique's capacity to improve upon the precision of `pwelch` furnished PAF estimates.

The SGF method produced very similar PAF and CoG estimates in moderate SNR condition, while the relation between these estimators was somewhat less consistent in the low SNR condition.
Interestingly, CoG estimates resulted in mild reductions in the average error, maximum error, and proportion of substantial estimate deviation registered by $PAF_{SG}$ in the moderate SNR condition.
This result provides encouraging evidence in favour of the SGF method's capacity to accurately localise the beginning and end of the alpha-band component (i.e. $f_1$ and $f_2$), at least when the peak is not too broadly distributed.
By contrast, the inferior performance of the CoG estimator relative to $PAF_{SG}$ across the most of low SNR simulations suggests the routine is less reliable when the alpha component is embedded amidst a relatively noisy signal.
This is unsurprising, since identifying the bounds of the alpha component is more challenging (and thus, more susceptible to error that biases CoG estimation) under noisy signal conditions.
Note, however, that while the CoG estimator was more prone to errors when SNR was low, the magnitude of estimate deviations from the target frequency did not reach the extremes observed in both methods of PAF estimation (≥ 1.75 Hz) when simulated components were broadly distributed.

## 4 Discussion
We have proposed a novel method for the automatic detection of the two most prevalent indices of individual alpha frequency in the literature.
This method pairs a common approach to the automated detection of local maxima (i.e. searching for first derivative zero crossings) with a well established method of extracting spectral peaks amidst noisy data (i.e. Savitzky-Golay filtering) to derive an estimate of peak alpha frequency (PAF).
It also extends the logic of the first derivative analysis to estimate the bounds the alpha peak component, thus enabling calculation of the alpha-band mean or centre of gravity frequency (CoG).
Like other automated peak detection (i.e. curve-fitting) algorithms in the literature, this method addresses key limitations of visual PSD analysis (e.g., proneness to subjective bias, inefficiency, replicability challenges), and improves upon simpler automated approaches that may be prone to various artifacts (e.g., failure to differentiate multiple peaks, assigning peaks to noisy spectral fluctuations).
Unlike these more sophisticated algorithms, however, our method is openly accessible and easy to integrate within existing EEGLAB and Python analysis pipelines.

Our results demonstrate that the SGF technique can extract a large proportion of IAF estimates from an empirical dataset, and that the general features of these estimates (intraindividual stability, interindividual variance, etc) are consistent with prior reports in the literature.
Furthermore, application of the technique on simulated data confirmed its capacity to render accurate estimates of peak location under relatively degraded SNR conditions.
Extended to more complex simulations, the SGF technique was shown to recover target IAF values with greater precision than were achieved by applying an automated peak detection method to non-smoothed PSDs.
We shall begin by discussing some of the important findings of the empirical analysis, before turning our attention to the simulation results.

### 4.1 Estimation of IAFs from an empirical EEG dataset
Savitzky-Golay filtering of `pwelch` generated PSD functions resulted in the extraction of a rather impressive number of IAF estimates from a moderate-sized dataset.
This outcome suggests the technique might offer meaningful benefits over traditional methods of analysis, which can be susceptible to substantial attrition if dominant peaks cannot be confidently distinguished from background noise [e.g. @bornkessel-schlesewsky2015].
We note also that our SGF method resulted in a higher proportion of PAF estimates than that produced by the Gaussian curve-fitting procedure implemented by Haegens and colleagues [-@haegens2014].
It may be the case then that our non-parametric approach, which attempts to smooth the PSD rather than fit a specified function to it, retains more data (i.e. excludes fewer individuals from IAF-related analyses) by virtue of its capacity to accommodate a broader range of alpha-band distributions.

While we cannot categorically rule out the possiblity that sampling error furnished us with a dataset comprised of an unusually high number of individuals manifesting singular, unambiguous alpha peaks, it is encouraging that the SGF technique was able to estimate such a high proportion of PAFs and CoGs.
By the same token, it is also reassuring that the two cases in which the technique was not able to extract IAF estimates did not demonstrate compelling evidence of any concerted alpha-band activity on visual inspection of the relevant PSD plots ([Fig_no_pafs](#no_pafs)).
It is perhaps also worth pointing out that the diverse age range of participants sampled within this study would most likely have posed a robust challenge to any automated peak frequency routine, given the typically reported changes in spectral power and alpha peak distribution associated with older adulthood [e.g., @dustman1999].
That the automated IAF estimation technique reported here was able to extract estimates for the vast majority of individuals, and that it did so using a fixed set of a priori-defined parameters assigned on the basis of preliminary testing in an independent dataset, speaks to its capacity to derive resting-state IAF estimates across a broad spectrum of the healthy population.

Although a high proportion of channel data derived across both pre- and post-experiment recordings were available for grand averaging in most cases, individuals for whom only a subset of channel estimates could be derived tended to demonstrate stronger alpha peak activity in the latter recording (see [Fig_stacked_chans](#stacked_chans)).
This trend was by no means universal, however, thus it seems prudent for researchers to include both pre- and post-experiment resting-state recordings within their protocols wherever possible.
As expected, pre- and post-experiment IAF estimates were highly intercorrelated; hence, it would seem reasonable to treat those estimates derived from singular recordings (i.e. cases in which only one recording yielded an estimate of PAF/CoG) as valid indicators of IAF.
Additional support for this view accrues from the observation that the inclusion of single-recording PAF estimates amongst grand-averaged PAFs did not result in any marked deviation from the expected (approximately) Gaussian distribution reported previously [@klimesch1996].

Comparison of PAF and CoG estimates also revealed a high degree of intercorrelation, despite showing some differences in the distribution of respective grand-averaged estimates.
Although this might prompt concerns of redundancy, we interpret this finding positively: the CoG seems to tap into a similar underlying neural process (or indeed, set of processes) as that indexed by PAF.
Although not necessary in the present analysis on account of the high proportion of PAFs that were extracted, this finding suggests that the CoG estimator might be substituted as an alternative marker of IAF in cases where the PAF cannot be determined (e.g., due to significant split-peaks).
In any case, given the dearth of research directly comparing these two measures (most IAF-related research involves some variant of PAF, perhaps on account of the additional complexities involved in calculating CoG), we suggest that it would be informative if future researchers were to report both of these indices in parallel.
Should it be the case that PAF and CoG track one another almost identically, then only one of these markers need be selected for the remainder of the analysis [as per @jann2010].
However, it could be the case that PAF and CoG diverge under certain circumstances, hence it may be hasty to dismiss the CoG as a redundant alternative to the PAF.
It is of course a notable advantage of the present method that it enables researchers to rapidly derive sample-wide estimates of both PAF and CoG. 
To the best of our knowledge, no other automated technique provides the functionality to simultaneously compute both of these IAF estimators.

### 4.2 Estimation of simulated IAFs
Our preliminary set of simulations showed that the SGF technique performed almost perfectly when 2 min synthetic signals contained at least 24 s of alpha-band oscillations (SNR = 0.20).
Indeed, the peak detection routine performed reasonably well when signals contained as little as 12 s of alpha-band activity, with fewer than 2% of simulated peaks being erroneously estimated.
This analysis served as a basic proof of concept insofar as it demonstrated that the SGF method is capable of extracting a high proportion of underlying peak frequencies without introducing systematic bias.
We acknowledge however that the estimation of sharply defined, single frequency alpha components is unlikely to be representative of genuine electrophysiological data in the vast majority of cases.
Indeed, if spectral data typically conformed to this pattern, distinguising the PAF would prove rather trivial (thus obviating the need for the SGF technique).
While it is encouraging then that the SGF technique performed well in these favourable conditions, it was necessary to demonstrate its capabilities when confronted with more complex, ecologically valid simulation data.

The multichannel simulations were designed to be more faithful to empirical resting-state EEG data, insofar as they comprised of a range of alpha components [thus emulating the multiple alpha generators supposed to underlie the dominant rhythm; @basar1997; @basar2012; @klimesch1999], and a variety of correlated (but nonidentical) estimates of these composite signals (i.e. multiple sources of information were available, emulating the scenario of multi-channel recordings).
This manipulation also enabled us to examine the accuracy of CoG estimates, which in the context of Gaussian-distributed alpha components should ideally converge with PAF.
The critical finding across all simulation conditions was that the SGF technique rendered PAF and CoG estimates that almost always improved upon PAF estimates derived from automated searches of averaged, unsmoothed channel spectra.
This finding held irrespective of whether estimator performance deficits were quantified in terms of the average error across simulated datasets, the magnitude of worst (i.e. most deviant) estimate errors, or the percentage of estimates in the dataset that deviated from the ground truth by more than 0.5 Hz.

IAF estimates from smoothed PSDs were on the whole considerably more accurate in all but the most broadly distributed of alpha-band components.
That the $\alpha = 1.0$ conditions precipitated a large increase in the percentage of estimates deviating from the target frequency by > 0.5 Hz might indicate that the predefined SGF settings were suboptimal (i.e. the filter frame width too narrow relative to the width of the target peak) in the context of such widely dispersed components.
In any case, the fact that error magnitude tended to be greater when PAF estimates were derived from unsmoothed PSDs suggests that the SGF approach is favourable even when the optimality of its analysis parameters is uncertain.

It is also interesting to note that, while CoG estimates were found to be more volatile in the low SNR condition (especially for broader peaks), it was less prone to extreme deviations from the target frequency than either the smoothed or unsmoothed PAF estimator.
This is perhaps unsurprising given that the CoG is sensitive to the overall dispersal of the alpha component, and as such may be skewed by background noise within the spectral region of interest.
On the other hand, the CoG is less likely to be unduly influenced by spurious outlying spectral activity on account of its central tendency-like characteristics.
This may explain then why the CoG estimator results in less extreme maximal error deviations in the broadest component conditions, even in the low SNR condition where it generates the highest percentage of deviant estimates as compared to the other two estimators.
This finding suggests that the CoG might constitute a preferable alternative to the PAF in cases where spectral data are particularly noisy or comprised of poorly-defined peaks (as may be the case in young children and certain clinical populations).
We again return to our earlier argument advocating the extraction and analysis of both PAF and CoG estimates in future IAF research.
This strategy would lead to a deeper understanding of the precise relationship between these two estimators (and the factors that predict their convergence/divergence), and might shed light on the circumstances in which one estimator may be preferred over the other.

### 4.3 Limitations and future developments
We aimed to create a straightforward automated routine that calculates reliable PAF and CoG estimates from posterior channel EEG data recorded during short periods of relaxed, eyes-closed wakefulness.
Although limited in its scope, we believe that the programme could easily be adapted for application across a broader range of empirical contexts (e.g., quantifying spectral peaks in different frequency bands during task-related activity; quantifying peak characteristics across different topographical regions).
It may prove more challenging, however, to accurately resolve estimates of IAF under conditions that are less conducive to the manifestation of a dominant alpha peak (or indeed, in populations that are known to manifest spectral characteristics distinct from those of neurotypical adults).
Further research is therefore needed to examine how well the present technique performs when its application is extended beyond the rather circumscribed conditions investigated herein.

One aspect of performance that was not investigated in the above analyses was whether the accuracy and precision of IAF estimates depend upon the method used to estimate the PSD. 
In its present implementation, the SGF programme relies upon Welch's modified periodogram method to first generate estimates of the PSD that are subsequently subjected to its filtering and differentiation operations.
It may be worthwhile then to investigate whether alternative methods of PSD estimation (e.g., continuous wavelet transform) significantly alter the performance of the SGF technique.

Another possible avenue for optimising the performance characteristics of the SGF would be the development of a function that automatically adapts the $F_w$ (filter width) and $k$ (polynomial degree) parameters in accordance with the approximate span of the dominant frequency component located within the initial search window.
This would involve implementing a two-stage process wherein the features of the alpha-band (or some other band, depending on the research question at stake) component are initially parameterised, and subsequently used to scale $F_w$ and $k$ in accordance with the heuristics described in section 2.2.2 (or some other predefined ratio preferred by the researcher).
Once these parameters have been determined according to the empirical qualities of the spectral data at hand, smoothing and component parameterisation would be performed as described above.

Finally, it would be desirable to create a package that incorporates the MATLAB implementation of the SGF routine within the EEGLAB graphical user interface.
Not only would this help to make the procedure accessible to the broadest possible range of EEGLAB users, it would also provide a convenient platform for integrating visualisations of the spectral analysis that, amongst other things, may assist in the diagnosis of suboptimal parameter settings and guide troubleshooting.
We intend to explore a number of these possibilities in future work.

## 5 Conclusion
We have developed a free, open source programme for automatically estimating individual alpha frequency in resting-state EEG data.
This programme has been shown to perform more accurately than a simpler automated peak detection routine, and may return a higher proportion of empirical IAF estimates than parametric curve-fitting techniques.
Furthermore, this method is not dependent on alpha-band reactivity differences, which may bias IAF estimation. 
In addition to its obvious advantages from the perspective of analytic replicability and efficiency, this method offers a convenient tool for IAF estimation that has the potential to improve the accuracy and precision of future IAF-related research.
The approach also promises to open up new lines of methodological inquiry, insofar as it facilitates direct comparison of two popular indices of IAF that have for the most part been studied independently of one another.

## References

<!-- Zhivomirov, H. (2013). pinknoise [Software]. Retrieved from https://au.mathworks.com/matlabcentral/fileexchange/42919-pink--red--blue-and-violet-noise-generation-with-matlab-implementation/content/pinknoise.m (Accessed 18/03/2017) -->
