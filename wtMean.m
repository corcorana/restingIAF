function [y1, y2, y3, y4] = wtMean(n, psd, nchan, peaks, cogs)
%% Estimate mean PSD over pre- & post-task resting-state recordings
% Weight PSDs according to proportion of included channels - extension of
% wtMeanPower to distinguish unsmoothed/smoothed spectra
% Also compute mean Centre of Gravity using same cross recording criteria
%
% Last modified by AC 26/02/2017 - add mean paf output to ensure grand
% mean PAF (averaged across recordings) is weighted by nchans. Also change
% conditional for CoG estimate to ensure NaNs properly handled.
%%
% Output:
%   y1 = vector of mean PSD estimates (unsmoothed)
%   y2 = vector of mean PSD estimates (smoothed)
%   y3 = mean peak alpha frequency
%   y4 = mean centre of gravity
%
% Inputs:
%   n = subject number
%   psd = structure containing PSD estimate (averaged across channels)
%   nchan = EEG struct containing number of channels included in analysis
%   peaks = mean peak across channels (from Q-weighted point estimates)
%   cogs = mean centre of gravity across channels for each condition
%%

wtP1 = nansum(psd(n, 1).selP)/nchan;
wtP2 = nansum(psd(n, 2).selP)/nchan;
    
if length(psd(n, 1).avSpec) == length(psd(n, 2).avSpec)

    totalQ = psd(n, 1).qWt + psd(n, 2).qWt;
    qWt1 = psd(n, 1).qWt/totalQ;
    qWt2 = psd(n, 2).qWt/totalQ;
    
    y1 = mean([psd(n, 1).avSpec .* wtP1, psd(n, 2).avSpec .* wtP2], 2);                   % NB: only weighted by number of selected channels (not Q)
    y2 = mean([psd(n, 1).avSmoo .* (wtP1+qWt1), psd(n, 2).avSmoo .* (wtP2 + qWt2)], 2);     % weighted by number of selected channels plus Q value

elseif isnan(psd(n, 1).avSpec)
	y1 = psd(n, 2).avSpec;
    y2 = psd(n, 2).avSmoo;
else
 	y1 = psd(n, 1).avSpec;
    y2 = psd(n, 1).avSmoo;
end

% average peaks (wt by number of chans)
if ~isnan(peaks(1, 1)) && ~isnan(peaks(1, 2))    % assumes comparison of 2 recordings
    y3 = sum([peaks(1, 1) * wtP1, peaks(1, 2) * wtP2])/(wtP1 + wtP2);
elseif isnan(peaks(1, 1)) && ~isnan(peaks(1,2))
    y3 = peaks(1, 2);
elseif isnan(peaks(1, 2)) && ~isnan(peaks(1,1))
    y3 = peaks(1, 1);
else y3 = NaN;
end

% average cogs (wt by number of chans)

wtG1 = nansum(psd(n, 1).selG)/nchan;
wtG2 = nansum(psd(n, 2).selG)/nchan;

if ~isnan(cogs(1, 1)) && ~isnan(cogs(1, 2))    % ? more efficient way of evaluating this
    y4 = sum([cogs(1, 1) * wtG1, cogs(1, 2) * wtG2])/(wtG1 + wtG2);
elseif isnan(cogs(1, 1)) && ~isnan(cogs(1,2))
    y4 = cogs(1, 2);
elseif isnan(cogs(1, 2)) && ~isnan(cogs(1,1))
    y4 = cogs(1, 1);
else y4 = NaN;
end

