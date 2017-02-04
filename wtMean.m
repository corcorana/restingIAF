function [y1, y2, y3] = wtMean(n, psd, nchan, cogs)
%% Estimate mean PSD over pre- & post-task resting-state recordings
% Weight PSDs according to proportion of included channels - extension of
% wtMeanPower to distinguish unsmoothed/smoothed spectra
% Also compute mean Centre of Gravity using same cross recording criteria
%%
% Output:
%   y1 = vector of mean PSD estimates (unsmoothed)
%   y2 = vector of mean PSD estimates (smoothed)
%   y3 = mean centre of gravity
%
% Inputs:
%   n = subject number
%   psd = structure containing PSD estimate (averaged across channels)
%   nchan = EEG struct containing number of channels included in analysis
%   cogs = average centre of gravity across channels for each condition
%%


if length(psd(n, 1).avSpec) == length(psd(n, 2).avSpec)
 	wtP1 = nansum(psd(n, 1).selP)/nchan;
    wtP2 = nansum(psd(n, 2).selP)/nchan;

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

% average cogs (wt by number of chans)

wtG1 = nansum(psd(n, 1).selG)/nchan;
wtG2 = nansum(psd(n, 2).selG)/nchan;

if sum([cogs(1, 1), cogs(1, 2)]) > 0    % ideally, should evaluate that both are not NaN
    y3 = sum([cogs(1, 1) * wtG1, cogs(1, 2) * wtG2])/(wtG1 + wtG2);
elseif isnan(cogs(1, 1))
    y3 = cogs(1, 2);
else y3 = cogs(1, 1);

end

