function [selP, qWt, avPeak, sdPeak, cog, avSpec, rmseSpec, avSmoo ] = meanChansSG(chanCogs, selG, peaks, qf, psd, i, j, nchan, cmin)
% Adaptation of chanMeans.m for PSD estimates post S-G filtering.
% Also provides information pertaining to number of viable channels which
% contribute to mean estimates. Produces averages of both peak point 
% estimates (with associated SDs) and full spectra (plus RMSE for non-
% smoothed spectra). Note, spectral means were included for comparison of
% estimation procedures (mean peak vs. unsmoothed spec vs. smoothed spec);
% currently only mean peak (avPeak) is taken as the PAF estimate.
%
% Last modified by AC 31/01/2017.
%%
% Outputs:
%   selP = number of channels contributing peak estimates to mean PAF
%   qWt = weight assigned to channel spectra on basis of peak quality (Qf)
%   avPeak = mean peak estimate across channels
%   sdPeak = std dev of mean peak estimate
%   cog = mean gravity across channels
%   avSpec = mean of channel spectra (unsmoothed)
%   rmseSpec = root mean squared error of avSpec
%   avSmoo = mean of channel spectra (smoothed)
%
% Inputs:
%   chanCogs = gravity estimates for each channel
%   selG = channels selected for mean gravity estimation
%   peaks = channel-wise PAF estimates
%   qf = peak quality (area bounded by inflection points)
%   psd = psd data structure
%   i = subj number
%   j = cond number
%   nchan = number of channels
%   cmin = min number of channels required to generate cross-channel mean
%%

% channel selection and weights
selP = isnan(peaks);            % evaluate whether channel provides estimate of PAF
selP = 1- selP;                 % invert such that 1 = channel available for analysis

qWt = nansum(qf)/sum(selP);     % average area under peak (Qf) across viable channels (use later for cross-recording comparisons)
chanWts = qf/max(qf);       	% channel weightings scaled in proportion to Qf value of channel with highest Qf

% average across peaks
if sum(selP) < cmin         	% if number of viable channels < minimal threshold, don't calculate cross-channel mean & sd
    avPeak = NaN;
    sdPeak = NaN;
else
    avPeak = nansum(bsxfun(@times, peaks, chanWts))/nansum(chanWts);        % compute (weighted) cross-channel average peak alpha frequency estimate
    sdPeak = nanstd(peaks);                                                 % calculate std dev of mean peak estimate (i.e. variance across channels within recording)          
end

% do the gravs

if sum(selG) < cmin
    cog = NaN;
else
    cog = nanmean(chanCogs, 2);      % as availability of Qf will vary depending on whether PAF was identified, do not use as a weighting factor here 
end

% average across spectra
select = zeros(length(psd(i,j,nchan).pxx), nchan, 2);       % f bin x channel x orig(1) vs. smooth (2)

% define alpha band - ?why ?do IAF band 
%[~, lower_alpha] = min(abs(f-w(1)));      % set lower bound for alpha band
%[~, upper_alpha] = min(abs(f-w(2)));      % set upper bound for alpha band
%k = lower_alpha:upper_alpha;                % restrict search to alpha band freq bins

for c = 1:nchan
    if selP(c) == 1
        select(:,c,1) = psd(i,j,c).pxx;
        select(:,c,2) = psd(i,j,c).d0;
    else
        select(:,c,1) = NaN;
        select(:,c,2) = NaN;
    end
end


% average across channels to estimate mean PSD for each recording,
% calculate sum of squared errors for alpha range estimates of PSD 
if sum(selP) < cmin
    avSpec = NaN;
    avSmoo = NaN;
    rmseSpec = NaN;    
else
    avSpec = nanmean(select(:,:,1), 2);         % NB: not weighted by Q values, as these are determined from smoothed spectra
    rmseSpec = sqrt(nansum(nansum(bsxfun(@minus, avSpec, select(:, :, 1)).^2)) /sum(selP));      % use mse (rather than sse) to account for number of included channels
    
    avSmoo = nanmean(select(:,:,2), 2);     
    %wtSmoo = bsxfun(@times, select(:,:,2), chanWts);
    %avWtSmoo = nanmean(wtSmoo, 2);
    %rmseSmoo = sqrt(nansum(nansum(bsxfun(@minus, select(k, :, 2), avSmoo(k,:)).^2)) /sum(sel));
    %rmseSmooWt = sqrt(nansum(nansum(bsxfun(@minus, avWtSmoo(k), wtSmoo(k,:)).^2)) /sum(sel));
    
    % mseSmoo is usually < mseSpec (thanks to decreased variation of
    % smoothed signal) & > mseSmooWt (thanks to disparities of mean distance introduced by unequal weights)
end


end
