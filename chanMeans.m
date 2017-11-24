function [selP, sums] = chanMeans(chanCogs, selG, peaks, specs, qf, cmin)
% Takes channel-wise estimates of peak alpha frequency (PAF) / centre of
% gravity (CoG) and calculates mean and standard deviation if cmin 
% condition satisfied. PAFs are weighted in accordance with qf, which aims
% to quantify the relative strength of each channel peak. 
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%% Outputs:
%   selP = channels contributing peak estimates to calculation of mean PAF
%   sums = structure containing summary estimates (m, std) for PAF and CoG
%
%% Required inputs:
%   chanCogs = vector of channel-wise CoG estimates
%   selG = vector (logical) of channels selected for individual alpha band estimation
%   peaks = vector of channel-wise PAF estimates
%   specs = matrix of smoothed spectra (helpful for plotting weighted estimates)
%   qf = vector of peak quality estimates (area bounded by inflection points)
%   cmin = min number of channels required to generate cross-channel mean

%% setup inputParser
p = inputParser;
p.addRequired('chanCogs',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector', 'positive'}));
p.addRequired('selG',...
                @(x) validateattributes(x, {'logical'}, ...
                {'vector'}));
p.addRequired('peaks',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector', 'positive'}));
p.addRequired('qf',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector', 'positive'}));
p.addRequired('cmin',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));          
p.parse(chanCogs, selG, peaks, qf, cmin)
%%

% channel selection and weights
selP = ~isnan(peaks);               % evaluate whether channel provides estimate of PAF

% qWt = nansum(qf)/sum(selP);       % average area under peak (Qf) across viable channels (depricated: was used when calculating cross-recording comparisons)
chanWts = qf/max(qf);               % channel weightings scaled in proportion to Qf value of channel manifesting highest Qf

% average across peaks
if sum(selP) < cmin                 % if number of viable channels < cmin threshold, don't calculate cross-channel mean & std
    sums.paf = NaN;
    sums.pafStd = NaN;
    sums.muSpec = NaN;
else                                % else compute (weighted) cross-channel average PAF estimate and corresponding std of channel PAFs
    sums.paf = nansum(bsxfun(@times, peaks, chanWts))/nansum(chanWts);
    sums.pafStd = nanstd(peaks);
    % estimate averaged spectra for plotting
    wtSpec = bsxfun(@times, specs, chanWts);
    sums.muSpec = nansum(wtSpec, 2)/nansum(chanWts);
end

% now for the gravs (no channel weighting, all channels included if cmin satisfied)
if sum(selG) < cmin
    sums.cog = NaN;
    sums.cogStd = NaN;
else
    sums.cog = nanmean(chanCogs, 2);
    sums.cogStd = nanstd(chanCogs);
end




end
