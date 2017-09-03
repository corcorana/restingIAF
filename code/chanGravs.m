function [cogs, sel, iaw] = chanGravs(d0, f, f1, f2)
% Takes smoothed channel spectra and associated estimates of individual 
% alpha bandwidth [f1:f2], calculate mean bandwidth, estimate CoG across
% all channels (as per Klimesch's group; e.g, 1990, 1993, & 1997 papers).
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%% Outputs:
%   cogs = centre of gravity derived from averaged f1:f2 frequency window
%   sel = channels contributing estimates of alpha window bandwidth
%   iaw = bounds of individual alpha window
%
%% Required inputs:
%   d0 = vector / matrix of smoothed PSDs
%   f = frequency bin vector
%   f1 = vector of f1 bin indices (lower bound, individual alpha window)
%   f2 = vector of f2 bin indices (upper bound, individual alpha window)
 
%% setup inputParser
p = inputParser;
p.addRequired('d0',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'2d'}));
p.addRequired('f',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector', 'nonnegative', 'increasing'}));
p.addRequired('f1',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.addRequired('f2',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'vector'}));
p.parse(d0, f, f1, f2)
%%

% trim off any NaNs in f1/f2 vectors
trim_f1 = f1(~isnan(f1));
trim_f2 = f2(~isnan(f2));

% derive average f1 & f2 values across chans, then look for nearest freq bin
mean_f1 = dsearchn(f, mean(f(trim_f1)));
mean_f2 = dsearchn(f, mean(f(trim_f2)));
iaw = [mean_f1, mean_f2];
        
% calculate CoG for each channel spectra on basis of averaged alpha window
cogs = zeros(1,size(d0,2));
for d = 1:length(cogs)
    if isempty(trim_f1) || isempty(trim_f2)
        cogs(d) = NaN;
    else
        cogs(d) = nansum(d0(mean_f1:mean_f2,d).*f(mean_f1:mean_f2)) / sum(d0(mean_f1:mean_f2,d));
    end
end

% report which channels contribute to averaged window
sel = ~isnan(f1);

end