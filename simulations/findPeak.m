function [p1, p2] = findPeak(f,y,w)
% Function used to derive local maximum estimates for comparison with SGF.
% Locates the maximal power spectral estimate within predefined alpha
% band and evaluates if power exceeds that of neighbouring bins.
% This latter information distinguishes true maxima from suprema.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%%
%   Outputs:
%       p1  = bin index of maximal PSD estimate within alpha window
%       p2  = logical, are neighbouring PSD estimates < p1
%
%   Inputs:
%       f   = frequency vector (centre freqs of each freq bin)
%       y   = pxx vector (power estimates for each freq bin)
%       w   = bounds of alpha-band search window (e.g., [7, 13 Hz])
%
%%
[~, lower_alpha] = min(abs(f-w(1)));      % set lower bound for alpha band
[~, upper_alpha] = min(abs(f-w(2)));      % set upper bound for alpha band

% find local maximum
k = lower_alpha:upper_alpha;        % restrict search as above
[~, bin] = max(y(k));               % find peak pxx estimate and bin index
ind = bin + lower_alpha-1;
p1 = f(ind);

% check neighbouring bins contain lower power estimates
if y(ind) > y(ind-1) && y(ind) > y(ind+1)
    p2 = true;
else
    p2 = false;
end

end