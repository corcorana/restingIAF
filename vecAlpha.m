function y = vecAlpha(x, m, w)
% Relatively old function which performs straightforward peak detection
% routine. Searches for max power estimate within specified spectral range
% (i.e. alpha-band interval). Checks peak power value > values in 
% neighbouring bins (protects against boundary artefacts).
%
% Last modified by AC 31/01/2017.
%%
% Output:
%   y = individual peak alpha frequency estimate vector
%
% Inputs:
%   x = frequency vector
%   m = mean pxx estimate vector
%   w = alpha bounds (lower & upper frequency values)
%%

% define alpha band 
[~, lower_alpha] = min(abs(x-w(1)));      % set lower bound for alpha band
[~, upper_alpha] = min(abs(x-w(2)));      % set upper bound for alpha band
k = lower_alpha:upper_alpha;                % restrict search to alpha band freq bins


[peak, bin] = max(m(k));             % find max (mean) pxx estimate (and corresponding bin index) within alpha band
f_bin = bin + (lower_alpha-1);              % correct indexing for full frequency band

if peak > m(f_bin-1)      % max power must be > that of preceding bin (prevent arbitrary max values at lower alpha boundary)
    y = x(f_bin);   
else
    y = NaN;
end