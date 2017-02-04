function [cogs, sel] = gravity(d0, f, f1, f2)
% Take smoothed spectra and associated estimates of individual alpha 
% bandwidth f1 - f2, calculate mean bandwidth, estimate CoG across all
% channels (as per Klimesch's group, e.g. 1990, 1993, 1997 papers).
%
% Developed from chanCogs, which only estimate CoG for channels rendering
% clear f1/f2 estimates.
%
% Last modified AC, 31/01/2017
%%
% Outputs:
%   cogs = centre of gravity derived from averaged f1/f2 frequency window
%   sel = number of channels contributing estimates of f1/f2
%
% Inputs:   
%   d0 = matrix of smoothed PSDs
%   f = freq vector
%   f1 = vector of channel f1s
%   f2 = vector of channel f2s
%%

% trim off NaNs
trim_f1 = f1(~isnan(f1));
trim_f2 = f2(~isnan(f2));

% derive average frequency value, then look for nearest bin centre freq
mean_f1 = dsearchn(f, mean(f(trim_f1)));
mean_f2 = dsearchn(f, mean(f(trim_f2)));

if isempty(trim_f1) || isempty(trim_f2)
    cogs = NaN;
else
    cogs = zeros(1,size(d0,2));
    for d = 1:length(cogs)
        cogs(d) = nansum(d0(mean_f1:mean_f2,d).*f(mean_f1:mean_f2)) / sum(d0(mean_f1:mean_f2,d));
    end
end

if isempty(trim_f1)
    sel = 0;
else sel = length(trim_f1);      % number of selected f1 values used to esitmate IAF window
end

end