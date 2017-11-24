function plotAvSpec(f, pSpec, ix, log)
% Simple function for plotting Q-weighted spectra from two sets of resting-
% state recordings. Recording must satisfy cmin for PAF estimation in order 
% to generate plot.
%
% Automatically selects all channels for plotting on a 2-25 Hz scale (can
% be modified from command line).
%
% This function and all custom-designed dependencies are part of the 
% `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
%
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%% Outputs:
%   figure displaying smoothed spectral estimates for analysed channels
%
%% Inputs:
%   f = frequency vector (for x-axis scaling)
%   pSpec = structure containing channel estimate data
%   ix = participant number (for indexing structure)
%   log = logical; plot linear (0) or log (1) scale

%%

figure
if log == 1         % log-scaled plots
    plot(f, 10*log10(pSpec(ix,1).sums.muSpec), 'Color', [0 .3 .7], 'LineWidth',2)
    hold on
    plot(f, 10*log10(pSpec(ix,2).sums.muSpec),  'Color', [1 .2 0], 'LineWidth',2)
elseif log == 0     % linear-scaled plots
    plot(f, pSpec(ix,1).sums.muSpec, 'Color', [0 .3 .7], 'LineWidth',2)
    hold on
    plot(f, pSpec(ix,2).sums.muSpec,  'Color', [1 .2 0], 'LineWidth',2)
end

xlim([2 25])
xlabel('Frequency (Hz)')
ylabel('Power (normalised)')
set(gca, 'FontSize', 14)