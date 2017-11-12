function plotSpec(f, pSpec, ix, jx, psd, log)
% Simple function for plotting PSDs from pSpec data structure.
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
%   figure displaying spectral estimates for analysed channels
%
%% Inputs:
%   f = frequency vector (for x-axis scaling)
%   pSpec = structure containing channel estimate data
%   ix = participant number (for indexing structure)
%   jx = recording number (for indexing structure)
%   psd = char; plot unsmoothed ('pxx') or smoothed ('d0') PSDs
%   log = logical; plot linear (0) or log (1) scale

%%
nchan = length(pSpec(ix, jx).chans);    % calculate number of included channels

figure
if log == 1         % log-scaled plots
    for kx = 1:nchan
        plot(f, 10*log10(eval(sprintf(['pSpec(' num2str(ix) ',' num2str(jx) ').chans(' num2str(kx) ').' psd]))))
        hold on
    end
elseif log == 0     % linear-scaled plots
    for kx = 1:nchan
        plot(f, eval(sprintf(['pSpec(' num2str(ix) ',' num2str(jx) ').chans(' num2str(kx) ').' psd])))
        hold on
    end
end

xlim([2 25])
xlabel('Frequency (Hz)')
ylabel('Power (normalised)')
set(gca, 'FontSize', 14)
