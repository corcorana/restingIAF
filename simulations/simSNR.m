function [y, rdm] = simSNR(ts, frex, snr, n)
% Create simulation data for PAF estimation under varied SNR conditions.
% Combines two vectors (target component signal and pink noise signal).
% Depends on `pinknoise` function to generate 1/f distributed noise signal.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%%
% Outputs:
%   y       = matrices of time series data (rows) for each simulation
%               - 1st dim = alpha signal centred at a randomly selected freq
%               - 2nd dim = randomly generated pink noise signal
%               - 3rd dim = combined signal (multiplied in time domain)
%   rdm     = vector of randomly sampled alpha freqs
%
% Inputs:
%   ts      = vector of time sample points
%   frex    = vector of alpha frequencies from which alpha signals sampled
%   snr     = desired signal-to-noise ratio (0 to 1)
%   n       = number of simulated signals
%%

% calculate number of alpha signal samples required for desired SNR
a_pnts = 1:round(length(ts)*snr);
unity = ones(1, length(ts)-length(a_pnts));     % remaining samples will be set to 1

% randomly sample frex vector (with replacement)
rng(812, 'twister');                            % rng set to replicate original analysis
rdm = frex(randsample(1:length(frex), n, true));      

% create mean centred / unity normalised alpha, pink, and combined signals
y = zeros(n, length(ts), 3);
for Fc = 1:n
    sig = sin(2*pi*rdm(Fc).*ts);
    sig = sig - mean(sig);
    rms = sqrt(mean(sig.^2));
    y(Fc,:,1) = sig/rms;
    y(Fc,:,2) = pinknoise(length(ts));
    y(Fc,:,3) = [y(Fc,a_pnts,1), unity].*y(Fc,:,2);
end


end