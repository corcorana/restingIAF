function sig = simParams(ts, frex, recip, snr, nchans, n)
% Simulate EEG signals comprised of a range of alpha frequency components
% and background pink noise distribution. Signals can be resampled to
% generate 'datasets' in which the same alpha band distribution is sampled
% across a number of channels (with random variation of background noise).
%
% Alpha signal constructed by concatenating segments of several sine wave
% signals, each of which oscillates at a unique frequency. The 
% proportion of each segment included in the composite series is determined
% by a weighted Gaussian function. The std dev of this function can be
% parametrically varied to give narrower or broader alpha distributions.
%
% The mixed alpha signal is multiplied by a randomly generated pink noise
% distribution to render the final time series. Signal-to-noise ratio (SNR)
% between alpha and noise can be varied between 0 (no alpha signal) and 1
% (equal contribution of alpha/noise signal).
%
% Depends on `pinknoise` function to generate 1/f distributed noise signal.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%%
%   Outputs:
%       sig     = structure comprising alpha signal & combined alpha/pink 
%                 noise signals (i.e. simulated channel data)
%
%   Inputs:
%       ts      = time series sample points
%       frex    = selected range of alpha frequencies
%       recip   = reciprocal of the standard deviation (alpha) of Gaussian window (higher = narrower distribution)
%       snr     = signal-to-noise ratio (alpha vs. pink noise signal)
%       nchans  = number of channels to simulate for each subject (randomly sampled freq)
%       n       = number of 'subject datasets' to simulate
%%

% calculate (extended) alpha component matrix
bound = (frex(2)-frex(1))*ceil(length(frex)/2);
fran = frex(1)-bound:frex(2)-frex(1):frex(end)+bound;

a_mat = zeros(length(fran),length(ts));
for Fc = 1:length(fran)
    alpha = sin(2*pi*fran(Fc).*ts);
    alpha = alpha - mean(alpha);
    rms = sqrt(mean(alpha.^2));
    a_mat(Fc,:) = alpha/rms;
end

% randomly sample frex vector (with replacement)
rng(713, 'twister');                    % rng set to replicate original analysis
rdm = frex(randsample(1:length(frex), n, true));   

% calculate number of alpha signal samples required for desired SNR
a_pnts = round(length(ts)*snr);
unity = ones(1, length(ts)-a_pnts);     % remaining samples will be set to 1

% sample from a_mat to build composite alpha signal
gw = gausswin(length(frex), recip);
gprop = gw/sum(gw);

% calculate alpha, pink noise, and combined signals
sig = struct('a',[],'y',[],'t',[]);
for Fc = 1:n
    k = dsearchn(fran', rdm(Fc));       % find index of centre frequency for gw
    f1 = k - floor(length(gw)/2);       % find index of first frequency for gw
    f2 = k + floor(length(gw)/2);       % find index of last frequency for gw
    fs = f1:f2;                         % vector of frequency row indices corresponding to window sample points
    tmp = struct('s', []);              % temporary struct holding varible length series across sampled frex
    for g = 1:length(gw)
        tmp(g).s = a_mat(fs(g), 1:round(a_pnts*gprop(g)));     % for each frequency vector included by the window, sample proportion of signal in relation to window function
    end
    sig(Fc,:).a = [horzcat(tmp(:).s), unity];       % output alpha signal

    x = zeros(nchans, length(sig(Fc,:).a));
    for chans = 1:nchans
        p = pinknoise(length(sig(Fc,:).a));
        x(chans, :) = sig(Fc,:).a.*p;
    end
    sig(Fc,:,:).y = x;
    sig(Fc).t = rdm(Fc);                            % output target peak
    
end

end
