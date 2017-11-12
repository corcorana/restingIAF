function [sig, gw] = simBimod(ts, frex, recip, prop, snr, nchans, n)
% Adapted from simParams. 
% The gw function used to construct the alpha band signal has been adapted
% to sample a bimodal gaussian distribution. This function is output for
% purposes of figure plotting.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%%
%   Outputs:
%       sig     = structure comprising alpha signal & combined alpha/pink 
%                 noise signals (i.e. simulated channel data)
%       gw      = gaussian window function
%
%   Inputs:
%       ts      = time series sample points
%       frex    = selected range of alpha frequencies
%       recip   = reciprocal of the standard deviation (alpha) of Gaussian window (higher = narrower distribution)
%       prop    = control dimensions of bimodal window
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
% set up bimodal window comprised of 2 Gaussians
gw = gausswin(35, recip);
bm = [gw(1:25);gw(10:35)];
% define relative proportion of area
if prop == 1                % equal gaussians (lower freq biased when combined with pink noise)
    bmprop = bm/sum(bm);

elseif prop == 2            % increase height of upper freq by 25%
    bmprop = [bm(1:25)/sum(bm); bm(26:end)/(sum(bm)/1.25)];
    bmprop = bmprop/sum(bmprop);    

elseif prop == 3            % increase height of upper freq by 50%
    bmprop = [bm(1:25)/sum(bm); bm(26:end)/(sum(bm)/1.5)];
    bmprop = bmprop/sum(bmprop);
end

% calculate alpha, pink noise, and combined signals
sig = struct('a',[],'y',[],'t',[]);
for Fc = 1:n
    k = dsearchn(fran', rdm(Fc));       % find index of centre frequency for bm
    f1 = k - floor(length(bm)/2);       % find index of first frequency for bm
    f2 = k + floor(length(bm)/2);       % find index of last frequency for bm
    fs = f1:f2;                         % vector of frequency row indices corresponding to window sample points
    tmp = struct('s', []);              % temporary struct holding varible length series across sampled frex
    for g = 1:length(bm)
        tmp(g).s = a_mat(fs(g), 1:round(a_pnts*bmprop(g)));     % for each frequency vector included by the window, sample proportion of signal in relation to window function
    end
    sig(Fc,:).a = [horzcat(tmp(:).s), unity];       % output alpha signal

    x = zeros(nchans, length(sig(Fc,:).a));
    for chans = 1:nchans
        p = pinknoise(length(sig(Fc,:).a));
        x(chans, :) = sig(Fc,:).a.*p;
    end
    sig(Fc,:,:).y = x;
    sig(Fc).t = rdm(Fc);                % output target (i.e. window centre freq)
    
end

end
