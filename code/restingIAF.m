function [pSum, pChans, f] = restingIAF(data, nchan, cmin, fRange, Fs, varargin)
% Primary function for running `restingIAF` analysis routine for estimating
% two indices of individual alpha frequency (IAF): Peak alpha frequency 
% (PAF) and the alpha mean frequency / centre of gravity (CoG).
%
% Calls on the Signal Processing Toolbox function `pwelch` to derive power 
% spectral density estimates of one or more vectors of EEG channel data, 
% which are subsequently smoothed by `sgfDiff` (depends on `sgolay`). 
% Alpha peak activity is parameterised by `peakBounds` and CoG estimated by 
% `chanGravs`. Channel-wise peak and gravity estimates are averaged across 
% channels by `chanMeans`.
%
% This function and all custom-designed dependencies are part of the 
% `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
%
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%% Outputs:
%   pSum = structure containing summary statistics
%   pChans = structure containing channel-wise spectral and alpha parameter data
%   f = trimmed vector of frequency bins resolved by `pwelch`
%
%% Required inputs:
%   data = vector or matrix containing continuous EEG channel data
%   nchan = number of channels in data
%   cmin = minimum number of channel estimtes that must be resolved in
%          order to calculate average PAF/CoG estimates
%   fRange = frequency range to be included in analysis (e.g., 1-40 Hz )
%   Fs = EEG sampling rate
%
%% Optional inputs:
%   w = bounds of alpha peak search window (default = [7 13])
%   Fw = frame width, Savitzky-Golay filter (default = 11, must be odd)
%   k = polynomial order, Savitzky-Golay filter (default = 5, must be < Fw)
%   mpow = error bound (standard deviation) used to determine threshold differentiating 
%          substantive peaks from background spectral noise (default = 1)
%   mdiff = minimal height difference distinguishing a primary peak from
%           competing peaks (default = 0.20)
%   tlen = length of taper window applied by `pwelch` (default = 4*Fs)
%   tover = length of taper window overlap in samples (default = 50% window length)
%   nfft = specify number of FFT points used to calculate PSD (default = next power of 2 above window length)

%% setup inputParser
p = inputParser;
p.addRequired('data',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'2d', 'nonempty'}));
p.addRequired('nchan',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
p.addRequired('cmin',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
p.addRequired('fRange',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'integer', 'nonnegative', 'increasing', 'size', [1,2]}));
p.addRequired('Fs',...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive', '>=', 2*fRange(2)}));

p.addOptional('w', [7 13],...
                @(x) validateattributes(x, {'numeric'}, ...
                {'nonnegative', 'increasing', 'size', [1,2], '>', fRange(1), '<', fRange(2)}));
p.addOptional('Fw', 11,...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive', 'odd'}));
p.addOptional('k', 5,...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive', '<', varargin{2} }));
p.addOptional('mpow', 1,...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'positive'}));
p.addOptional('mdiff', .20,...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', '>=', 0, '<=', 1}));
p.addOptional('tlen', (Fs*4),...
                @(x) validateattributes(x, {'numeric'}, ...
                {'scalar', 'integer', 'positive'}));
p.addOptional('tover', [],...
                @(x) validateattributes(x, {'numeric'}, ...
                {'integer', 'positive'}));
p.addOptional('nfft', [],...
                @(x) validateattributes(x, {'numeric'}, ...
                {'integer', 'positive'}));           
p.parse(data, nchan, cmin, fRange, Fs, varargin{:})

w = p.Results.w;
Fw = p.Results.Fw;
k = p.Results.k;
mpow = p.Results.mpow;
mdiff = p.Results.mdiff;
tlen = p.Results.tlen;
tover = p.Results.tover;
nfft = p.Results.nfft;
%%

pChans = struct('pxx', [], 'minPow', [], 'd0', [], 'd1', [], 'd2', [], 'peaks', [], 'pos1', [], 'pos2', [], 'f1', [], 'f2', [], 'inf1', [], 'inf2', [], 'Q', [],'Qf', [], 'gravs', [], 'selP', [], 'selG', [] );      % struct for channel data (PSD estimates & derivatives, some additional info)

for kx = 1:nchan
	% perform pwelch routine to extract PSD estimates by channel
    [pxx, f] = pwelch(data(kx,:), hamming(tlen), tover, nfft, Fs, 'psd');

  	% delimit range of freq bins to be included in analysis
 	frex = dsearchn(f, fRange(1)):dsearchn(f, fRange(2));      
   	f = f(frex);
	pxx = pxx(frex);      % truncate PSD to frex range

   	% normalise truncated PSD
   	pChans(kx).pxx = pxx / mean(pxx);            
            
  	% calculate minPower vector
  	[pfit, sig] = polyfit(f, log10(pChans(kx).pxx), 1);     % fit 1st order poly (regression line) to normalised spectra (log-scaled)
    [yval, del] = polyval(pfit, f, sig);                    % derive yval coefficients of fitted polynomial and delta (std dev) error estimate
  	pChans(kx).minPow = yval + (mpow * del);                % takes [minPowThresh * Std dev] as upper error bound on background spectral noise
                
  	% apply Savitzky-Golay filter to fit curves to spectra & estimate 1st and 2nd derivatives
    [pChans(kx).d0, pChans(kx).d1, pChans(kx).d2] = sgfDiff(pChans(kx).pxx, Fw, k, Fs, tlen);

  	% calculate frequency resolution
  	exp_tlen = nextpow2(tlen);
  	fres = Fs/2.^exp_tlen;
        
   	% take derivatives, find peak(s) and boundaries of alpha band
   	[pChans(kx).peaks, pChans(kx).pos1, pChans(kx).pos2, pChans(kx).f1, pChans(kx).f2, pChans(kx).inf1, pChans(kx).inf2, pChans(kx).Q, pChans(kx).Qf] = peakBounds(pChans(kx).d0, pChans(kx).d1, pChans(kx).d2, f, w, pChans(kx).minPow, mdiff, fres);

end
        
% estimate gravities for smoothed spectra (average IAF window across channels)
[ gravs, selG ] = chanGravs([pChans(:).d0], f, [pChans(:).f1], [pChans(:).f2] );

% calculate average pt estimates/spectra across k-th channels for each j-th recording
[ selP, pSum ] = chanMeans(gravs, selG, [pChans(:).peaks], [pChans(:).Qf], cmin); 

% retain gravity estimates and selected channels in channel data struct
for ix = 1:nchan
    pChans(ix).gravs = gravs(ix);
    pChans(ix).selP = selP(ix);
    pChans(ix).selG = selG(ix); 
end
    
% get total number of chans that contributed to PAF/CoG estimation
pSum.pSel = sum(selP);
pSum.gSel = sum(selG);


end