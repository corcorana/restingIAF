%% Savitsky-Golay smoothing & differentiation of power spectral density
% Script developed for extraction of peak alpha frequency (PAF) & alpha 
% centre of gravity (CoG) estimates from resting-state EEG. 
%
% Power spectral density (PSD) estimated using the pwelch.m routine and 
% subjected to a Savitzky-Golay filter in order to smooth and differentiate 
% the PSD function. Peaks are inferred from downward-going zero-crossings 
% in the 1st derivative of the smoothed spectrum, troughs (limits of 
% individual alpha bandwidth) from upward-going zero-crossings. 
%
% Several parameters must be defined to specify the S-G filter operation, 
% criteria for peak extraction and averaging. 
% Manuscript describing method currently in preparation - contact Andrew 
% Corcoran (andrew.corcoran1@monash.edu) for further details or assistance 
% with debugging, etc.
%
% Script currently organised to loop through 1:n participants, reading in
% eyes-closed resting-state EEG from 2 conditions (pre/post experiment).
%
% WARNING - this programme (and its associated dependencies) are still in 
% development. More output is generated than necessary, however this
% information is sometimes useful for sanity checking / comparing the
% convergence of various estimators (and hence has been retained).
%
% Uncomment options at end of script to save IAF/PSD data in local path.
%
% AC, UniSA Magill, 31/01/2017. Programme developed in MATLAB 2014b/2015a.
% Last modified by AC 26/02/2017 to ensure PAFs are weighted by number of
% channels contributing estimates across recordings. Also re-assigned gravs
% to psd struct to avoid strange problem with matrix not correctly
% initialising (leading to subscripted dimension error on first run)
%% Required functions:
%   alphaParams
%   gravity
%   meanChansSG
%   pwelch (signal processing toolbox)
%   savGolDiff
%   sgolay (signal processing toolbox)
%   vecAlpha
%   wtMean
%%

addpath('eeglab path goes here');                % require EEGlab to read in and filter EEG data
dataPath = 'path to your data goes here';

%% set file name conventions in line 96 %%

% fire up eeglab
if ~exist('EEG','var')
    eeglab;
end

% ensure double precision switched on
pop_editoptions('option_single', 0);

%% set up some initial params 
n = 18;             % number of subjects to include
w = [7 13];         % alpha band limits
poly = 5;           % polynomial order for S_G curve fitting procedure
Fw = 11;            % frame width for curve fitting procedure
minPeakDiff = .20;  % minimum proportion of peak height by which highest peak must surpass second highest peak (if >1 candidate peak identified)
t = 4;
Fs = 512;           % sampling rate (BioSemi files first)
h = Fs*4;           % length (samples) pwelch hamming window
cmin = 3;           % minimum number of channels with clear peaks required to compute cross-channel average estimate of alpha peak

freqU = 30;

% to get around unequal number of freq bands due to diff sampling rate,
% zero pad spectra to enable matrix dimension matching (add input wtMean)
numFrex = 250;

%% initialise structures
psd = struct('pxx', [], 'd0', [], 'd1', [], 'd2', [], 'selP', [], 'grav', [], 'selG', [], 'sdPeak', [], 'avSpec', [], 'rmseSpec', [], 'avSmoo', [], 'CoG', []);      % this will collect PSD estimates & derivatives
peaks = zeros(n, EEG.nbchan, 2);
pos1 = zeros(n, EEG.nbchan, 2);
pos2 = zeros(n, EEG.nbchan, 2);
f1 = zeros(n, EEG.nbchan, 2);
f2 = zeros(n, EEG.nbchan, 2);
inf1 = zeros(n, EEG.nbchan, 2);
inf2 = zeros(n, EEG.nbchan, 2);
Q = zeros(n, EEG.nbchan, 2);            % handy for scaling plots
Qf = zeros(n, EEG.nbchan, 2);
cog = zeros(n, 2);
paf = zeros(n, 2);

mu_spec = zeros(n, numFrex);                         % this will collect mean PSD values computed across recordings
mu_smoo = zeros(n, numFrex);                         % this will collect mean PSD values computed across recordings
mu_peak = zeros(n, 1);
mu_cog = zeros(n, 1); 

iaf = zeros(n, 12);                                  % this will collect peak alpha frequency estimates

%% run analysis loops
for idx_i = 1:n     % for each i-th subject
    
    for idx_j = 1:2       % for each j-th resting-state recording
        
        % load pre-processed data file
        EEG = pop_loadset('filename', sprintf(['joad_%02d_', num2str(idx_j), '.set'], idx_i), 'filepath',dataPath);
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
            
        % calculate PSD estimates for each k-th channel
        Fs = EEG.srate;                  % makes sure Fs is correct (read in from current data structure)
        for idx_k = 1:EEG.nbchan
            % perform pwelch routine to extract PSD estimates by channel
        	[psd(idx_i, idx_j, idx_k).pxx, f] = pwelch(EEG.data(idx_k,:), hamming(Fs*4), 50, [], Fs, 'psd');

            % delimit range of freq bins to be included in analysis
            frex = 1:dsearchn(f, freqU);      
            f = f(frex);
            psd(idx_i, idx_j, idx_k).pxx = psd(idx_i, idx_j, idx_k).pxx(frex);      % truncate PSD to frex range
            
            % normalise truncated PSD
            psd(idx_i, idx_j, idx_k).pxx = psd(idx_i, idx_j, idx_k).pxx / mean(psd(idx_i, idx_j, idx_k).pxx);            
            
            % calculate minPower vector
            [p, s] = polyfit(f, psd(idx_i, idx_j, idx_k).pxx, 1);           % fit 1st order poly (regression line) to normalised spectra
            [y, delta] = polyval(p, f, s);                                  % derive y coefficients of fitted polynomial and delta (std dev) error estimate
            minPow = y+delta;                      % takes upper error bound as threshold
            
            % consider a robust alternative to above
            %[b, stats] = robustfit(f, psd(idx_i, idx_j, idx_k).pxx);
            %psd(idx_i, idx_j, idx_k).minPow = (b(1)+b(2)*f)+2*stats.s;
            
            % apply Savitzky-Golay filter to fit curves to spectra & estimate 1st and 2nd derivatives
            [psd(idx_i, idx_j, idx_k).d0, psd(idx_i, idx_j, idx_k).d1, psd(idx_i, idx_j, idx_k).d2] = savGolDiff(psd(idx_i, idx_j, idx_k).pxx, poly, Fw, Fs, Fs*4);
            
            % take derivatives, find peak(s) and boundaries of alpha band
            [peaks(idx_i, idx_k, idx_j), pos1(idx_i, idx_k, idx_j), pos2(idx_i, idx_k, idx_j), f1(idx_i, idx_k, idx_j), f2(idx_i, idx_k, idx_j), inf1(idx_i, idx_k, idx_j), inf2(idx_i, idx_k, idx_j), Q(idx_i, idx_k, idx_j), Qf(idx_i, idx_k, idx_j)] = alphaParams(psd(idx_i, idx_j, idx_k).d0, psd(idx_i, idx_j, idx_k).d1, psd(idx_i, idx_j, idx_k).d2, f, w, minPow, minPeakDiff, t);
            
        end
        
        % estimate gravities for smoothed spectra (average IAF window across channels)
        [ psd(idx_i, idx_j).grav, psd(idx_i, idx_j).selG ] = gravity([psd(idx_i, idx_j, :).d0], f, f1(idx_i, :, idx_j), f2(idx_i, :, idx_j));
                          
        % calculate average pt estimates/spectra across k-th channels for each j-th recording
        [ psd(idx_i, idx_j).selP, psd(idx_i, idx_j).qWt, psd(idx_i, idx_j).sdPeak, paf(idx_i, idx_j), cog(idx_i, idx_j), psd(idx_i, idx_j).avSpec, psd(idx_i, idx_j).rmseSpec, psd(idx_i, idx_j).avSmoo ] = meanChansSG(psd(idx_i, idx_j).grav, psd(idx_i, idx_j).selG, peaks(idx_i, :, idx_j), Qf(idx_i, :, idx_j), psd, idx_i, idx_j, EEG.nbchan, cmin); 
                    
    end         % end j-th recording loop

    % average mean PSD estimates across (j-th) recordings
    [mspec, msmoo, mu_peak(idx_i,:), mu_cog(idx_i, :)] = wtMean(idx_i, psd, EEG.nbchan, paf(idx_i,:), cog(idx_i,:));
    
    pad = zeros(1, numFrex - length(mspec));        % for zero-padding
    mu_spec(idx_i, :) = [mspec' pad];
    mu_smoo(idx_i, :) = [msmoo' pad];
    
    % estimate individual peak alpha frequency
    iaf(idx_i,1) = vecAlpha(f, mu_spec(idx_i,:), w);
    iaf(idx_i,2) = psd(idx_i, 1).rmseSpec;
    iaf(idx_i,3) = psd(idx_i, 2).rmseSpec;
    iaf(idx_i,4) = vecAlpha(f, mu_smoo(idx_i,:), w);
    iaf(idx_i,5) = sum(psd(idx_i, 1).selP);
    iaf(idx_i,6) = sum(psd(idx_i, 2).selP);
    iaf(idx_i,7) = mu_peak(idx_i);
    iaf(idx_i,8) = psd(idx_i, 1).sdPeak;
    iaf(idx_i,9) = psd(idx_i, 2).sdPeak;
    iaf(idx_i,10) = mu_cog(idx_i);
    iaf(idx_i,11) = sum(psd(idx_i, 1).selG);
    iaf(idx_i,12) = sum(psd(idx_i, 2).selG);

end     % end i-th subject loop

%% print summary of IAF estimates to console
S_num = (1:n)';
colnames = { 'S', 'nChans1', 'nChans2', 'muPeaks', 'sd1', 'sd2', 'CoG', 'nGravs1', 'nGravs2' };
IAF = array2table([S_num, iaf(:, 5:12)], 'VariableNames', colnames)
fprintf('Estimates derived via pwelch method, initial alpha window = %.1f - %.1f Hz.', w(1), w(2));

%% save output data
% save('psd_struct.mat', 'psd')
% save('iaf.mat', 'IAF')