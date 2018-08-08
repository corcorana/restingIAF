%% `restingIAF` tutorial script
% This script is designed to provide a brief example of how the various
% functions belonging to the `restingIAF` package can be combined within
% EEGLAB to enable replicable analysis of multiple data files.
%
% The script analyses pre- and post-experiment eyes-closed resting-state
% data from 3 subjects, and outputs the grand-averaged PAF and CoG
% estimates and associated summary data (std dev of channel-wise estimates,
% number of channels included in analysis) to the console.
% 
% The data files called by this script have been preprocessed to retain 6
% parietal channels of interest and reduce recordings to 2 min duration. 
% They have been bandpass filtered (1 - 40 Hz) and downsampled to 250 Hz.
%
% For more information, including access to the original (i.e. minimally-
% processed) data files, visit the figshare repository for this dataset:
%
%   https://figshare.com/articles/Muspelheim_data/3412312
% 
% Use the following key to trace data files back to their original version:
%
%   tute_01_1 == ali0039_azv
%   tute_01_2 == ali0039_azn
%   tute_02_1 == ali0038_azv
%   tute_02_2 == ali0038_azn
%   tute_03_1 == ali0037_azv
%   tute_03_2 == ali0037_azn
%
% This script (and all custom-designed dependencies) are part of the 
% `restingIAF` package, (c) Andrew W. Corcoran, 2016-2018.
%
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.

%% Preliminary setup
 
% this command ensures restingIAF functions and tutorial files are
% accessible to my path (NB: may not work with all operating systems)
addpath(genpath('../'))    

% fire up eeglab (may also need to add to path)
if ~exist('EEG','var')
    eeglab;
end

% ensure double precision is switched on
pop_editoptions('option_single', 0);

% define path to the data files -- useful if files not stored in current
% path; redundant here (note syntax conventions may differ across systems)
dataPath = './datasets/';

% define initial parameters
ns = 3;             % number of subjects for analysis
nr = 2;             % number of recordings per subject

cmin = 3;           % minimum number of channel estimates required for 
                    % cross-channel averages (for tutorial data: min == 1, max == 6)
fRange = [1 40];    % spectral range (set to filter passband)
w = [7 13];         % alpha peak search window (Hz)
Fw = 11;            % SGF frame width (11 for ~0.24 Hz resolution)
k = 5;              % SGF polynomial order

% initialise data matrices / structures
pSpec = struct('chans', [], 'sums', []);
nchan = nan(1, 2);
muPaf = nan(ns, 1);
muCog = nan(ns, 1); 

%% Analysis loops
for ix = 1:ns     % for each i-th subject
    
    for jx = 1:nr       % for each j-th resting-state recording
        
        % setup filename / path
        fileName = sprintf(['tute_%02d_', num2str(jx), '.set'], ix);
        filePattern = fullfile(dataPath, fileName);
    
        if exist(filePattern, 'file')       % check if filename in dir
            
            % load pre-processed data file
            EEG = pop_loadset('filename', fileName, 'filepath', dataPath);
            [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );

            % select data, set remaining params
            data = EEG.data;
            nchan(jx) = EEG.nbchan;
            Fs = EEG.srate;
            % if data are epoched, concatenate
            if length(size(data)) == 3
                data = reshape(data, nchan(jx), []);
            end

            % run `restingIAF`
            % NOTE: only required inputs specified, further optional inputs 
            % are also available (see `restingIAF` help) 
            [pSpec(ix, jx).sums, pSpec(ix, jx).chans, f]...
                = restingIAF(data, nchan(jx), cmin, fRange, Fs, w, Fw, k);

        else
            sprintf(['unable to find ', fileName, ', skipping file'])

        end

    end         % end j-th recording loop

    % weighted average of mean IAF estimates across (j-th) recordings
    [muPaf(ix, :), muCog(ix, :)] = meanIAF([pSpec(ix, :).sums], nchan, cmin);
    
end

%% Print summary of IAF estimates to console (tabular structure assumes nr == 2)
% find & fill in empties to enable struct fields to be concatenated
emptySums = arrayfun(@(pSpec) isempty(pSpec.sums), pSpec);
[pSpec(emptySums).sums] = deal(struct('paf', NaN, 'pafStd', NaN, 'cog', NaN,...
    'cogStd', NaN, 'muSpec', NaN, 'pSel', NaN, 'gSel', NaN, 'iaw', [1 1]));

% pull out 1st and 2nd recording estimates
sum_1 = [pSpec(:,1).sums];
sum_2 = [pSpec(:,2).sums];

% tabulate summary statistics
S_num = (1:ns)';
colnames = { 'S', 'PAF', 'sd1', 'sd2', 'nPafs1', 'nPafs2', 'CoG',...
    'sd_1', 'sd_2', 'nGravs1', 'nGravs2' };
IAF_estimates = array2table([S_num, muPaf, [sum_1(:).pafStd]', [sum_2(:).pafStd]',...
    [sum_1(:).pSel]', [sum_2(:).pSel]', muCog, [sum_1(:).cogStd]', [sum_2(:).cogStd]',...
    [sum_1(:).gSel]', [sum_2(:).gSel]' ], 'VariableNames', colnames)
fprintf('Estimates derived via pwelch method, initial alpha window = %.1f - %.1f Hz.', w(1), w(2));
