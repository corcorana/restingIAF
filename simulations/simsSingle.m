%% Simulation analysis script 1 - Single component SNR
% Script for reproducing 1st set of simulations reported in Corcoran et al 
% (2017). Code for figures and tables included.
%
% Ensure `restingIAF` and dependencies are in path.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%% Simulation parameters
% generic signal construction
T = 120;                % time-series duration (sec)
Fs = 250;               % srate (samples/sec)
ts = 0 : 1/Fs : T;      % sample vector
fc_vec = 7.5:0.1:12.5;  % target freq vector
n = 1000;               % number of simulated signals

% restingIAF
fRange = [1 40];        % freq range (same as empirical dataset)
cMin = 1;               % in this context cMin param is irrelevant (averaging across these independent 'channels' is nonsensical)
w = [7 13];             % alpha peak search window
Fw = 11;                % SGF frame width
k = 5;                  % SGF polynomial degree

%% simSNR - extract single alpha frequency component under varied SNR conditions
% null simulation (no target alpha component in signal)
simN = struct('chans', [], 'sums', [], 'maxi', []);

simNull = simSNR(ts,fc_vec,0,n);

data = simNull(:,:,3);                  % select combined alpha/pinknoise matrix for analysis
nbchan = size(data, 1);                 % treat each simulated signal as a channel
[simN.sums, simN.chans, f] = restingIAF(data, nbchan, cMin, fRange, Fs, w, Fw, k);  
for kx = 1:nbchan
 	[simN.maxi(kx).a, simN.maxi(kx).b] = findPeak(f, simN.chans(kx).pxx, w);
end

% number of (false) alpha component estimates extracted
n_pafsN = zeros(2, 1);
n_pafsN(2) = sum(([simN.maxi(:).b]==1));
n_pafsN(1) = simN.sums.pSel;

% eliminate spurious lower bound estimates
[simN.maxi([simN.maxi(:).b]==0).a] = deal(NaN);  
fpos = [simN.maxi(:).a];
hfpos = hist(fpos, 7:.25:13);


%% graded SNR simulations
[sim1, rdm] = simSNR(ts,fc_vec,0.05,n);
sim2 = simSNR(ts,fc_vec,0.1,n);
sim3 = simSNR(ts,fc_vec,0.15,n);
sim4 = simSNR(ts,fc_vec,0.2,n);
sim5 = simSNR(ts,fc_vec,0.25,n);
sim6 = simSNR(ts,fc_vec,0.3,n);
sim7 = simSNR(ts,fc_vec,0.4,n);
sim8 = simSNR(ts,fc_vec,0.5,n);

% apply `restingIAF` to simulated data & analyse how well SGF estimates do
simSnr = struct('chans', [], 'sums', [], 'maxi', []);

for snr = 1:8
    sim = eval(sprintf('sim%01d', snr));    % each SNR condition treated as one dataset
    data = sim(:,:,3);                      % select combined alpha/pinknoise matrix for analysis
    nbchan = size(data, 1);                 % treat each simulated signal as a channel
    [simSnr(snr).sums, simSnr(snr).chans, f] = restingIAF(data, nbchan, cMin, fRange, Fs, w, Fw, k);  

    % also collect peak estimates using local maximum (LM) technique
    for kx = 1:nbchan
        [simSnr(snr).maxi(kx).a, simSnr(snr).maxi(kx).b] = findPeak(f, simSnr(snr).chans(kx).pxx, w);
    end
end


%% table 1
% number of PAF component estimates extracted (of 1000 sims)
n_pafs = zeros(2, 8);
for ix = 1:8
    n_pafs(2, ix) = sum(([simSnr(ix).maxi(:).b]==1));
    n_pafs(1, ix) = simSnr(ix).sums.pSel;
end
% eliminate spurious lower bound estimates
[simSnr(1).maxi([simSnr(1).maxi(:).b]==0).a] = deal(NaN);  

peaks = zeros(n, 8, 2);
for ix = 1:8
    peaks(:, ix, 1) = [simSnr(ix).chans(:).peaks];
    peaks(:, ix, 2) = [simSnr(ix).maxi(:).a];
end

failLm1 = sum(isnan(peaks(:,:,2)));
failSg1 = sum(isnan(peaks(:,:,1)));

% root mean squared error (Hz) of PAF estimates
rmseLm1 = sqrt(nanmean(bsxfun(@minus, peaks(:,:,2), rdm').^2));
rmseSg1 = sqrt(nanmean(bsxfun(@minus, peaks(:,:,1), rdm').^2)); 

% maximum difference (Hz) between target and estimated PAF
mDiffLm1 = max(abs(bsxfun(@minus, peaks(:,:,2), rdm')));
mDiffSg1 = max(abs(bsxfun(@minus, peaks(:,:,1), rdm')));

% number of estimates that deviated from the target by >1 freq bin (~.24 Hz)
errLm1 = [[1:n]' bsxfun(@minus, peaks(:,:,2), rdm')];
devBinLm = sum(abs(errLm1(:,2:9))>.24);
errSg1 = [[1:n]' bsxfun(@minus, peaks(:,:,1), rdm')];
devBinSg = sum(abs(errSg1(:,2:9))>.24);

% print to console
tab1a = [n-failLm1;1000-failSg1]
tab1b = [rmseLm1;rmseSg1; mDiffLm1;mDiffSg1; devBinLm;devBinSg]

% breakdown LM SNR = 0.05 devBin by magnitude of deviation
sum(abs(errLm1(:,2:9))>=1 & abs(errLm1(:,2:9))<2.6);
sum(abs(errLm1(:,2:9))>2.5);


% boxplots comparing LM/SGF estimate error across SNR levels (unpublished)
figure
titvec = {'SNR = 0.05'; 'SNR = 0.10'; 'SNR = 0.15'; 'SNR = 0.20'; 'SNR = 0.25'; 'SNR = 0.30'; 'SNR = 0.40'; 'SNR = 0.50'};
for bp = 1:2
    subplot(1,2,bp)
    boxplot([peaks(:,bp,2)-rdm', peaks(:,bp,1)-rdm'], 'colors','k', 'Symbol','ko', 'jitter',0.5)
    hold on
    line([0 3], [0 0], 'Color', [0.3 0.3 0.3], 'LineStyle', '--')
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),'r','FaceAlpha',.3);
    patch(get(h(2),'XData'),get(h(2),'YData'),'b','FaceAlpha',.3);
    set(findobj('-regexp','Tag','\w*Whisker'),'LineStyle','-')
    title(titvec{bp})    
    if bp == 1
        ylabel(' ')
        ylim([-6 1.5])
    else
        ylim([-2 1.2])
    end
    set(gca, 'xticklabel', {' ' ' '})
    set(gca, 'FontSize', 14)
end

figure
for bp = 3:8
    subplot(2,3,bp-2)
    boxplot([peaks(:,bp,2)-rdm', peaks(:,bp,1)-rdm'], 'colors','k', 'Symbol','ko', 'jitter',0.5)
    hold on
    line([0 3], [0 0], 'Color', [0.3 0.3 0.3], 'LineStyle', '--')
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),'r','FaceAlpha',.3);
    patch(get(h(2),'XData'),get(h(2),'YData'),'b','FaceAlpha',.3);
    set(findobj('-regexp','Tag','\w*Whisker'),'LineStyle','-')
    title(titvec{bp})    
    if bp == 3
        ylabel(' ')
    end
    if bp == 3 || bp == 4 || bp == 5
        ylim([-1 1])
    elseif bp == 6 || bp == 7 || bp == 8
        ylim([-0.5 0.5])
    end
    if bp == 6
        set(gca, 'xticklabel', {'LM' 'SG'})
        ylabel('Estimate error (Hz)')
    else
        set(gca, 'xticklabel', {' ' ' '})
    end
    set(gca, 'FontSize', 14)
end


%% figure 5: visualise e.g. smoothed/unsmoothed PSD estimates across SNR levels
figure
indvec = [614 346 988 359 711 110 44 282];      % randomly sampled
for ix = 1:4
    subplot(2,4,ix)
    plot(f, simSnr(ix).chans(indvec(ix)).pxx, 'Color', [0 .3 .7], 'LineWidth',2)
    hold on
    plot(f, simSnr(ix).chans(indvec(ix)).d0,  'Color', [1 .2 0], 'LineWidth',2)
    xlim([3 17])
    ylim([0 15])
    title(titvec{ix})
    text(11.3,13.5,['F\alpha = ', num2str(rdm(indvec(ix)))], 'FontSize',12)
    text(11.3,12.2,['PAF_{SG} = ',num2str(peaks(indvec(ix),ix,1),3)], 'Color', [1 .2 0], 'FontSize',12)
    text(11.3,11,['PAF_{LM} = ', num2str(peaks(indvec(ix),ix,2),3)], 'Color', [0 .3 .7], 'FontSize',12)
    set(gca, 'FontSize', 14)
end
for ix = 5:8
    subplot(2,4,ix)
    plot(f, simSnr(ix).chans(indvec(ix)).pxx, 'Color', [0 .3 .7], 'LineWidth',2)
    hold on
    plot(f, simSnr(ix).chans(indvec(ix)).d0, 'Color', [1 .2 0], 'LineWidth',2)
    xlim([3 17])
    ylim([0 25])
    title(titvec{ix})
    text(11.3,22.8,['F\alpha = ', num2str(rdm(indvec(ix)))], 'FontSize',12)
    text(11.3,20.5,['PAF_{SG} = ',num2str(peaks(indvec(ix),ix,1),3)], 'Color', [1 .2 0], 'FontSize',12)
    text(11.3,18.5,['PAF_{LM} = ', num2str(peaks(indvec(ix),ix,2),3)], 'Color', [0 .3 .7], 'FontSize',12)
    set(gca, 'FontSize', 14)
end
subplot(2,4,5)
xlabel('Frequency (Hz)')
ylabel('Power (a.u.)')
