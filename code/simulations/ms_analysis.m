%% Simulation analysis script
% Script for reproducing simulations reported in Corcoran et al (2017). 
% Code for figures and tables included.
%
% Each of the 3 sets of simulations (simSNR, simParams, simBimod) can run
% independently, however params for signal construction must be initialised.
% Ensure `restingIAF` and dependencies are in path.
%
% Part of the `restingIAF` package, (c) Andrew W. Corcoran, 2016-2017.
% Visit github.com/corcorana/restingIAF for licence, updates, and further
% information about package development and testing.
%
%% params for signal construction (required by all sets of simulations)
T = 120;                % time-series duration (sec)
Fs = 250;               % srate (samples/sec)
ts = 0 : 1/Fs : T;      % sample vector
fc_vec = 7.5:0.1:12.5;  % target freq vector
fRange = [1 40];        % freq range (for restingIAF)

%% simSNR - extract single alpha frequency component under varied SNR conditions
[sim1, rdm] = simSNR(ts,fc_vec,0.05,1000);
sim2 = simSNR(ts,fc_vec,0.1,1000);
sim3 = simSNR(ts,fc_vec,0.15,1000);
sim4 = simSNR(ts,fc_vec,0.2,1000);
sim5 = simSNR(ts,fc_vec,0.25,1000);
sim6 = simSNR(ts,fc_vec,0.3,1000);
sim7 = simSNR(ts,fc_vec,0.4,1000);
sim8 = simSNR(ts,fc_vec,0.5,1000);

% apply `restingIAF` to simulated data & analyse how well SGF estimates do
simSnr = struct('chans', [], 'sums', [], 'maxi', []);
cMin = 1;           % in this context cMin param is irrelevant (averaging across these independent 'channels' is nonsensical)

for snr = 1:8
    sim = eval(sprintf('sim%01d', snr));    % each SNR condition treated as one dataset
    data = sim(:,:,3);                      % select combined alpha/pinknoise matrix for analysis
    nbchan = size(data, 1);                 % treat each simulated signal as a channel
    [simSnr(snr).sums, simSnr(snr).chans, f] = restingIAF(data, nbchan, cMin, fRange, Fs);  

    % also collect peak estimates using local maximum (LM) technique
    for kx = 1:nbchan
        [simSnr(snr).maxi(kx).a, simSnr(snr).maxi(kx).b] = findPeak(f, simSnr(snr).chans(kx).pxx, [7 13]);
    end
end


% table 1
% number of PAF component estimates extracted (of 1000 sims)
n_pafs = zeros(2, 8);
for ix = 1:8
    n_pafs(2, ix) = sum(([simSnr(ix).maxi(:).b]==1));
    n_pafs(1, ix) = simSnr(ix).sums.pSel;
end
% eliminate spurious lower bound estimates
[simSnr(1).maxi([simSnr(1).maxi(:).b]==0).a] = deal(NaN);  

peaks = zeros(1000, 8, 2);
for ix = 1:8
    peaks(:, ix, 1) = [simSnr(ix).chans(:).peaks];
    peaks(:, ix, 2) = [simSnr(ix).maxi(:).a];
end

% root mean squared error (Hz) of PAF estimates
rmseLm1 = sqrt(nanmean(bsxfun(@minus, peaks(:,:,2), rdm').^2));
rmseSg1 = sqrt(nanmean(bsxfun(@minus, peaks(:,:,1), rdm').^2)); 

% maximum difference (Hz) between target and estimated PAF
mDiffLm1 = max(abs(bsxfun(@minus, peaks(:,:,2), rdm')));
mDiffSg1 = max(abs(bsxfun(@minus, peaks(:,:,1), rdm')));

% number of estimates that deviated from the target by >1 freq bin (~.24 Hz)
errLm1 = [[1:1000]' bsxfun(@minus, peaks(:,:,2), rdm')];
devBinLm = sum(abs(errLm1(:,2:9))>.24);
errSg1 = [[1:1000]' bsxfun(@minus, peaks(:,:,1), rdm')];
devBinSg = sum(abs(errSg1(:,2:9))>.24);


% print to console
tab1 = [peaks(:,:,1);peaks(:,:,2); rmseLm1;rmseSg1; mDiffLm1;mDiffSg1;mDiffCg2; devBinLm;devBinSg]


% breakdown LM SNR = 0.05 devBin by magnitude of deviation
sum(abs(errLm1(:,2:9))>=1 & abs(errLm1(:,2:9))<2.6)
sum(abs(errLm1(:,2:9))>2.5)


% figure 16: boxplots comparing LM/SGF estimate error across SNR levels
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


% figure 17: visualise smoothed/unsmoothed PSD estimates across SNR levels
figure
indvec = [614 346 988 359 711 110 44 282 906 263];      % randomly sampled
for ix = 1:4
    subplot(2,4,ix)
    plot(f, simSnr(ix).chans(indvec(ix)).pxx, 'Color', [0 .3 .7], 'LineWidth',2)
    hold on
    plot(f, simSnr(ix).chans(indvec(ix)).d0, 'LineWidth',2)
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
    plot(f, simSnr(ix).chans(indvec(ix)).d0, 'LineWidth',2)
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

%% simParams - parametric variation of composite alpha signals (embedded in random channel noise) 
n = 100;

% simulate channel data from varied distributions of alpha components at
% low (0.15) and mod (0.4) SNR
simP1 = simParams(ts, fc_vec, 1, .15, 9, n);
simP2 = simParams(ts, fc_vec, 2.5, .15, 9, n);
simP3 = simParams(ts, fc_vec, 4, .15, 9, n);
simP4 = simParams(ts, fc_vec, 1, .4, 9, n);
simP5 = simParams(ts, fc_vec, 2.5, .4, 9, n);
simP6 = simParams(ts, fc_vec, 4, .4, 9, n);

% subject to pwelch / SGF
simPar = struct('sums', [], 'chans', [], 'f', [], 'maxi', []);
nbchan = 9;
cMin = 3;

for sx = 1:6
    for ix = 1:100
        data = eval([sprintf('simP%01d', sx), sprintf('(%01d).y', ix)]);
        [simPar(ix, sx).sums, simPar(ix, sx).chans, simPar(ix, sx).f] = restingIAF(data, nbchan, cMin, fRange, Fs, [7 13]); 
     
        % average across channel PSDs and find local maximum
        muPSD = nanmean([simPar(ix, sx).chans(:).pxx], 2);
        [simPar(ix, sx).maxi.a, simPar(ix, sx).maxi.b] = findPeak(simPar(ix, sx).f, muPSD, [7 13]);
    end
end

% extract estimate arrays
pafLm = arrayfun(@(simParam) (simParam.maxi.a), simPar);
pafSg = arrayfun(@(simParam) (simParam.sums.paf), simPar);
cog = arrayfun(@(simParam) (simParam.sums.cog), simPar);

% calculate error matrices
ta2 = [simP1(:).t];      % extract target alpha freq vector (same across all conditions)
errLm2 = bsxfun(@minus, pafLm, ta2');
errSg2 = bsxfun(@minus, pafSg, ta2');
errCg2 = bsxfun(@minus, cog, ta2');

% figure 18: boxplots of estimate error (SNR x alpha)
figure
for bp = 1:3
    subplot(1,3,bp)
    boxplot([errLm2(:,bp), errSg2(:,bp), errCg2(:,bp)], 'colors','k', 'Symbol','ko')
    hold on
    line([0 4], [0 0], 'Color', [0.3 0.3 0.3], 'LineStyle', '--')
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),'y','FaceAlpha',.3);
    patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.3);
    patch(get(h(3),'XData'),get(h(3),'YData'),'b','FaceAlpha',.3);
    set(findobj('-regexp','Tag','\w*Whisker'),'LineStyle','-')
    if bp ==1
        title('\alpha = 1.0')
        ylabel('Estimate error (Hz)')
        ylim([-2 2])
    elseif bp ==2
        title('\alpha = 2.5')
        ylim([-1.5 1.5])
    elseif bp ==3
        title('\alpha = 4.0')
        ylim([-1 1])
    end
    set(gca, 'xticklabel', {'LM' 'SG' 'CoG'})
    set(gca, 'FontSize', 14)
end
suptitle('SNR = 0.15')

figure      
for bp = 4:6
    subplot(1,3,bp-3)
    boxplot([errLm2(:,bp), errSg2(:,bp), errCg2(:,bp)], 'colors','k', 'Symbol','ko')
    hold on
    line([0 4], [0 0], 'Color', [0.3 0.3 0.3], 'LineStyle', '--')
    h = findobj(gca,'Tag','Box');
    patch(get(h(1),'XData'),get(h(1),'YData'),'y','FaceAlpha',.3);
    patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.3);
    patch(get(h(3),'XData'),get(h(3),'YData'),'b','FaceAlpha',.3);
    set(findobj('-regexp','Tag','\w*Whisker'),'LineStyle','-')
    if bp ==4
        title('\alpha = 1.0')
        ylabel('Estimate error (Hz)') 
        ylim([-2 2])
    elseif bp ==5
        title('\alpha = 2.5')
        ylim([-1.5 1.5])
    elseif bp ==6
        title('\alpha = 4.0')
        ylim([-1 1])
    end
    set(gca, 'xticklabel', {'LM' 'SG' 'CoG'})
    set(gca, 'FontSize', 14)
end
suptitle('SNR = 0.40')


% number of simulations that failed to generate estimates
failLm2 = sum(arrayfun(@(simParam) isnan(simParam.maxi.a), simPar));
failSg2 = sum(arrayfun(@(simParam) isnan(simParam.sums.paf), simPar));
failCg2 = sum(arrayfun(@(simParam) isnan(simParam.sums.cog), simPar));

% table 2
rmseLm2 = sqrt(nanmean(errLm2.^2));
rmseSg2 = sqrt(nanmean(errSg2.^2));
rmseCg2 = sqrt(nanmean(errCg2.^2));

mDiffLm2 = max(abs(errLm2));
mDiffSg2 = max(abs(errSg2));
mDiffCg2 = max(abs(errCg2));

% percentage of esimates deviating from target freq by > 0.5 Hz
pDevLm = nansum(abs(errLm2)>.5);
pDevSg = round((nansum(abs(errSg2)>.5)./(100-failSg2))*100);
pDevCg = nansum(abs(errCg2)>.5);

% median (std dev) number of channels (per set of 9) contributing to
% estimate
pChans2 = arrayfun(@(simParam) (simParam.sums.pSel), simPar);
gChans2 = arrayfun(@(simParam) (simParam.sums.gSel), simPar);

% print to console
tab2 = [rmseLm2;rmseSg2;rmseCg2; mDiffLm2;mDiffSg2;mDiffCg2; pDevLm;pDevSg;pDevCg; median(pChans2);std(pChans2);median(gChans2);std(gChans2)]

%% simBimod - same general scheme as simParams, but with bimodal window function (hold alpha @ 2.5)
n = 100;

[simB1, gw] = simBimod(ts, fc_vec, 2.5, 1, .15, 9, n);      % output gw for plotting
simB2 = simBimod(ts, fc_vec, 2.5, 2, .15, 9, n);
simB3 = simBimod(ts, fc_vec, 2.5, 3, .15, 9, n);
simB4 = simBimod(ts, fc_vec, 2.5, 1, .4, 9, n);
simB5 = simBimod(ts, fc_vec, 2.5, 2, .4, 9, n);
simB6 = simBimod(ts, fc_vec, 2.5, 3, .4, 9, n);

% subject to pwelch / SGF
simBim = struct('sums', [], 'chans', [], 'f', [], 'maxi', []);
nbchan = 9;
cMin = 3;

% NB: extended alpha window due to potential for higher frequency PAFs when
% window centred about 12-12.5
for sx = 1:6
    for ix = 1:100
        data = eval([sprintf('simB%01d', sx), sprintf('(%01d).y', ix)]);
        [simBim(ix, sx).sums, simBim(ix, sx).chans, simBim(ix, sx).f] = restingIAF(data, nbchan, cMin, fRange, Fs, [7 15]); 
     
        % average across channel PSDs and find local maximum
        muPSD = nanmean([simBim(ix, sx).chans(:).pxx], 2);
        [simBim(ix, sx).maxi.a, simBim(ix, sx).maxi.b] = findPeak(simBim(ix, sx).f, muPSD, [7 15]);
    end
end

% check for ineligible LM ests
for xi = 1:n
    for xj = 1:6
        if simBim(xi,xj).maxi.b==0
            simBim(xi,xj).maxi.a=NaN;
        end
    end
end

% extract estimate arrays
pafLm = arrayfun(@(simBim) (simBim.maxi.a), simBim);
pafSg = arrayfun(@(simBim) (simBim.sums.paf), simBim);
cog = arrayfun(@(simBim) (simBim.sums.cog), simBim);

% calculate error matrices
ta3 = [simB1(:).t];      % extract target alpha freq vector (same across all conditions - should be identical to simParams)
errLm3 = bsxfun(@minus, pafLm, ta3');
errSg3 = bsxfun(@minus, pafSg, ta3');
errCg3 = bsxfun(@minus, cog, ta3');


% figure 19: boxplot estimate spread (centred on sampling window)
% get deviations (relative to centre of window function)
figure
bp = 1;
for sp = [2 5 8 3 6 9]
	subplot(3,3,sp)
	boxplot([errLm3(:,bp), errSg3(:,bp), errCg3(:,bp)], 'colors','k', 'Symbol','ko')
	hold on
	line([0 4], [0 0], 'Color', [0.3 0.3 0.3], 'LineStyle', '--')
    line([0 4], [.8 .8], 'Color', [0.3 0.3 0.3], 'LineStyle', ':')
	line([0 4], [-.8 -.8], 'Color', [0.3 0.3 0.3], 'LineStyle', ':')
	h = findobj(gca,'Tag','Box');
	patch(get(h(1),'XData'),get(h(1),'YData'),'y','FaceAlpha',.3);
	patch(get(h(2),'XData'),get(h(2),'YData'),'r','FaceAlpha',.3);
	patch(get(h(3),'XData'),get(h(3),'YData'),'b','FaceAlpha',.3);
	ylim([-2.2 2.2])
	set(findobj('-regexp','Tag','\w*Whisker'),'LineStyle','-')
    if sp ==2
        title('SNR = 0.15')
    elseif sp ==3
        title('SNR = 0.40')
    end
    set(gca, 'xticklabel', {'LM' 'SG' 'CoG'})
    set(gca, 'FontSize', 14)
    bp = bp+1;
end

subplot(3,3,1)
plot(gw)
hold on
plot(17:51,gw)
camroll(90)
title('Sampling Window')
xlim([0 52])
ylim([0 1.7])
xlabel('Low $f$ -- High $f$', 'Interpreter', 'latex')
ylabel('Scale factor = 0')
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'FontSize', 14)

subplot(3,3,4)
plot(gw)
hold on
plot(17:51,gw.*1.25)
camroll(90)
xlim([0 52])
ylim ([0 1.7])
xlabel('Low $f$ -- High $f$', 'Interpreter', 'latex')
ylabel('Scale factor = +.25')
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'FontSize', 14)

subplot(3,3,7)
plot(gw)
hold on
plot(17:51,gw.*1.5)
camroll(90)
xlim([0 52])
ylim ([0 1.7])
xlabel('Low $f$ -- High $f$', 'Interpreter', 'latex')
ylabel('Scale factor = +.50')
set(gca, 'xtick', [], 'ytick', [])
set(gca, 'FontSize', 14)

% number of failed cases
failLm3 = arrayfun(@(simBimod) (simBimod.maxi.a), simBim);
failSg3 = arrayfun(@(simBimod) (simBimod.sums.paf), simBim);
failCg3 = arrayfun(@(simBimod) (simBimod.sums.cog), simBim);

% table 3
rmseLm3 = sqrt(nanmean(errLm3.^2));
rmseSg3 = sqrt(nanmean(errSg3.^2));
rmseCg3 = sqrt(nanmean(errCg3.^2));

mDiffLm3 = max(abs(errLm3));
mDiffSg3 = max(abs(errSg3));
mDiffCg3 = max(abs(errCg3));

pChans3 = arrayfun(@(simBimod) (simBimod.sums.pSel), simBim);
gChans3 = arrayfun(@(simBimod) (simBimod.sums.gSel), simBim);

% print to console
tab3 = [rmseLm3;rmseSg3;rmseCg3; mDiffLm3;mDiffSg3;mDiffCg3; median(pChans3);std(pChans3);median(gChans3);std(gChans3)]
