%% Simulation analysis script 2 - Multiple channels, mixed components
% Script for reproducing 2nd (single-peak) and 3rd (split-peak) sets of 
% simulations reported in Corcoran et al (2017).
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
n = 100;                % number of simulated datasets

% restingIAF
fRange = [1 40];        % freq range (same as empirical dataset)
nbchan = 9;             % number of channels per dataset
cMin = 3;               % in this context cMin param is irrelevant (averaging across these independent 'channels' is nonsensical)
w = [7, 13];            % alpha peak search window
Fw = 11;                % SGF frame width
k = 5;                  % SGF polynomial degree

%% simParams - parametric variation of composite alpha signals (embedded in random channel noise) 
simPar = struct('sums', [], 'chans', [], 'maxi', []);

% simulate channel data from varied distributions of alpha components 
% (inverse sd = 1.0/2.5/4.0) at low (0.15) and mod (0.40) SNR
simP1 = simParams(ts, fc_vec, 1, .15, 9, n);
simP2 = simParams(ts, fc_vec, 2.5, .15, 9, n);
simP3 = simParams(ts, fc_vec, 4, .15, 9, n);
simP4 = simParams(ts, fc_vec, 1, .4, 9, n);
simP5 = simParams(ts, fc_vec, 2.5, .4, 9, n);
simP6 = simParams(ts, fc_vec, 4, .4, 9, n);

for sx = 1:6
    for ix = 1:n
        data = eval([sprintf('simP%01d', sx), sprintf('(%01d).y', ix)]);
        [simPar(ix, sx).sums, simPar(ix, sx).chans, f] = restingIAF(data, nbchan, cMin, fRange, Fs, w, Fw, k); 
     
        % average across channel PSDs and find local maximum
        muPSD = nanmean([simPar(ix, sx).chans(:).pxx], 2);
        [simPar(ix, sx).maxi.a, simPar(ix, sx).maxi.b] = findPeak(f, muPSD, w);
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

%% Figure 6: boxplots of estimate error (SNR x alpha)
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

%% table 2
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
simBim = struct('sums', [], 'chans', [], 'maxi', []);
wX = [6 14];             % extended window (capture extreme peaks)

% simulate recovery of split-peaks (inverse sd = 2.5) when right-most peak 
% amplitude is equal (1), +0.25 (2), or +0.50 (3) > left
[simB1, gw] = simBimod(ts, fc_vec, 2.5, 1, .15, 9, n);      % output gw for plotting
simB2 = simBimod(ts, fc_vec, 2.5, 2, .15, 9, n);
simB3 = simBimod(ts, fc_vec, 2.5, 3, .15, 9, n);
simB4 = simBimod(ts, fc_vec, 2.5, 1, .4, 9, n);
simB5 = simBimod(ts, fc_vec, 2.5, 2, .4, 9, n);
simB6 = simBimod(ts, fc_vec, 2.5, 3, .4, 9, n);


for sx = 1:6
    for ix = 1:100
        data = eval([sprintf('simB%01d', sx), sprintf('(%01d).y', ix)]);
        [simBim(ix, sx).sums, simBim(ix, sx).chans, f] = restingIAF(data, nbchan, cMin, fRange, Fs, wX, Fw, k); 
     
        % average across channel PSDs and find local maximum
        muPSD = nanmean([simBim(ix, sx).chans(:).pxx], 2);
        [simBim(ix, sx).maxi.a, simBim(ix, sx).maxi.b] = findPeak(f, muPSD, wX);
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


%% figure 7: boxplot estimate spread (centred on sampling window)
% show estimate deviations (relative to centre of window function)
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
failLm3 = sum(arrayfun(@(simBimod) isnan(simBimod.maxi.a), simBim));
failSg3 = sum(arrayfun(@(simBimod) isnan(simBimod.sums.paf), simBim));
failCg3 = sum(arrayfun(@(simBimod) isnan(simBimod.sums.cog), simBim));

%% table 3
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
