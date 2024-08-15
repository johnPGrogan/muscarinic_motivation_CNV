function TrioEEGAnalyseByTrials(timelock);
% Analyse EEG data that can be split by the result.data trial information.
% This requires that the same number of trials are found in the EEG
% recordings and in the result.data structure, which may not be the case
% due to missing triggers or aborted trials (use TrioEEGAnalyseByBins if
% this is the case).
% This script will load up the individual EEG data, split it by the trial
% types (or other factors e.g. location, direction, captured), and will run
% mean-wise and trial-wise analyses. 
% This can allow for comparison of EEG based on the responses (e.g. by RT,
% accuracy, capture, etc), as well as individual trial analyses.

%% set up stuff
% clc; clear all; eeglab; close all;
% hh = findobj('tag', 'EEGLAB');
% if isempty(hh); eeglab; end % call this to setup path and initialise

%% options

if ~exist('timelock','var') || isempty(timelock)
    timelock = 'readyBS';
end

if regexpi(timelock, 'sacc|targ|fb')
    plotDistr = 1; % include distr in plots?
else
    plotDistr = 0;
end

loadFile = 1;

flipDrugs = 1;
fileType = 'erpTrials'; % or 'eegTrials' if erp averaging not run


%%
set(0, 'DefaultAxesFontSize',14)
load('chanlocs.mat');

triFolder = 'D:/groganj/Work/Trihexy/'; % base folder for data
if ~exist(triFolder,'file')
    triFolder = 'D:/Backups/Trihexy/';
end

analysisFolder = [triFolder 'EEGAnalysis/']; % where to load EEG data from

switch timelock
    case 'iti'
        suffix = '_processed_filt resampled_be_iti_bs100_artflag_trigflag_erp.erp'; % unbaselined, for baselining
        timesToPlot = [-100 500];
        toi = [80 150; 220 320]; % times of interest
        doBaseline = 0; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'visCue'
        suffix = '_processed_filt resampled_be_visCue_bs100_artflag_trigflag_erp.erp'; % unbaselined, for baselining
        timesToPlot = [-100 500];
        toi = [200 300; 300 max(xTimes)]; % times of interest
        doBaseline = 0; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'audCue'
        suffix = '_processed_filt resampled_be_audCue_none_artflag_trigflag_erp.erp'; % unbaselined, for baselining
        timesToPlot = [-200 000];
        toi = [-100 0; -200 0]; % times of interest
        doBaseline = 0; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'audCue_200'
        suffix = '_processed_filt resampled_be_audCue_200_artflag_trigflag_erp.erp';
        timesToPlot = [-200 1100];
        toi = [200 280; 900 1100]; % times of interest
        doBaseline = 0; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'ready'
        suffix = '_processed_filt resampled_be_readyCue_bs_audCue200_artflag_trigflag_erp.erp';
        timesToPlot = [-100 1500];
        toi = [200 280; 1200 1500]; % times of interest
        doBaseline = 0; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'readyBS' % ready, not baselined, will manually rebaseline
        suffix = '_processed_filt resampled_be_readyCue_none_artflag_trigflag_erp.erp';
        timesToPlot = [-100 1500];
        toi = [200 280; 1200 1500]; % times of interest
        doBaseline = 1; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'ready500'
        suffix = '_processed_filt resampled_be_readyCue500_bs_audCue200_artflag_trigflag_erp.erp';
        timesToPlot = [-100 1500];
        toi = [200 280; 1200 1500]; % times of interest
        doBaseline = 1; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'targ'
        suffix = '_processed_filt resampled_be_targ_bs_audCue200_artflag_trigflag_erp.erp';
        timesToPlot = [-100 400];
        toi = [100 160; 200 300]; % times of interest
        doBaseline = 1; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'sacc'
        suffix = '_processed_filt resampled_be_firstSacc_bs_audCue200_artflag_trigflag_erp.erp';
        timesToPlot = [-500 200];
        toi = [-200 -100; -50 50]; % times of interest
        doBaseline = 2; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    case 'fb'
        suffix = '_processed_filt resampled_be_fb_bs_100_artflag_trigflag_erp.erp';
        timesToPlot = [-100 400];
        toi = [170 220; 250 320]; % times of interest
        doBaseline = 1; % 0=none, 1= [-100 0] here, 2 = 100ms pre-target (done later)
    otherwise
        error('timelock not found');
end


%% 
loadedFile = sprintf('TrioEEGAnalyseByTrials_%s_output.mat', timelock);
if loadFile && exist(loadedFile,'file')
    old = load(loadedFile);
    new = workspace2struct('old');
    % check this loaded data matches this file
    fn = {'suffix','flipDrugs','fileType'};
    for i = 1:length(fn)
        if ~equals(old.(fn{i}), new.(fn{i}))
            error('loaded data does not match options set in this file')
        end
    end
    output = old.output;
    clearvars new old;
else
    output = LoadPipelineFiles(analysisFolder, suffix, flipDrugs, fileType);
    save(loadedFile,'-v7.3'); % save this output here
end

struct2workspace(output, false); % unload struct into workspace, overwriting EEG
xTimes = times;

eegRawData = single(eegRawData); % need to make double or rmanova doesn't work
%% re-baseline?

if doBaseline == 1
    baselineWindow = [-100 0]; % xTimes
    baselineInds = isBetween(xTimes, baselineWindow);

    baselineMean = nanmean(eegRawData(:,baselineInds,:,:,:),2); % get mean

    eegRawData = eegRawData - baselineMean; % subtract
end


%% pick your factors

% make all result.data structs have the same fields
behData = ensureAllStructsAssignable(behData, 0);

% un-cell them
behData = permute(reshape(nancat(1, behData{:}), nPP, 2, []),[3,1,2]); % [tr pp ses]

% make a function to extract factors easily
extractVar = @(behData, x) reshape(nancat(1,behData.(x)), nTrials, nPP, nSes, [] ); 
% x = fieldname, extra dims become 4th dim

isDist = extractVar(behData, 'distractor');
isRew = extractVar(behData, 'rewardMax');

trialType = 2*isDist + (isRew==50)+1;
% trial type
% tt = sq(nancat(tt)); % [ tr pp ses]
% check they match
% all(tt == trialType | isnan(trialType),'all') % they do

% also get location + direction
fixLoc = extractVar(behData, 'lastTarget'); % 1=left,2=right,3=bottom
targLoc = extractVar(behData, 'targetPos'); % 1=left,2=right,3=bottom
targDir = extractVar(behData, 'targetDir'); % 1 = clockwise, 2 = anti-clockwise


% get reward received
outcome = extractVar(behData, 'reward');
fbSound = outcome;
fbSound(outcome > 30) = 2;
fbSound(outcome > 9 & outcome <= 30) = 1;
fbSound(outcome <= 9) = 0;
fbSound2 = permute(groupMeans(fbSound,1,trialType,'dim'),[2,1,3,4]); %[npp conds drug trial]

%% load beh analysis also, can use those data to split
%  e.g. whether they were captured


%% check EEG and result.data match up 

nEEGTrials = sq(sum(~all(isnan(eegRawData),[1,2]),3));
nBehTrials = sq(sum(~isnan(trialType)));

if any(nEEGTrials ~= nBehTrials,'all')
    error('the number of EEG and behavioural trials does not match');
end

%% set up filtering

behFlagNames = {'isPaused','driftFlag1','driftFlag2'}; % conditions for flagging
behRTLimits = [100 1500]; % ms
useR = 1; % filter if result.data.R < 0
blinkTrigs = { {'31|32'}, {'B','54'} }; % filter blinks between B and 54
% saccTrigs = { {'52|53'}, {'B','54'} }; % filter if no saccades found between B and 54

isTrial = ~isnan(trialType); % trial was attempted, not missing 
[filters, filtersSummary] = SetupFilteringByData(behData, isTrial, ...
    'behFlagNames', behFlagNames,'behRTLimits', behRTLimits, 'useR', useR,...
    'rejInds', rejInds, 'trigs', trigs, 'blinkTrigs', blinkTrigs);%, 'saccTrigs', saccTrigs);


%% actually filter them

% change shape to allow indexing
eegRawData(repmat(permute(filters.toExclude,[4,5,1,2,3]),nChans,nTimes)) = NaN; % NaN is excluded


%% split EEG data by factors

eegByCond = groupMeans(eegRawData,3, repmat(permute(trialType,[4,5,1,2,3]),nChans, nTimes),'dim');
[nChans, nTimes, nBins, nPP, nDrugs, nTr] = size(eegByCond); 


%% re-baseline to 100ms before target onset
if doBaseline == 2
    load('TrioEEGBehAnalysis.mat','vars');
    srt = vars.srt;
    for iPP = 1:nPP
        for iC = 1:4
            for iDr = 1:2
                for iTr = 1:nTr
                    baselineMean(:,1,iC,iPP,iDr,iTr) = nanmean(eegByCond(:,isBetween(xTimes, vars.srt(iPP,iC,iDr,iTr) - [100 0]),iC,iPP,iDr,iTr),2);
                end
            end
        end
    end

    eegByCond = eegByCond - baselineMean;
end
%%
[nChans, nTimes, nBins, nPP, nDrugs, nTr] = size(eegByCond); 

% split into [rew dist]
eegByCond2 = permute(reshape(eegByCond, nChans, nTimes, 2, 2, nPP, nDrugs, nTr), [5 2 3 4 6 1 7]);
% [pp times rew dist drug chan tr]

%% split the filtering also

filtersByCond = structfun(@(x) permute(groupMeans(x, 1, trialType,'dim'),[2,1,3,4]), filters,'UniformOutput',0); 
%[nPP nBins nDrugs nTr]

filtersByCondSummary = structfun(@(x) sq(nansum(x,4)), filtersByCond,'UniformOutput',0); % summarise
filtersByCondSummary.toInclude = filtersByCondSummary.isTrial - filtersByCondSummary.toExclude; % number included

% do exclusions differ by condition?
exclAnova = rmanova(reshape(filtersByCondSummary.toExclude,nPP,2,2,2),{'pp','rew','distr','drug'},'categorical',2:4)

%% numbers excluded
% sq(nanmean(filtersSummary.toInclude,'all')) % mean num included
% sq(nanstd(filtersSummary.toInclude,[],'all')) % SD
% sq(nanmean(filtersSummary.toInclude ./ filtersSummary.isTrial,'all')) * 100 % %included
% sq(nanstd(filtersSummary.toInclude ./ filtersSummary.isTrial,[],'all')) * 100 %SD
% 
% % per cond
% sq(nanmean(filtersByCondSummary.toInclude,'all')) % mean num per cond
% sq(nanstd(filtersByCondSummary.toInclude,[],'all')) %SD
% sq(nanmean(filtersByCondSummary.toInclude ./ filtersByCondSummary.isTrial,'all'))*100 %mean % per cond
% sq(nanstd(filtersByCondSummary.toInclude ./ filtersByCondSummary.isTrial,[],'all'))*100 % SD


%% get grand averages

avDims = [6 5]; % dims to average over, 6=trials, 5=drugs, 4 = distr
permDims = [4,2,3,1]; % dims to permute with, length = ndims(eegByCond) - length(avDims)
eegByCondMean = permute(nanmean(eegByCond, avDims),permDims); % [pp time bin chan]

%% can also take var/sd over trials or do windowing here
eegByCondStd = nanstd(eegByCond2, [], 7); % SD over trials [npp time rew distr drug chan

% moving std across time
eegByCondMovStd = permute(movstd(eegByCond, 5, [], 2), [4,2,3,1,5,6]);

% moving std of diffs
eegByCondMovStdDiff = [NaN(nPP,1,nBins, nChans, nDrugs, nTr), permute(movstd(diff(eegByCond,[],2), 5, [], 2), [4,2,3,1,5,6])];


% take SD within windows
windows = timesToPlot(1):100:timesToPlot(2);
midWindows = windows(1:end-1) + diff(windows)/2; % midpoints of windows

% get SD across window within each trial
eegByCondWindowStd = EEGSplitByWindow(eegByCond, xTimes, windows);
% eegWindowMean = nanmean(eegByCondWindowStd,7); % mean voltage in window

% std within window, mean over trials
eegByCondWindowStd = permute(nanmean(nanstd(eegByCondWindowStd,[],7),avDims), permDims); 



%% average across conditions to pick a window
load('chanlocs.mat');

%%% pick channels to plot
% chansToPlot = chanNames;
chansToPlot = {'Cz'};
% chansToPlot = {'FPz','Fz','FCz','Cz','CPz','Pz','POz'};%,'elX', 'elY'}; % all midlines
% chansToPlot = {'FC1','FCz','FC2','C3','Cz','C4','CP1','CPz','CP2'};
% chanInds = find(ismember(chanNames, chansToPlot));



condNames = {'targLo','targHi','distLo','distHi'};
drugNames = {'Placebo','THP'};
condNames2 = {'lowPlacebo','highPlacebo','lowTHP','highTHP'};
drugRewNames = {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'};

grandAverage = nanmean(eegByCondMean,3); %[pp time chan]

f = figure();
% subplot(2,2,1);
h = PlotEEGErrorBars(grandAverage, 'chanNames',chanNames, 'chansToPlot', chansToPlot,...
    'XTimes', xTimes, 'timesToPlot', timesToPlot, 'condNames', [],...
    'superTitle','', 'yLabels', {'\muV'},'xLabels',{'time from ready cue (ms)'},'subplotInds',GetSubPlotShape(length(chansToPlot)));

% % show TOI
hold on;
fill(toi(1,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');
fill(toi(2,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');


%% topoplot grand average

% animated
% AnimateTopoplots(sq(nanmean(grandAverage))', xTimes, 181:641, chanLocs, {''}, [], 'readyCue_voltage.gif',crameri('vik'),0);

% specific times
% cmap = othercolor('PRGn11');
cmap = 'jet';
mapLimits = [-6 6];
% toi = [100 160; 200 300];
for i = 1:size(toi,1)
    figure();
    
    topoplot(sq(nanmean(grandAverage(:,isBetween(xTimes, toi(i,:)),:,:),[1 2]))', chanLocs(1:61),...
            'mapLimits',[-1 1]*6, 'colormap',cmap,...
            'style', 'both', 'plotrad',.55, 'headrad', .5,'emarker', {'.','k',[],1},...
            'numcontour', 6,'electrodes', 'on', 'nosedir', '+X');
%     title([num2str(round(xTimes(find(xTimes>=toi(i),1,'first')))) ' ms']);
    title(sprintf('%d:%dms', toi(i,1), toi(i,2)));
    colorbar;
end

%% plot grand average by cond


figure();
h = PlotEEGErrorBars(eegByCondMean, 'chanNames',chanNames, 'chansToPlot', chansToPlot,...
    'XTimes', xTimes(1:end-1), 'timesToPlot', timesToPlot, 'condNames', condNames,...
    'yLabels', {'\muV'},'yLine',1);
box off;
% % show TOI
hold on;
fill(toi(1,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');
fill(toi(2,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');
legend([h{:,1}], condNames,'Location','NorthWest'); % overwrite to remove grey boxes

%% plot average over targ/dist, but with drugs

figure();
eegByRewDrug = sq(nanmean(eegByCond2, [7 4]));%
cols = [0 0 1; 0 .5 .9;1 0 0;  .9 .5 0];
% subplot(2,2,[1 2])
h = PlotEEGErrorBars(reshape(eegByRewDrug,nPP, nTimes, 4, nChans), 'chanNames',chanNames, 'chansToPlot', {'Cz'},...
    'XTimes', xTimes, 'timesToPlot', timesToPlot, 'condNames', drugRewNames,...
    'yLabels', {'\muV'},'lineCols',cols,'xLabels',{'time from preparation cue (ms)'},'subplotInds',0);
box off
% show TOI
hold on;
fill(toi(1,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');
fill(toi(2,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');
legend([h{:,1}], drugRewNames,'Location','Best'); % overwrite to remove grey boxes
%% plot STD windowed

figure();
h = PlotEEGErrorBars(eegByCondWindowStd, 'chanNames',chanNames, 'chansToPlot', chansToPlot,...
    'XTimes', midWindows, 'timesToPlot', minMax(midWindows,2), 'condNames', condNames,...
    'superTitle','SD within window', 'yLabels', {'SD \muV'},'yLine',0);

%% trial-variance

figure();
h = PlotEEGErrorBars(reshape(nanmean(eegByCondStd,4),nPP,nTimes,4,nChans), 'chanNames',chanNames, 'chansToPlot', chansToPlot,...
    'XTimes', xTimes, 'timesToPlot', timesToPlot, 'condNames', condNames2,...
    'superTitle','voltage', 'yLabels', {'\muV SD'},'yLine',0);

% take mean SD in window
meanSDWindow = nanmean(eegByCondStd(:,isBetween(xTimes, toi(1,:)),:,:,:,:),2);
meanSDWindowAnova = rmanova(sq(meanSDWindow(:,:,:,:,:,30)), {'pp','rew','distr','drug'});

tWindows = timesToPlot(1):100:timesToPlot(2);
for i = 1:length(tWindows)-1
    disp(tWindows(i));

    % use this if doing multiple channels
%     meanWindow = nanmean(permute(GetMeanAmpl(eegByCond2, times, [tWindows(i) tWindows(i+1)]),[1 3 4 5 6 7 2]),6); % mean within window per trial
    
    % quicker for one channel
    meanWindowSD = sq(nanmean(eegByCondStd(:,isBetween(xTimes, [tWindows(i) tWindows(i+1)]),:,:,:,30),2));
    
    windowSDAnova{i} = rmanova(meanWindowSD,{'pp','rew','distr','drug'},'categorical',2:4);
    
    windowSDPVals(:,i) = windowSDAnova{i}.pValue(2:end);

end

figure();
imagep(windowSDPVals, windowSDAnova{1}.Term(2:end));
set(gca,'XTick',1:2:length(tWindows)-1, 'XTickLabels', tWindows(1:2:end))

%% do diff waves for rew, drug, rew*drug

tInds = isBetween(xTimes, timesToPlot);

chanInd = ismember(chanNames, 'Cz');
% average over tr+distr together 
eegMeans = sq(nanmean(eegByCond2(:,tInds,:,:,:,chanInd,:),[7 4])); %[pp T rew drug]

dims = [4;3]; % leave rew, drug
eegDiff = NaN(nPP, sum(tInds), 2);
for i = 1:2
    eegDiff(:,:,i) = diff(sq(nanmean(eegMeans,dims(i))),[],3); % [pp T ch]
end

%%
effNames = {'Incentive Effect', 'THP Effect'};
figure();
c = get(gca,'ColorOrder');
set(gca,'ColorOrder', c([3,4,5],:),'nextplot','replacechildren');
h = errorBarPlot(eegDiff, 'area',1,'xaxisvalues',xTimes(tInds),'alpha',.2);
hold on;
% fill windows too?
fill(toi(1,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .2, 'LineStyle','none');
fill(toi(2,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .2, 'LineStyle','none');
xlabel('time from preparation cue (ms)');
ylabel('difference wave \Delta \muV');
xline(0, '--k');
yline(0, '--k');
legend([h{:,1}], effNames,'Location','Best');

%% now plot incentive effect in each drug

eegIncentiveEffects = sq(diff(eegMeans,[],3)); % [pp T drug]
figure();
% c = get(gca,'ColorOrder');
% set(gca,'ColorOrder', c([1 2],:),'nextplot','replacechildren');
h = errorBarPlot(eegIncentiveEffects, 'area',1,'xaxisvalues',xTimes(tInds),'alpha',.2);
hold on;
% fill windows too?
fill(toi(1,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .2, 'LineStyle','none');
fill(toi(2,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .2, 'LineStyle','none');
xlabel('time from preparation cue (ms)');
ylabel('difference wave \Delta \muV');
xline(0, '--k');
yline(0, '--k');
legend([h{:,1}], drugNames,'Location','Best');


%% use pop_ploterps

% figure();
% myERPTopoPlot(EEG, grandAverage, drugNames, timesToPlot, 0)


%% take mean amplitudes and analyse grand averages

iChan = 30;
% look before the target appears
amplWindow = toi(1,:); % get amplitude around response time
% [pp times rew dist drug chan tr]
meanAmplTrials = permute(GetMeanAmpl(eegByCond2, xTimes, amplWindow),[1 3 4 5 6 7 2]); % mean within window per trial
meanAmpl = nanmean(meanAmplTrials,6); % mean over trials [pp rew distr drug chan]


%% bar plot means

figure();
if plotDistr
%     subplot(2,2,3);
%     include targ/distr
    h = errorBarPlot(reshape(meanAmpl(:,:,:,:,iChan),nPP,4,2),'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
    JoinUpErrorPoints(h, [1 2; 3 4]);
    box off
    set(gca,'XTick',1:4,'XTickLabel',{'0p','50p'});
    ylabel('mean \muV'); 
    xlabel(['no distractor', '           ', 'with distractor']);
    xlim([.5 4.5])
    legend(h, drugNames, 'Location','Best');
else
    h = errorBarPlot(sq(nanmean(meanAmpl(:,:,:,:,iChan),3)),'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
    JoinUpErrorPoints(h, [1 2]);
    box off
    set(gca,'XTick',1:2,'XTickLabel',{'0p','50p'});
    % ylabel(sprintf('mean uV from %d:%dms', amplWindow(1), amplWindow(2))); 
    ylabel('mean \muV');
    xlabel('incentive');
    xlim([.5 2.5])
    legend(h, drugNames, 'Location','Best');
end
%% and just after (although this should probably be filtered a bit better

amplWindowPost = toi(2,:); % get amplitude around response time
meanAmplTrialsPost = permute(GetMeanAmpl(eegByCond2, xTimes, amplWindowPost),[1 3 4 5 6 7 2]); % mean within window per trial
meanAmplPost = nanmean(meanAmplTrialsPost,6); % mean over trials [pp rew distr drug chan]



%% bar plot means


figure();
if plotDistr
%     subplot(2,2,4);
    % % include targ/distr
    h = errorBarPlot(reshape(meanAmplPost(:,:,:,:,30),nPP,4,2),'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
    JoinUpErrorPoints(h, [1 2; 3 4]);
    box off
    set(gca,'XTick',1:4,'XTickLabel',{'0p','50p'});
    ylabel('mean \muV'); 
    xlabel(['no distractor', '           ', 'with distractor']);
    xlim([.5 4.5])
%     legend(h, drugNames, 'Location','Best');
else
    h = errorBarPlot(sq(nanmean(meanAmplPost(:,:,:,:,iChan),3)),'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
    JoinUpErrorPoints(h, [1 2]);
    box off
    set(gca,'XTick',1:2,'XTickLabel',{'0p','50p'});
    % ylabel(sprintf('mean uV from %d:%dms', amplWindow(1), amplWindow(2))); 
    ylabel('mean \muV');
    xlabel('incentive');
    xlim([.5 2.5])
    % legend(h, drugNames, 'Location','Best');
end

%% do rmanova at each time point

tWindows = timesToPlot(1):100:timesToPlot(2);
for i = 1:length(tWindows)-1
    disp(tWindows(i));

    % use this if doing multiple channels
%     meanWindow = nanmean(permute(GetMeanAmpl(eegByCond2, times, [tWindows(i) tWindows(i+1)]),[1 3 4 5 6 7 2]),6); % mean within window per trial
    
    % quicker for one channel
    meanWindow = sq(nanmean(nanmean(eegByCond2(:,isBetween(xTimes, [tWindows(i) tWindows(i+1)]),:,:,:,iChan,:),2),7));
    
    windowsAnova{i} = rmanova(meanWindow,{'pp','rew','distr','drug'},'categorical',2:4);
    
    windowPVals(:,i) = windowsAnova{i}.pValue(2:end);

end

figure();
imagep(windowPVals, windowsAnova{1}.Term(2:end));
set(gca,'XTick',1:2:length(tWindows)-1, 'XTickLabels', tWindows(1:2:end))

%% mean-wise GLM
% 
% % PSN ~ drug * rew + (1 | pp)
% 
% t = dePivot(meanAmpl(:,:,:,:,iChan)); % get factors
% 
% % zscore, apart from DV
% t(:,2:end-1) = nanzscore(t(:,2:end-1)); % zscore IVs, not DV
% 
% tNames = {'pp','rew','distr','drug','y'}; % y is DV
% t = array2table(t,'VariableNames',tNames);
% 
% % can add random effect of trial too: (1 | trial)
% formula = 'y ~ 1 + rew*drug*distr + (1 | pp)'; % full model plus random intercept by participant
% psnReg = fitglme(t, formula,'DummyVarCoding','effects'); % run with effects (so coefficients within predictor sum to zero)
% % psnReg(:,i) = psnReg{i}.Coefficients.pValue(2:end);
% 

%% do trial wise GLM


t = double(dePivot(sq(meanAmplTrials(:,:,:,:,iChan,:)),'KeepNaN',1)); % get factors

% zscore, apart from DV
t(:,2:end) = nanzscore(t(:,2:end)); % zscore IVs, not DV

tNames = {'pp','rew','distr','drug','trial','y'}; % y is DV
t = array2table(t,'VariableNames',tNames);

% can add random effect of trial too: (1 | trial)
formula = 'y ~ 1 + rew*drug*distr + (1 | pp)'; % full model plus random intercept by participant
psnTrialReg = fitglme(t, formula,'DummyVarCoding','effects') % run with effects (so coefficients within predictor sum to zero)


t2 =  dePivot(sq(meanAmplTrialsPost(:,:,:,:,iChan,:)),'KeepNaN',1); % get psn

t.y2 = double(nanzscore(t2(:,end)));
formula = 'y2 ~ 1 + rew*drug*distr + (1 | pp)'; % full model plus random intercept by participant
psnTrialReg2 = fitglme(t, formula,'DummyVarCoding','effects') % run with effects (so coefficients within predictor sum to zero)

% % add in the fixation location
% fixLocTT = permute(groupMeans(fixLoc,1,trialType,'dim'),[2,1,3,4]);
% fixLocTT = reshape(fixLocTT,nPP,2,2,2,[]); % [pp rew dist drug trial]
% tFixLoc = dePivot(fixLocTT,'KeepNaN',1);
% 
% tFixLoc = nanzscore(tFixLoc(:,end)); % zscore the fixation location
% 
% t = horzcat(t, array2table(tFixLoc, 'VariableNames', {'fixLoc'})); % add into table
% 
% formula = 'y ~ 1 + rew*drug*fixLoc + (1 | pp)'; % full model plus random intercept by participant
% psnTrialRegFixLoc = fitglme(t, formula,'DummyVarCoding','effects'); % run with effects (so coefficients within predictor sum to zero)
% 

% could add in baseline here like that GLM baselining method?


%% trialwise in bins

tWindows = timesToPlot(1):100:timesToPlot(2);

clearvars psnTrialWindowsB psnTrialWindowsP;
formula = 'y3 ~ 1 + rew*drug*distr + (1 | pp)'; % full model plus random intercept by participant
for i = 1:length(tWindows)-1
    disp(tWindows(i));

    meanWindow = sq(nanmean(eegByCond2(:,isBetween(xTimes, [tWindows(i) tWindows(i+1)]),:,:,:,iChan,:),2));
    
    t.y3 = nanzscore(col(meanWindow));
    psnTrialWindows{i} = fitglme(t, formula);
    
    psnTrialWindowsP(:,i) = psnTrialWindows{i}.Coefficients.pValue(2:end);
    psnTrialWindowsB(:,i) = psnTrialWindows{i}.Coefficients.Estimate(2:end);
end

figure();
imagep(psnTrialWindowsP, psnTrialWindows{1}.Coefficients.Name(2:end));
set(gca,'XTick',1:2:length(tWindows)-1, 'XTickLabels', tWindows(1:2:end))


%% topoplot each effect in CNV window for all electrodes


for iCh = 1:61
    t.chan =  nanzscore(col(meanAmplTrialsPost(:,:,:,:,iCh,:))); % get psn
    
    psnTrialRegChans = fitglme(t, 'chan ~ 1 + rew*distr*drug + (1 | pp)');
    
    psnTrialChansP(:,iCh) = psnTrialRegChans.Coefficients.pValue(2:end);
    psnTrialChansB(:,iCh) = psnTrialRegChans.Coefficients.Estimate(2:end);
end

%%
cmap = othercolor('PRGn11');
mapLimits = [-.05 .05];
alpha = .05 /61;
figure();
for i = 1:7
    subplot(2,4,i);
    topoplot(psnTrialChansB(i,:), chanLocs(1:61),...
            'mapLimits',mapLimits, 'colormap',cmap,...
            'style', 'both', 'plotrad',.55, 'headrad', .5,'emarker', {'.','k',[],1},...
            'numcontour', 6,'electrodes', 'on', 'nosedir', '+X',...
            'emarker2', {find(psnTrialChansP(i,:)<alpha),'o','y',2});
    title(psnTrialRegChans.CoefficientNames(i+1));
    colorbar;

end

drawnow;
%% do permutation and mass univar tests on grand averages
doUnivars = 1;
GNDEffects = {};
if doUnivars
    
    plotUnivars = 0; % don't plot univariate stuff
    
%     eegByCond2Trimmed = reshape(permute(eegWindowMean,[4,2,3,5,1,6]), nPP,[],2,2,2,nChans,nTr);
%     EEG.times = midWindows;
    eegByCond2Trimmed = eegByCond2(:,isBetween(xTimes, timesToPlot),:,:,:,:,:);
    
%     % smooth it?
%     eegByCond2Trimmed = movmean(eegByCond2Trimmed,5, 2);
    
    EEG.times = xTimes(isBetween(xTimes, timesToPlot));
    if ~isfield(EEG,'chanlocs') || isempty(EEG.chanlocs)
        EEG.chanlocs = chanLocs;
    end

    avDims = {  [7 4 5];    % tr, dist, drug. to get reward effect
                [7 3 5];    % tr, rew, drug. dist effect
                [7 3 4];    % tr, rew + dist. drug effect
                [7 5];      % tr, drug. rew*dist [targLo targHi distLo distHi]
                [7 4];      % tr, dist. rew*drug [lowPlac highPlac lowDrug highDrug]
                [7 3];      % tr, rew. dist*drug [targPlac distPlac targDrug distDrg]
                };

    drugDiffWaveInds = {[2 1]; % rew effect: hi - low
                        [2 1]; % dist effect: dist-targ
                        [2 1]; % drug effect: drug-placebo
                        [2 1; 4 3; 6 5]; % rew*dist: [targHi-targLow, distHi-distLow, 6-5]
                        [2 1; 4 3; 6 5]; % rew*drug: [placHi-low, drugHi-Low, 6-5]
                        [2 1; 4 3; 6 5]; % distr*drug: [placDistr-Targ, drugDistr-targ, 6-5]
                        };

    nPermsEff = 2500;
    nPerms = {  nPermsEff;
                nPermsEff;
                nPermsEff;
                [0 0 nPermsEff];
                [0 0 nPermsEff];
                [0 0 nPermsEff];
                };
    drugPermNames = {   {'low','high'};
                        {'targ','dist'};
                        {'Placebo','Drug'};
                        {'targLow','targHi','distLow','distHi'};
                        {'lowDrug', 'lowPlac', 'highDrug', 'highPlac'};
                        {'targDrug', 'targPlac', 'distDrug', 'distPlac'};
                        };

    GNDEffects = cell(size(avDims,1),1);
    if plotDistr; inds = 1:6; % inlcude distr
    else inds = [1 3 5]; end % skip distr
    for i = inds%1:size(avDims,1) % for each comparison
        % average, combining back over factors
        m = reshape(sq(nanmean(eegByCond2Trimmed, avDims{i})), nPP, length(EEG.times), [], nChans);

        % also mean number of good and mean
        gTT = reshape(nanmean(reshape(filtersByCondSummary.toInclude,nPP,2,2,2), avDims{i}(2:end)-1),nPP,[]);
        bTT = reshape(nanmean(reshape(filtersByCondSummary.toExclude,nPP,2,2,2), avDims{i}(2:end)-1),nPP,[]);
    %     disp(size(gTT));

        GNDEffects{i} = myGNDTesting(m, 1:61, EEG, drugPermNames{i},...
                        drugDiffWaveInds{i},gTT, bTT, nPerms{i}, suffix, [], plotUnivars);

    end


end

%% ceiling effects?

cnv = reshape(permute(meanAmplTrialsPost(:,:,:,:,30,:),[2,4,1,6,3,5]),2,2,nPP,[]); 
%[rew drug nPP [distr/trials]]

figure();

for j = 1:2
    for i = 1:2
        subplot(2,2,(i-1)*2+j);
    %     plot(reshape(cnv(:,i,:,1),2,[]),'-x');
    %     set(gca,'XTick',1:2, 'XTickLabel', {'0p','50p'});
        hist(col(cnv(j,i,:,:)),100);
        xlabel('CNV');
        
        title(condNames2( (i-1)*2+j) );
        xlim([-50 50]);
        xline(0, ':k');
    end
end

makeSubplotScalesEqual(2,2);

%% predict beh?

load('TrioEEGBehAnalysis.mat','trialTab','varNames','ylabels');
vNames = {'residual velocity','RT','distractor pull'};
nVars = 3;

trialTab.CNV = double(nanzscore(col(meanAmplTrialsPost(:,:,:,:,30,:))));
trialTab.P3a = double(nanzscore(col(meanAmplTrials(:,:,:,:,30,:))));
toiNames = {'P3a','CNV'};

extractVar = @(x,inds) [x.Coefficients.Estimate(inds), x.Coefficients.SE(inds), x.Coefficients.tStat(inds), x.Coefficients.pValue(inds), x.Coefficients.Lower(inds), x.Coefficients.Upper(inds)];

formula = '%s ~ 1 + rew*distr*drug + %s + (1 | pp)';
toiRegStats = zeros(3,2,8,6);
for i = 1:2
    for j = 1:3
        reg = fitglme(trialTab, sprintf(formula, varNames{j}, toiNames{i}));
        
        toiRegStats(j,i,:,:) = extractVar(reg, 2:9);
    end
end


termNames = reg.CoefficientNames(2:end);
%%

alpha = .05 / 6;

f = figure();
set(f, 'DefaultAxesFontSize',16)
h = barwitherr(toiRegStats(:,:,4,2)*1.96, toiRegStats(:,:,4,1));
h(1).FaceColor = [0.4660    0.6740    0.1880];
h(2).FaceColor = [0.3010    0.7450    0.9330];
[h.LineWidth] = deal(1.2);
xticklabels(vNames(1:nVars));
ylabel('\beta coefficient')
box off

for i = 1:2

    
    stars = col(toiRegStats(:,i,4,4)') < alpha;
%     y = repmat(max([h.YData])*1.1, 1,3);    
    y = repmat(-.04, 1,3);
    y(~stars) = NaN;
    hold on;
    plot((1:3) + .3*i - .45, y, '*','MarkerSize',10, 'Color', h(i).FaceColor, 'LineWidth',1.5);
    
end

legend(h, toiNames, 'Location','Best');
%% save

save(sprintf('trioEEGAnalyseByTrials_%s.mat', timelock), '-v7.3',...
    'eegRawData','trialType','filters','xTimes','GNDEffects','amplWindow',...
    'suffix','doBaseline','chanNames','EEG','eegByCond','condNames','nPP',...
    'nBins','timelock','meanAmplTrialsPost','meanAmplTrials','timesToPlot',...
    'toi','plotDistr')

%% save a smaller version

% trim times?
if strcmp(timelock, 'audCue_200')
    toi = [900 1100]; % just this
    tInds = isBetween(xTimes, toi);
    eegByCond = eegByCond(:,tInds,:,:,:,:,:);
    xTimes = xTimes(tInds);
end

save(sprintf('trioEEGAnalyseByTrials_%s_small.mat', timelock), '-v7.3',...
    'xTimes','eegByCond','timesToPlot','toi','chanNames','GNDEffects')

end