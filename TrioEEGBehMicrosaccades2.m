% TrioEEGBehMicrosaccades2
% check for corrleations between microsaccade number / drift speed and the
% 3 behavioural metrics
% and maybe also with CNV

clc;clear all; close all;
%% correlate with behaviour

% load up beh
load('TrioEEGBehAnalysis.mat', 'vars','varNames'); %[pp cond4 drug tr]

%% load up ready2, which has distr
load(sprintf('TrioEEGBehMicrosaccades_%s.mat','ready2Stats'))

% take mean speed over period
meanSpeed = sq(nanmean(ready2Stats.meanSpeedCond,2))*1000; % [pp cond4 drug tr]
% meanDens = ready2Stats.meanDens; % pp rew distr drug

[nPP, nConds, nDrugs, nTr] = size(meanSpeed);

%% load each pp's #microsaccades

load('TrioEEGBehReady2MicroFixation.mat','O')

%% extract
nMicro = NaN(nPP,nConds,nDrugs,nTr);
for i = 1:20
    for j = 1:2
        n = nancat(O(i,j).nMicro); %[tr 1 cond]
        nMicro(i,:,j,1:size(n,1)) = sq(n)';
    end
end

%% filter out excluded trials here

load('./trioEEGAnalyseByTrials_readyBS.mat', 'filters','trialType');

% filters.toExclude is [nTr, nPP, nDrugs]
filtersByCond = structfun(@(x) permute(groupMeans(x, 1, trialType,'dim'),[2,1,3,4]), filters,'UniformOutput',0); 
%[pp cond drug tr]

meanSpeed(filtersByCond.toExclude==1) = NaN;
nMicro(filtersByCond.toExclude==1) = NaN;

for iV = 1:3
    vars.(varNames{iV})(filtersByCond.toExclude==1) = NaN;
end



%% run one regression to check it matches previous reg

regTab = struct2table(structfun(@col, keepFields(vars, varNames(1:3)), 'Uni',0)); 
behTab = array2table(nanzscore(dePivot(reshape(vars.velres,nPP,2,2,2,nTr), 'KeepNaN',1)), 'VariableNames', {'pp1','rew','distr','drug','trial','velres'});
regTab = horzcat(behTab(:, 1:end-1), regTab); % drop velres here
regTab.meanSpeed = col(meanSpeed);
regTab.nMicro = col(nMicro); % replace

% % zscore DVs?
regTab.meanSpeed_z = nanzscore(regTab.meanSpeed);
regTab.nMicro_z = nanzscore(regTab.meanSpeed);
regTab.velres = nanzscore(regTab.velres);
regTab.srt = nanzscore(regTab.srt);
regTab.depAngle= nanzscore(regTab.depAngle);
% 
% %

fitglme(regTab, 'velres ~ 1 + rew*distr*drug + (1 | pp1)')% check matches, it does!

% %% loop and check var ~ speed
% 
% % first test mean s
% formula = '%s ~ 1 + meanSpeed_z + (1 | pp1)';
% for iV = 1:3
%     f{iV} = fitglme(regTab, sprintf(formula, varNames{iV}));
% end
% 
% %% now while controlling
% 
% formula = '%s ~ 1 + meanSpeed_z + rew*distr*drug + (1 | pp1)';
% for iV = 1:3
%     f2{iV} = fitglme(regTab, sprintf(formula, varNames{iV}));
% end
% 

%% load CNV

load('trioEEGAnalyseByTrials_readyBS.mat','meanAmplTrialsPost')

%% get CNV channel

regTab.cnv = col(meanAmplTrialsPost(:,:,:,:,30,:));

fitglme(regTab, 'cnv ~ 1 + meanSpeed_z +  (1 | pp1)') % no
fitglme(regTab, 'cnv ~ 1 + meanSpeed_z + rew*distr*drug + (1 | pp1)') % no
fitglme(regTab, 'cnv ~ 1 +  rew*distr*drug + (1 | pp1)') % no


fitglme(regTab, 'meanSpeed ~ 1 + cnv + rew*distr*drug + (1 | pp1)')


%% load all times+channels
% this is the same as Figure 4, but controlling for ocSpeed/nMicro (without
% p-value perm testing as it takes ages)

timelock = 'readyBS';

load(sprintf('trioEEGAnalyseByTrials_%s.mat',timelock), 'eegByCond','xTimes')
[nChans, nTimes, nBins, nPP, nDrugs, nTr] = size(eegByCond);

timesToPlot = [-200 1500];
xTimes2 = xTimes(isBetween(xTimes, timesToPlot));

eeg = reshape(permute(eegByCond,[4,3,5,6,2,1]),[],nTimes,nChans);
eeg = eeg(:,isBetween(xTimes, timesToPlot),:); %[nPP*nConds*nDrugs*nTr, nT, nChans]

% % % this bit will get the absolute movement of eye pos (from baseline)
% eeg = eeg(:,isBetween(xTimes, timesToPlot),62:65); % just h/v-eog + el-X/y
% xy = complex(eeg(:,:,3), eeg(:,:,4)); 
% eeg(:,:,end+1) = abs(xy); % magnitude of movement
% eeg(:,:,end+1) = angle(xy); % angle
% chanNames2 = [chanNames(62:65), 'el-abs', 'el-angle'];

% get absolute eye-pos (from baseline) from eye-tracker x+y coords
eyePos = abs(complex(eeg(:,:,64), eeg(:,:,65))); % abs distance
eeg = eeg(:,:,1:61); % remove eog + eyetracker channels
eeg = nanzscore(eeg);


% get beh data
load('TrioEEGBehAnalysis.mat', 'trialTab','varNames','nVars');
trialTab.meanSpeed = nanzscore(col(meanSpeed)); % zscore
trialTab.nMicro = nanzscore(col(nMicro)); % zscore

% remove nan rows
toRemove = all(isnan(eeg),[2 3]);
eeg(toRemove,:,:) = [];
trialTab(toRemove,:) = [];
eyePos(toRemove,:) = [];
eyePos = nanzscore(eyePos); % zscore, [nTr, nT] % save for later


%% run one regression as a test
[~,nT,nCh] = size(eeg);
iV = 1;
iM = 1;

measures = {'meanSpeed', 'nMicro','eyePos'};
formula = sprintf('%s ~ 1 + %s + v + rew*distr*drug + (1 | pp)', varNames{iV}, measures{iM});

% run one to find ind
trialTab.v = eeg(:,1,1);
reg = fitglme(trialTab, formula);
disp(reg);
inds = ismember(reg.CoefficientNames, {'v', measures{iM}}); % store both of these
nInds = sum(inds);

%% now run all channels/times, per variable and control

% runOneSet is a subfunction 

for iV = 1:3
    for iM = 1:2 % not eyepos
        if exist(sprintf('%s_regress_eeg_%s_%s.mat', timelock, varNames{iV}, measures{iM}),'file'); continue; end
        formula = sprintf('%s ~ 1 + %s + v + rew*distr*drug + (1 | pp)', varNames{iV}, measures{iM});
        disp(formula);

        [betaVals, pVals, tVals] = runOneSet(trialTab, eeg, formula, nT, nCh, nInds, inds); % 7 mins per

        save(sprintf('%s_regress_eeg_%s_%s.mat', timelock, varNames{iV}, measures{iM}), ...
             'betaVals','pVals','tVals','formula')

    end
end

%% also do it with eyePos at each time-sample - runs a little differently

[~,nT,nCh] = size(eeg);
iV = 1;
iM = 3; % eyePos

% measures = {'eyePos'};
formula = sprintf('%s ~ 1 + %s + v + rew*distr*drug + (1 | pp)', varNames{iV}, measures{iM});

% run one to check
trialTab.v = eeg(:,1,1);
trialTab.eyePos = eyePos(:,1);
reg = fitglme(trialTab, formula); 
disp(reg);
inds = ismember(reg.CoefficientNames, {'v', measures{iM}}); % store both these terms
nInds = sum(inds);

%%

for iV = 1:3
    for iM = 3
        if exist(sprintf('%s_regress_eeg_%s_%s.mat', timelock, varNames{iV}, measures{iM}),'file'); continue; end
        formula = sprintf('%s ~ 1 + %s + v + rew*distr*drug + (1 | pp)', varNames{iV}, measures{iM});
        disp(formula);

%         [betaVals, pVals, tVals] = runOneSet(trialTab, eeg, formula, nT, nCh, nInds, inds); % 7 mins per
        [betaVals, pVals, tVals] = deal(single(NaN(nT, nCh, nInds)));
        tic;
        parfor iCh = 1:nCh
            disp(iCh);
            for iT = 1:nT
                trialTab1 = trialTab;
                trialTab1.v = eeg(:,iT,iCh);
                trialTab1.eyePos = eyePos(:,iT);% insert this time's eyePos
                reg = fitglme(trialTab1, formula); % ~300ms per
        
                %         ind = strcmpi(reg.CoefficientNames,'rew:drug');
        
                betaVals(iT,iCh,:) = reg.Coefficients.Estimate(inds);
                pVals(iT,iCh,:) = reg.Coefficients.pValue(inds);
                tVals(iT,iCh,:) = reg.Coefficients.tStat(inds);
            end
        end
        toc

        save(sprintf('%s_regress_eeg_%s_%s.mat', timelock, varNames{iV}, measures{iM}), ...
             'betaVals','pVals','tVals','formula')

    end
end

%% also run without covars - same data as in Figure 4
% can load up if already done
% 
% formula = sprintf('%s ~ 1 + v + rew*distr*drug + (1 | pp)', varNames{iV});
% 
% % run one to find ind
% trialTab.v = eeg(:,1,1);
% reg = fitglme(trialTab, formula);
% disp(reg);
% inds = ismember(reg.CoefficientNames, {'v', measures{iM}});
% nInds = sum(inds);
% 
% 
% for iV = 1:3
%     
%     if exist(sprintf('%s_regress_eog_%s.mat', timelock, varNames{iV}),'file'); continue; end
%     formula = sprintf('%s ~ 1 + v + rew*distr*drug + (1 | pp)', varNames{iV});
%     disp(formula);
%     
%     [betaVals, pVals, tVals] = runOneSet(trialTab, eeg, formula, nT, nCh, nInds, inds); % 7 mins per
%     
%     save(sprintf('%s_regress_eog_%s.mat', timelock, varNames{iV}), ...
%          'betaVals','pVals','tVals','formula')
% end

%% load up all
measures = {'meanSpeed', 'nMicro','eyePos'};

[allBVals, allPVals] = deal(single(NaN(nT,nCh,nInds,3,3)));
for iV = 1:3
    for iM = 1:length(measures)
        d = load(sprintf('%s_regress_eeg_%s_%s.mat', timelock, varNames{iV}, measures{iM}), ...
             'betaVals','pVals');
        allBVals(:,:,:,iV,iM) = d.betaVals;
        allPVals(:,:,:,iV,iM) = d.pVals;
    end
end
ind = 2; % 1=speed, 2=v
alpha = .05 ./ numel(allPVals); % bonferroni

% terms are wrong way round for eyePos, flip them
allBVals(:,:,[1 2],:,3) = allBVals(:,:,[2 1],:,3);
allPVals(:,:,[1 2],:,3) = allPVals(:,:,[2 1],:,3);


%% plot

topoTimes = [-100 300 650 1000 1350];
timelockName = 'preparation cue onset'; % time from...

if ~exist('othercolor','file')
    addpath '../../../../Documents/TCD/GeneralScripts/MatlabPackages/othercolor/';
end
cmap = othercolor('PRGn11');


mapLims = [-.15 .15] ;
chansToPlot = {'FPz','Cz','POz'};
yLabs = {'residual peak velocity','RT','distractor pull','endptVarCorr'};


% figure();
% close all;
set(0, 'DefaultAxesFontSize',14);
for iV = 1:3
    for iM = 1:length(measures)
        if iV==1
            cbarArgs = {'Position',[0.0625 0.1952 0.0232 0.1905]};
        else
            cbarArgs = {};
        end
        f = figure();
        set(f, 'DefaultAxesFontSize',20);
    
        
    
        PlotTopoAndTrace(allBVals(:,:,ind,iV,iM)', allPVals(:,:,ind,iV,iM)',xTimes2,...
            'topoTimes', topoTimes, 'chansToPlot',chansToPlot,...
            'alpha', -1, 'mapLimits',mapLims, 'cmap', cmap,...
            'yLab', '\beta coefficient', 'title', sprintf('%s', yLabs{iV}),...
            'doLegend',iV==1, 'cbarArgs',cbarArgs, 'subplotInds',[2 length(topoTimes) 1],...
            'xLab', sprintf('time from %s (ms)', timelockName), 'sepFigs',0);
        box off;
    
        % save figures
        for j = 1:length(f)
%             saveas((f(j)), sprintf('./Figs/PlotRegsFE_covars_%s_%s_%s.svg', timelock, varNames{iV}, measures{iM}))
        end
    end
end


%% diff from true regressions

timesToPlot = [-100 1500]; % this was the window used for this
tInds = isBetween(xTimes,timesToPlot);
trueBetas = single(NaN(sum(tInds),61,3));
for iV = 1:3
    a = load(sprintf('%s_regressCluster3_%s.mat',timelock, varNames{iV}),'betaVals');
    trueBetas(:,:,iV) = a.betaVals(1:410,:);
end

tInds2 = isBetween(xTimes2, timesToPlot); % trim to common times
diffBetas = sq(allBVals(tInds2,:,ind,:,:)) - trueBetas; %[T ch vars meas]
mapLims = [-.02 .02] ;

% set(0, 'DefaultAxesFontSize',14);
for iV = 1:3
    for iM = 1:length(measures)
        if iV==1
            cbarArgs = {'Position',[0.0625 0.1952 0.0232 0.1905]};
        else
            cbarArgs = {};
        end
        f = figure();
        set(f, 'DefaultAxesFontSize',20);
    
        
    
        PlotTopoAndTrace(diffBetas(:,:,iV,iM)', allPVals(tInds2,:,ind,iV,iM)',xTimes2(tInds2),...
            'topoTimes', topoTimes, 'chansToPlot',chansToPlot,...
            'alpha', -1, 'mapLimits',mapLims, 'cmap', cmap,...
            'yLab', '\beta coefficient', 'title', sprintf('%s', yLabs{iV}),...
            'doLegend',iV==1, 'cbarArgs',cbarArgs, 'subplotInds',[2 length(topoTimes) 1],...
            'xLab', sprintf('time from %s (ms)', timelockName), 'sepFigs',0);
        box off;
    
        % save figures
        for j = 1:length(f)
%             saveas((f(j)), sprintf('./Figs/PlotRegsFE_covarsDiff_%s_%s_%s.svg', timelock, varNames{iV}, measures{iM}))
        end
    end
end


%% functions



function [betaVals, pVals, tVals] = runOneSet(trialTab, eeg, formula, nT, nCh, nInds, inds)

    [betaVals, pVals, tVals] = deal(single(NaN(nT, nCh, nInds)));
    tic;
    parfor iCh = 1:nCh
        disp(iCh);
        for iT = 1:nT
            trialTab1 = trialTab;
            trialTab1.v = eeg(:,iT,iCh);
            reg = fitglme(trialTab1, formula); % ~300ms per
            
    %         ind = strcmpi(reg.CoefficientNames,'rew:drug');
            
            betaVals(iT,iCh,:) = reg.Coefficients.Estimate(inds);
            pVals(iT,iCh,:) = reg.Coefficients.pValue(inds);
            tVals(iT,iCh,:) = reg.Coefficients.tStat(inds);
        end
    end
    toc

end