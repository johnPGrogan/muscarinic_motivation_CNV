function output = TrioEEGBehAnalysis(runProcess)
% function output = TrioEEGBehAnalysis()
% Process each datafile if not already done. Combine group data, process
% and do quick plots and stats.
% Inputs: runProces: 1 to reprocess datafiles, 0 to load (default)
% 

if ~exist('runProcess','var') 
    runProcess = 0;
end
close all

%% find files
behDataFolder = 'D:/groganj/Work/Trihexy/BehData/'; % change this to where you store your data

files = what(behDataFolder); % get files

% function to find file names
cellregexpi = @(cellArray, pattern) ~cellfun(@isempty, regexpi(cellArray, pattern));

files.mat = files.mat(cellregexpi(files.mat, '^TrioRewardCue3_EEG_TRI')); % file name pattern
files.mat = files.mat(~cellregexpi(files.mat, '_0.mat')); % remove practice versions

files.mat = files.mat(~cellregexpi(files.mat, 'TRI_015')); % exclude TRI_015
files.mat = files.mat(~cellregexpi(files.mat, 'TRI_024|TRI_027|TRI_031|TRI_033|TRI_036')); % remove incompletes

nFiles = length(files.mat);

ppInfo = cellfun(@(x) x(end-12:end-4), files.mat, 'UniformOutput',0); % full info
ppID = cellfun(@(x) x(1:7), ppInfo,'UniformOutput',0); % just ppID
sesNo = cell2mat(cellfun(@(x) str2num(x(end)), ppInfo, 'UniformOutput',0)); % 1 or 2

%% process each file to extract data

if runProcess==0 && exist( 'TrioEEGBehProcessFiles.mat', 'file')
    disp('loading previous processed file')
    load('TrioEEGBehProcessFiles.mat');
else
    disp('previous processed file not found. loading and processing each mat/edf file now...')
    
    useParFor = 1; % use parfor or for?
    result = TrioEEGBehProcessFiles(files.mat, behDataFolder, useParFor);
end

%% get trial types
 
tt = nancat(2, result.trialType);
exclude = nancat(2, result.exclude);

tt(exclude) = NaN; % remove ones where R < 0 i.e. incomplete, paused, aborted etc

%% split by pp x sess

ttSes = groupMeans(tt,2,sesNo,'dim'); % split into ses 1 or 2
% [pp, session, trial]

%% recode into A/B
% load drug conds
[~,~,raw] = xlsread( '../../Experiments/Forms/ACh/TrihexyDrugConds.xlsx');

%% remove NaN rows
raw(cell2mat(cellfun(@(x) any(isnan(x)), raw(:,1),'UniformOutput',0)),:) = [];

xlsID = raw(:,1); % get IDs
uniqID = unique(ppID); % unique only

% match
[~,xlsInds] = ismember(uniqID, xlsID);

drugConds = raw(xlsInds, 3:4);
drugA = double(strcmp(drugConds(:,1), 'A'));
% NaNs
drugA(cell2mat(cellfun(@(x) any(isnan(x)), drugConds(:,1),'UniformOutput',0))) = NaN;


% flip if drug was first: [placebo drug]
ttCond = ttSes;
ttCond(drugA==1,:,:) = ttSes(drugA==1,[2 1],:); % [pp, plac/drug, trial]

% NaNs
ttCond(isnan(drugA),:,:) = NaN;

%% extract vars and split by trial type

vars = struct();%preset
varNames = result(1).varNames;
varNames = [varNames, {'endpt','correct'}]; % also get these
for i = 1:length(varNames)
    v = nancat(2, result.(varNames{i})); % extract into nTr x nFiles
    v = groupMeans(v, 2, sesNo,'dim'); % [pp ses tr]
    
    v(drugA==1,:,:) = v(drugA==1,[2 1],:); % switch their conds round
    
    vars.(varNames{i}) = permute(groupMeans(v, 3, ttCond, 'dim'), [1,3,2,4]); %[pp tt ses tr]
end

% struct2workspace(vars);
[nPP, ~, nSes, nTr]  = size(vars.rt);

%% get residual velocity
% for each pp, ses, separately

vel = groupMeans(nancat(2,result.vel),2,sesNo,'dim');
vel(drugA==1,:,:) = vel(drugA==1,[2 1],:); %[placebo drug]
ampl = groupMeans(nancat(2,result.ampl),2,sesNo,'dim');
ampl(drugA==1,:,:) = ampl(drugA==1,[2 1],:); % [placebo drug]
nTrAll = size(ampl,3);

velres = NaN(nPP,nSes,nTrAll);
for iPP = 1:nPP
    for iSes = 1:2
        [~,~,velres(iPP,iSes,:)] = regress( sq(vel(iPP,iSes,:)),  [ ones(nTrAll,1)  sq(ampl(iPP,iSes,:)) ] );
    end
end

vars.velres = permute(groupMeans(velres, 3, ttCond, 'dim'), [1,3,2,4]); %[pp tt ses tr]
varNames = ['velres', varNames];

%% get endptVar for correct trials only

endpts = vars.endpt;
endpts(vars.correct == 0) = NaN; % remove incorrect
vars.endptVarCorr = abs(endpts - nanmean(endpts, [2 4])); % abs diff from mean within exp session

%% log/transform variables

vars.srt= log(vars.srt);
vars.rt = log(vars.rt);

vars.endptVar = log(vars.endptVar);
vars.endptVarCorr = log(vars.endptVarCorr);

%% set vars to analyse
varNames = {'velres','srt','depAngle','endptVarCorr'};%,'endptVar','vel', 'ampl','rt', 'distrCapture', 'maxDev','correct'};

%% plot means
figure();
lineColours = get(gca,'ColorOrder');
labels = {'0p','50p'; 'no distractor', 'with distractor'; 'Placebo','THP'};
ylabels = {'residual peak velocity (deg/s)','saccadic RT (ms)','distractor pull (deg)','endpt var corr log(deg)',...
    'endpt var log(deg)','velocity (deg/s)', 'amplitude (deg)', 'RT (ms)',...
     'prop capture', 'max deviation (deg)', 'prop correct'};

 
plotIndivs = 0;

nVars = length(varNames);
for i = 1:nVars

    % repeated measures anova
    meanVar = reshape(nanmean(vars.(varNames{i}),4), [],2,2,2);
    [stats{i}, regStats{i}] = rmanova(meanVar, {'pp','rew','distr','drug'},'DummyVarCoding','effects','categorical',2:4);
    pVals(:,i) = stats{i}.pValue(2:end);
    
    
    figure();
%     subplot(ceil(sqrt(nVars)),ceil(nVars/ceil(sqrt(nVars))),i);
    h = errorBarPlot(nanmean(vars.(varNames{i}),4), 'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
    hold on;
    for j = 1:2
        plot(h(j).XData(1:2), h(j).YData(1:2), 'LineWidth', 2, 'Color', lineColours(j,:));
        plot(h(j).XData(3:4), h(j).YData(3:4), 'LineWidth', 2, 'Color', lineColours(j,:));
    end
    
    box off
    set(gca,'XTick',1:4,'XTickLabel',labels(1,:));
    ylabel(ylabels{i}) 
    xlabel([labels{2,1} '       ' labels{2,2}]);
    xlim([.5 4.5])
    
    
    if plotIndivs
        for j = 1:2
            x = (1:4) + (rand(20,4)-.5).*.1 + (3 - 2*j)/20;
            plot(x, nanmean(vars.(varNames{i})(:,:,j,:),4),'x','Color',h(j).Color,'MarkerSize',10);
        end
    end
%     if pVals(1,i)<=.05; text(0.5, 0.9, 'r','Units','Normalized','FontSize',14);end
%     if pVals(2,i)<=.05; text(0.5, 0.8, 'd','Units','Normalized','FontSize',14);end
%     if pVals(3,i)<=.05; text(0.5, 0.7, 'o','Units','Normalized','FontSize',14);end
end
legend(labels(3,:),'Location','Best')

statsTab = horzcat( table(stats{1}.Term(2:end)), array2table(pVals));
statsTab.Properties.VariableNames(2:end) = varNames;

disp(statsTab);

%% trial wise glme

isLogistic = strcmp(varNames, 'distrCapture') | strcmp(varNames,'correct'); % which vars are logistic
isZScored = [false true(1,5) ~isLogistic]; % don't zscore them or pp

% make a table with col per factor and variable
t = structfun(@(x) dePivot(reshape(x,nPP,2,2,2,[]),'KeepNaN',1), vars,'UniformOutput',0); % reshape into col
disCap = t.distrCapture(:,end); % do not z score
corr = t.correct(:,end); % do not z score
t = structfun(@nanzscore, t, 'UniformOutput',0); % zscore all
% % or 
% t.velres = nanzscore(t.velres); % just the factors - will overwrite velres later

trialTab = array2table(t.velres, 'VariableNames', {'pp','rew','distr','drug','trial','y'});
t = structfun(@(x) x(:,end), t, 'UniformOutput',0);
trialTab = horzcat( trialTab(:,1:end-1), struct2table(t));

trialTab.distrCapture = disCap; % this cannot be z scored
trialTab.correct = corr; % this cannot be z scored

formula = '%s ~ 1 + rew*distr*drug + (1 | pp)'; % full model plus random intercept by participant
  
for i = 1:nVars
   
    if isLogistic(i)
        args = {'link','logit','distribution','binomial'};
    else
        args = {'link','identity','distribution','normal'};
    end
    
    regTrialStats{i} = fitglme(trialTab, sprintf(formula, varNames{i}),'DummyVarCoding','effects', args{:}); % run with effects (so coefficients within predictor sum to zero)
    regTrialPVals(:,i) = regTrialStats{i}.Coefficients.pValue(2:end);
    
end

regTrialStatsTab = horzcat( table(regTrialStats{1}.CoefficientNames(2:end)'), array2table(regTrialPVals));
regTrialStatsTab.Properties.VariableNames(2:end) = varNames;

disp(regTrialStatsTab);

%% hist sRT by correct/incorrect

for i = 1:2
    [y(i,:),x(i,:)] = ksdensity(vars.srt(vars.correct==(i-1)));
end
figure();
h = plot(x',y', '-');
xlabel(ylabels(2));
ylabel('density');
legend(h, {'Incorrect', 'Correct'}, 'Location','NorthEast');

%%

output = workspace2struct('result'); % don't save result as it is in TrioEEGBehProcessFiles.mat

save('TrioEEGBehAnalysis.mat','-v7.3','-struct','output') % save after processing



end

