function TrioEEGBehMicrosaccades(fixResult, cuesToAnalyse)
% analyse microsaccades in the ready period and early target period.
% microsaccade density, ocular drift speed, FFT of movements in early and
% late period between ready cue and target onset.
% 
% requires edft

%% load up

if ~exist('cuesToAnalyse','var') || isempty(cuesToAnalyse) 
    cuesToAnalyse = {'ready','ready2'}; % do all
end

if ~exist('fixResult','var') || isempty(fixResult)
%     load('TrioEEGBehProcessFiles.mat','result');
    load('TrioEEGBehProcessFilesFix.mat','fixResult');
end

load('TrioEEGBehAnalysis.mat','drugA','ttCond','sesNo');

saves = {}; % will add var names to save in here 

%% startviscue - startaudcue
if any(strcmpi(cuesToAnalyse, 'vis'))
    isRew = mod(ttCond,2) == 0; % 0=low, 1=high
    [infoVis, trajVis] = ExtractInfo(fixResult, 'infoVis', 'normTrajVis', drugA, sesNo);

    trajVis(:,:,:,501:end) = []; % trim to same length

    % run analyse each fixation period
    OVis = LoopFix(infoVis, trajVis, isRew, 1, [1 500], 'TrioEEGBehVisMicroFixation.mat');

    % plot
    N = 250;
    vnames = {'pp','rew','drug'}; % names of factors for rmanova
    rewNames = {'low','high'}; % no targ/dist in ready period
    t = 0:250:500; % time inds for xticks
    visStats = PlotAndStats(OVis, N, rewNames, t, vnames, infoVis);

    [visStats.nMacroCond, visStats.nMacroMeanStats] = AnalyseMacros(infoVis, isRew, t, vnames, rewNames);

    saves = [saves {'visStats'}];
end
%% startaudcue - startreadycue
if any(strcmpi(cuesToAnalyse, 'aud'))
    isRew = mod(ttCond,2) == 0; % 0=low, 1=high
    [infoAud, trajAud] = ExtractInfo(fixResult, 'infoAud', 'normTrajAud', drugA, sesNo);

    trajAud(:,:,:,1501:end) = NaN; % trim to same length

    % run analyse each fixation period
    OAud = LoopFix(infoAud, trajAud, isRew, 1, [1 1500], 'TrioEEGBehAudMicroFixation.mat');

    % plot
    N = 500;
    vnames = {'pp','rew','drug'}; % names of factors for rmanova
    rewNames = {'low','high'}; % no targ/dist in ready period
    t = 0:750:1500; % time inds for xticks
    audStats = PlotAndStats(OAud, N, rewNames, t, vnames, infoAud);

    [audStats.nMacroCond, audStats.nMacroMeanStats] = AnalyseMacros(infoAud, isRew, t, vnames, rewNames);
    
    saves = [saves {'audStats'}];
end
%% do ready period
if any(strcmpi(cuesToAnalyse, 'ready'))
    % isRew = nancat(2,result.isRew); % is high reward
    % isRew = groupMeans(isRew,2,sesNo,'dim');  % split
    % isRew(drugA==1,:,:) = isRew(drugA==1,[2 1],:); % flip
    isRew = mod(ttCond,2) == 0; % 0=low, 1=high
    [infoReady, trajReady] = ExtractInfo(fixResult, 'infoReady', 'normTrajReady', drugA, sesNo);

    trajReady(:,:,:,1501:end) = NaN; % trim to same length

    % run analyse each fixation period
    OReady = LoopFix(infoReady, trajReady, isRew, 1, [1 1500], 'TrioEEGBehReadyMicroFixation.mat');

    % plot
    N = 500;
    vnames = {'pp','rew','drug'}; % names of factors for rmanova
    rewNames = {'low','high'}; % no targ/dist in ready period
    t = 0:750:1500; % time inds for xticks
    readyStats = PlotAndStats(OReady, N, rewNames, t, vnames, infoReady);

    [readyStats.nMacroCond, readyStats.nMacroMeanStats] = AnalyseMacros(infoReady, isRew, t, vnames, rewNames);
    
    saves = [saves {'readyStats'}];
end

%% ready period incl distr?

if any(strcmpi(cuesToAnalyse, 'ready2'))
    % isRew = nancat(2,result.isRew); % is high reward
    % isRew = groupMeans(isRew,2,sesNo,'dim');  % split
    % isRew(drugA==1,:,:) = isRew(drugA==1,[2 1],:); % flip
%     isRew = mod(ttCond,2) == 0; % 0=low, 1=high
    [infoReady2, trajReady2] = ExtractInfo(fixResult, 'infoReady', 'normTrajReady', drugA, sesNo);

    trajReady2(:,:,:,1501:end) = NaN; % trim to same length

    % run analyse each fixation period
    OReady2 = LoopFix(infoReady2, trajReady2, ttCond, 1, [1 1500], 'TrioEEGBehReady2MicroFixation.mat');

    % plot
    N = 500;
    vnames = {'pp','rew','dist','drug'}; % names of factors for rmanova
    condNames = {'targLo','targHi','distLo','distHi'};
    t = 0:750:1500; % time inds for xticks
    ready2Stats = PlotAndStats(OReady2, N, condNames, t, vnames, infoReady2);

    [ready2Stats.nMacroCond, ready2Stats.nMacroMeanStats] = AnalyseMacros(infoReady2, ttCond, t, vnames, condNames);
    
    saves = [saves {'ready2Stats'}];
end


%% do pre-target period
if any(strcmpi(cuesToAnalyse, 'targ'))
    [infoTarg, trajTarg] = ExtractInfo(fixResult, 'infoMicro', 'normTrajMicro', drugA, sesNo);

    excludeMacros = 0; % 1=exclude all macrosaccades, which allows us to use the data up until first macrosaccade
    % 0 = don't, so will just exclude any trial with a macrosaccade within timerange
    if excludeMacros
        % trajTarg is already set to NaN from sRT of first macro saccade onwards
        % (or entire trial if none found)
        % do the same for the saccade info (we can just set all macros to NaN and
        % stripnan
        fn = fieldnames(infoTarg(1));
        for i = 1:numel(infoTarg)
            isMacro = infoTarg(i).sAmpl >= 1; % set all macro saccades to NaN

            % set all to NaN
            for j = 1:length(fn)
                infoTarg(i).(fn{j})(isMacro) = NaN;
            end

        end
    end

    timeRange = [1 200]; % time period to look in 

    OTarg = LoopFix(infoTarg, trajTarg, ttCond, 0,timeRange, 'TrioEEGBehTargMicroFixation.mat');

    % plot
    N = 100;
    vnames = {'pp','rew','dist','drug'}; % names of factors for rmanova
    condNames = {'targLo','targHi','distLo','distHi'};
    t = 0:100:200; % time inds for xticks
    targStats = PlotAndStats(OTarg, N, condNames, t, vnames);

    [targStats.nMacroCond, targStats.nMacroMeanStats] = AnalyseMacros(infoTarg, ttCond, t, vnames, condNames);
    
    saves = [saves {'targStats'}];
end

%% save

for i = 1:length(saves)
    save(sprintf('TrioEEGBehMicrosaccades_%s.mat',saves{i}), saves{i});
end
end


function [info, traj] = ExtractInfo(result, infoName, trajName, drugA, sesNo)
% get info and traj in size [pp drug tr [sample]]



info = nancat(2, result.(infoName)); % [1 nPP*2]

info = groupMeans(info', 1, sesNo, 'dim')'; % [nPP nSes]
info(drugA==1,:) = info(drugA==1, [2 1]); % flip [placebo drug]

% normalised trajectory
traj = nancat(3, result.(trajName)); %[nTr nSamp nPP*2]
traj = permute(groupMeans(traj,3, sesNo,'dim'), [2,3,4,1]);  % [pp ses tr sampl]
traj(drugA==1,:,:,:) = traj(drugA==1,[2 1],:,:);% flip drugs

% [nPP, nDrug, nTr, ~] = size(traj);


end


function O = LoopFix(info, traj, condition, doFFTHalves, timeRange, loadName)

[nPP,nDrug,~,~] = size(traj);


%% analyse fixation


ppd = 1; % don't need to convert as already in degrees
% doFFTHalves = 1;
if exist(loadName,'file') % load if done already
    o1 = load(loadName,'O');
    O = o1.O;
else % process them
    for i = 1:nPP
        for j = 1:nDrug
            disp([i j]);
            nTr = size(info(i,j).sAmpl,1); % n trials for this session
%             nSamples = max(apply(2, @(x) find(isnan(x),1)-1, (sq(traj(i,j,1:nTr,:))))); % max number of samples for this session
            O1 = analyseFixationPeriod(sq(traj(i,j,1:nTr,:)), info(i,j), sq(condition(i,j,1:nTr)), timeRange, 0, ppd);
            if doFFTHalves
                % FFT on first and second halves
                fft1 = analyseFixationPeriod(sq(traj(i,j,1:nTr,:)), info(i,j), sq(condition(i,j,1:nTr)), [timeRange(1) floor(timeRange(2)/2)], doFFTHalves, ppd);
                fft2 = analyseFixationPeriod(sq(traj(i,j,1:nTr,:)), info(i,j), sq(condition(i,j,1:nTr)), [ceil(timeRange(2)/2) timeRange(2)], doFFTHalves, ppd);
                O1.fft1 = fft1.fft;
                O1.fft2 = fft2.fft;
            end
            O(i,j) = O1; % store here
        end
    end
    save(loadName, 'O', '-v7.3'); % save this here as it takes ages to run
end

end


function output = PlotAndStats(O, N, condNames, t, vnames, info)

drugNames = {'Placebo','Drug'};


% extract

microDens = nancat(4, nancat(3, O(:,1).microDensity), nancat(3, O(:,2).microDensity)); % density over time
meanSpeedCond = nancat(5, nancat(4, O(:,1).meanSpeedCond), nancat(4, O(:,2).meanSpeedCond)); % mean ocular drift speed

if all(isfield(O, {'fft1', 'fft2'})) % is fft split into two halves?
    fft1 = nancat(5, nancat(4, O(:,1).fft1), nancat(4, O(:,2).fft1)); % mean ocular drift speed
    fft2 = nancat(5, nancat(4, O(:,1).fft2), nancat(4, O(:,2).fft2)); % mean ocular drift speed
    nMeas = 4;
else
    fft1 = nancat(5, nancat(4, O(:,1).fft), nancat(4, O(:,2).fft)); % mean ocular drift speed
    fft2 = NaN(size(fft1));
    nMeas = 3;
end


% change sizes
microDens = permute(microDens, [3,1,2,4]); %[pp 100 conds drugs]
meanSpeedCond = permute(meanSpeedCond, [4,1,3,5,2]); %[pp time cond drug tr]
fft1 = permute(fft1, [4,1,3,5,2]); %[pp freq cond drug tr]
fft2 = permute(fft2, [4,1,3,5,2]); %[pp freq cond drug tr]

[nPP, ~, nConds, ~] = size(microDens); % get these
%%
srate = 1000; %Hz
nyquist = srate/2;
freqs = linspace(0,nyquist,floor(N/2)+1);
plotfreq = 1:min([125 length(freqs)]);

%% stats - use permOLS

for i = 1:nConds
    [~,fixP(i,:)] = permutationOLS(diff(microDens(:,:,i,1:2),[],4), [],[],[],'cluster',true,'clustermethod','mean','two_tailed',true);
end
    

%% plot
n=2;%2=PD ON vs OFF, 3=HC too
figure();
doPerm = 1;
useClust = 1;
nPerms = 5000;

% t = 0:750:1500;
for i = 1:nConds
    subplot(nMeas,nConds,i)
    errorBarPlot(sq(microDens(:,:,i,1:n)),'area',1, 'doStats', 0);%,'xaxisvalues',xi);
    title(condNames{i});
    if i==1; ylabel ('microsaccades'); xlabel('time during fixation period (ms)'); end
    set(gca,'XTick',0:50:100,'XTickLabel',t)
    if doPerm
        [~,p]=permutationOLS( diff(microDens(:,:,i,1:2),[],4), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p);
    end
    
    subplot(nMeas,nConds,nConds+i);
    errorBarPlot(sq(nanmean(meanSpeedCond(:,:,i,1:n,:),5)),'area',1, 'doStats', 0);
    if i==1; ylabel('drift speed');xlabel('time during fixation period (ms)'); end
    xlim([min(t) max(t)]);
    set(gca,'XTick',t,'YTick',0:.1:1);
    if doPerm
        [~,p]=permutationOLS( diff(nanmean(meanSpeedCond(:,:,i,1:2,:),5),[],4), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p,'yVal',0);
    end
    
    subplot(nMeas,nConds,(2*nConds)+i); 
%     plotfreq=50:250; %  don't report ultra-high or low frequencies
    errorBarPlot(sq(nanmean(log(fft1(:,plotfreq,i,1:n,:)),5)),'area',1, 'xaxisvalues',freqs(plotfreq), 'doStats', 0);
    if i==1; ylabel 'early log power'; xlabel 'frequency (Hz)'; end
%     set(gca,'XTick',50:100:250)
    if doPerm
        [~,p]=permutationOLS( diff(nanmean(log(fft1(:,plotfreq,i,1:2,:)),5),[],4), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p);
    end
    
    if nMeas == 4
        subplot(nMeas,nConds,(3*nConds)+i);
    %     plotfreq=50:250; %  don't report ultra-high or low frequencies
        errorBarPlot(sq(nanmean(log(fft2(:,plotfreq,i,1:n,:)),5)),'area',1, 'xaxisvalues',freqs(plotfreq), 'doStats', 0);
        if i==1; ylabel 'late log power'; xlabel 'frequency (Hz)'; end
    %     set(gca,'XTick',50:100:250)
        if doPerm
            [~,p]=permutationOLS( diff(nanmean(log(fft2(:,plotfreq,i,1:2,:)),5),[],4), [],[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
            hold on;
            pbar(p);
        end
    end
end
drawnow

for i = 1:nMeas, makeSubplotScalesEqual(nMeas,nConds,i*nConds-(nConds-1):i*nConds); end

subplot(nMeas,nConds,nMeas*nConds); 
h = findobj(gca,'Type','Line');
legend(flipud(h),drugNames,'Location','Best');

%% plot rew effects for each condition

s = 10; % movmean window size
s2 = 10; % window size for fft
n=2;
figure();
doPerm = 1;
useClust = 1;
nPerms = 5000;
% plotfreq = 1:length(freqs);
for i = 1:(nConds/2)
    subplot(nMeas,nConds/2,i)
    y = movmean(sq(diff(microDens(:,:,i*2-1:i*2,1:n),[],3)),s,2);
    errorBarPlot(y,'area',1, 'doStats', 0);%,'xaxisvalues',xi);
    
    title([condNames{i*2} ' - ' condNames{i*2-1}]);
    if i==1; ylabel ('microsaccades'); xlabel('time during fixation period (ms)'); end
    set(gca,'XTick',0:50:100,'XTickLabel',t)
    yline(0,':k');
    if doPerm
        y1 = diff(y,[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1(:,:,1), x,[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    
    
    subplot(nMeas,nConds/2,i+nConds/2); 
    y = movmean(sq(diff(nanmean(meanSpeedCond(:,:,i*2-1:i*2,1:n,:),5),[],3)),s,2);
    errorBarPlot(y,'area',1, 'doStats', 0);%,'xaxisvalues',xi);
    if i==1; ylabel('drift speed');xlabel('time during fixation period (ms)'); end
    yline(0,':k');
    xlim([0 max(t)]);
    set(gca,'YTick',-1:.1:1,'XTick',t);
    if doPerm
        y1 = diff(y,[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1(:,:,1), x,[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    
    subplot(nMeas,nConds/2,i+2*nConds/2); 
    y = movmean(sq(nanmean(diff(log(fft1(:,plotfreq,i*2-1:i*2,1:n,:)),[],3),5)),s2,2);
    errorBarPlot(y,'area',1, 'xaxisvalues',freqs(plotfreq), 'doStats', 0);
    if i==1; ylabel 'early log power'; xlabel 'frequency (Hz)'; end
    yline(0,':k');
    set(gca,'XTick',0:250:500)
    if doPerm
        y1 = diff(y,[],3);
        x = ones(nPP,1);
        [~,p]=permutationOLS( y1(:,:,1), x,[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
        hold on;
        pbar(p, 'yVal', min(ylim));
    end
    
    if nMeas==4
        subplot(nMeas,nConds/2,i+3*nConds/2);
        y = movmean(sq(nanmean(diff(log(fft2(:,plotfreq,i*2-1:i*2,1:n,:)),[],3),5)),s2,2);
        errorBarPlot(y,'area',1, 'xaxisvalues',freqs(plotfreq), 'doStats', 0);
        if i==1; ylabel 'late log power'; xlabel 'frequency (Hz)'; end
        yline(0,':k');
        if doPerm
            y1 = diff(y,[],3);
            x = ones(nPP,1);
            [~,p]=permutationOLS( y1(:,:,1), x,[],[],'cluster',useClust,'clustermethod','mean','two_tailed',true,'nperms',nPerms);
            hold on;
            pbar(p, 'yVal', min(ylim));
        end
    end
    set(gca,'XTick',0:250:500)
end
drawnow

if nConds==4
    for i = 1:nMeas, makeSubplotScalesEqual(nMeas,nConds/2,i*2-1:i*2); end
end

subplot(nMeas,nConds/2,nMeas*nConds/2); 
h = findobj(gca,'Type','Line');
legend(flipud(h),drugNames,'Location', [0.7772 0.2584 0.1268 0.0726]);

%% mean microsaccade density across 1500ms
% dens = permute(reshape(microDens,[nPP,100,2,2,2]),[1,3,4,5,2]);

meanDens = sq(nanmean(microDens, 2)); % mean across time

% vnames = {'pp','rew','drug'};
ls = '-';
if nConds==4
    meanDens = reshape(meanDens, nPP, 2,2,2);
    ls = 'none'; % line style
end
meanDensStats = rmanova(meanDens, vnames);


figure();
h = errorBarPlot(sq(nanmean(microDens,2)),'plotargs',{'LineWidth',2,'LineStyle',ls}, 'doStats', 0);
if nConds==4
    JoinUpErrorPoints(h, [1 2; 3 4]);
end
set(gca,'XTick',1:nConds, 'XTickLabel',condNames);
xlim([.5 nConds+.5])
ylabel('# microsaccades')
xlabel('condition')
legend(drugNames, 'Location','Best')
box off

%% mean drift speed across 1500ms

meanDrift = sq(nanmean(meanSpeedCond,[5,2]));
if nConds==4
    meanDrift = reshape(meanDrift, nPP, 2,2,2);
end
meanDriftStats = rmanova(meanDrift, vnames);


figure();
h = errorBarPlot(sq(nanmean(meanSpeedCond,[5,2])),'plotargs',{'LineWidth',2,'LineStyle',ls}, 'doStats', 0);
if nConds==4
    JoinUpErrorPoints(h, [1 2; 3 4]);
end

set(gca,'XTick',1:nConds, 'XTickLabel',condNames);
xlim([.5 nConds+.5])
ylabel('mean ocular drift speed (deg/s)')
xlabel('condition')
legend(drugNames,'Location','Best')
box off



%% 

output.meanDens = meanDens;
output.meanDrift = meanDrift;
output.meanSpeedCond = meanSpeedCond;
output.fft1 = fft1;
output.fft2 = fft2;
output.meanDriftStats = meanDriftStats;
output.meanDensStats = meanDensStats;

end


function [nMacroCond, nMacroMeanStats] = AnalyseMacros(info, condition, t, vnames, condNames)

nConds = length(condNames);
drugNames = {'Placebo','Drug'};

%% count number of macrosaccades too

sAmpl = nancat(4,nancat(3,info(:,1).sAmpl),nancat(3,info(:,2).sAmpl));
sRT = nancat(4,nancat(3,info(:,1).sRT),nancat(3,info(:,2).sRT));

isMacro = sAmpl >= 1 & sRT >= min(t) & sRT <= max(t);

nMacro = sq(sum(isMacro,2)); %[tr pp drug]

% split by cond

nMacroCond = permute(groupMeans(nMacro,1,permute(condition,[3,1,2]),'dim'),[2,1,3,4]); %[pp cond drug trial]

%% rmanova and plot

nMacroMean = nanmean(nMacroCond,4); % mean across trials

ls = '-';
if nConds==4
    nMacroMean = reshape(nMacroMean, size(nMacroMean,1), 2,2,2);
    ls = 'none'; % line style
end
nMacroMeanStats = rmanova(nMacroMean, vnames);


figure();
h = errorBarPlot(sq(nanmean(nMacroCond,4)),'plotargs',{'LineWidth',2,'LineStyle',ls}, 'doStats', 0);
if nConds==4
    JoinUpErrorPoints(h, [1 2; 3 4]);
end
set(gca,'XTick',1:nConds, 'XTickLabel',condNames);
xlim([.5 nConds+.5])
ylabel('mean # macrosaccades')
xlabel('condition')
legend(drugNames, 'Location','Best')
box off

end
