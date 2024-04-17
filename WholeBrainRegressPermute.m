function WholeBrainRegressPermute2(varName, nPerms, saveFolder, eeg1, trialTabTrim)
% function WholeBrainRegressPermute2(varName, nPerms, saveFolder, eeg1, trialTabTrim)
% regress EEG at each time point and channel against behaviour, use
% permutation testing (and clustering to adjust p values afterwards)
% shuffle eeg once per permutation, and pass that to a regression loop for
% each chan*time. should be quicker than previous way. Saves tValsPerms
% only
% 
% Inputs: 
%   varName = string: 'velres' ,'srt', or 'depAngle'
%   nPerms = scalar: number of permutations to run
%   saveFolder = string: path to folder where to save data
%   eeg1: [nPP nConds nDrugs nTr nTimes nChans] matrix of eeg voltages
%   trialTabTrim: table  with behavioural data + conditions, with the rows
%       corresponding to the NaN voltages missing (eeg will be trimmed
%       before added to the table)
% 
% Outputs:
%   None. It saves the tValsPerm [nT nChans nPerms] in the savefolder,
%   timestamped.
% 

%% load

if ~exist('eeg1','var') || isempty(eeg1)
    timelock = 'readyBS';
    file1 = sprintf('RegressPermuteStuff2_%s.mat',timelock);
    load(file1','trialTabTrim','eeg1');

end

[~, ~, ~, ~, nT, nChans] = size(eeg1);

%% whole brain regressions with permutation clusters?

formula = sprintf('%s ~ 1 + rew*distr*drug + v + (1 | pp)', varName);

tValsPerm = zeros(nT, nChans, nPerms);

for iP = 1:nPerms
    
    fprintf('\nperm: %d:', iP);
%     trialTab1 = trialTab; % copy

    % shuffle trials within conds etc
    eegShuffled = shuffleWithNaNs(eeg1);

    % make into [19520, nT, nCh] and nanzscore - quicker in matrix
    eegShuffled = nanzscore(reshape(eegShuffled,[],nT,nChans));
    
    % can also removed NaN at this point?
    toRemove = all(isnan(eegShuffled),[2 3]);
    eegShuffled(toRemove,:,:) = [];   

    tic;
    % regress each t+ch in a loop, using this shuffled data. about 1 hr per
    tValsPerm(:,:,iP) = regressEachChanTime(trialTabTrim, formula, eegShuffled); 
    toc
end

% now get the real values
% [tVals, betaVals] = regressEachChanTime(trialTab, formula, eeg1); % ~ 35 mins

saveName = fullfile(saveFolder, ['WholeBrainRegressPermute_' varName '_' datestr(now, 'yy.mm.dd HH.MM.SS') '.mat']);
save(saveName, 'tValsPerm');

end

function [tVals, betaVals] = regressEachChanTime(trialTab, formula, eeg)
% run regressions on each channel and time point, within once shuffled
% permutation. return tVals which will be cluster-corrected later

[~,nT,nCh] = size(eeg);
[tVals, betaVals] = deal(NaN(nT,nCh));
trialTab1=trialTab; % if using parfor
parfor iCh = 1:nCh
    trialTab1=trialTab; % if using parfor
    fprintf('%d,',iCh);
    for iT = 1:nT
%         trialTab.v = nanzscore(col(eeg(:,:,:,:,iT,iCh)));
        trialTab1.v = eeg(:,iT,iCh);
        reg = fitglme(trialTab1, formula);
        
        betaVals(iT,iCh) = reg.Coefficients.Estimate(5);
%         pVals(iT,iCh) = reg.Coefficients.pValue(5);
        tVals(iT,iCh) = reg.Coefficients.tStat(5);
   
    end
end
end


function eegShuffled = shuffleWithNaNs(eeg)
% shuffle trials within pp/bin/drug, while skipping NaNs
% takes ~35 secs

empty = isnan(eeg(:,:,:,:,1,1)); % find empties

[nPP,nBins,nDrugs,~,~,~] = size(eeg); % get sizes
eegShuffled = NaN(size(eeg));
for iPP = 1:nPP
    for iB = 1:nBins
        for iDr = 1:nDrugs
            eegShuffled(iPP,iB,iDr,~empty(iPP,iB,iDr,:),:,:) = Shuffle(sq(eeg(iPP,iB,iDr,~empty(iPP,iB,iDr,:),:,:))); % shuffle valid trials
        end
    end
end


end