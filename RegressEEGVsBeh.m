function RegressEEGVsBeh(timelock, timesToPlot)
% function RegressEEGVsBeh(timelock, timesToPlot)
% regress each time point and electrode separately against main 3
% behavioural variables
% Inputs:
%   timelock is a string of the epoch name to load e.g. "readyBS" or "sacc"
%   timesToPlot are the [min max] times to plot in the epoch


%% regress each time point

if nargin < 2
    timelock = 'readyBS';
    timesToPlot = [-200 1500];
end

load(sprintf('trioEEGAnalyseByTrials_%s.mat',timelock), 'eegByCond','xTimes')
[nChans, nTimes, nBins, nPP, nDrugs, nTr] = size(eegByCond); 


load('TrioEEGBehAnalysis.mat', 'trialTab','varNames','nVars');


eeg = reshape(permute(eegByCond,[4,3,5,6,2,1]),[],nTimes,nChans);
eeg = eeg(:,isBetween(xTimes, timesToPlot),1:61);
eeg = nanzscore(eeg);
for iV = 1:3
    formula = sprintf('%s ~ 1 + rew*distr*drug + v + (1 | pp)', varNames{iV});
    % include other two as covars?
%     formula = sprintf('%s ~ 1 + rew*distr*drug + v + %s + %s + (1 | pp)', varNames{iV}, varNames{mod(iV,3)+1}, varNames{mod(iV+1,3)+1});
    [~,nT,nCh] = size(eeg);
    [betaVals, pVals, tVals] = deal(zeros(nT, nCh));
    parfor iCh = 1:nCh
        disp(iCh);
        for iT = 1:nT
            trialTab1 = trialTab;
            trialTab1.v = eeg(:,iT,iCh);
            reg = fitglme(trialTab1, formula);

            ind = strcmpi(reg.CoefficientNames,'v');

            betaVals(iT,iCh) = reg.Coefficients.Estimate(ind);
            pVals(iT,iCh) = reg.Coefficients.pValue(ind);
            tVals(iT,iCh) = reg.Coefficients.tStat(ind);
        end
    end

    save(sprintf('%s_regress_%s.mat', timelock, varNames{iV}), ...
        'betaVals','pVals','tVals','formula','timesToPlot')
    beep;
end