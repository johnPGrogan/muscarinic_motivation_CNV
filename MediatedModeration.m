% MediatedModeration

timelock = 'readyBS';
load(sprintf('TrioEEGAnalyseByTrials_%s.mat',timelock), 'meanAmplTrials','meanAmplTrialsPost');

load('TrioEEGBehAnalysis.mat','trialTab');
%%

trialTab.erp1 = col(meanAmplTrials(:,:,:,:,30,:)); % early ERP -Cz
trialTab.erp2 = col(meanAmplTrialsPost(:,:,:,:,30,:)); % late - Cz
trialTab.erp2F = col(meanAmplTrialsPost(:,:,:,:,2,:)); % FPz
%% moderated mediation

clc;
beh = 'srt'; % or velres
eeg = 'erp2'; % or erp2F
main = 'rew';
moderator = 'drug';
formulas = {sprintf('%s ~ 1 + %s + (1 | pp)',beh, main);
            sprintf('%s ~ 1 + %s + (1 | pp)', eeg, main);
            sprintf('%s ~ 1 + %s + %s + (1 | pp)', beh, eeg, main);
            sprintf('%s ~ 1 + %s*%s + (1 | pp)', beh, main, moderator);
            sprintf('%s ~ 1 + %s*%s + (1 | pp)', eeg, main, moderator);
            sprintf('%s ~ 1 + %s*%s + %s*%s + (1 | pp)', beh, main, moderator, eeg, moderator);
            };
       
nF = length(formulas);
regs = cell(nF,1);
for iF = 1:nF
    regs{iF} = fitglme(trialTab, formulas{iF});

end


%% tests:

% rew effects velres
isEff1 = regs{1}.Coefficients.pValue(2) <= .05
regs{1}.Coefficients.Estimate(2)

% rew effects pns
isEff2 = regs{2}.Coefficients.pValue(2) <= .05
regs{2}.Coefficients.Estimate(2)

% psn mediates rew effect on velres
isMed = regs{3}.Coefficients.pValue(3) <= .05
regs{3}.Coefficients.Estimate(3)

effIsReduced = abs(regs{3}.Coefficients.Estimate(2)) < abs(regs{1}.Coefficients.Estimate(2)) 

mediationSize = abs(regs{1}.Coefficients.Estimate(2)) - abs(regs{3}.Coefficients.Estimate(2)) 
%% mediated moderation

isMod = regs{4}.Coefficients.pValue(4) <= .05 %b_43: velres ~ rew:drug
regs{4}.Coefficients.Estimate(4)

isMod2 = regs{5}.Coefficients.pValue(4) <= .05 && regs{6}.Coefficients.pValue(4) <= .05 %b_53 (psn ~ rew:drug) & b_64  (velres ~ psn)
regs{5}.Coefficients.Estimate(4)%b_53 (psn ~ rew:drug)
regs{6}.Coefficients.Estimate(4)% b_64  (velres ~ psn)

isMod3 = regs{5}.Coefficients.pValue(2) <= .05 && regs{6}.Coefficients.pValue(6) <= .05 %b_51 (psn ~ rew) & b_65 (velres ~ drug:psn)
regs{5}.Coefficients.Estimate(2)%b_51 (psn ~ rew) 
regs{6}.Coefficients.Estimate(6)%b_65 (velres ~ drug:psn)

% b_43 (velres ~ rew:drug) & (b_63 (velres ~ rew:drug) < b_43 (velres ~ rew:drug)) 
isMediatedMod = regs{4}.Coefficients.pValue(4) <= .05 && abs(regs{6}.Coefficients.Estimate(5)) < abs(regs{4}.Coefficients.Estimate(4))
regs{4}.Coefficients.Estimate(4)%b_43 (velres ~ rew:drug)
regs{6}.Coefficients.Estimate(5)%(b_63 (velres ~ rew:drug) 
regs{4}.Coefficients.Estimate(4)%b_43 (velres ~ rew:drug)) 

% b_43 - b_63
mediatedModeration = abs(regs{4}.Coefficients.Estimate(4)) - abs(regs{6}.Coefficients.Estimate(5))


%% moderated mediation - not fonud

isRew = regs{4}.Coefficients.pValue(2) < .05 % b4_1 is sig
isNotMod = regs{4}.Coefficients.pValue(4) > .05 % b4_3 is not sig

%either/or:
% b5_3 and b6_4 are sig - i.e. THP moderates rew eff on psn, and psn effects velres
isMed1 = regs{5}.Coefficients.pValue(4) <= .05 && regs{6}.Coefficients.pValue(4) <= .05

% or B6_5 and b5_1 are sig: i.e. mediation depends on moderator and
% reward has effect on psn
isMed2 = regs{6}.Coefficients.pValue(6) <= .05 && regs{5}.Coefficients.pValue(2) <= .05



% does 

