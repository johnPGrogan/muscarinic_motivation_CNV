% MediatedModeration
clc; clear all; close all;
timelock = 'readyBS';
load(sprintf('TrioEEGAnalyseByTrials_%s.mat',timelock), 'meanAmplTrials','meanAmplTrialsPost');
%%
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

%% drop nan trials from all

toRemove = isnan(trialTab.(eeg)) | isnan(trialTab.(beh));
trialTab = trialTab(~toRemove,:);

%% run regs
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

mediationSize = abs(regs{3}.Coefficients.Estimate(2)) - abs(regs{1}.Coefficients.Estimate(2)) % negative, one-tailed


% sobel's test
a=regs{1}.Coefficients.Estimate(2); sa=regs{1}.Coefficients.SE(2);
b=regs{3}.Coefficients.Estimate(3); sb=regs{3}.Coefficients.SE(3);
t = (a*b) / (sqrt(b^2 * sa^2 + a^2*sb^2)) % t-stat
pU = 1-tcdf(abs(t), regs{3}.DFE) % directional one-tailed p-value
% pL = tcdf(-abs(t), regs{3}.DFE)

% p2 = (pU + pL)% two-tailed p-value. not used here though as we are directional




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
mediatedModeration = abs(regs{6}.Coefficients.Estimate(5)) - abs(regs{4}.Coefficients.Estimate(4)) % negative, one-tailed


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

%% bootstrapping?

% draw n samples with replacement from data
% recalc 4 regressions:
% a,b,c,c'
% repeat nPerms times to build approximate distribution of true indirect
% effects, and give conf-int, e.g. 2.5% and 97.5% for 0.5 alpha

% can just run some?
formulas2 = {sprintf('%s ~ 1 + %s + (1 | pp)',beh, main); % c
            sprintf('%s ~ 1 + %s + (1 | pp)', eeg, main); % a
            sprintf('%s ~ 1 + %s + %s + (1 | pp)', beh, eeg, main); %c2
            sprintf('%s ~ 1 + %s + %s + (1 | pp)', beh, eeg, main); %b
            sprintf('%s ~ 1 + %s*%s + (1 | pp)', beh, main, moderator); %b43
            sprintf('%s ~ 1 + %s*%s + (1 | pp)', eeg, main, moderator); %b53
            sprintf('%s ~ 1 + %s*%s + %s*%s + (1 | pp)', beh, main, moderator, eeg, moderator); %b63
            };
nF = length(formulas2);
fInds = [2 2 2 3 4 4 5] % which index to store per formula

terms = {'c', 'a', 'c2', 'b', 'b43', 'b53', 'b63', 'mediation', 'mediatedModeration'};

% remove all nans
toRemove = isnan(trialTab.(eeg)) | isnan(trialTab.(beh));
trialTab2 = trialTab(~toRemove,:);
% trialTab2=trialTab;

% make trial-type
trialTab2.tt = trialTab2.rew + trialTab2.distr/10 + trialTab2.drug/100;
uTT = unique(trialTab2.tt);
% sub-sample n per pp + cond
uPP = unique(trialTab2.pp);
nPerms = 2500;

% vn = trialTab2.Properties.VariableNames; % store
% % n = 50; % per pp per cond - or should this be same as n trials per pp/tt?
% b = NaN(nF, nPerms); 
% 
% parfor iP = 1:nPerms
%     if mod(iP, nPerms/10)==0; disp(iP); end
%     
%     trialTabBoot = trialTab2(1,:); % copy, will delete first row later
%     for iPP = 1:20
%         for iTT = 1:8
%             inds = find(trialTab2.pp == uPP(iPP) & trialTab2.tt==uTT(iTT));
%     
%             indsToKeep = inds(randi(length(inds),length(inds),1)); % draw with replacement
%     
%             trialTabBoot = vertcat(trialTabBoot, trialTab2(indsToKeep,:));
%         end
%     end
%     trialTabBoot(1,:) = []; % remove this row
%     
%     trialTabBoot = varfun(@nanzscore, trialTabBoot); % zscore again, to deal with diff means etc
%     trialTabBoot.Properties.VariableNames = vn; % replace
%     
%     for iF = 1:nF
%         fit = fitglme(trialTabBoot, formulas2{iF});
%         b(iF,iP) = fit.Coefficients.Estimate(fInds(iF));

%     end
% 
% end
% 
% %% estims
% 
% mediation = abs(b(3,:)) - abs(b(1,:));
% % mediationSize = abs(regs{3}.Coefficients.Estimate(2)) - abs(regs{1}.Coefficients.Estimate(2)) % negative, one-tailed
% b(nF+1,:) = mediation; % add here
% % terms{nF+1} = 'c2 - c';
% 
% medMod = abs(b(7,:)) - abs(b(5,:));
% b(nF+2,:) = medMod;
% % b(4,iP) = abs(regs{6}.Coefficients.Estimate(5)) - abs(regs{4}.Coefficients.Estimate(4)) % negative, one-tailed
% 
% 
% ci = prctile(b, [2.5 97.5], 2); % take lower and upper 95% CI
% meanBetas = nanmean(b,2);
% 
% disp([meanBetas, ci])
% 
% %% plot
% 
% figure();
% barwitherr(ci, meanBetas);
% xticklabels(terms);
% ylabel('\beta coefficient');
% 
% %% ksdensity
% 
% 
% [f,x] = ksdensities(b');
% % % f = f ./ sum(f,2); % norm
% % % f = f + (1:nF+1)'; % offset on yscale
% 
% figure();
% h = plot(x', f');
% legend(h, terms);


%% permutation approach?

% on each perm
% permute EEG randomly, across pps/conditions
% store betas for each
% these are the null distributions
% so now find proportion of null >= true

n = height(trialTab2);
bPerm = single(NaN(nF, nPerms)); 
parfor iP = 1:nPerms
    trialTabPerm = trialTab2; % copy

    if mod(iP, nPerms/10)==0; disp(iP); end

    % randomly permute just EEG? within each pp
    for iPP = 1:20
%         inds = find(trialTab2.pp == uPP(iPP));
        % within trials?
        for iTT = 1:8
            inds = find(trialTab2.pp == uPP(iPP) & trialTab2.tt==uTT(iTT));

            trialTabPerm.(eeg)(inds,:) = trialTabPerm.(eeg)(inds(randperm(length(inds))),:);
        end
    end
        
    for iF = 1:nF
        fit = fitglme(trialTabPerm, formulas2{iF});
        bPerm(iF,iP) = fit.Coefficients.Estimate(fInds(iF));

    end
end

bPerm(nF+1,:) = abs(bPerm(3,:))  - abs(bPerm(1,:));
bPerm(nF+2,:) = abs(bPerm(7,:)) - abs(bPerm(5,:));

%% get true

for iF = 1:nF
    fit = fitglme(trialTab2, formulas2{iF});
    bTrue(iF,1) = fit.Coefficients.Estimate(fInds(iF));
end
bTrue(nF+1,:) = abs(bTrue(3,:)) - abs(bTrue(1,:));
bTrue(nF+2,:) = abs(bTrue(7,:)) - abs(bTrue(5,:));

%% ksdensity, and plot true

figure();
t = tiledlayout('flow'); ax=[];
for i = 1:nF+1

    ax(i) = nexttile(t);

    ksdensity(bPerm(i,:));
%     hist(bPerm(i,:),100);
    xline(bTrue(i));
    title(terms{i});
end

%% calc prop

pLower = mean(bPerm <= bTrue,2);
% pUpper = mean(bPerm >= abs(bTrue),2);

disp(pLower)
% p = (pUpper + pLower)/2;
% disp(p);

%% format the regressions into easy to read 

fInds1 = [1 2; 1 2; 3 3; 3 2];
terms1 = {'c','a','b','c2'};

extractVar = @(x,inds) [x.Coefficients.Estimate(inds), x.Coefficients.SE(inds), x.Coefficients.tStat(inds), x.Coefficients.pValue(inds)];

mediationTable = NaN(length(fInds1),4);
for i=1:size(fInds1,1)
    mediationTable(i,:) = extractVar(regs{fInds1(i,1)}, fInds1(i,2));
end

mediationTable = array2table(mediationTable, 'VariableNames', {'beta','SE','t','p'});
% mediation effects
mediationTable.beta(5) = abs(mediationTable.beta(4)) - abs(mediationTable.beta(1));
mediationTable.SE(5) = NaN; % none
mediationTable.t(5) = NaN; % none
mediationTable.p(5) = pLower(end-1); % permutation-test p-value
mediationTable.Properties.RowNames = [terms1'; 'mediationEffect'];

mediationTable.isSignificant = mediationTable.p <= .05;

%% medmod table

fInds2 = [4 4; 5 4; 6 6; 6 5;];
terms2 = {'c','a','b','c2'};

extractVar = @(x,inds) [x.Coefficients.Estimate(inds), x.Coefficients.SE(inds), x.Coefficients.tStat(inds), x.Coefficients.pValue(inds)];

medModTable = NaN(length(fInds2),4);
for i=1:size(fInds2,1)
    medModTable(i,:) = extractVar(regs{fInds2(i,1)}, fInds2(i,2));
end

medModTable = array2table(medModTable, 'VariableNames', {'beta','SE','t','p'});
% mediation effects
medModTable.beta(5) = abs(medModTable.beta(4)) - abs(medModTable.beta(1));
medModTable.SE(5) = NaN; % none
medModTable.t(5) = NaN; % none
medModTable.p(5) = pLower(end); % permutation-test p-value
medModTable.Properties.RowNames = [terms1'; 'medModEffect'];

medModTable.isSignificant = medModTable.p <= .05;

%% save

save(sprintf('TrioEEGMediations_%s_%s.mat',timelock, beh),...
    'formulas', 'mediationTable','medModTable');
