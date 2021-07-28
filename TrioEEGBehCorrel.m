% TrioEEGBehCorrel
clc; clear all; close all;
load('TrioEEGBehAnalysis.mat','varNames','trialTab','vars' );
nVars=3;
[nPP, nConds, nDrugs, nTr] = size(vars.velres);
%% do all vars?
% varNames = {'vel','ampl','srt','endptVarCorr','depAngle','velres'};
% nVars = length(varNames);

%% regress each variable against each other (control for factors)

[r,p] = deal(NaN(nVars));
for i = 1:nVars
    for j = 1:nVars
        if i == j
            r(i,j) = 0;
            p(i,j) = NaN;
        else
            f = fitglme(trialTab, sprintf('%s ~ 1 + %s + rew*distr*drug + (1 | pp)', varNames{i}, varNames{j}));
            r(i,j) = f.Coefficients.Estimate(5);
            p(i,j) = f.Coefficients.pValue(5);
        end
    end
end

%% plot a correlation matrix


figure(1);clf
colormap(flipud(crameri('vik')));
imagesc(r, [-.5 .5]);
set(gca,'XTick',1:nVars, 'XTickLabel', varNames, 'YTick',1:nVars, 'YTickLabel', varNames,'YDir','reverse')
colorbar

hold on
pLims = [.0000001];
for ii = length(pLims)
    [i,j] = find(p <= pLims(ii));
    text(j,i,repmat('*',1,ii),'FontSize',20);
end


%% correl rew effects

% zscore across each pp
vars1 = vars; % no zscoring
vars1.srt = exp(vars.srt);
% vars1 = structfun(@(x) reshape(nanzscore(reshape(x,nPP,[])')',nPP,4,2,nTr), vars1,'UniformOutput',0);

effNames = {'incentive','distractor','THP'};
figure();
scatterArgs = {'pearson',0,'plot_ci',1,'text',1,'plotline',1,'showzero',1,'alpha', .05};
varsPairs = [2 1; 2 3; 1 3];
nEff = 1;
nVs = 3;
for iEff = 1:nEff
    varsEff = structfun(@(x) sq(diff(nanmean(reshape(x,nPP,2,2,2,nTr),5),[],iEff+1)), vars1, 'UniformOutput',0); 

    % plot rew eff against each others

    for iVs = 1:nVs
        subplot(2,2,iVs);
%         subplot(nVs,nEff+1,(iVs-1)*(nEff+1) + iEff)

%     plot(nanmean(varsEff.srt,[3 2]), nanmean(varsEff.velres, [3 2]), 'x');
%     lsline;
        scatterRegress(nanmean(varsEff.(varNames{varsPairs(iVs,1)}),[3 2]), nanmean(varsEff.(varNames{varsPairs(iVs,2)}), [3 2]), scatterArgs{:});
    

        xlabel([' \Delta ' (varNames{varsPairs(iVs,1)})]);
        ylabel(['\Delta ' (varNames{varsPairs(iVs,2)})]);
        
        if iVs
            title([effNames{iEff} ' effect']);
        end

    end
end
%%
% also do the THP*incentive effect
varsInt = structfun(@(x) sq(diff(diff(nanmean(reshape(x,nPP,2,2,2,nTr),5),[],2),[],4)), vars1, 'UniformOutput',0);
for iVs = 1:nVs
    subplot(nVs,nEff+1, iVs*(nEff+1));
    scatterRegress(nanmean(varsInt.(varNames{varsPairs(iVs,1)}),[2]), nanmean(varsInt.(varNames{varsPairs(iVs,2)}), [2]), scatterArgs{:});
    xlabel(['\Delta ' (varNames{varsPairs(iVs,1)})]);
    ylabel(['\Delta ' (varNames{varsPairs(iVs,2)})]);
    if iVs==1
        title('THP*incentive effect');
    end
end
%%

% [nPP distr drug]

% rew eff table
d = nanzscore(dePivot(varsEff.velres,'KeepNaN',1)); % get factors
% into table
rewTab = array2table(d(:,1:end-1),'VariableNames',{'pp','distr','drug'});
% get DVs - rew eff
d = struct2table(structfun(@(x) nanzscore(col(x)), varsEff, 'UniformOutput',0));
rewTab = horzcat(rewTab, d);


[r2,p2] = deal(NaN(nVars));
for i = 1:nVars
    for j = 1:nVars
        if i == j
            r2(i,j) = 0;
            p2(i,j) = NaN;
        else
%             f = fitglme(rewTab, sprintf('%s ~ 1 + %s + distr*drug + (1 | pp)', varNames{i}, varNames{j}));
%             r2(i,j) = f.Coefficients.Estimate(5);
%             p2(i,j) = f.Coefficients.pValue(5);
            [r2(i,j),p2(i,j)] = corr(rewTab.(varNames{i}), rewTab.(varNames{j}), 'rows','pairwise','type','Spearman');
        end
    end
end

%%

figure(2);clf
colormap(flipud(crameri('vik')));
imagesc(r2, [-.5 .5]);
set(gca,'XTick',1:nVars, 'XTickLabel', varNames, 'YTick',1:nVars, 'YTickLabel', varNames,'YDir','reverse')
colorbar

hold on
pLims = [.05];
for ii = length(pLims)
    [i,j] = find(p2 <= pLims(ii));
    text(j,i,repmat('*',1,ii),'FontSize',20);
end


%% trialwise measure of reward effects?
% the random effect of reward on srt correlates with that on velres

% regress reward out, and then within each condition correlate the
% residuals? NO! this is the random variance, not the reward effect

f = fitglme(trialTab, 'velres ~ 1 + rew*distr*drug + (1 | pp)');
velresRewResid = f.residuals;
f = fitglme(trialTab, 'srt ~ 1 + rew*distr*drug + (1 | pp)');
srtRewResid = f.residuals;


%% 
distLabels = {'no distractor' ,'with distractor'};
figure();
% cols = [0 0 1; 0 .5 .9; 1 0 0;  .9 .5 0];
cols = get(gca,'ColorOrder');
cols = cols([1 1 2 2],:);

for i = 2
%     subplot(1,2,i)
    set(gca,'ColorOrder', cols);
    [~,~,~,~,h] = conditionalPlot(reshape(permute((vars.srt(:,[i*2-1 i*2],:,:)),[4,1,2,3]),nTr,nPP,[]), reshape(permute(vars.depAngle(:,[i*2-1 i*2],:,:),[4,1,2,3]),nTr,nPP,[]));
    [h(2:2:end).Visible] = deal('off');
    [h(1:2:end).Marker] = deal('none');
    [h(1:4:end).LineStyle] = deal('--');
    [h(1:2:end).LineWidth] = deal(2);
    ylabel('distractor pull (deg)');
    xlabel('saccadic RT (ms)');
%     if i==1
        legend(h(1:2:end), {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'}, 'Location','Best');
%     end
    title(distLabels(i));
    axis([180 403 -11 20]);
end

