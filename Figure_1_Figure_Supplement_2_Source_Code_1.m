function PlotFigure1_S1(vars,varNames)
% plot sRT vs velres, incentive effects
% then separated by drug condition

    if nargin==0
        load('Figure 1 - Source Data 1.mat', 'vars', 'varNames')
        % vars is struct with fields of [nPP, nConds, nDrug] trial means
    end
    
    [nPP, nConds, nDrugs] = size(vars.srt); % mean-level data
    
    effNames = {'incentive','distractor','THP'};
    scatterArgs = {'pearson',0,'plot_ci',2,'text',1,'plotline',2,'showzero',1,'alpha', .05};
    varsPairs = [2 1;];
    iEff = 1; % incentive-effects
    
    % take difference in var across incentive levels
    varsEff = structfun(@(x) sq(diff(reshape(x,nPP,2,2,2),[],iEff+1)), vars, 'UniformOutput',0); 
    
    % plot rew eff against each others
    iVs = 1;
    
    figure();
    subplot(1,2,1);
    scatterRegress(nanmean(varsEff.(varNames{varsPairs(iVs,1)}),[3 2]), nanmean(varsEff.(varNames{varsPairs(iVs,2)}), [3 2]), scatterArgs{:});
    xlabel(['\Delta ' (varNames{varsPairs(iVs,1)})]);
    ylabel(['\Delta ' (varNames{varsPairs(iVs,2)})]);
    title([effNames{iEff} ' effect']);
    
    % now plot per drug cond
    varsInt = structfun(@(x) sq(diff(diff(reshape(x,nPP,2,2,2),[],2),[],4)), vars, 'UniformOutput',0);
    subplot(1,2,2);
    scatterRegress(nanmean(varsInt.(varNames{varsPairs(iVs,1)}),[2]), nanmean(varsInt.(varNames{varsPairs(iVs,2)}), [2]), scatterArgs{:});
    xlabel(['\Delta ' (varNames{varsPairs(iVs,1)})]);
    ylabel(['\Delta ' (varNames{varsPairs(iVs,2)})]);
    title('THP*incentive effect');
