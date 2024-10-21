function PlotFigure3(timelock)
% plot mean ERP, difference waves, diff-difference waves
% mean ERPs, indiv-pp mean ERPs, and regressions of means vs beh, a

    if nargin==0 
        timelock = 'Figure 3 - Source Data 1.mat';
    end
    if regexp(timelock, '.mat')
        % try loading this as a file name
        data = load(timelock, 'audCue_200','readyBS'); % both epochs
        
        % plot incentive cue and p3a/cnv separately, due to non-contiguous
        % times

        timelocks = {'audCue_200','readyBS'};
        for i = 1:2
            PlotIt(data.(timelocks{i}));
        end
    
    else
        % if not using source data, load up data, store as struct, then
        % call as above
        if strcmp(timelock, 'readyBS')
        
            load('trioEEGAnalyseByTrials_readyBS_small.mat',...
                'eegByCond','chanNames','timesToPlot','xTimes','toi');
            erpNames = {'P3a','CNV'};
    
            
        elseif strcmp(timelock, 'audCue_200')
            load('trioEEGAnalyseByTrials_audCue_200_small.mat',...
                'eegByCond','chanNames','timesToPlot','xTimes','toi');
            erpNames = {'incentiveCueERP'};
            
        
        else
            error('timelock is not filename or cue-name');
            
        end
    
        % take means and reformat
    
        [nChans,nTimes,nConds,nPP,nDrugs,nTr]  = size(eegByCond);
        iCh = 30;  % Cz channel to keep
    
        % keep only one channel, change dimension order
        % [pp, time, rew, distr, drugs, tr, ch=30]
        eegByCond2 = permute(reshape(eegByCond(iCh,:,:,:,:,:), 1, nTimes, 2, 2, nPP, nDrugs,nTr), [5 2 3 4 6 7 1]);
    
        [nPP, nTimes, nRew,nDistr,nDrugs,nTr] = size(eegByCond2);
    
        nToi = size(toi,1); % number of erp windows
    
        % take means within window, per trial
        meanAmpls = NaN(nPP,nRew,nDistr,nDrugs,nTr,nToi);
        for i = 1:nToi
            meanAmpls(:,:,:,:,:,i) = sq(nanmean(eegByCond2(:, isBetween(xTimes, toi(i,:)),:,:,:,:),2)); %[pp rew distr dr tr]
        end
    
        % take mean over trials + distractor-conditions
        eegByCond2 = sq(nanmean(eegByCond2,[4 6])); %[pp, time, rew, drugs]

    
        %% load regression table and run regs
        load('TrioEEGBehAnalysis.mat','trialTab','varNames');
        nVars = 3;
    
    
        % put means into table
        for i = 1:nToi
            trialTab.(erpNames{i}) = double(nanzscore(col(meanAmpls(:,:,:,:,:,i)))); % make sure double, for fitglme
        end
           
        % setup this extraction func [b, se, t, p]
        extractVar = @(x,inds) [x.Coefficients.Estimate(inds), x.Coefficients.SE(inds), x.Coefficients.tStat(inds), x.Coefficients.pValue(inds), x.Coefficients.Lower(inds), x.Coefficients.Upper(inds)];
        
        % run regs
        formula = '%s ~ 1 + rew*distr*drug + %s + (1 | pp)';
        toiRegStats = zeros(3,nToi,8,6); % [vars, erp, terms, measures]
        for i = 1:nToi
            for j = 1:nVars
                reg = fitglme(trialTab, sprintf(formula, varNames{j}, erpNames{i}));
                
                toiRegStats(j,i,:,:) = extractVar(reg, 2:9);
            end
        end

        %% store into struct
        data = struct();
        data.meanAmpls = meanAmpls;
        data.eegByCond2 = eegByCond2;
        data.timesToPlot = timesToPlot;
        data.xTimes = xTimes;
        data.toi = toi;
        data.erpNames=erpNames;
        data.toiRegStats = toiRegStats;
    
        %% call plotting function on that 

        PlotIt(data);
    
    end
end

function PlotIt(data)
%     load('chanlocs.mat'); % channel data

% unpack data
meanAmpls = data.meanAmpls;
eegByCond2 = data.eegByCond2;
timesToPlot = data.timesToPlot;
xTimes = data.xTimes;
toi = data.toi;
erpNames = data.erpNames;
toiRegStats = data.toiRegStats;
clear data;


    [nPP, nTimes, nRew, nDrugs] = size(eegByCond2);
    nToi = size(toi,1);
    iCh = 30; % Cz

    vNames = {'residual velocity','RT','distractor pull'};

    
    %% plot ERP figure

    drugRewNames = {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'};
    cols = [0 0 1; 0 .5 .9;1 0 0;  .9 .5 0];
    
    figure();
    h = PlotEEGErrorBars(reshape(eegByCond2,nPP, nTimes, 4), 'chanNames',{'Cz'}, 'chansToPlot', {'Cz'},...
        'XTimes', xTimes, 'timesToPlot', timesToPlot, 'condNames', drugRewNames,...
        'yLabels', {'\muV'},'lineCols',cols,'xLabels',{'time from preparation cue (ms)'},'subplotInds',0);
    box off
    % show TOI
    hold on;
    for i = 1:nToi
        fill(toi(i,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');
    end
    
    legend([h{:,1}], drugRewNames,'Location','Best'); % overwrite to remove grey boxes


    %% Fig3h - mean CNV difference waves for Incentive and THP
    
    tInds = isBetween(xTimes, timesToPlot); % smaller window?
    
    % make difference waves for incentive/drug effect, by first taking mean
    % over the other one, then the difference
    dims = [4;3]; % leave rew, drug
    eegDiff = NaN(nPP, sum(tInds), 2);
    for i = 1:2
        eegDiff(:,:,i) = diff(sq(nanmean(eegByCond2(:,tInds,:,:),dims(i))),[],3); % [pp T ch]
    end
    
    effNames = {'Incentive Effect', 'THP Effect'};
    figure();
    c = get(gca,'ColorOrder');
    set(gca,'ColorOrder', c([3,4,5],:),'nextplot','replacechildren');
    h = errorBarPlot(eegDiff, 'area',1,'xaxisvalues',xTimes(tInds),'alpha',.2);
    hold on;
    ylim([-3 2]);
    % fill windows too?
    for i = 1:nToi
        fill(toi(i,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');
    end
    xlabel('time from preparation cue (ms)');
    ylabel('difference wave \Delta \muV');
    xline(0, '--k');
    yline(0, '--k');
    xlim(timesToPlot)
    legend([h{:,1}], effNames,'Location','Best');


    %% Fig 3i - mean CNV difference waves for Incentive, for THP+Placebo separately

    drugNames = {'Placebo','THP'};
    % take incentive difference in each drug cond
    eegIncentiveEffects = sq(diff(eegByCond2(:,tInds,:,:),[],3)); % [pp T drug]

    figure();
    h = errorBarPlot(eegIncentiveEffects, 'area',1,'xaxisvalues',xTimes(tInds),'alpha',.2);
    hold on;
    % fill windows too?
    for i = 1:nToi
        fill(toi(i,[1 2 2 1]), col(ylim .* [1;1])', [.5 .5 .5], 'FaceAlpha', .5, 'LineStyle','none');
    end
    xlabel('time from preparation cue (ms)');
    ylabel('difference wave \Delta \muV');
    xline(0, '--k');
    yline(0, '--k');
    xlim(timesToPlot)
    legend([h{:,1}], drugNames,'Location','Best');


    %% Fig 3b,c,d - mean window-ampl per ERP for Incentive*THP
    
    for i = 1:nToi
        % average over distractor/trials dimensions together
        meanAmpl = sq(nanmean(meanAmpls(:,:,:,:,:,i),[3 5])); % [pp, rew, drug]

        figure();    
        h = errorBarPlot(meanAmpl,'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
        JoinUpErrorPoints(h, [1 2]);
        box off
        set(gca,'XTick',1:2,'XTickLabel',{'0p','50p'});
        ylabel('mean \muV');
        xlabel('incentive');
        xlim([.5 2.5])
        title(erpNames{i});
    end

    %% plot these with individual-pp means

    for i = 1:nToi

        % average over distractor/trials dimensions together
        meanAmpl = sq(nanmean(meanAmpls(:,:,:,:,:,i),[3 5])); % [pp, rew, drug]

        figure();    
        h = errorBarPlot(meanAmpl,'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
        JoinUpErrorPoints(h, [1 2]);
        
        hold on; % plot indivs
        for j = 1:2
            plot( (1:2)' + (rand(nPP,2)-.5)'*.1, meanAmpl(:,:,j)', '-x', 'LineWidth', 1, 'Color', [h(j).Color .2]);
        end


        box off
        set(gca,'XTick',1:2,'XTickLabel',{'0p','50p'});
        ylabel('mean \muV');
        xlabel('incentive');
        xlim([.5 2.5])
        title(erpNames{i});
    end

    %% Fig 3e,f,g - regression coeffs for mean ERP vs behaviour

    alpha = .05 / 9; % bonferroni, across 3 erps and 3 vars
    
    % now plot beta +SE
    for j = 1:nToi
        f = figure();
        set(f, 'DefaultAxesFontSize',16)
        h = barwitherr(toiRegStats(:,j,4,2)*1.96, toiRegStats(:,j,4,1));
        [h.LineWidth] = deal(1.2);
        ylim([-.05, .15]);
        xticklabels(vNames);
        ylabel('\beta coefficient')
        box off
        
        stars = col(toiRegStats(:,j,4,4)') < alpha;
    %     y = repmat(max([h.YData])*1.1, 1,3);    
        y = repmat(-.04, 1,3);
        y(~stars) = NaN;
        hold on;
        plot((1:3) , y, '*','MarkerSize',10, 'Color', h.FaceColor, 'LineWidth',1.5);
                
        title(erpNames{j});
    end



end