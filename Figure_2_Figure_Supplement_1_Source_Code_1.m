function PlotFigure2_S1(vars)
% function PlotFigure2_S1()
% Figure 2 - Figure Supplement 1
% sRT vs distractor pull, per incentive*thp level
% only distractor-present trials

%% load if not given
    if nargin==0
        load('Figure 2 - Source Data 1.mat', 'vars')
        % vars is struct with fields of [nPP, nConds, nDrug, nTr] 
        % srt is raw RT, not log
    end

%% plot - only distractor present trials

    [nPP,nCond,nDr,nTr] = size(vars.srt);
    distLabels = {'no distractor' ,'with distractor'};
    
    figure();
    cols = get(gca,'ColorOrder');
    cols = cols([1 1 2 2],:);
    set(gca,'ColorOrder', cols);
    
    [~,~,~,~,h] = conditionalPlot(reshape(permute(vars.srt(:,3:4,:,:),[4,1,2,3]), nTr,nPP,4),...
        reshape(permute(vars.depAngle(:,3:4,:,:),[4,1,2,3]), nTr,nPP,4),[],'dostats',0);
    [h(2:2:end).Visible] = deal('off');
    [h(1:2:end).Marker] = deal('none');
    [h(1:4:end).LineStyle] = deal('--');
    [h(1:2:end).LineWidth] = deal(2);
    ylabel('distractor pull (deg)');
    xlabel('saccadic RT (ms)');
    legend(h(1:2:end), {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'}, 'Location','Best');
    axis([180 403 -11 20]);
    box off
