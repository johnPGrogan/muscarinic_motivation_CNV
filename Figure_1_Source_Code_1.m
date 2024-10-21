function PlotFigure1(vars,varNames)
% Plot Figure 1c-d
% if no inputs, will load figure source data file
    
    if nargin==0
        load('Figure 1 - Source Data 1.mat', 'vars', 'varNames')
        % vars is struct with fields of [nPP, nConds, nDrug] trial means
    end
    
        
    lineColours = get(gca,'ColorOrder');
    labels = {'0p','50p'; 'no distractor', 'with distractor'; 'Placebo','THP'};
    ylabels.velres = 'residual peak velocity (deg/s)';
    ylabels.srt = 'saccadic RT (ms)';
    ylabels.depAngle = 'distractor pull (deg)';
    
    
    plotIndivs = 1;
    
    for i = 1:length(varNames)
           
        figure();
    %     subplot(ceil(sqrt(nVars)),ceil(nVars/ceil(sqrt(nVars))),i);
        h = errorBarPlot(vars.(varNames{i}), 'type','line','plotargs',{'LineWidth',2,'LineStyle','none'});
        hold on;
        for j = 1:2
            plot(h(j).XData(1:2), h(j).YData(1:2), 'LineWidth', 2, 'Color', lineColours(j,:));
            plot(h(j).XData(3:4), h(j).YData(3:4), 'LineWidth', 2, 'Color', lineColours(j,:));
        end
        
        box off
        set(gca,'XTick',1:4,'XTickLabel',labels(1,:));
        ylabel(ylabels.(varNames{i})) 
        xlabel([labels{2,1} '       ' labels{2,2}]);
        xlim([.5 4.5])
        
        
        if plotIndivs
            for j = 1:2
                x = (1:4) + (rand(20,4)-.5).*.1 + (3 - 2*j)/20;
                plot(x, nanmean(vars.(varNames{i})(:,:,j,:),4),'x','Color',h(j).Color,'MarkerSize',10);
                % join up pairs?
                for k = 1:2
                    plot(x(:,k*2-1:k*2)', nanmean(vars.(varNames{i})(:,k*2-1:k*2,j,:),4)', '-', 'Color', [h(j).Color .2]);
                end
            end
        end
    end
    legend(labels(3,:),'Location','Best')
