% CompareBetaPerToi

%% get cnv

load('TrioEEGBehAnalysis.mat','trialTab','varNames');
nVars=3;
vNames = {'velocity','RT','distractor pull'};

load('trioEEGAnalyseByTrials_readyBS.mat','meanAmplTrialsPost','meanAmplTrials','eegByCond','xTimes')
trialTab.CNV = nanzscore(col(meanAmplTrialsPost(:,:,:,:,30,:)));
trialTab.P3a = nanzscore(col(meanAmplTrials(:,:,:,:,30,:)));


a = load('trioEEGAnalyseByTrials_audCue_200.mat','meanAmplTrialsPost','eegByCond','xTimes');
trialTab.incentive = nanzscore(col(a.meanAmplTrialsPost(:,:,:,:,30,:)));


toiNames = {'incentive','P3a','CNV'};
nTois = length(toiNames);

extractVar = @(x,inds) [x.Coefficients.Estimate(inds), x.Coefficients.SE(inds), x.Coefficients.tStat(inds), x.Coefficients.pValue(inds)];

formula = '%s ~ 1 + rew*distr*drug + %s + (1 | pp)';
toiRegStats = zeros(nVars,nTois,8,4);
for i = 1:nTois
    for j = 1:nVars
        reg = fitglme(trialTab, sprintf(formula, varNames{j}, toiNames{i}));
        
        toiRegStats(j,i,:,:) = extractVar(reg, 2:9);
    end
end


termNames = reg.CoefficientNames(2:end);
drugNames = {'Placebo','THP'};

%% plot trace

x1 = reshape(permute(nanmean(reshape(a.eegByCond(30,isBetween(a.xTimes,[900 1100]),:,:,:,:),[],2,2,20,2,122),[6,3]),[4,1,2,5,3]),20,[],4);
x2 = reshape(permute(nanmean(reshape(eegByCond(30,:,:,:,:,:),[],2,2,20,2,122),[6,3]),[4,1,2,5,3]),20,[],4);

drugRewNames = {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'};
toi = [200 280; 1200 1500]; % times of interest

toiCols = [.5 1 .5; .5 .5 1; 1 .5 .5];
cols = [0 0 1; 0 .5 .9;1 0 0;  .9 .5 0];

fs = 24; % font size
set(0, 'DefaultAxesFontSize',fs);
figure();
% subplot(4,1,1:2)
set(gca, 'DefaultAxesFontSize',fs,'ColorOrder',cols)

xline(0,'--','Color',[1 1 1].*.2, 'LineWidth',2);
yline(0,'--','Color',[1 1 1].*.2, 'LineWidth',2);
hold on;
fill([-500 -300 -300 -500], [-10 -10 6 6], toiCols(1,:), 'FaceAlpha', .5, 'LineStyle','none');
fill(toi(1,[1 2 2 1]), [-10 -10 6 6], toiCols(2,:), 'FaceAlpha', .5, 'LineStyle','none');
fill(toi(2,[1 2 2 1]), [-10 -10 6 6], toiCols(3,:), 'FaceAlpha', .5, 'LineStyle','none');

h = errorBarPlot(x2, 'area',1,'xaxisvalues',xTimes, 'plotargs', {'LineWidth',3});
h1 = errorBarPlot(x1, 'area',1,'xaxisvalues',a.xTimes(isBetween(a.xTimes, [900 1100]))-1400, 'plotargs', {'LineWidth',3});
xlim([-500 1500]);
legend([h{:,1}], drugRewNames,'Location','South'); % overwrite to remove grey boxes
title('Cz');
ylabel('\muV');
xlabel('time from preparation cue (ms)');

axis([-500 1500 -10 6]);

%% plot each mean

x = nancat(2, sq(nanmean(nanmean(a.meanAmplTrialsPost(:,:,:,:,30,:),6),3)), ...
    sq(nanmean(nanmean(meanAmplTrials(:,:,:,:,30,:),6),3)),...
    sq(nanmean(nanmean(meanAmplTrialsPost(:,:,:,:,30,:),6),3)));

figure();
for i = 1:3
    subplot(3,3,i+3);cla
    set(gca, 'DefaultAxesFontSize',fs)
    h = errorBarPlot(x(:,i*2-1:i*2,:), 'plotargs',{'LineWidth',3});
%     JoinUpErrorPoints(h, [1 2; 3 4; 5 6]);
    box off
    set(gca,'XTick',1:2,'XTickLabel',{'0p','50p'});
    % ylabel(sprintf('mean uV from %d:%dms', amplWindow(1), amplWindow(2))); 
    
    xlabel('incentive');
    xlim([.75 2.25])
    if i==1
        legend(h, drugNames, 'Location','Best');
        ylabel('mean \muV');
    end
    title(toiNames{i});
    yl = mean(ylim);
    ylim(yl + [-2 2]);
    %     yt = yticks;
%     yticks(linspace(min(yt),max(yt),3));
    yticks(-10:1:10);
end
% makeSubplotScalesEqual(3,3,4:6);

% plot

alpha = .05 /9;

% figure();
for j = 1:3
    subplot(3,3,6+j);cla
    set(gca, 'DefaultAxesFontSize',fs)
    h = barwitherr(toiRegStats(:,j,4,2)*1.96, toiRegStats(:,j,4,1));
    h.FaceColor = toiCols(j,:);
    [h.LineWidth] = deal(1.2);
    xticklabels(vNames(1:nVars));
    if j==1
        ylabel('\beta coefficient')
    end
    box off

%     for i = 1:nTois


        stars = col(toiRegStats(:,j,4,4)) < alpha;
    %     y = repmat(max([h.YData])*1.1, 1,3);    
        y = repmat(.15, 1,3);
        y(~stars) = NaN;
        hold on;
        plot(1:3, y, '*','MarkerSize',14, 'Color', h.FaceColor, 'LineWidth',1.5);

%     end

%     legend(h, toiNames, 'Location','Best');
    ylim([ (min(ylim)) .16])
    yticks([0 .1]);
end

makeSubplotScalesEqual(3,3,7:9)

%% mean ERPs + indivs

x = nancat(2, sq(nanmean(nanmean(a.meanAmplTrialsPost(:,:,:,:,30,:),6),3)), ...
    sq(nanmean(nanmean(meanAmplTrials(:,:,:,:,30,:),6),3)),...
    sq(nanmean(nanmean(meanAmplTrialsPost(:,:,:,:,30,:),6),3)));

figure();
for i = 1:3
    subplot(2,2,i);cla
    set(gca, 'DefaultAxesFontSize',fs)
    h = errorBarPlot(x(:,i*2-1:i*2,:), 'plotargs',{'LineWidth',2});
    
    hold on;
    for j = 1:2
        plot(rand(20,2)'/20 + [.9 2.05]', x(:,i*2-1:i*2,j)', '-x', 'Color', [h(j).Color .5]);
    end

    box off
    set(gca,'XTick',1:2,'XTickLabel',{'0p','50p'});
    % ylabel(sprintf('mean uV from %d:%dms', amplWindow(1), amplWindow(2))); 
    
    xlabel('incentive');
    xlim([.75 2.25])
    ylabel('mean \muV');
    if i==1
        legend(h, drugNames, 'Location','Best');
    end
    title(toiNames{i});
end
% makeSubplotScalesEqual(3,3,4:6);
