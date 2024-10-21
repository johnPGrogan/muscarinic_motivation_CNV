function PlotFigure2(vars, varNames)
% Plot figure 2

%% load up if not given

    if nargin==0
        load('Figure 2 - Source Data 1.mat', 'vars')
        % vars is struct with fields of [nPP, nConds, nDrug, nTr] 
        varNames = {'depAngle'}; % just this
    end


%% first plot fig2b - calling PlotFigure1 on means

    varMeans = struct();
    varMeans.depAngle = nanmean(vars.depAngle, 4); % mean over trials
    PlotFigure1(varMeans, varNames);


%% Figure 2c

    [nPP,nConds,nDr,nTr] = size(vars.depAngle);
    
    % calculate varFreqs per pp/cond/drug
    angles = linspace(-90, 90, 100);
    varFreqs = zeros(nPP,length(angles),nConds,nDr);
    for iPP = 1:nPP
        for iC = 1:nConds
            for iDr = 1:nDr
                varFreqs(iPP,:,iC,iDr) = ksdensity(sq(vars.depAngle(iPP,iC,iDr,:)), angles);
                %[pp angle cond drug]
            end
        end
    end
    % reshape
    varFreqs2 = reshape(varFreqs, nPP, [],2,2,2); %[pp angles rew distr drug var]
    
    depFreqs = sq(varFreqs2(:,:,:,2,:)); %[pp angles rew thp] - distractor-present trials only
    
    
    % set(0, 'DefaultAxesFontSize',16)
    figure();
    set(gca,'ColorOrder', [0 0 0]);
    h = errorBarPlot(nanmean(depFreqs,[4 3]),'area',1,'xaxisvalues',angles,'plotargs',{'LineWidth',3});
    xlabel('Distractor pull (deg)');
    ylabel('Density');
    xlim(minMax(angles,2));
    xline(60,'--k','LineWidth',3);
    xline(0,'-k','LineWidth',3);
    legend([h{:,1}], {'Grand Mean'}, 'Location','Best');
    set(gca,'XTick',-90:30:90)
    box off
    
    
    %% Figure 2d - Incentive and THP effects

    clusterArgs = {'cluster',1,'clustermethod','mass','two_tailed',true,'nperms',5000};


    % one with effects of rew,thp,rew*thp
    depFreqEffs = nancat(4, diff(nanmean(depFreqs,4),[],3), diff(nanmean(depFreqs,3),[],4)); %[pp angle eff1 eff2]
    
    figure();
    %     subplot(3,length(inds),j+(length(inds)));
    c = get(gca,'ColorOrder');
    set(gca,'ColorOrder', c(3:4,:));
    h = errorBarPlot(sq(depFreqEffs(:,:,1,:)),'area',1,'xaxisvalues',angles,'plotargs',{'LineWidth',3});
    xlabel('Distractor pull (deg)');
    ylabel('Density'); 
    xlim(minMax(angles,2));
    yline(0,':k','LineWidth',3);
    xline(60,'--k','LineWidth',3);
    xline(0,'-k','LineWidth',3);
    set(gca,'XTick',-90:30:90)
    ylim([-.004 .0035]);
    for i = 1:2
        [~,p] = permutationOLS(depFreqEffs(:,:,1,i),[],[],[],clusterArgs{:});
        hold on;
        pbar(p, 'xVals', angles, 'yVal', -.004,'plotargs',{'LineWidth',5,'Color',h{i,1}.Color});
    end
    legend([h{:,1}], {'Incentive effect', 'THP effect'}, 'Location','Best');
    box off

