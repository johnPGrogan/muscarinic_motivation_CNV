% TrioEEGBehHists
% histogram depangle + others
clc; clear all; close all;
%% hist vars

load('TrioEEGBehAnalysis.mat','vars','varNames','ylabels','labels');
varNames = varNames([1 2 3]);
nVars = 3; % remove endptvar
vars.srt = exp(vars.srt);
drugNames = {'placebo','drug'};
condNames = {'targLo','targHi','distLo','distHi'};
[nPP,nConds,nDr,nTr] = size(vars.velres);
%%
figure();
for iDr = 1:2
    for iC = 1:4
        subplot(2,4,(iDr-1)*4+iC);

        histogram(col(vars.depAngle(:,iC,iDr,:)),100,'Normalization','probability');
        if iDr==1
            title(condNames{iC});
        else
            xlabel('departure angle (deg)')
        end
        if iC==1
            ylabel(drugNames{iDr});
        end
    end
end

makeSubplotScalesEqual(2,4)

%%

% close all
% varNames = fieldnames(vars); % to do all vars
inds = [1 3 2 4];
ksArgs = {};
clear h
for iV = 1:length(varNames)
    figure();
    for iC = 1:4
        subplot(2,2,inds(iC));
        hold on;
        for iDr = 1:2
            [y,x] = ksdensity(col(vars.(varNames{iV})(:,iC,iDr,:)),ksArgs{:});
            h(iDr) = plot(x,y,'LineWidth',1.2);
        end
        title(condNames{iC});
%         xlabel('departure angle (deg)')
%         xlim([-180 180]);
%         xticks(-180:60:180);
    end

    legend(h, drugNames,'Location','Best');
    makeSubplotScalesEqual(2,2);
    SuperTitle(varNames{iV});
end


%% drug eff - per person

angles = linspace([-200; exp(4.5); -90; -6;], [200; exp(7); 90; 2.3], 100);

varFreqs = zeros(nPP,length(angles),nConds,nDr,nVars);

for iV = 1:nVars
    for iPP = 1:nPP
        for iC = 1:nConds
            for iDr = 1:nDr
                varFreqs(iPP,:,iC,iDr,iV) = ksdensity(sq(vars.(varNames{iV})(iPP,iC,iDr,:)), angles(iV,:));
                %[pp angle cond drug]
            end
        end
    end
end
%% plot - error bar across

figure();
for iV = 1:nVars
    for iDr = 1:2
        subplot(2,nVars,(iDr-1)*nVars+iV);
        h = errorBarPlot(varFreqs(:,:,:,iDr,iV), 'area',1, 'xaxisvalues',angles(iV,:));

        xlabel(varNames{iV});
        ylabel('average density');
        title(drugNames{iDr});
    end
end
legend([h{:,1}], condNames, 'Location','Best');

%% diffs

varFreqs2 = reshape(varFreqs, nPP, [],2,2,2,nVars); %[pp angles rew distr drug var]
varFreqsEff = nancat(5, sq(nanmean(varFreqs2,[4 5])),... %[pp angle level var eff
    sq(nanmean(varFreqs2,[3 5])),...
    sq(nanmean(varFreqs2,[3 4])));
% labels = {'0p', '50p'; 'target', 'distractor'; 'placebo', 'drug'};



diffNames = {'reward','distractor','drug'};

figure();
c = get(gca,'ColorOrder');
c = c([3 4 5 6 1 2 7],:);
for iV = 1:nVars
    for iEff = 1:3
        subplot(3,nVars,(iEff-1)*nVars+iV);
        set(gca,'ColorOrder',c([iEff*2-1 iEff*2 7],:));
        h = errorBarPlot(cat(3, varFreqsEff(:,:,:,iV,iEff), diff(varFreqsEff(:,:,:,iV,iEff),[],3)),'area',1,'xaxisvalues', angles(iV,:));
        yline(0,':k');
        xlim(minMax(angles(iV,:),2));
        if iV==1
            ylabel(diffNames{iEff});
        elseif iV==4
            legend([h{:,1}], [labels{iEff,:} {'effect'}], 'Location','Best');
        end
        if iEff==1
            title(varNames{iV});
        elseif iEff==3
            xlabel(varNames{iV});
        end
    end
%     makeSubplotScalesEqual(3,nVars,[1 1+nVars 1+nVars*2] + iV - 1)
end

%% plot grand average hist, then one plot of effects, and one of interactions

% make interactions - e.g. rew eff in each drug, distr eff in each drug,
% rew eff in each distr

varInteractions = nancat(5, sq(diff(nanmean(varFreqs2,4),[],3)),... % av over dist, rew eff
                            sq(diff(nanmean(varFreqs2,3),[],4)),... % av over rew, distr eff
                            sq(diff(nanmean(varFreqs2,5),[],3)));   % av over drugs, rew eff
% [pp angles levels inter/eff var]                        
interNames = {  'placebo: rew eff' ,'drug: rew eff';
                'placebo: distr eff', 'drug: distr eff';
                'targ: rew eff', 'distr: rew eff'};
interacPerV = [1 1 2];

figure();
for iV = 1:3
%     figure();

    subplot(3,3,(iV));
    errorBarPlot(nanmean(varFreqs(:,:,:,:,iV),[4 3]),'area',1,'xaxisvalues',angles(iV,:));
%     xlabel(varNames{iV});
    ylabel('Density');
    xlim(minMax(angles(iV,:),2));
    title(varNames{iV});

    subplot(3,3,(iV)+3);
    h = errorBarPlot(sq(diff(varFreqsEff(:,:,:,iV,:),[],3)),'area',1,'xaxisvalues',angles(iV,:));
    yline(0,':k');
    xlim(minMax(angles(iV,:),2));
%     xlabel(varNames{iV});
    ylabel('effect');
    if iV==1
        legend([h{:,1}], diffNames, 'Location','Best');
    end


    subplot(3,3,(iV)+6);
    h = errorBarPlot(varInteractions(:,:,:,interacPerV(iV),iV),'area',1,'xaxisvalues',angles(iV,:));
    yline(0,':k');
    xlim(minMax(angles(iV,:),2));
    xlabel(ylabels{iV});
    ylabel('effect');
    if iV==1
        legend([h{:,1}], interNames(interacPerV(iV),:), 'Location','Best');
    end

end


%% polar plots

figure();
iV=3;

for iF = 1:3 % drug level
    subplot(2,2,iF);
    h = polarplot(deg2rad(-angles(iV,:))',sq(nanmean(varFreqsEff(:,:,:,iV,iF),[1]))*10 , '-');
    title(diffNames(iF));
end

legend(h, labels(2,:),'Location','Best');

%% depAngle only 
% distr trials only

clusterArgs = {'cluster',1,'clustermethod','mass','two_tailed',true,'nperms',5000};

set(0, 'DefaultAxesFontSize',16)
figure();
inds = [2];

for j = 1:length(inds)
    depFreqs = sq(varFreqs2(:,:,:,inds(j),:,3)); %[pp angles rew thp] - distr only
    % depFreqs = sq(nanmean(varFreqs2(:,:,:,:,:,3),4)); %[pp angles rew thp] - all trials

    % one figure of grand mean 

%     subplot(3,length(inds),j)
    set(gca,'ColorOrder', [0 0 0]);
    h = errorBarPlot(nanmean(depFreqs,[4 3]),'area',1,'xaxisvalues',angles(3,:),'plotargs',{'LineWidth',3});
    xlabel('Distractor pull (deg)');
    if j==1; ylabel('Density');end
    xlim(minMax(angles(3,:),2));
%     yline(0,':k','LineWidth',3);
    xline(60,'--k','LineWidth',3);
    xline(0,'-k','LineWidth',3);
    legend([h{:,1}], {'Grand Mean'}, 'Location','Best');
    set(gca,'XTick',-90:30:90)
    box off
    
%     title(labels(2,j));
    
    
    % one with effects of rew,thp,rew*thp
    depFreqEffs = nancat(4, diff(nanmean(depFreqs,4),[],3), diff(nanmean(depFreqs,3),[],4),...
        sq(diff(depFreqs,[],3))); %[pp angle eff1 eff2]

    figure();
%     subplot(3,length(inds),j+(length(inds)));
    c = get(gca,'ColorOrder');
    set(gca,'ColorOrder', c(3:4,:));
    h = errorBarPlot(sq(depFreqEffs(:,:,1,1:2)),'area',1,'xaxisvalues',angles(3,:),'plotargs',{'LineWidth',3});
    xlabel('Distractor pull (deg)');
    if j==1; ylabel('Density'); end
    xlim(minMax(angles(3,:),2));
    yline(0,':k','LineWidth',3);
    xline(60,'--k','LineWidth',3);
    xline(0,'-k','LineWidth',3);
    set(gca,'XTick',-90:30:90)
    ylim([-.004 .0035]);
    for i = 1:2
        [~,p] = permutationOLS(depFreqEffs(:,:,1,i),[],[],[],clusterArgs{:});
        hold on;
        pbar(p, 'xVals', angles(3,:), 'yVal', -.004,'plotargs',{'LineWidth',5,'Color',h{i,1}.Color});
    end
    legend([h{:,1}], {'Incentive effect', 'THP effect'}, 'Location','Best');
    box off
    % thp*rew

    figure();
%     subplot(3,length(inds),j+length(inds)*2)
    h = errorBarPlot(depFreqEffs(:,:,:,3),'area',1,'xaxisvalues',angles(3,:),'plotargs',{'LineWidth',3});
    xlabel('Distractor pull (deg)');
    if j==1; ylabel('Incentive effect: Density'); end
    xlim(minMax(angles(3,:),2));
    yline(0,':k','LineWidth',3);
    xline(60,'--k','LineWidth',3);
    xline(0,'-k','LineWidth',3);
    set(gca,'XTick',-90:30:90)
    for i = 1:2
        [~,p] = permutationOLS(depFreqEffs(:,:,i,3),[],[],[],clusterArgs{:});
        hold on;
        pbar(p, 'xVals', angles(3,:), 'yVal', -.004,'plotargs',{'LineWidth',5,'Color',h{i,1}.Color});
    end
    ylim([-.004 .0035]);
    legend([h{:,1}], {'Placebo', 'THP'}, 'Location','Best');
end

% makeSubplotScalesEqual(3,2,3:6)