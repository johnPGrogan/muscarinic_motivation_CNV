% SplitDepAngleByRT

load('TrioEEGBehAnalysis.mat','vars')
[nPP,nConds,nDr,nTr] = size(vars.srt);
%% plot rew eff for depangle for slow/fast RTs
% only distr trials

depAngle = vars.depAngle(:,[3 4],:,:);
srt = exp(vars.srt(:,[3 4],:,:));

% medianRT = nanmedian(x, [4 3 2]);
% rtCutOffs = [medianRT, medianRT];
rtCutOffs = prctile(srt, [50 50], [4 3 2]); % lower and upper thirds

srtByRT = nancat(5, srt, srt);
depAngleByRT = nancat(5, depAngle, depAngle);

toNaN = nancat(5, srt >= rtCutOffs(:,1), srt <= rtCutOffs(:,2));
srtByRT(toNaN) = NaN;
depAngleByRT(toNaN) = NaN;

%% do rt vs depAngle

figure();
cols = get(gca,'ColorOrder');
cols = cols([1 1 2 2],:);
set(gca,'ColorOrder', cols);

[~,~,~,~,h] = conditionalPlot(reshape(permute(srt,[4,1,2,3]), nTr,nPP,4), reshape(permute(depAngle,[4,1,2,3]), nTr,nPP,4));
[h(2:2:end).Visible] = deal('off');
[h(1:2:end).Marker] = deal('none');
[h(1:4:end).LineStyle] = deal('--');
[h(1:2:end).LineWidth] = deal(2);
ylabel('distractor pull (deg)');
xlabel('saccadic RT (ms)');
legend(h(1:2:end), {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'}, 'Location','Best');
axis([180 403 -11 20]);
box off

%% now do diffs on the YDATA from previous

% these are the percentiles
y = reshape(nancat(1,h(1:2:end).YData),2,2,[]);  %[rew drug prctile]
x = reshape(nancat(1,h(1:2:end).XData),2,2,[]);

effNames = {'incentive','drug'};
effLabels = {'Placebo','THP';'0p','50p'};

lnStyles = {'-','-'; '--','-'};

%% also average over each?

figure();

for i = 1:2
    subplot(2,2,i);
    h1 = plot(1:80, nanmean(sq(diff(y,[],i)),1), '-k','LineWidth',2);
    title(effNames{i});
    if i==1
        ylabel('effect on distractor pull');
    end
    xlabel('percentile of RT');
    yline(0, ':k');
    xticklabels(0:25:100);
end

for i = 1:2 % dim
    subplot(2,2,i+2);
    h1 = plot(1:80, sq(diff(y,[],i)), 'LineWidth',2);
    if i==2
        [h1.Color] = deal([0 0 0]);
    end
    h1(1).LineStyle = lnStyles{i,1};
    if i==1
        ylabel('effect on distractor pull');
    end
    xlabel('percentile of RT');
    yline(0, ':k');
    legend(effLabels(i,:), 'Location','Best')
    xticklabels(0:25:100);
end

makeSubplotScalesEqual(2,2);




%%

speeds = {'fast' ,'slow'};
figure();
% cols = [0 0 1; 0 .5 .9; 1 0 0;  .9 .5 0];
cols = get(gca,'ColorOrder');
cols = cols([1 1 2 2],:);

for i = 1:2
    subplot(1,2,i)
    set(gca,'ColorOrder', cols);
    [~,~,~,~,h] = conditionalPlot(reshape(permute(srtByRT(:,:,:,:,i),[4,1,2,3]),nTr,nPP,[]), reshape(permute(depAngle,[4,1,2,3]),nTr,nPP,[]));
    [h(2:2:end).Visible] = deal('off');
    [h(1:2:end).Marker] = deal('none');
    [h(1:4:end).LineStyle] = deal('--');
    [h(1:2:end).LineWidth] = deal(2);
    ylabel('distractor pull (deg)');
    xlabel('saccadic RT (ms)');
%     if i==1
        legend(h(1:2:end), {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'}, 'Location','Best');
%     end
    title(speeds(i));
    axis([180 403 -11 20]);
end


%% 

figure();

for i = 1:2
    subplot(1,2,i);
    errorBarPlot(nanmean(depAngleByRT(:,:,:,:,i),4), 'plotargs', {'LineWidth',2});
    set(gca,'XTick',1:2, 'XTickLabel', {'0p','50p'});
    xlabel('with distractor');
    ylabel('distractor pull (deg)');
    title(speeds{i});
    xlim([.5 2.5]);
    
    
end
makeSubplotScalesEqual(1,2);

%% lme

t = dePivot(permute(depAngleByRT,[1 2 3 5 4]),'KeepNaN',1);
t = nanzscore(t);
t = array2table(t, 'VariableNames', {'pp','rew','drug','late','trial','depAngle'});

lme{1} = fitglme(t, 'depAngle ~ 1 + rew*drug*late + (1 | pp)');

% in each separately
lme{2} = fitglme(t(t.late<0,:), 'depAngle ~ 1 + rew*drug + (1 | pp)');
lme{3} = fitglme(t(t.late>0,:), 'depAngle ~ 1 + rew*drug + (1 | pp)');


%% shift function on depAngle by RT
% get RT deciles per person in each reward condition + drug
% take the reward effect in each decile
% plot average of these for THP + PLACEBO
p = 0:10:100;
rtDeciles = prctile(srt, p, 4); % [pp rew drug]
rtDecileInds = zeros(size(srt));
for i = 1:length(p)-1
    rtDecileInds(srt > rtDeciles(:,:,:,i) & srt <= rtDeciles(:,:,:,i+1)) = i;
end
    
% split 
srtByDeciles = groupMeans(srt, 4, rtDecileInds, 'dim');

% remove 0th 
srtByDeciles(:,:,:,1,:) = []; 

figure();
h = errorBarPlot(reshape(permute(nanmean(srtByDeciles,5),[1,4,2,3]),nPP,[],4),'area',1,'xaxisvalues',p(2:end)-5);
ylabel('sRT (ms)');
xlabel('Decile');
legend([h{:,1}], {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'}, 'Location','Best');

%% split depAngle



% split 
depAngleByDeciles = groupMeans(depAngle, 4, rtDecileInds, 'dim');

% remove 0th 
depAngleByDeciles(:,:,:,1,:) = []; 

figure();
h = errorBarPlot(reshape(permute(nanmean(depAngleByDeciles,5),[1,4,2,3]),nPP,[],4),'area',1,'xaxisvalues',p(2:end)-5);
ylabel('distractor pull (deg)');
xlabel('Decile');

legend([h{:,1}], {'Placebo 0p','Placebo 50p','THP 0p','THP 50p'}, 'Location','Best');


%% rew eff

depAngleRewEff = sq(diff(depAngleByDeciles,[],2)); %[pp drug decile tr]

figure();
h = errorBarPlot(permute(nanmean(depAngleRewEff,4),[1,3,2]),'area',1,'xaxisvalues',p(2:end)-5);
ylabel('rew eff: distractor pull (deg)');
xlabel('Decile');
legend([h{:,1}], {'Placebo','THP'}, 'Location','Best');
