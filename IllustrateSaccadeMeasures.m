% IllustrateSaccadeMeasures
% draw some saccades (times vs pos) and illustrate sRT, endptvar, velres,
% depAngle

set(0, 'DefaultAxesFontSize',14)
%%

r = load('TrioEEGBehProcessFiles.mat');

%%
result = r.result(10);

%% get trials starting from bottom, going right

trialInds = result.d.previousTargetPos==2 & result.d.targetDir==1 & result.d.distractor==1;

isCaptured = result.distrCapture(trialInds);

% get some trials
normTraj = result.normTraj(trialInds,:);

% rotate by -60deg (right is positive now)
normTraj = normTraj .* exp(-1/3 * pi * 1i); 

% plot(permute(normTraj,[2,1]))


%% after sRT, add final value until end - looks like fixed position
x = real(normTraj);
% t1 = [10 90];
sRT = result.srt(trialInds);
sDur = result.info.sDur(trialInds);
sEndT = sRT + sDur;
for i = 1:length(sRT)
%     x1 = nanmean(x(i,srt(i)+t1(1):srt(i)+t1(2)));
    x(i,sRT(i)+sDur(i):end) = NaN;
    endpts(i,1) = x(i,sRT(i)+sDur(i)-1);
end

%% % plot t vs real pos

nTrials = 5;

goodTr = find(~isCaptured,100);
% hand pick these
goodTr = goodTr([1 2 3 4 8 12]);
badTr = find(isCaptured,nTrials);

t = 1:350;

lw = {'LineWidth',2};

figure();
ax = [1 t(end)+50 -7 9];
% fill([ax(1) ax(2) ax(2) ax(1)], [ax(3) ax(3) ax(4) ax(4)], [0 0 0]);
axis(ax);
hold on

h2 = plot(x(badTr,t)','g-',lw{:});
h1 = plot(x(goodTr,t)','b-',lw{:});
xlabel('time from target onset (ms)');
ylabel('horizontal eye position (deg)');

box off
xlim([t(1)  t(end) + 50]);

% super endpts
plot(sEndT(badTr), endpts(badTr), 'go','MarkerSize',6,lw{:});
h = plot(sEndT(goodTr), endpts(goodTr), 'bo','MarkerSize',6,lw{:});


% % super ksdens
[yy, xx] = ksdensity(endpts(goodTr));%, 4:.1:8,'Support',[4 8]);
fill(t(end) + 10 + [yy, fliplr(yy)]*30, [zeros(size(xx)), xx], 'b','FaceAlpha',1);

% % means
% plot(t(end) + 20, nanmean(endpts(badTr)), 'rx','MarkerSize',10,lw{:});
plot(t(end) + 10, nanmean(endpts(goodTr)), 'kx','MarkerSize',10,lw{:});

% y2 = interp1(xx,yy,nanmean(endpts(goodTr))) .*30;
plot(t(end)+ [10 40], repmat(nanmean(endpts(goodTr)),1,2), ':k','MarkerSize',10,lw{:});
% 
% % join one to mean
[~,i] = max([h.YData]);
plot(repmat(t(end) + 40,2,1), [h.YData(i) nanmean(endpts(goodTr))], '-k',lw{:});
% % plot(t(end) + 20, h(i).YData , 'ow',lw{:});
plot([sEndT(goodTr(i)) t(end)+40], [h.YData(i) h.YData(i)], ':k',lw{:});


legend([h1(1) h2(1)], {'Correct','Distractor Capture'}, 'Location','SouthWest');


%% without endpt ks

figure();
ax = [1 t(end) -7 9];
% fill([ax(1) ax(2) ax(2) ax(1)], [ax(3) ax(3) ax(4) ax(4)], [0 0 0]);
axis(ax);
hold on

h2 = plot(x(badTr,t)','g-',lw{:});
h1 = plot(x(goodTr,t)','b-',lw{:});
xlabel('time from target onset (ms)');
ylabel('eye position (deg)');

box off
xlim([t(1)  t(end)]);
yticks([-6 0 6]);
yticklabels([-11 0 11]);

% super endpts
plot(sEndT(badTr), endpts(badTr), 'go','MarkerSize',6,lw{:});
h = plot(sEndT(goodTr), endpts(goodTr), 'bo','MarkerSize',6,lw{:});


legend([h1(1) h2(1)], {'Correct','Distractor Capture'}, 'Location','SouthWest');

%% plot full pos, show 

normTraj2 = normTraj .* exp(1/3 * pi * 1i); 

locs = (result.locs - result.locs(1,:));
cols = [.5 .5 .5; .2 .2 .2; .8 .8 .8];

f = figure();
hold on;
fill([-4 14 14 -4], [-14 -14 4 4], [0 0 0]);
for i = 1:3
    DrawCircles(locs(i,:), 2.5, 1, cols(i,:),'LineWidth',1);
end

rowPlot(locs(1,:)', 'w+','LineWidth',2,'MarkerSize',10)
rowPlot(locs(2:3, :)', 'k+','LineWidth',2,'MarkerSize',10)

plot(permute(normTraj2(goodTr(1:3),:),[2,1]), '-b',lw{:});
plot(permute(normTraj2(badTr(1:3),:),[2,1]), '-g',lw{:});

rowPlot(locs(1:2,:)','--w');

axis off
%% resid vel

figure();
inds = 1:120;
col = 'k';%[0 0.4470 0.7410];

ax = [7 14 300 550];
% fill([ax(1) ax(2) ax(2) ax(1)], [ax(3) ax(3) ax(4) ax(4)], [0 0 0]);
hold on;
h = plot(r.result(3).ampl(inds), r.result(3).vel(inds)*1000,'x','LineWidth',1.5,'Color', col); 
hold on;
% main seq lines
l = lsline; 
l.LineWidth = 1.5;
l.LineStyle = '--';
l.Color = 'k';


% show residuals
PlotResidLines(r.result(3).ampl(inds), r.result(3).vel(inds)*1000, '-', 'Color',col);

% % plot invis thing for legend
% h(4) = plot(NaN, '--', 'Color', 'k','LineWidth',1.5);
% h(5) = plot(NaN, '-', 'Color', col);

% legend(h([4 5 1]), {'main sequence', 'residual velocity', 'single trial'}, 'Location','NorthWest');

xlabel('saccade amplitude (deg)')
ylabel('peak saccade velocity (deg/s)')
box off
yticks(300:100:600);

%% plot thp + placebo

figure();
inds = 1:60;
cols = {'b','r'};%[0 0.4470 0.7410];

ax = [7 14 300 710];

for pp = 18
%     subplot(4,5,pp);
    for i = 1:2
        iPP = pp*2 - 1*(i==1);
        col = cols{i};
% fill([ax(1) ax(2) ax(2) ax(1)], [ax(3) ax(3) ax(4) ax(4)], [0 0 0]);
hold on;
h(i) = plot(r.result(iPP).ampl(inds), r.result(iPP).vel(inds)*1000,'x','LineWidth',1.5,'Color', col); 
hold on;

% show residuals
PlotResidLines(r.result(iPP).ampl(inds), r.result(iPP).vel(inds)*1000, '-', 'Color',col);

% % plot invis thing for legend
% h(4) = plot(NaN, '--', 'Color', 'k','LineWidth',1.5);
% h(5) = plot(NaN, '-', 'Color', col);

% legend(h([4 5 1]), {'main sequence', 'residual velocity', 'single trial'}, 'Location','NorthWest');
    end

    % main seq lines
l = flip(lsline); 
[l.LineWidth] = deal(1.5);
[l.LineStyle] = deal('--');
[l.Color] = deal(cols{:});

    axis(ax);
xlabel('saccade amplitude (deg)')
ylabel('peak saccade velocity (deg/s)')
box off
yticks(300:100:600);
end
legend(h, {'Placebo','Drug'});