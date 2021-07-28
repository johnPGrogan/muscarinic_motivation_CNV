function  TrihexyVasAnalysis()
% TrihexyVasAnalysis()
% Load up VAS and FinalVAS data, extract the ratings for each question, and
% plot means and individuals
% Would be good to split by condition A or B later

%%

% func to do regexpi on cell and convert to matrix of logicals
cellRegexpi = @(cellArray, pattern) ~cellfun(@isempty, regexpi(cellArray, pattern));

behDataFolder = 'D:/groganj/Work/Trihexy/BehData';
files = what(behDataFolder);

% keep only files with 'tri_'
isPP = cellRegexpi(files.mat, 'VAS_Tri_');
files.mat(~isPP) = [];

% remove any sesNo = 0
is0 = cellRegexpi(files.mat, '0.mat');
files.mat(is0) = [];

% remove pps
toRemove = cellRegexpi(files.mat, 'TRI_015|TRI_036|TRI_033|TRI_031|TRI_027|TRI_024');
files.mat(toRemove) = [];

vasFiles = files.mat;

isPP = cellRegexpi(vasFiles, '^VAS_');
vasFiles(~isPP) = [];


finalVasFiles = files.mat;
isPP = cellRegexpi(finalVasFiles, 'FinalVAS_');
finalVasFiles(~isPP) = [];



%% questions

vasText  = {{'Alert','Drowsy'}
               {'Calm','Excited'}
               {'Strong','Feeble'}
               {'Muzzy','Clear-headed'}
               {'Well-coordinated','Clumsy'}
               {'Lethargic','Energetic'}
               {'Troubled','Tranquil'}
               {'Mentally slow','Quick-witted'}
               {'Tense','Relaxed'}
               {'Attentive','Dreamy'}
               {'Incompetent','Proficient'}
               {'Happy','Sad'}
               {'Antagonistic','Friendly'}
               {'Interested','Bored'}
               {'Withdrawn','Sociable'}
               {'Depressed','Elated'}
               {'Self-centred','Outward-going'}
               };

vasQs = cellfun(@(x) x{2}, vasText, 'UniformOutput',0); % 

finalVasText  = {{'I think I''m on the drug','I think I''m on placebo'}
               {'No headache','Headache'}
               {'No stomach ache','Stomach ache'}
               {'No nausea','Nauseous'}
               {'No dizziness','Dizzy'}
               {'Normal vision','Blurred vision'}
               {'No muscle pain','Muscle pain'}
               {'No muscle twitches','Muscle twitches'}
               {'Good mood','Bad mood'}
               };
finalVasQs = {'placebo', 'headache', 'stomachAche', 'nauseous', 'dizzy', 'blurredVision', 'musclePain', 'muscleTwitches', 'badMood'}';

allVasQs = {vasQs, finalVasQs};

taskNames = {'initial VAS', 'final VAS'};
%% extract scores

vas = VasExtract(vasFiles, behDataFolder);
vas = vas ./ 2; % normalise to -100:100 (from (-200:200))
vas(vas > 100) = 100; % could be outside response bar, set to limits
vas(vas < -100) = -100;

[finalVas,ppIDs] = VasExtract(finalVasFiles, behDataFolder); % second VAS
finalVas = finalVas ./ 2.5; % normalise to -100:100
finalVas(finalVas > 100) = 100; % could be outside response bar, set to limits
finalVas(finalVas < -100) = -100;


allVas = nancat(3, vas, finalVas); % combine

%% figure

fig1 = figure();

for i = 1:2
    subplot(2,1,i)
    errorBarPlot(allVas(:,:,i), 'type', 'bar', 'plotIndividuals', 1);
    set(gca, 'XTick', 1:length(allVasQs{i}), 'XTickLabel', allVasQs{i}, 'XTickLabelRotation', 45);
    title(taskNames{i});
    xlim([0 length(allVasQs{i})+1])
    ylim([-100 100]);
    ylabel('agreement')
    
    % 
    [~,p] = ttest(allVas(:,:,i));
    
    p = double(p<.05);
    p(p==0) = NaN;
    hold on;
    plot(p* max(ylim), '*k','MarkerSize',6);
    
    

end
xlabel('question')

%% split by condition



ppIDs = cellfun(@(x) upper(x(5:11)), vasFiles,'UniformOutput',0);

[~,txt] = xlsread('../../Experiments/Forms/ACh/TrihexyDrugConds.xlsx');

[~,ppInd] = ismember(unique(ppIDs), txt(:,1));
condOrder = txt(ppInd,3:4);
condOrder = col(condOrder'); % into column

condOrderNum = strcmp(condOrder, 'A')+1; % 1=B, 2=A

% isA = [2,1,2,2,1,1,2,1,1,2,2];
% a = [1;2]; b = [2;1];
% condOrder = col([b,a,b,b,a,a,b,a,a,b,b]);
allVasCond = permute( groupMeans(allVas,1, condOrderNum,'dim'), [4,1,2,3] ); % [pp q time drug]
fig2 = figure();

clear p
for i = 1:2
%     for j = 1:2
        subplot(2,1,i)%(i-1)*2+j)
        errorBarPlot(sq(allVasCond(:,:,i,:)), 'type', 'bar');%, 'plotIndividuals', 0);
        set(gca, 'XTick', 1:length(allVasQs{i}), 'XTickLabel', allVasQs{i}, 'XTickLabelRotation', 45);
        title(taskNames{i});
        xlim([0 length(allVasQs{i})+1])
        ylim([-100 100]);
        ylabel('agreement')

        % 
        [~,p{i}] = ttest(allVasCond(:,:,i,1), allVasCond(:,:,i,2));

        p2 = double(p{i}<.05);
        p2(p2==0) = NaN;
        hold on;
        plot(p2 * max(ylim)*.9, '*k','MarkerSize',6);
    
%     end
end
xlabel('question')

legend({'Placebo','Drug'}, 'Location','Best'); % A is drug


%% raincloud plots

cols2 = get(gca, 'ColorOrder');

alpha = .5;
limits2 = [-100 100];

for k = 1:2
    figs{k} = figure();
    clear h;
    
    n = sum(~isnan(allVasCond(1,:,k,1)));
    hold on;
    for j = 1:n
        for i = 1:2
            data = allVasCond(:,j,k,i);
            if any(data,'all')

                [y,x] = ksdensity(data, linspace(limits2(1), limits2(2), 100), 'BandWidth', diff(limits2)/20);

                x1 = [x, fliplr(x)]; % mirror it
                y1 = [y, zeros(size(y))]; % put zeros
                y1 = y1 ./ max(y1,[],'all'); % normalise
                y1 = y1./2 + j;
                h(i) = fill(x1, y1, cols2(i,1:3),'FaceAlpha',alpha);

        %             plot(nanmean(data), j+1,'o','Color',cols2(i,:),'LineWidth',1);

                % scatter
                plot(data, j + randn(size(data))*.1, '.', 'Color', cols2(i,:));
            end
        end
    end


    set(gca,'YTick',1:n, 'YTickLabel',allVasQs{k});

    ylim([.5 n+1]);
    xlim(limits2);

    xline(0,'-k');
    %     xlabel(parNames{iP});
    ylabel('question');

    legend((h), {'Placebo','Drug'},'Location','NorthEast');
    
end

%% dichotomise placebo

guessedPlacebo = sq(allVasCond(:,1,2,:)) > 0;

% get counts: rows = answers, cols = conditions
obsData = [sum(guessedPlacebo); sum(~guessedPlacebo)];

[chiP,chiSq,chiDf] = ChiSquareTest(obsData);




%% save

% saveas(fig2, 'TrihexyVAS.fig');
saveas(fig2, 'TrihexyVAS.jpg');

save('TrihexyVasAnalysis.mat', 'allVas','allVasCond',...
    'vasText','vasQs','finalVasQs','finalVasText',...
    'ppIDs','vasFiles','finalVasFiles','p');

end

function [vas, ppIDs, resStruct] = VasExtract(files, behDataFolder)

%% get ppIDs

ppIDs = cellfun(@(x) x(end-8:end-4), files, 'UniformOutput',0);


sesNo = str2num(cellfun(@(x) x(end), ppIDs));
ppNums = cellfun(@str2num, cellfun(@(x) x(1:3), ppIDs,'UniformOutput',0));

%% load

results = cellfun(@(x) x.result.data, cellfun(@load, fullfile(behDataFolder, files), 'UniformOutput', 0), 'UniformOutput',0);

resStruct = sq(nancat(results)); % into struct array

resQ = [resStruct.question];
vasRating = [resStruct.vas];
% results = groupMeans(results, 1, sesNo, 'dim');

vas = groupMeans(vasRating, 2, resQ, 'dim');

end
