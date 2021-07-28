% TrihexyQuestionnaireAnalysis
% Analyse demographics + questionnaires
clc; clear all; close all;


% file name
qFile = 'C:\Users\groganj\OneDrive - Nexus365\Comp_Neurology_Lab\Data\Trihexy\Questionnaires\TrihexyQuestionnaires.xlsx';
[num,txt,raw] = xlsread(qFile,'Sheet1','A1:BW26');

%% get into a table

colNames = raw(1,:);
t = cell2table(raw(2:end,:),'VariableNames',colNames);

xlsID = t.ppID;

% exclude pps
load('TrioEEGBehAnalysis.mat','uniqID')

toKeep = ismember(xlsID, uniqID);
t(~toKeep,:) = []; % remove them

%% get means etc

isDrugFirst = strcmp(t.session1, 'A'); % is drug first

ageMean = nanmean(t.age);
ageSD = nanstd(t.age);


allDelay = [t.delay_s1, t.delay_s2] .* 1440; % convert to minutes

allDelayMean = nanmean(allDelay(:));
allDelaySD = nanstd(allDelay(:));

% get [placebo drug]
delayDrug = allDelay;
delayDrug(isDrugFirst,:) = delayDrug(isDrugFirst,[2 1]); % [plac drug]

delayDrugMeans = nanmean(delayDrug);
delayDrugSD = nanstd(delayDrug);
[~,p] = ttest(delayDrug(:,1), delayDrug(:,2));

%% questionnaires

% total is first column

das = [t.dasTotal, t.dasExec, t.dasEmot, t.dasBCI];
bis = [t.bisTotal, t.bisNonPlan, t.bisMotor, t.bisAtten];
dass = [t.dassTotal, t.dassDep, t.dassAnx, t.dassStress];

qAll = nancat(3, das, bis, dass);

qMeans = nanmean(qAll);
qSD = nanstd(qAll);


%% save

save('TrihexyQuestionnaireAnalysis.mat');