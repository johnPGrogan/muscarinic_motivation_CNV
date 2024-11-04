% TrioEEG_PlotAllFigures
% 
% For each figure in the paper, load up the saved/processed data (source
% files), and plot the figures
% 

clc; clear; close all; 
set(0, 'DefaultAxesFontSize', 14)
%% code required
% matlab stats:
%   ksdensity.m
% ERPLab:
%   topoplot.m
% matlib:
%   scatterRegress, nancat, permutationOLS, conditionalPlot, errorBarPlot, sq
% my functions:
%   minMax, isBetween, PlotTopoAndTrace, pbar

%% data files required
% TrioEEGBehAnalysis.mat, 'vars','varNames','ylabels','labels'
% chanlocs.mat
% TrioEEGAnalyseByTrials_readyBS.mat, 'eegByCond','chanNames','timesToPlot','xTimes','toi', 'GNDEffects'
% trioEEGAnalyseByTrials_audCue200.mat, 'eegByCond','chanNames','timesToPlot','xTimes','toi', 'GNDEffects'
% WholeBrainRegressions_perm_Fig4.mat, 'allBetas','allPVals','xTimes2'
% WholeBrainRegressions3_perm_Fig4_S1.mat, 'allBetas','allPVals','xTimes2'
% WholeBrainRegressions_Eye_perm_Fig4_S2.mat, 'allBetas','allPVals','xTimes2','ind'
% TrioEEGMediations_readyBS_srt.mat, 'mediationTable','medModTable'
% TrioEEGMediations_readyBS_velres.mat, 'mediationTable','medModTable'


%% Figure 1

    % Figure 1 - Source Data 1.mat
    load('TrioEEGBehAnalysis.mat','vars');
    varNames = {'velres','srt'}; % just velres + srt
    vars.srt = exp(vars.srt); % log -> raw RT
    vars = structfun(@(x) nanmean(x,4), vars,'Uni',0); % take means over trials
    
    PlotFigure1(vars, varNames); %load('TrioEEG_elife_Fig1_data.mat', 'vars', 'varNames')


%% Figure 1 - Figure Supplement 1
       
    % Figure 1 - Source Data 1.mat
    load('TrioEEGBehAnalysis.mat','vars');
    varNames = {'velres','srt'}; % just velres + srt
    vars.srt = exp(vars.srt); % log -> raw RT
    vars = structfun(@(x) nanmean(x,4), vars,'Uni',0); % take means over trials

    PlotFigure1_S1(vars, varNames); %load('TrioEEG_elife_Fig1_data.mat', 'vars', 'varNames')


%% Figure 2
    % distractor effects
    %% fig 2b,c,d

    % Figure 2 - Source Data 1.mat
    load('TrioEEGBehAnalysis.mat','vars');
    PlotFigure2(vars, {'depAngle'}); % this calls PlotFigure1 to make fig2b


%% Figure 2 - Figure Supplement 1
   
    % Figure 2 - Source Data 1.mat
    load('TrioEEGBehAnalysis.mat','vars'); % srt + depAngle
    vars.srt = exp(vars.srt); % log -> raw RT
    PlotFigure2_S1(vars);


%% Figure 3

    %% P3a and CNV figures
    % mean ERPs by Incentive*THP
    % difference waves: Incentive effect overall, and per THP/placebo
    % mean ERP within P3a and CNV window
    % regression coeffs of means against each behavioural variable
    
    % Figure 3 - Source Data 1.mat
    PlotFigure3('readyBS');


    %% same figures for incentive-cue
    % mean window-ampl for incentive-cue ERP
    % ERP waveform is discontiguous part of Figure3a
    % (difference waves are not in the paper)

    % Figure 3 - Source Data 1.mat
    PlotFigure3('audCue_200');



%% Figure 3 - Figure Supplement 1
    % pp-level mean ERPs (windowed): incentive*THP, prep,p3a,cnv

    % is done in PlotFigure3 function


    
    
%% Figure 3 - Figure Supplement 2
% time*channel difference waves: GND outputs (t-stat and p-values)
% incentive, THP, incentive*THP effects

    % Figure 3 - Figure Supplement 4 - Source Data 1.mat
    timelock = 'readyBS';
    load(sprintf('TrioEEGAnalyseByTrials_%s.mat',timelock), 'GNDEffects');

    PlotFigure3_S2(GNDEffects);   



%% Figure 4
% time*channel regressions coefficients of var~voltage, for each variable
% + p-values

    % Figure 4 - Source Data 1.mat
    load('WholeBrainRegressions_perm_Fig4.mat', 'allBetas','allPVals','xTimes2')
    PlotFigure4(allBetas,allPVals,xTimes2);

    


%% Figure 4 - Figure Supplement 1
% time*channel reg coeffs, when controlling for other 2 vars
% beta+p-values

    % Figure 4 - Figure Supplement 1 - Source Data 1.mat
    load('WholeBrainRegressions3_perm_Fig4_S1.mat', 'allBetas','allPVals','xTimes2')
    PlotFigure4(allBetas,allPVals,xTimes2);
    

%% Figure 4 - Figure Supplement 2
% time*channel reg coeffs, when controlling for microsaccades or ocular-drift
% beta+p-values

    measures = {'meanSpeed', 'nMicro'};

    % Figure 4 - Figure Supplement 2 - Source Data 1.mat
    load('WholeBrainRegressions_Eye_perm_Fig4_S2.mat','allBetas','allPVals','xTimes2','ind')
    allBetas = sq(allBetas(:,:,ind,:,:)); % only the v term
    allPVals= sq(allPVals(:,:,ind,:,:)); % these are all 1, no perm testing done

    for iM = 1:2
        PlotFigure4(allBetas(:,:,:,iM), allPVals(:,:,:,iM), xTimes2);
    end



%% Figure 5 
    % mediation analyses stats, velres+srt, cnv
    % 4 tests - b, t, df, p
    % + medEff
    
    % plus medMod, srt:CNV
    
    % load up the files and just print the tables
    timelock = 'readyBS';
    
    beh = 'velres'; % velres ~ CNV
    % Figure 5 - Source Data 1.mat
    v = load(sprintf('TrioEEGMediations_%s_%s.mat',timelock, beh), 'mediationTable');
    fprintf('velres ~ Incentive is not mediated by CNV\n');
    disp(v.mediationTable);
    
    beh = 'srt'; % srt ~ CNV
    % Figure 5 - Source Data 2.mat
    s = load(sprintf('TrioEEGMediations_%s_%s.mat',timelock, beh), 'mediationTable','medModTable');
    
    fprintf('srt ~ Incentive is mediated by CNV\n');
    disp(s.mediationTable);
    
    fprintf('srt ~ THP*Incentive moderation is mediated by CNV\n');
    disp(s.medModTable);
    
    
%% end




