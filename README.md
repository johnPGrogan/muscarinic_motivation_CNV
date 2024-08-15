%%%% TrihexyTrioAllAnalysis

% This contains details on how to replicate the analysis from the paper:
% Grogan et al., 2024. Muscarinic receptors mediate motivation via 
% preparatory neural activity in humans.

% Anonymous data are available here (https://osf.io/zuq5c/)
% The paths in the scripts given here will need to be updated to point to
% where this data is saved.

% Requirements:
% Matlab r2018b
% EEGlab v2019.1
% ERPlab v8.0.0
% DMGroppe Mass Univariate ERP Toolbox
% Matlib (www.github.com/sgmanohar/matlib)
% Othercolor & Crameri for colormaps
% 


%% Task code
% The task was run using matlab r2018b and PsychToolbox v3.
% The code is at https://osf.io/dqzmf/


%% Pre-processing details

% The anonymous data have already been pre-processed, using the script
% TrioEEGPipelineCall.m.

%% Behavioural analysis

% The eye-tracking data was processed using TrioEEGBehProcess.m
% It can be analysed by calling:
TrioEEGBehAnalysis

% Correlations can be run via
TrioEEGBehCorrel

% RT vs distractor pull conditional plots via:
SplitDepAngleByRT

% Distractor pull ksdensities via:
TrioEEGBehHists

%% EEG analysis

% TrioEEGAnalyseByTrials.m takes the name of the epoch, and runs main
% analyses. The preparation epoch is "readyBS" and the incentive-cue epoch
% is "audCue_200". The saved files are large (>1.5GB). This will plot
% traces of the voltages for some channels, and analyse the mean voltages
% for Cz. It will also do Mass Univariate Testing (cluster-based) for
% effects and interactions, and regress the mean voltages against
% behaviour.

% Analyse the incentive cue period
TrioEEGAnalyseByTrials('audCue_200');

% Analyse the preparation cue period
TrioEEGAnalyseByTrials('readyBS');

%% Whole-brain regressions

% This will run the whole-brain regressions, regressing each time-point and
% channel against the three behavioural variables
RegressEEGVsBeh('readyBS', [-200 1500]);

% this was actually run on a cluster, to speed up, using this:
% WholeBrainRegressPermuteParallel

% and this will plot it, along with the Mass Univariate testing from above
PlotTopoAndEffect

%% Mediation analysis

% This will run the mediation analysis and mediated moderation
% Change the variables 'beh' and 'eeg' to look at the other variables
MediatedModeration

%% supplemental analyses

% whole-brain regressions while controlling for other 2 variables (on
% cluster)
% Fig 4 - supplement 1
WholeBrainRegressPermuteParallel2


% whole-brain regressions while controlling for eye-movements
% (microsaccades, drift, eye-pos). no permutation testing for speed
% Fig 4 - supplement 2
TrioEEGBehMicrosaccades2


