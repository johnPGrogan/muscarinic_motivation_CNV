%%%% TrihexyTrioAllAnalysis

% This contains details on how to replicate the analysis from the paper:
% Grogan et al., 2021. Muscarinic receptors mediate motivation via 
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

%% plot figures from final data
TrioEEG_PlotAllFigures; 
% uses these files from the OSF
% TrioEEGBehAnalysis.mat, 
% chanlocs.mat
% TrioEEGAnalyseByTrials_readyBS.mat
% trioEEGAnalyseByTrials_audCue200.mat
% WholeBrainRegressions_perm_Fig4.mat
% WholeBrainRegressions3_perm_Fig4_S1.mat
% WholeBrainRegressions_Eye_perm_Fig4_S2.mat
% TrioEEGMediations_readyBS_srt.mat
% TrioEEGMediations_readyBS_velres.mat


%% Pre-processing details

% The anonymous data have already been pre-processed, using the script
% TrioEEGPipelineCall.m.

%% Behavioural analysis

% The eye-tracking data was processed using TrioEEGBehProcess.m
% It can be analysed by calling:
TrioEEGBehAnalysis % uses TrioEEGBehProcessFiles.mat

% Correlations can be run via
TrioEEGBehCorrel % uses TrioEEGBehAnalysis.mat

% RT vs distractor pull conditional plots via:
SplitDepAngleByRT % uses TrioEEGBehAnalysis.mat

% Distractor pull ksdensities via:
TrioEEGBehHists % uses TrioEEGBehAnalysis.mat

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
% uses chanlocs.mat, TrioEEGAnalyseByTrials_audCue_200_output.mat, TrioEEGBehAnalysis.mat

% Analyse the preparation cue period
TrioEEGAnalyseByTrials('readyBS');
% uses chanlocs.mat, TrioEEGAnalyseByTrials_readyBS_output.mat, TrioEEGBehAnalysis.mat

%% Whole-brain regressions

% This will run the whole-brain regressions, regressing each time-point and
% channel against the three behavioural variables
RegressEEGVsBeh('readyBS', [-200 1500]); % uses trioEEGAnalyseByTrials_readyBS.mat, TrioEEGBehAnalysis.mat

% this was actually run on a cluster, to speed up, using this:
% WholeBrainRegressPermuteParallel

% and this will plot it, along with the Mass Univariate testing from above
PlotTopoAndEffect
% uses TrioEEGAnalyseByTrials_readyBS.mat, chanlocs.mat, and the files
% created by those regression scripts above

%% supplemental analyses

% whole-brain regressions while controlling for other 2 variables (on
% cluster)
% Fig 4 - supplement 1
WholeBrainRegressPermuteParallel2


% whole-brain regressions while controlling for eye-movements
% (microsaccades, drift, eye-pos). no permutation testing for speed
% Fig 4 - supplement 2
TrioEEGBehMicrosaccades2
% uses TrioEEGBehAnalysis.mat, TrioEEGBehMicrosaccades_ready2Stats,
% TrioEEGBehReady2MicroFixation.mat, trioEEGAnalyseByTrials_readyBS,

%% Mediation analysis

% This will run the mediation analysis and mediated moderation
% Change the variables 'beh' and 'eeg' to look at the other variables
MediatedModeration
% uses 'TrioEEGAnalyseByTrials_readyBS.mat' and 'TrioEEGBehAnalysis.mat'



