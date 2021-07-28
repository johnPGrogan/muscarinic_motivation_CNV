% TrioEEGPipelineCall.m
% 
% call trioEEGPipeline with the epochInfo and flagInfo saved here as a
% record.

clc; clear all; close all;


% %% ITI
% sectionsToRun = struct( 'convert',  0,... % convert ov2mat
%                         'load',     0,... % trioEEGLoad (incl trigger stuff)
%                         'merge',    0,... % whether to remerge split files
%                         'delTrigs', 0,... % delete trigs from deleted trials
%                         'saccTrig', 0,... % insert sacc/blink trigs
%                         'check',    0,... % check info, times, trials match
%                         'filter',   0,... % notch and bandpass
%                         'resample', 0,... % resample
%                         'ica',      0,... % run ICA
%                         'epoch',    1,... % epoch 
%                         'flag',     1,... % flag artefacts
%                         'flagBeh',  0,... % flag trials by beh data
%                         'flagTrig', 1,... % flag beh by triggers (works if diff trial numbers)
%                         'average',  1,... % ERP average by bins, save
%                         'plot',     0  ); % plot ERPs
%                         
% epochInfo.bdfFile           = './trio_EEG_iti_BDF.txt';
% epochInfo.timeWindow        = [-500 500];
% epochInfo.baselineMethod    = 'window'; % 'none','window','withinEpoch','separateEpoch'
% epochInfo.baselineWindow    = [-100 0]; % or 'none' for no baselining
% epochInfo.baselineBdfFile   = [];%'./trio_EEG_audCue_BDF.txt'; % or [] to use main epoch file
% epochInfo.baselineTriggers  = []; % [], or regexpi string
% epochInfo.timeWindow2       = [-500 500]; % [] or window smaller than timewindow
% epochInfo.suffix            = '_be_iti_bs100'; % to append to EEG.setname
% epochInfo.save              = 0; % whether to save
% epochInfo.fieldsToCheck     = {'ppID', 'dates', 'trialTypes'}; % for CheckEEGMatData
% 
%     
%     
% flagInfo.artWindow  = [-100 300]; % time window for flagging
% flagInfo.save       = 0; % whether to save output of flagging    
% % flag based on behaviours
% flagInfo.behArtefacts = {'paused', 1,'driftFlag1', [1 1 1], 'driftFlag2', [1 1 1]};
% % flag if these triggers occur within the trial limits
% flagInfo.flagTriggers = {'triggers', '24|27|30|45|72|75',... % fixdebug/dco/driftt, lostfix, pause start/end
%                         'trialLimTrigs', {'^6$', '^6$|66|75|78|84'},... % itiStart, startrest, endpause, break, endexpt
%                         'flag', [1 9]};
% 
% 
% 
% [processedFiles, me] = trioEEGPipeline(sectionsToRun, epochInfo, flagInfo);
% disp(me);
% 
% %% viscue
% sectionsToRun = struct( 'convert',  0,... % convert ov2mat
%                         'load',     0,... % trioEEGLoad (incl trigger stuff)
%                         'merge',    0,... % whether to remerge split files
%                         'delTrigs', 0,... % delete trigs from deleted trials
%                         'saccTrig', 0,... % insert sacc/blink trigs
%                         'check',    0,... % check info, times, trials match
%                         'filter',   0,... % notch and bandpass
%                         'resample', 0,... % resample
%                         'ica',      0,... % run ICA
%                         'epoch',    1,... % epoch 
%                         'flag',     1,... % flag artefacts
%                         'flagBeh',  0,... % flag trials by beh data
%                         'flagTrig', 1,... % flag beh by triggers (works if diff trial numbers)
%                         'average',  1,... % ERP average by bins, save
%                         'plot',     0  ); % plot ERPs
%                         
% epochInfo.bdfFile           = './trio_EEG_visCue_BDF.txt';
% epochInfo.timeWindow        = [-100 500];
% epochInfo.baselineMethod    = 'window'; % 'none','window','withinEpoch','separateEpoch'
% epochInfo.baselineWindow    = [-100 0]; % or 'none' for no baselining
% epochInfo.baselineBdfFile   = [];%'./trio_EEG_audCue_BDF.txt'; % or [] to use main epoch file
% epochInfo.baselineTriggers  = []; % [], or regexpi string
% epochInfo.timeWindow2       = [-100 500]; % [] or window smaller than timewindow
% epochInfo.suffix            = '_be_visCue_bs100'; % to append to EEG.setname
% epochInfo.save              = 0; % whether to save
% epochInfo.fieldsToCheck     = {'ppID', 'dates', 'trialTypes'}; % for CheckEEGMatData
% 
%     
%     
% flagInfo.artWindow  = [-100 500]; % time window for flagging
% flagInfo.save       = 0; % whether to save output of flagging    
% % flag based on behaviours
% flagInfo.behArtefacts = {'paused', 1,'driftFlag1', [1 1 1], 'driftFlag2', [1 1 1]};
% % flag if these triggers occur within the trial limits
% flagInfo.flagTriggers = {'triggers', '24|27|30|45|72|75',... % fixdebug/dco/driftt, lostfix, pause start/end
%                         'trialLimTrigs', {'^6$', '^6$|66|75|78|84'},... % itiStart, startrest, endpause, break, endexpt
%                         'flag', [1 9]};
% 
% 
% 
% [processedFiles, me] = trioEEGPipeline(sectionsToRun, epochInfo, flagInfo);
% disp(me);
% 
% %% pre-audcue (for baseline purposes) - not baselined itself
% 
% % change which beh stuff is checked - no RT
% 
% sectionsToRun = struct( 'convert',  0,... % convert ov2mat
%                         'load',     0,... % trioEEGLoad (incl trigger stuff)
%                         'merge',    0,... % whether to remerge split files
%                         'delTrigs', 0,... % delete trigs from deleted trials
%                         'saccTrig', 0,... % insert sacc/blink trigs
%                         'check',    0,... % check info, times, trials match
%                         'filter',   0,... % notch and bandpass
%                         'resample', 0,... % resample
%                         'ica',      0,... % run ICA
%                         'epoch',    1,... % epoch 
%                         'flag',     1,... % flag artefacts
%                         'flagBeh',  0,... % flag trials by beh data
%                         'flagTrig', 1,... % flag beh by triggers (works if diff trial numbers)
%                         'average',  1,... % ERP average by bins, save
%                         'plot',     0  ); % plot ERPs
%                         
% epochInfo.bdfFile           = './trio_EEG_audCue_BDF.txt';
% epochInfo.timeWindow        = [-200 0];
% epochInfo.baselineMethod    = 'none'; % 'none','window','withinEpoch','separateEpoch'
% epochInfo.baselineWindow    = 'none'; % or 'none' for no baselining
% epochInfo.baselineBdfFile   = [];%'./trio_EEG_audCue_BDF.txt'; % or [] to use main epoch file
% epochInfo.baselineTriggers  = []; % [], or regexpi string
% epochInfo.timeWindow2       = [-200 0]; % [] or window smaller than timewindow
% epochInfo.suffix            = '_be_audCue_none'; % to append to EEG.setname
% epochInfo.save              = 0; % whether to save
% epochInfo.fieldsToCheck     = {'ppID', 'dates', 'trialTypes'}; % for CheckEEGMatData
% 
%     
%     
% flagInfo.artWindow  = [-200 0]; % time window for flagging
% flagInfo.save       = 0; % whether to save output of flagging    
% % flag based on behaviours
% flagInfo.behArtefacts = {'paused', 1,'driftFlag1', [1 1 1], 'driftFlag2', [1 1 1]};
% % flag if these triggers occur within the trial limits
% flagInfo.flagTriggers = {'triggers', '24|27|30|45|72|75',... % fixdebug/dco/driftt, lostfix, pause start/end
%                         'trialLimTrigs', {'^6$', '^6$|66|75|78|84'},... % itiStart, startrest, endpause, break, endexpt
%                         'flag', [1 9]};
% 
% 
% 
% [processedFiles, me] = trioEEGPipeline(sectionsToRun, epochInfo, flagInfo);
% disp(me);

%% post-audcue (for rew resp e.g. P3b)

sectionsToRun = struct( 'convert',  0,... % convert ov2mat
                        'load',     0,... % trioEEGLoad (incl trigger stuff)
                        'merge',    0,... % whether to remerge split files
                        'delTrigs', 0,... % delete trigs from deleted trials
                        'saccTrig', 0,... % insert sacc/blink trigs
                        'check',    0,... % check info, times, trials match
                        'filter',   0,... % notch and bandpass
                        'resample', 0,... % resample
                        'ica',      0,... % run ICA
                        'epoch',    1,... % epoch 
                        'flag',     1,... % flag artefacts
                        'flagBeh',  0,... % flag trials by beh data
                        'flagTrig', 1,... % flag beh by triggers (works if diff trial numbers)
                        'average',  1,... % ERP average by bins, save
                        'plot',     0  ); % plot ERPs
                        
epochInfo.bdfFile           = './trio_EEG_audCue_BDF.txt';
epochInfo.timeWindow        = [-1600 1400];
epochInfo.baselineMethod    = 'window'; % 'none','window','withinEpoch','separateEpoch'
epochInfo.baselineWindow    = [-200 0]; % or 'none' for no baselining
epochInfo.baselineBdfFile   = [];%'./trio_EEG_audCue_BDF.txt'; % or [] to use main epoch file
epochInfo.baselineTriggers  = []; % [], or regexpi string
epochInfo.timeWindow2       = [-1600 1400]; % [] or window smaller than timewindow
epochInfo.suffix            = '_be_audCue_200'; % to append to EEG.setname
epochInfo.save              = 0; % whether to save
epochInfo.fieldsToCheck     = {'ppID', 'dates', 'trialTypes'}; % for CheckEEGMatData
    
    
flagInfo.artWindow  = [-200 550]; % time window for flagging
flagInfo.save       = 0; % whether to save output of flagging    
% flag based on behaviours
flagInfo.behArtefacts = {'paused', 1,'driftFlag1', [1 1 1], 'driftFlag2', [1 1 1]};
% flag if these triggers occur within the trial limits
flagInfo.flagTriggers = {'triggers', '24|27|30|45|72|75',... % fixdebug/dco/driftt, lostfix, pause start/end
                        'trialLimTrigs', {'^6$', '^6$|66|75|78|84'},... % itiStart, startrest, endpause, break, endexpt
                        'flag', [1 9]};



[processedFiles, me] = trioEEGPipeline(sectionsToRun, epochInfo, flagInfo);
disp(me);

%% post-readyCue

sectionsToRun = struct( 'convert',  0,... % convert ov2mat
                        'load',     0,... % trioEEGLoad (incl trigger stuff)
                        'merge',    0,... % whether to remerge split files
                        'delTrigs', 0,... % delete trigs from deleted trials
                        'saccTrig', 0,... % insert sacc/blink trigs
                        'check',    0,... % check info, times, trials match
                        'filter',   0,... % notch and bandpass
                        'resample', 0,... % resample
                        'ica',      0,... % run ICA
                        'epoch',    1,... % epoch 
                        'flag',     1,... % flag artefacts
                        'flagBeh',  0,... % flag trials by beh data
                        'flagTrig', 1,... % flag beh by triggers (works if diff trial numbers)
                        'average',  1,... % ERP average by bins, save
                        'plot',     0  ); % plot ERPs
                        
epochInfo.bdfFile           = './trio_EEG_readyCue_BDF.txt';
epochInfo.timeWindow        = [-200 1700];
epochInfo.baselineMethod    = 'none'; % 'none','window','withinEpoch','separateEpoch'
epochInfo.baselineWindow    = 'none'; % or 'none' for no baselining
epochInfo.baselineBdfFile   = [];%'./trio_EEG_audCue_BDF.txt'; % or [] to use main epoch file
epochInfo.baselineTriggers  = '36|39'; % [], or regexpi string
epochInfo.timeWindow2       = [-200 1700]; % [] or window smaller than timewindow
epochInfo.suffix            = '_be_readyCue_none'; % to append to EEG.setname
epochInfo.save              = 0; % whether to save
epochInfo.fieldsToCheck     = {'ppID', 'dates', 'trialTypes'}; % for CheckEEGMatData
    
    
flagInfo.artWindow  = [-200 1500]; % time window for flagging
flagInfo.save       = 0; % whether to save output of flagging    
% flag based on behaviours
flagInfo.behArtefacts = {'paused', 1,'driftFlag1', [1 1 1], 'driftFlag2', [1 1 1]};
% flag if these triggers occur within the trial limits
flagInfo.flagTriggers = {'triggers', '24|27|30|45|72|75',... % fixdebug/dco/driftt, lostfix, pause start/end
                        'trialLimTrigs', {'^6$', '^6$|66|75|78|84'},... % itiStart, startrest, endpause, break, endexpt
                        'flag', [1 9]};


[processedFiles, me] = trioEEGPipeline(sectionsToRun, epochInfo, flagInfo);
disp(me);



% %% pre-target
% 
% 
% 
% sectionsToRun = struct( 'convert',  0,... % convert ov2mat
%                         'load',     0,... % trioEEGLoad (incl trigger stuff)
%                         'merge',    0,... % whether to remerge split files
%                         'delTrigs', 0,... % delete trigs from deleted trials
%                         'saccTrig', 0,... % insert sacc/blink trigs
%                         'check',    0,... % check info, times, trials match
%                         'filter',   0,... % notch and bandpass
%                         'resample', 0,... % resample
%                         'ica',      0,... % run ICA
%                         'epoch',    1,... % epoch 
%                         'flag',     1,... % flag artefacts
%                         'flagBeh',  0,... % flag trials by beh data
%                         'flagTrig', 1,... % flag beh by triggers (works if diff trial numbers)
%                         'average',  1,... % ERP average by bins, save
%                         'plot',     0  ); % plot ERPs
%                         
% epochInfo.bdfFile           = './trio_EEG_targ_BDF.txt';
% epochInfo.timeWindow        = [-8000 1000];
% epochInfo.baselineMethod    = 'withinEpoch'; % 'none','window','withinEpoch','separateEpoch'
% epochInfo.baselineWindow    = [-200 0]; % or 'none' for no baselining
% epochInfo.baselineBdfFile   = [];%'./trio_EEG_audCue_BDF.txt'; % or [] to use main epoch file
% epochInfo.baselineTriggers  = '36|39'; % [], or regexpi string
% epochInfo.timeWindow2       = [-1600 1000]; % [] or window smaller than timewindow
% epochInfo.suffix            = '_be_targ_bs_audCue200'; % to append to EEG.setname
% epochInfo.save              = 0; % whether to save
% epochInfo.fieldsToCheck     = {'ppID', 'dates', 'trialTypes', 'RT','triggerRT'}; % for CheckEEGMatData
%     
%     
% flagInfo.artWindow  = [-1500 0]; % time window for flagging
% flagInfo.save       = 0; % whether to save output of flagging    
% % flag based on behaviours
% flagInfo.behArtefacts = {'paused', 1,'driftFlag1', [1 1 1], 'driftFlag2', [1 1 1],...
%                          'saccRT', 10, 'behRT', [100 2000]};
% % flag if these triggers occur within the trial limits
% flagInfo.flagTriggers = {'triggers', '24|27|30|45|72|75',... % fixdebug/dco/driftt, lostfix, pause start/end
%                         'trialLimTrigs', {'^6$', '^6$|66|75|78|84'},... % itiStart, startrest, endpause, break, endexpt
%                         'flag', [1 9]};
% 
% 
% 
% [processedFiles, me] = trioEEGPipeline(sectionsToRun, epochInfo, flagInfo);
% disp(me);
% 
% %% pre-saccade
% 
% 
% 
% sectionsToRun = struct( 'convert',  0,... % convert ov2mat
%                         'load',     0,... % trioEEGLoad (incl trigger stuff)
%                         'merge',    0,... % whether to remerge split files
%                         'delTrigs', 0,... % delete trigs from deleted trials
%                         'saccTrig', 0,... % insert sacc/blink trigs
%                         'check',    0,... % check info, times, trials match
%                         'filter',   0,... % notch and bandpass
%                         'resample', 0,... % resample
%                         'ica',      0,... % run ICA
%                         'epoch',    1,... % epoch 
%                         'flag',     1,... % flag artefacts
%                         'flagBeh',  0,... % flag trials by beh data
%                         'flagTrig', 1,... % flag beh by triggers (works if diff trial numbers)
%                         'average',  1,... % ERP average by bins, save
%                         'plot',     0  ); % plot ERPs
%                         
% epochInfo.bdfFile           = './trio_EEG_firstSacc_BDF.txt';
% epochInfo.timeWindow        = [-6000 500];
% epochInfo.baselineMethod    = 'withinEpoch'; % 'none','window','withinEpoch','separateEpoch'
% epochInfo.baselineWindow    = [-200 0]; % or 'none' for no baselining
% epochInfo.baselineBdfFile   = [];%'./trio_EEG_audCue_BDF.txt'; % or [] to use main epoch file
% epochInfo.baselineTriggers  = '36|39'; % [], or regexpi string
% epochInfo.timeWindow2       = [-1600 500]; % [] or window smaller than timewindow
% epochInfo.suffix            = '_be_firstSacc_bs_audCue200'; % to append to EEG.setname
% epochInfo.save              = 0; % whether to save
% epochInfo.fieldsToCheck     = {'ppID', 'dates', 'trialTypes', 'RT','triggerRT'}; % for CheckEEGMatData
%     
% flagInfo.artWindow  = [-1500 -200]; % time window for flagging muscle/eye artefacts
% flagInfo.save       = 0; % whether to save output of flagging    
% % flag based on behaviours
% flagInfo.behArtefacts = {'paused', 1,'driftFlag1', [1 1 1], 'driftFlag2', [1 1 1],...
%                          'saccRT', 10, 'behRT', [100 2000]};
% % flag if these triggers occur within the trial limits
% flagInfo.flagTriggers = {'triggers', '24|27|30|45|72|75',... % fixdebug/dco/driftt, lostfix, pause start/end
%                         'trialLimTrigs', {'^6$', '^6$|66|75|78|84'},... % itiStart, startrest, endpause, break, endexpt
%                         'flag', [1 9]};
% 
% 
% [processedFiles, me] = trioEEGPipeline(sectionsToRun, epochInfo, flagInfo);
% disp(me);
% 
% 
% %% feedback
% 
% sectionsToRun = struct( 'convert',  0,... % convert ov2mat
%                         'load',     0,... % trioEEGLoad (incl trigger stuff)
%                         'merge',    0,... % whether to remerge split files
%                         'delTrigs', 0,... % delete trigs from deleted trials
%                         'saccTrig', 0,... % insert sacc/blink trigs
%                         'check',    0,... % check info, times, trials match
%                         'filter',   0,... % notch and bandpass
%                         'resample', 0,... % resample
%                         'ica',      0,... % run ICA
%                         'epoch',    1,... % epoch 
%                         'flag',     1,... % flag artefacts
%                         'flagBeh',  0,... % flag trials by beh data
%                         'flagTrig', 1,... % flag beh by triggers (works if diff trial numbers)
%                         'average',  1,... % ERP average by bins, save
%                         'plot',     0  ); % plot ERPs
%                         
% epochInfo.bdfFile           = './trio_EEG_fb_BDF.txt';
% epochInfo.timeWindow        = [-400 1000];
% epochInfo.baselineMethod    = 'window'; % 'none','window','withinEpoch','separateEpoch'
% epochInfo.baselineWindow    = [-100 0]; % or 'none' for no baselining
% epochInfo.baselineBdfFile   = [];%'./trio_EEG_audCue_BDF.txt'; % or [] to use main epoch file
% epochInfo.baselineTriggers  = '36|39'; % [], or regexpi string
% epochInfo.timeWindow2       = [-400 1000]; % [] or window smaller than timewindow
% epochInfo.suffix            = '_be_fb_bs_100'; % to append to EEG.setname
% epochInfo.save              = 0; % whether to save
% epochInfo.fieldsToCheck     = {'ppID', 'dates', 'trialTypes', 'RT'}; % for CheckEEGMatData
%     
% flagInfo.artWindow  = [0 500]; % time window for flagging muscle/eye artefacts
% flagInfo.save       = 0; % whether to save output of flagging    
% % flag based on behaviours
% flagInfo.behArtefacts = {'paused', 1,'driftFlag1', [1 1 1], 'driftFlag2', [1 1 1],...
%                          'saccRT', 10, 'behRT', [100 2000]};
% % flag if these triggers occur within the trial limits
% flagInfo.flagTriggers = {'triggers', '24|27|30|45|72|75',... % fixdebug/dco/driftt, lostfix, pause start/end
%                         'trialLimTrigs', {'^6$', '^6$|66|75|78|84'},... % itiStart, startrest, endpause, break, endexpt
%                         'flag', [1 9]};
% 
% 
% [processedFiles, me] = trioEEGPipeline(sectionsToRun, epochInfo, flagInfo);
% disp(me);
% 
% 
