function fixResult = TrioEEGBehGetFixations(result, useParFor)
% function result = TrioEEGBehGetFixations(result, useParFor)
% Loops through array of result structures, getting fixation data for each,
% saving a  file for each
% 
% Inputs:   result = result structure array (from TrioEEGBehProcess
%           useParFor: 1 = parfor loop (requires parallel programming
%                       toolbox), 0 = normal for loop
% Outputs: result = structure array with data for each session
% 

doFix=1;
useParFor=1;
%%

if exist('useParFor','var') && useParFor

    parfor i = 1:length(result)
%         disp(fileNames{i});
%         result(i) = ProcessOneFile(fileNames{i}, behDataFolder); % process targ-accept period
        
        disp(i);
            % also process previous periods, incl viscue, audcue, readycue, and
            % targ period with microsaccades included
        fixResult(i) = GetFixPeriod(result(i).result, result(i).d, result(i).dpp, result(i).normTraj);
    end
    
else
        
    for i = 1:length(result)
%         disp(fileNames{i});
%         result(i) = ProcessOneFile(fileNames{i}, behDataFolder);
        disp(i);
        fixResult(i) = GetFixPeriod(result(i).result, result(i).d, result(i).dpp, result(i).normTraj);
    end

end

%%
% save('TrioEEGBehProcessFiles.mat','-v7.3','result') % save after processing

if doFix
    save('TrioEEGBehProcessFilesFix.mat','-v7.3','fixResult');
end

end



function output = GetFixPeriod(result, d, dpp, normTraj)

%% viscue period

[trajVis, infoVis]=snipSaccades(result.s,'startviscue_t','startaudcue_t',  ...
    'centre', result.params.screenSize/2 - [0 100], ...
    'rotate', -2*pi* (d.previousTargetPos-1)/3, ...
    'flip', d.targetDir==-1,  'flipangle', 5*-pi/6, ...
    'minsize',0,... % include microsaccades
    'verbose',1, 'browse',0, 'traces',1);


infoVis= ConvertInfo(infoVis, dpp); % convert pixels into degrees. velocity into deg/sec

[trajVis,normTrajVis] = invertAndNormTraj(trajVis,dpp); % convert to degrees, flip Y, start at 0

%% audcue period

[trajAud, infoAud]=snipSaccades(result.s,'startaudcue_t','startreadycue_t',  ...
    'centre', result.params.screenSize/2 - [0 100], ...
    'rotate', -2*pi* (d.previousTargetPos-1)/3, ...
    'flip', d.targetDir==-1,  'flipangle', 5*-pi/6, ...
    'minsize',0,... % include microsaccades
    'verbose',1, 'browse',0, 'traces',1);


infoAud = ConvertInfo(infoAud, dpp); % convert pixels into degrees. velocity into deg/sec
[trajAud, normTrajAud] = invertAndNormTraj(trajAud,dpp);

%% fixation period


[trajReady, infoReady]=snipSaccades(result.s,'startreadycue_t','starttarget_t',  ...
    'centre', result.params.screenSize/2 - [0 100], ...
    'rotate', -2*pi* (d.previousTargetPos-1)/3, ...
    'flip', d.targetDir==-1,  'flipangle', 5*-pi/6, ...
    'minsize',0,... % include microsaccades
    'verbose',1, 'browse',0, 'traces',1);


infoReady = ConvertInfo(infoReady, dpp); % convert pixels into degrees. velocity into deg/sec

% don't filter these saccades

% sort traj
[trajReady,normTrajReady] = invertAndNormTraj(trajReady,dpp);



%% count microsaccades

[~, infoMicro]=snipSaccades(result.s,'starttarget_t','saccadeaccepted_t',  ...
    'centre', result.params.screenSize/2 - [0 100], ...
    'rotate', -2*pi* (d.previousTargetPos-1)/3, ...
    'flip', d.targetDir==-1,  'flipangle', 5*-pi/6, ...
    'minsize',0,...
    'verbose',1, 'browse',0, 'traces',1);

infoMicro = ConvertInfo(infoMicro, dpp); % convert pixels into degrees. velocity into deg/sec

% don't filter these ones

% trim normTraj to only have samples before the first macro
normTrajMicro = normTraj;
for i = 1:size(normTraj,1)
    macroInd = find(infoMicro.sAmpl(i,:) > 1, 1, 'first');
    if macroInd % trim from first macro onwards
        n = infoMicro.sRT(i, macroInd);
    else % if no macro, nan entire trial
        n = 1;
    end
    normTrajMicro(i, n:size(normTraj,2) ) = NaN; % from first saccade or entire if none detected
end


%%

output = workspace2struct('result|^d|^i$|normTraj$'); % don't save these as already saved

end

function info = ConvertInfo(info, dpp)
    % convert pixels into degrees
    degNames = {'sEndpt','sAmpl','sVec','sCurvS','sSpd'};
    for i = 1:length(degNames)
        info.(degNames{i}) = info.(degNames{i}) * dpp;
    end

    % move speed into deg/sec
    info.sSpd = info.sSpd * 1000;
    
    

end

function [traj,normTraj] = invertAndNormTraj(traj,dpp)
% convert to degrees, invert Y, normalise to starting point

% need to invert Y
traj = complex(real(traj), -imag(traj))  * dpp;

% remove NaNs
traj(isnan(traj)) = complex(NaN,NaN);

% normalise starting point
normTraj = traj - traj(:,1);
end


