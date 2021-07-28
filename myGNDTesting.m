function GND = myGNDTesting(data, chanInds, EEG, condNames, diffBins, nGoodTrials, nBadTrials, nPerms, expDesc, fNames, doUnivars)
% function GND = myGNDTesting(erpFiles, data, chanInds, EEG, condNames, diffBins, doPerm, expDesc, fNames)
% Put data into an empty GND structure, then run a set of tests on
% difference waves
% Inputs:
%   data = [nPP nTimes nBins nChans] matrix of EEG data (means over trials)
%   chanInds = indices of channels to keep
%   EEG = eeglab data structure with chanlocs, srate, times 
%   condNames = cell array of names of conds/bins
%   diffBins = [N 2] matrix of bin indices to make difference waves from,
%       it will do a separate comparison for each row
%   nGoodTrials = [nPP nBin] number of trials included
%   nBadTrials = [nPP nBin] number of trials excluded
%   nPerms = scalar or vector (to match size(diffBins,1)) of number of
%       permutations to do for each comparison. 0 = don't do them
%   expDesc = [optional] experiment descriptor (e.g. file suffix)
%   fNames = [optional] list of subject/file names to put into GND
%   doUnivars = scalar or vector (to match size(diffBins,1)) of boolean for
%       whether to do univariate tests. default: 0 = don't do them
% 
% Outputs:
%   GND = GND structure containing results of tests

%%
if ~exist('expDesc','var') || isempty(expDesc)
    expDesc = '';
end
if ~exist('fNames','var') || isempty(fNames)
    fNames = [];
end
if numel(nPerms)==1
    nPerms = repmat(nPerms, size(diffBins,1),1);
end
if ~exist('doUnivars','var') || isempty(doUnivars)
    doUnivars = zeros(size(nPerms));
elseif numel(doUnivars) == 1
    doUnivars = repmat(doUnivars, size(diffBins,1),1);
end
%% mass univariate toolbox for permutation tests - pool A and B

% get files to reformat
% erpFiles = erpFiles(~cellfun(@isempty, erpFiles)); % remove empties

% % this will average across all pps/sessions
% % convert .erps into GND 
% GND = erplab2GND(fullfile(analysisFolder,erpFiles),...
%         'bsln', NaN,...                                         % no baselining
%         'exclude_chans', {'HEOG','VEOG','elX','elY'},...        % don't keep eye stuff
%         'out_fname', 'no save');                                % don't save


% fill it up
GND = FillEmptyGND(data, chanInds, EEG, condNames, nGoodTrials, nBadTrials, expDesc, fNames);

% replace mean(mean) with mean(SD)
% GND = ReplaceGNDWithStd(GND, eegStd);


% INVERT  IPSI/CONTRA?


timeWindow = [min(GND.time_pts) max(GND.time_pts)]; % times to do permutation across

nBinsOrig = size(GND.grands,3);
for i = 1:size(diffBins,1)
    % create difference wave: bin1 - bin2
    GND = bin_dif(GND, diffBins(i,1), diffBins(i,2), [GND.bin_info(diffBins(i,1)).bindesc '--' GND.bin_info(diffBins(i,2)).bindesc]);
    % gui_erp(GND, 'bin', 3); % plot

    
    binInd = nBinsOrig + i;

    % FDR controlled univariate t-tests
    if doUnivars(i)
        GND = tfdrGND(GND, binInd, 'method', 'by', 'time_wind', timeWindow, 'save_GND', 'no','plot_raster','no');
    end

    drawnow;
    
    if nPerms(i)
        % cluster permutation test
        GND2 = clustGND(GND, binInd, 'time_wind', timeWindow, 'save_GND', 'no','plot_raster','no', 'n_perm',nPerms(i));
        GND.perm_t_tests(i) = GND2.t_tests(end); % stick in in a different field so aligns with bins
    end
    
end
