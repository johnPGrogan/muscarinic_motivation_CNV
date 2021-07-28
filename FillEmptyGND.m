function GND = FillEmptyGND(data, chanInds, EEG, binNames, nGoodTrials, nBadTrials, expDesc, fNames)
% function GND = FillEmptyGND(data, chanInds, EEG, condNames, expDesc, fNames)
% Fill up an empty GND structure with data and other fields (will call
% MakeEmptyGND subfunction)
% 
% Inputs:
%   data = [nPP nTimes nBins nChans] data matrix
%   chanInds = indices of channels to include
%   EEG = eeglab structure, containing chanlocs, times, srate
%   binNames = cell array of names of bins/conds
%   nGoodTrials = [nPP nBin] number of trials included
%   nBadTrials = [nPP nBin] number of trials excluded
%   expDesc = [optional] descriptor to add
%   fNames = [optional] cell array of individual file/subject names
% 
% Outputs:
%   GND = GND structure 

%%

if ~exist('expDesc','var') || isempty(expDesc)
    expDesc = '';
end
if ~exist('fNames','var') || isempty(fNames)
    fNames = [];
end

[nPP, ~, nBins, ~] = size(data);


%% get an empty GND structure

GND = MakeEmptyGND(); % get an empty struct

%%
GND.exp_desc = expDesc;

GND.grands = permute(nanmean(data(:,:,:,chanInds),1),[4,2,3,1]); % [chans time bins] - mean
GND.grands_stder = permute(nanstd(data(:,:,:,chanInds)),[4,2,3,1]); % st err
GND.grands_t = GND.grands ./ GND.grands_stder; % t val = SNR

GND.sub_ct = repmat(size(data,1),1,nBins); % number of subjects per bin
GND.chanlocs = EEG.chanlocs(chanInds); % channel info

for i = 1:nBins 
    GND.bin_info(i).bindesc = binNames{i};% bin name
    GND.bin_info(i).condcode = 1;
end

GND.time_pts = EEG.times; % time points
GND.bsln_wind = NaN; % no baselining within this
GND.srate = EEG.srate;

GND.indiv_fnames = fNames; % pp names
GND.indiv_subnames = GND.indiv_fnames;

GND.indiv_bin_ct = nGoodTrials; % num trials included per [pp bin]
GND.indiv_bin_raw_ct = nGoodTrials + nBadTrials; % num trials included + excluded

GND.indiv_erps = permute(data(:,:,:,chanInds),[4,2,3,1]); %[chan, time, bin, pp]

GND.indiv_art_ics = repmat({NaN},1,nPP);

% diff stuff

end

function GND = MakeEmptyGND()
% set up an empty GND structure to match the field structure generated from
% erplab2GND

%% set names to fill

nameVals = {'exp_desc', 'An Experiment';
            'filename', [];
            'filepath', [];
            'saved',    'no';
            'grands',   [];
            'grands_stder', [];
            'grands_t', [];
            'sub_ct', [];
            'chanlocs', struct();
            'bin_info', struct('bindesc','','condcode',[]);   
            'condesc', {'Experiment (not cal pulses)'};
            'time_pts', [];
            'bsln_wind', [];   
            'odelay', [];
            'srate', [];     
            'indiv_fnames', {};
            'indiv_subnames', {};
            'indiv_traits', [];
            'indiv_bin_ct', [];
            'indiv_bin_raw_ct', [];
            'indiv_erps', {};
            'indiv_art_ics', {}; 
            'cals', [];
            'history',[];       
            't_tests',[];
            };
        
        
%% fill them

GND = cell2struct(nameVals(:,2), nameVals(:,1),1);
        
end
