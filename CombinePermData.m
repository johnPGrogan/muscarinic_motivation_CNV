function n = CombinePermData(saveFolder, varNames)
% load up each perm file, combine 3rd dim, save per var

if ~exist('varNames','var') || isempty(varNames)
    varNames = {'velres', 'srt','depAngle'};
end
nVars = length(varNames);

if ~exist('saveFolder','var') || isempty(saveFolder)
    saveFolder = './Data/';%'C:\Users\groganj1\OneDrive - Nexus365\AnalysisOD\Trio_EEG';
end
    
for iV = 1:nVars
    f = what(saveFolder);
    files = f.mat(cellRegexpi(f.mat, sprintf('WholeBrainRegressPermuteParallel3_%s',varNames{iV}))>0);
    if isempty(files)
        continue; 
    end
    % this should now just cat them, across nperm sets
    d.tValsPerm = [];
    for i = 1:length(files)
        data = load(fullfile(saveFolder, files{i}));
        isEmpty = all(data.tValsPerm==0,[1 2]);
        d.tValsPerm = cat(3, d.tValsPerm, data.tValsPerm(:,:,~isEmpty));
    end

    %% check for duplicates - due to not resetting random seeds

    nUniq = length(unique(d.tValsPerm(1,1,:)));
    if nUniq < size(d.tValsPerm,3)
        error('duplications in permutations detected');
    end

    fprintf('%s has %d permutations\n', varNames{iV}, size(d.tValsPerm,3))
    n(iV) = size(d.tValsPerm,3);
    if n(iV) > 2500
        % trim
        d.tValsPerm = d.tValsPerm(:,:,1:2500);
        n(iV) = size(d.tValsPerm,3);
    end
    save(sprintf('%s_combined_%s.mat', fullfile(saveFolder, 'WholeBrainRegressPermuteParallel3'), varNames{iV}),'-struct','d');

    
end

