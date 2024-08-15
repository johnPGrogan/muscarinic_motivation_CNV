function pval = FindClustersLikeGND(trueTVals, permTVals, chan_hood, tail, df)
% Function to find significant clusters after having already run a set of
% permutation tests, using the resulting null t-value distribution
% The code below is copied from the Mass Univariate Toolbox (DMGroppe)
% functions clustGND.m & clust_perm1, and requires that toolbox be on the path.
%
% Inputs:
%   trueTVals = the t-values for the 'true' effect, i.e. without any
%       permutation, in [nTimes, nChans] size
%   permTVals = the t-values from the permutation tests [nTimes, nChans, nPerms]
%   chan_hood = [nCh * nCh] symmetric binary matrix indicating which channels are
%               neighbors. If chan_hood(a,b)=1, then Channel A and Channel
%               B are nieghbors: spatial_neighbors(chanLocs, 0.61, [])
%   tail = [-1, 0, or 1], the tail of the distribution/test to use. -1 is
%       the lower tail (alt hypothesis that the effect is below the null),
%       +1 is the upper tail (e.g. effect > 0), and 0 is two-tailed (i.e.
%       that there is a difference)
%   df = degrees of freedom to use
%
%
% Outputs:
%   pval = [chans times] matrix of clustered permutation pvals
% 
% John Grogan, 2022.


%% get info

[nTimes, nChans, nPerms] = size(permTVals);

%% get channel neighbours

% chan_dist = 0.61; % default from the toolboxes,
% head_radius = [];
% 
% chan_hood = spatial_neighbors(chanLocs(1:nChans), chan_dist, head_radius);

%% set desired t_thresh

% tail = 0;

thresh_p = .05; % desired alpha
fwer = .05;

% df = 29;%18576; % df in fitglme (small differences don't matter)

if tail
    %one tailed test
    thresh_t=tinv(thresh_p,df); %note unless thresh_p is greater than .5, thresh_t will be negative
else
    %two tailed test
    thresh_t=tinv(thresh_p/2,df);
end


%% on each perm, get cluster

% nPerms = 100;
mn_clust_mass = zeros(1,nPerms);

parfor iPerm = 1:nPerms

    t = permTVals(:,:,iPerm)'; %[chans times] % get the permuted t-vals

    [clust_ids, n_clust]=find_clusters(t, thresh_t, chan_hood, -1); % find the clusters that pass thresh


    %get most extremely negative t-score (sign isn't important since we asumme
    %symmetric distribution of null hypothesis for one sample test)
    mn_clust_mass(iPerm) = find_mn_mass(clust_ids, t, n_clust);

end

% Estimate true FWER of test
if tail==0
    %two-tailed
    tmx_ptile=prctile(mn_clust_mass, 100*fwer/2);
    est_alpha=mean(mn_clust_mass<=tmx_ptile)*2;
else
    %one tailed
    tmx_ptile=prctile(mn_clust_mass,100*fwer);
    est_alpha=mean(mn_clust_mass<=tmx_ptile);
end
fprintf('Desired family-wise error rate: %f\n',fwer);
fprintf('Estimated actual family-wise error rate: %f\n',est_alpha);

%% get pvals

trueTVals = trueTVals'; % make into [chans times] for this bit

pval=ones(nChans,nTimes); % will rotate back afterwards

if tail==0 % pos + neg clusters
    %positive clusters
    [clust_ids, n_clust]=find_clusters(trueTVals,-thresh_t,chan_hood,1); %note thresh_t is negative by default

    clust_info.pos_clust_pval=ones(1,n_clust);
    clust_info.pos_clust_mass=zeros(1,n_clust);
    clust_info.pos_clust_ids=clust_ids;
    for a=1:n_clust
        use_ids=find(clust_ids==a);
        clust_mass=sum(trueTVals(use_ids));
        clust_p=mean(mn_clust_mass<=(-clust_mass))*2; %multiply by 2 since we're effectively doing Bonferroni correcting for doing two tests (an upper tail and lower tail)
        pval(use_ids)=clust_p;
        clust_info.pos_clust_pval(a)=clust_p;
        clust_info.pos_clust_mass(a)=clust_mass;
    end

    %negative clusters
    [clust_ids, n_clust]=find_clusters(trueTVals,thresh_t,chan_hood,-1); %note thresh_t is negative by default
    clust_info.neg_clust_pval=ones(1,n_clust);
    clust_info.neg_clust_mass=zeros(1,n_clust);
    clust_info.neg_clust_ids=clust_ids;
    for a=1:n_clust
        use_ids=find(clust_ids==a);
        clust_mass=sum(trueTVals(use_ids));
        clust_p=mean(mn_clust_mass<=clust_mass)*2; %multiply by 2 since we're effectively doing Bonferroni correcting for doing two tests (an upper tail and lower tail)
        pval(use_ids)=clust_p;
        clust_info.neg_clust_pval(a)=clust_p;
        clust_info.neg_clust_mass(a)=clust_mass;
    end
elseif tail==1
    %upper tailed
    [clust_ids, n_clust]=find_clusters(trueTVals,-thresh_t,chan_hood,1); %note thresh_t is negative by default
    clust_info.pos_clust_pval=ones(1,n_clust);
    clust_info.pos_clust_mass=zeros(1,n_clust);
    clust_info.pos_clust_ids=clust_ids;
    for a=1:n_clust
        use_ids=find(clust_ids==a);
        clust_mass=sum(trueTVals(use_ids));
        clust_p=mean(mn_clust_mass<=(-clust_mass));
        pval(use_ids)=clust_p;
        clust_info.pos_clust_pval(a)=clust_p;
        clust_info.pos_clust_mass(a)=clust_mass;
    end
else
    %lower tailed
    [clust_ids, n_clust]=find_clusters(trueTVals,thresh_t,chan_hood,-1); %note thresh_t is negative by default
    clust_info.neg_clust_pval=ones(1,n_clust);
    clust_info.neg_clust_mass=zeros(1,n_clust);
    clust_info.neg_clust_ids=clust_ids;
    for a=1:n_clust
        use_ids=find(clust_ids==a);
        clust_mass=sum(trueTVals(use_ids));
        clust_p=mean(mn_clust_mass<=clust_mass);
        pval(use_ids)=clust_p;
        clust_info.neg_clust_pval(a)=clust_p;
        clust_info.neg_clust_mass(a)=clust_mass;
    end
end


pval = pval';% rotate back into [nTimes nChans]


end



%% subfunc

function mn_clust_mass=find_mn_mass(clust_ids,data_t,n_clust)
% copied from clust_perm1.m in Mass_Univariate_Toolbox

mn_clust_mass=0;

%looking for most negative cluster mass
for z=1:n_clust,
    use_ids=(clust_ids==z);
    use_mass=sum(data_t(use_ids));
    if use_mass<mn_clust_mass,
        mn_clust_mass=use_mass;
    end
end

end