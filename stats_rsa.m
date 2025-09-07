%% init
addpath('~/CoSMoMVPA/mvpa/')
addpath('~/Repository/CommonFunctions/matplotlib/')
stats = {};

%% load adult data
fprintf('loading adult data\n');
adult_data=load('adultRDMs.mat');
stats.adult_timevec = adult_data.timevect;

%% load infant data
fprintf('loading infant data\n');
fns = dir('../derivatives/results/sub-*_rdm.mat');
T = readtable('../participants.tsv','FileType','text');
infant_data={};cc=clock();mm='';
for f=1:numel(fns)
    fn = fullfile(fns(f).folder,fns(f).name);
    x=load(fn);
    infant_data.grouprdm_5hz(f,:,:,:) = x.rdm;
    stats.infant_timevec = x.res_rdm.a.fdim.values{1};
    stats.sub_id{f} = strrep(fns(f).name,'_rdm.mat','');
    stats.age(f) = T.age_months(find(contains(T.participant_id,stats.sub_id(f))));
    stats.demographics(f,:) = T(find(contains(T.participant_id,stats.sub_id(f))),:);
    mm=cosmo_show_progress(cc,f/numel(fns),sprintf('%i/%i loading ../results/%s\n',f,numel(fns),fns(f).name),mm);
end
fprintf('finished\n');

%% split by age group
stats.infant_group_young = stats.age<=6;
stats.infant_group_old = stats.age>6;

%% time time correlations and permutation test
loweridx = find(tril(ones(200),-1));
nboot = 1000;
clusterformingthreshold = .05;
clustermeasure = 'maxsum';
fprintf('time time correlation young\n');
X = squeeze(nanmean(infant_data.grouprdm_5hz(stats.infant_group_young,:,:,:)));
Y = squeeze(mean(adult_data.grouprdm_5hz(:,:,loweridx)))';
stats.time_time_young = cluster_permutation_test(X,Y,nboot,clusterformingthreshold,clustermeasure);
fprintf('time time correlation old\n');
X = squeeze(nanmean(infant_data.grouprdm_5hz(stats.infant_group_old,:,:,:)));
Y = squeeze(mean(adult_data.grouprdm_5hz(:,:,loweridx)))';
stats.time_time_old = cluster_permutation_test(X,Y,nboot,clusterformingthreshold,clustermeasure);

%% correlations with models and permutation test
loweridx = find(tril(ones(200),-1));
nboot = 1000;
clusterformingthreshold = .05;
clustermeasure = 'max';
stats.models = [];
stats.models(:,1) = pdist(ceil((1:200)/100)','jaccard'); %animacy
stats.models(:,2) = pdist(ceil((1:200)/20)','jaccard'); %category
stats.models(:,3) = pdist(ceil((1:200)/4)','jaccard'); %object
stats.modelnames = {'animacy','category','object'};
fprintf('model correlation young\n');
X = squeeze(nanmean(infant_data.grouprdm_5hz(stats.infant_group_young,:,:,:)));
for m=1:3
    stats.(['infant_young_model_',stats.modelnames{m}]) = cluster_permutation_test(X,stats.models(:,m),nboot,clusterformingthreshold,clustermeasure);
end
fprintf('model correlation old\n');
X = squeeze(nanmean(infant_data.grouprdm_5hz(stats.infant_group_old,:,:,:)));
for m=1:3
    stats.(['infant_old_model_',stats.modelnames{m}]) = cluster_permutation_test(X,stats.models(:,m),nboot,clusterformingthreshold,clustermeasure);
end
fprintf('model correlation adult\n');
X = squeeze(nanmean(adult_data.grouprdm_5hz(:,:,:,:)));
for m=1:3
    stats.(['adult_model_',stats.modelnames{m}]) = cluster_permutation_test(X,stats.models(:,m),nboot,clusterformingthreshold,clustermeasure);
end

%%
save('../derivatives/results/rsa_stats.mat','stats','-v7.3')

