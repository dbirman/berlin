function run_berlin( cfolder )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
load(fullfile('~/data/berlin_experiment',sprintf('%s_data.mat',cfolder)));

%% Parse stimvols
% csv dir1 dir2 mask correct match

udir1 = unique(data.design(:,2));
umask = unique(data.design(:,4));
for di = 1:length(udir1)
    for mi = 1:length(umask)
        svs = data.design(:,1);
        sv{mi,di} = svs(logical((data.design(:,2)==udir1(di)).*(data.design(:,4)==umask(mi))));
    end
end
%% Deconvolve
% just for fun?
concatInfo.runTransition = data.runtrans;
curd = constructD(mean(data.rtSeries{2}),sv(4,1:2),0.5,18,concatInfo,'none','deconv',0);
decon = getr2timecourse(curd.timecourse,curd.nhdr,curd.hdrlenTR,curd.scm,curd.framePeriod,curd.verbose);
decon = rmfield(decon,'scm');
decon = rmfield(decon,'covar');
plot(decon.ehdr');

%% Get instances
inst = cell(size(sv));
for di = 1:size(sv,1)
    for mi = 1:size(sv,2)
        svs = sv{di,mi};
        for si = 1:length(svs)
            inst{di,mi}(end+1,:) = mean(ts(:,svs(si)+8:svs(si)+12)');
        end
    end
end

%% PCA?
ts = ts
%% Compute permutation test (balanced)
group1 = [inst{4,1}];
group2 = [inst{4,2}];
insts = [group1;group2];
groups = [ones(size(group1,1),1); 2*ones(size(group2,1),1)];
insts = zscore(insts);
out = classify(insts,insts,groups);
% overfit score
overfit = sum(out==groups)/length(groups);
% permute
for i = 1:1000
    pgroups = groups(randperm(length(groups)));
    out = classify(insts,insts,pgroups);
    pfit(i) = sum(out==groups)/length(groups);
end
pv = find(overfit<=sort(pfit),1)/length(pfit); if isempty(pv), pv = 1/length(pfit); end

%% Compute actual performance (7-fold CV)
fold  = 1;
aps = nchoosek(1:length(groups),fold);
for i = 1:length(aps)
    testidx = aps(i);
    trainidx = 1:length(aps);
    trainidx = trainidx(trainidx~=testidx);
%     cvidx = randperm(length(groups));
%     testidx = cvidx(1:fold);
%     trainidx = cvidx(fold+1:end);
    testdata = insts(testidx,:);
    grouptest = groups(testidx,:);
    traindata = insts(trainidx,:);
    groupdata = groups(trainidx,:);
    out = classify(testdata,traindata,groupdata);
    cvfit(i) = sum(out==grouptest)/length(grouptest);
end





