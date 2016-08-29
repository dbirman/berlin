function savedata_berlin(cfolder,version)

% Generic function to save data from berlin runs

% (1) Extract the left and right flat map voxel timecourses from the
% concatenation timeseries
% (2) Generate the design by loading the contrast, coherence, and timing
% responses, using the basecon/basecoh of the run (default timing=5)
% This is the ROI_data structure and gets saved in ~/data/cohcon_localizer

%% Move to Data WD
mrQuit
cd(fullfile('~/data/berlin_experiment/',cfolder));
folders = dir(pwd);
skip = 1;
for fi = 1:length(folders)
    if ~isempty(strfind(folders(fi).name,'Concatenation')), skip = 0; end
end
if skip
    disp(sprintf('Data folder %s has not been prepared for analysis',cfolder));
    return
end
%% Setup a view + Load Concatenation
view = newView();
view = viewSet(view,'curGroup','Concatenation');
view = viewSet(view,'curScan',1); % make sure scan # is correct

%% Get the mean timeseries using the reversed pRF
view = loadAnalysis(view,sprintf('erAnal/%s','all')); % check analysis name!
allROIs = {'V1','MT'};
% allROIs = {'V1','V2','V3','V4','V3a','V3b','V7','LO1','LO2','MT'};
% ROIs = {'V1'};
% % % pfxs = {'l','r'};
% % % allROIs = {};
% % % for ri = 1:length(ROIs)
% % %     for pi = 1:length(pfxs)
% % %         allROIs{end+1} = strcat(pfxs{pi},ROIs{ri});
% % %     end
% % % end
if ~isfield(view,'analyses') || isempty(view.analyses)
    disp(sprintf('Data folder %s has not been prepared for analysis',cfolder));
    return
end    
analysis = view.analyses{1};
rois = loadROITSeries(view,allROIs,view.curScan,view.curGroup,'keepNAN=true');
if isempty(rois)
    disp(sprintf('Data folder %s has not been prepared for analysis',cfolder));
    return
end    

%% Run leaveOneOut

% [stimvol, conds,~] = getStimvol(view,'dir2','taskNum=1','phaseNum=2','segmentNum=5');
[stimvol,~,~] = getStimvol(view,'dir1','taskNum=1','phaseNum=2','segmentNum=2');
inst = getInstances(view,rois,stimvol,'startLag=8','blockLen=6');
inst_ = leaveOneOut(inst,'balancByRemovI');

%% perm
nb = 10;
bs = 8;
v1 = zeros(1,nb*bs); mt = zeros(1,nb*bs);
disppercent(-1/nb,'Permuting!');
for j=1:nb
    parfor i = (j-1)*bs+1:(j-1)*bs+bs
        inst_p = leaveOneOut(inst,'permutationBal=1');
        v1(i) = inst_p{1}.classify.leaveOneOut.correct;
        mt(i) = inst_p{2}.classify.leaveOneOut.correct;
    end
    disppercent(j/nb,'Block finished');
end
disppercent(inf);

%%
tv1 = inst_{1}.classify.leaveOneOut.correct;
tmt = inst_{2}.classify.leaveOneOut.correct;
pv1 = 1-find(tv1<=sort(v1),1)/length(v1);
pmt = 1-find(tmt<=sort(mt),1)/length(mt);

%% 

tSeries = cell(1,length(rois)); % meaned tSeries with prf
rtSeries = cell(1,length(rois));
r2s = cell(1,length(rois));
rights = cell(1,length(rois));
lefts = cell(1,length(rois));
corrs = cell(1,length(rois));

scanDims = viewGet(view,'scanDims');

left = []; right = []; r2 = [];
% check which overlay is which to make sure they get sorted properly
for i = 1:length(analysis.overlays)
    cOverlay = analysis.overlays(i);
    if strfind(cOverlay.name,'r2')
        r2 = cOverlay.data{1};
        disp(sprintf('(sd_loc) Setting overlay %i to r2',i));
    elseif strfind(cOverlay.name,'left')
        left = cOverlay.data{1};
        disp(sprintf('(sd_loc) Setting overlay %i to left',i));
    elseif strfind(cOverlay.name,'right')
        right = cOverlay.data{1};
        disp(sprintf('(sd_loc) Setting overlay %i to right',i));
    else
        warning('failure');
        keyboard
    end
end

if isempty(left) || isempty(right)
    disp('(sd_loc) Failed to find the overlays, averaging across entire ROI');
end

% just incase we want this at some point?
rois = getSortIndex(view,rois,r2);

for ri = 1:length(rois)
    r = rois{ri};
    r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));
    
    cr2 = r2(r.linearScanCoords); 
    
    rtSeries{ri} = r.tSeries(r.sortindex,:);
    
    %%
%     if ~isempty(left) && ~isempty(right) 
%         rightO = right(r.linearScanCoords);
%         leftO = left(r.linearScanCoords);
%         rightO(isnan(rightO))=0;
%         leftO(isnan(leftO))=0;
% 
%         idxs = ~any(isnan(r.tSeries),2);
%         if any(~idxs)
%             warning('Failure');
%             keyboard
%         end
%     %     tSeriesnoNaN = r.tSeries;
%     %     tSeriesnoNaN(isnan(tSeriesnoNaN)) = 0;
% 
%         r2s{ri} = cr2;
%         rights{ri} = rightO;
%         lefts{ri} = leftO;
% 
%         if strcmp(r.name(1),'l')
%             tSeries{ri} = (rightO*r.tSeries)/sum(rightO);
%             corrs{ri} = corr(cr2',rightO');
%         else
%             tSeries{ri} = (leftO*r.tSeries)/sum(leftO);
%             corrs{ri} = corr(cr2',leftO');
%         end
%     else
%         tSeries{ri} = mean(r.tSeries);
%         corrs{ri} = 0;
%         rights{ri} = [];
%         lefts{ri} = [];
%         r2s{ri} = cr2;
%     end
% % % % % % %     r2cutoff = 0.5;
% % % % % % %     while sum(cr2>r2cutoff)<25
% % % % % % %         r2cutoff = r2cutoff-0.001;
% % % % % % %     end
% % % % % % %     rtSeries{ri} = r.tSeries(cr2>r2cutoff,:);
end

%% testing
% r = rois{19};    r.linearScanCoords = sub2ind(scanDims,r.scanCoords(1,:),r.scanCoords(2,:),r.scanCoords(3,:));
% cr2 = r2(r.linearScanCoords);
% clf, hold on
% r2cutoff = 0.5;
% while sum(cr2>r2cutoff)<25
%     r2cutoff = r2cutoff-0.001;
% end
% point5 = mean(r.tSeries(cr2>r2cutoff,:));
% rightO = left(r.linearScanCoords);
% rightO = rightO;
% rightO(isnan(rightO))=0;
% rightf = (rightO*r.tSeries)/sum(rightO);
% plot(point5(1000:2200));
% plot(rightf(1000:2200),'r');
% a = axis;
% axis([a(1) a(2) .95 1.05]);


%% Pull Stimvols
pull = struct;
[pull.dir1.vol, pull.dir1.conds,~] = getStimvol(view,'dir1','taskNum=1','phaseNum=2','segmentNum=2');
[pull.dir2.vol, pull.dir2.conds,~] = getStimvol(view,'dir2','taskNum=1','phaseNum=2','segmentNum=2');
[pull.mask.vol, pull.mask.conds,~] = getStimvol(view,'contrast','taskNum=1','phaseNum=2','segmentNum=2');
[pull.correct.vol, pull.correct.conds, ~] = getStimvol(view,'correct','taskNum=1','phaseNum=2','segmentNum=2');
[pull.match.vol,pull.match.conds,~] = getStimvol(view,'match','taskNum=1','phaseNum=2','segmentNum=2');
pull.dir1.str = 'dir1'; pull.dir2.str = 'dir2'; pull.mask.str = 'contrast';
pull.correct.str = 'correct'; pull.match.str = 'match';
%% Parse stimvols
pull = parseConds(pull);

%% Convert to Long Form
pfields = fields(pull);
allsv = [];
for pi = 1:length(pfields)
    pull.(pfields{pi}).data = volval2long(pull.(pfields{pi}).vol,pull.(pfields{pi}).val);
    allsv = [allsv pull.(pfields{pi}).data(:,1)'];
end
allsv = unique(allsv);

%% Generate design
design = zeros(10000,6);
count = 1;

for ci = 1:length(allsv)
    csv = allsv(ci);
    
    % if we find the stimvol in norespvols (i.e. resp was NaN) we will just
    % ignore this trial. This means that we won't model the response at all
    % (this isn't ideal as a solution... better might be to blank out the data for the duration of the trial
    % but it's reasonable since these trials are rare).
        
    dat = struct;
    for pi = 1:length(pfields)
        local = pull.(pfields{pi});
        values = local.data(:,2);
        svs = local.data(:,1);
        dat.(pfields{pi}) = values(find(svs==csv,1));
    end

    if ~isempty(dat.correct)
        design(count,:) = [csv dat.dir1 dat.dir2 dat.mask dat.correct dat.match];
    % sv dir1 dir2 mask correct match
    count = count+1;
    end
end
design = design(1:count-1,:);

if size(design,1)~=length(allsv)
    warning('%i trials were dropped due to no responses',length(allsv)-size(design,1));
end

%% Save data
savefolder = '~/data/berlin_experiment';
fname = sprintf('%s_data.mat',cfolder);
fname = fullfile(savefolder,fname);

% data.tSeries = tSeries;
data.rtSeries = rtSeries;
data.ROIs = allROIs;
% data.roi.r2 = r2s;
% data.roi.right = rights;
% data.roi.left = lefts;
% data.roi.corrs = corrs;
data.design = design;
concatInfo = viewGet(view,'concatInfo');
data.runtrans = concatInfo.runTransition;
data.concatInfo = concatInfo;
data.TR = 0.5;
data.version = version;

save(fname,'data');

clear data

disp(sprintf('Saving data for %s completed successfully',cfolder));