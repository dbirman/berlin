%% Power Analysis Script

% we're going to simulate the DMS experiment and try to figure out how much
% data is necessary to get the effects we're interested in.

% we will assume that fixations are sampled from a 2-d gaussian using
% mvnrnd, with pdf mvnpdf.

% note that sigma has to be diagonal, which can be passed in just as [s0
% s1]

%% SET PARAMS

% NULL PARAMS
x0 = 0;
x1 = 0;

s0 = 1;
s1 = 1;

mu0 = [x0 x1];
sigma0 = [s0 s1];
% TEST PARAMS
x0 = 0;
x1 = 0; % under the assumption that the center of fixation is always constant

s0 = 1;
s1 = 2;
% effect size calculation?

mu1 = [x0 x1];
sigma1 = [s0 s1];
%% Display true fixation density
rng = -3:.1:3;
[X,Y] = meshgrid(rng,rng);
subplot(121)
map = mvnpdf([X(:) Y(:)],mu0,sigma0);
map = reshape(map,size(X));
imagesc(rng,rng,map);
colormap(gray(200))
axis square
subplot(122)
map = mvnpdf([X(:) Y(:)],mu1,sigma1);
map = reshape(map,size(X));
imagesc(rng,rng,map);
colormap(gray(200))
axis square

%% For X samples from Y subjects, calculate the optimal s0/s1 values, then
% bootstrap across subjects to determine probability of getting true effect
samplesizes = [200 400];
subjectsizes = [5 10];
repeats = 20;
runcount = length(samplesizes)*length(subjectsizes);

disppercent(-1/repeats,'LETS DO THIS');
% Save the values obtained for comparisons (p-values!)
data = cell(repeats,length(samplesizes),length(subjectsizes));
for ri = 1:repeats
    for sai = 1:length(samplesizes)
        samples = samplesizes(sai);
        for sui = 1:length(subjectsizes)
            subjects = subjectsizes(sui);
            count = (sai-1)*length(subjectsizes)+sui;
            
            % save data form this round
            data{ri,sai,sui}.fit0 = zeros(subjects,2);
            data{ri,sai,sui}.fit1 = zeros(subjects,2);
            data{ri,sai,sui}.boot0 = zeros(subjects,2,2);
            data{ri,sai,sui}.boot1 = zeros(subjects,2,2);
            for si = 1:subjects
                % get the samples
                samples0 = mvnrnd(mu0,sigma0,samples);
                samples1 = mvnrnd(mu1,sigma1,samples);
                % bootstrap each distributions mean values
                for val = 1:2
                    data{ri,sai,sui}.boot0(si,val,:) = bootci(1000,@mean,samples0(:,val));
                    data{ri,sai,sui}.boot1(si,val,:) = bootci(1000,@mean,samples1(:,val));
                end
                % fit an MVN to each set of samples
                fit0 = fitmvn(samples0);
                data{ri,sai,sui}.fit0(si,:) = fit0.params;
                fit1 = fitmvn(samples1);
                data{ri,sai,sui}.fit1(si,:) = fit1.params;
            end
        end
    end
    disppercent(ri/repeats);
end
disppercent(inf);

%% Analysis time--reorganize to get all repeats
for sai = 1:length(samplesizes)
    samples = samplesizes(sai);
    for sui = 1:length(subjectsizes)
        subjects = subjectsizes(sui);
        % collect fit0 data
        fit0 = zeros(2,subjects,repeats);
        fit1 = zeros(2,subjects,repeats);
        
        for lj = 1:subjects
            for lk = 1:repeats
                fit0(:,lj,lk) = data{ri,sai,sui}.fit0(lj,:);
                fit1(:,lj,lk) = data{ri,sai,sui}.fit1(lj,:);
            end
        end
        % compute bootstraps
        fit0_ci = zeros(2,repeats,2);
        fit1_ci = zeros(2,repeats,2);
        for ri = 1:repeats
            dat = squeeze(fit0(:,:,repeats));
            fit0_ci(1,ri,:) = bootci(500,@mean,dat(1,:));
           	fit0_ci(2,ri,:) = bootci(500,@mean,dat(2,:));
            dat = squeeze(fit1(:,:,repeats));
            fit1_ci(1,ri,:) = bootci(500,@mean,dat(1,:));
           	fit1_ci(2,ri,:) = bootci(500,@mean,dat(2,:));
        end
        % compute statistics
    end
end