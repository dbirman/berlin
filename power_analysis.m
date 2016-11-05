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
samplesizes = 100;
subjectsizes = 3;
runcount = length(samplesizes)*length(subjectsizes);

% Save the values obtained for comparisons (p-values!) 

for sai = 1:length(samplesizes)
    samples = samplesizes(sai);
    for sui = 1:length(subjectsizes)
        subjects = subjectsizes(sui);
        count = (sai-1)*length(subjectsizes)+sui;
        
        %
        boots = zeros(
        for si = 1:subjects
            % get the samples
            samples0 = mvnrnd(mu0,sigma0,samples);
            samples1 = mvnrnd(mu1,sigma1,samples);
            % bootstrap each distributions mean values
            for val = 1:2
                boot_ci(1,val,:,count) = bootci(1000,@mean,samples0(:,val));
                boot_ci(2,val,:,count) = bootci(1000,@mean,samples1(:,val));
            end
            % fit an MVN to each set of samples
            fit0 = fitmvn(samples0);
            fit1 = fitmvn(samples1);
        end
    end
end