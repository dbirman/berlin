function fit = fitmvn(data)
% FITMVN - Takes x/y fixation coordinates and fits the optimal 2-d gaussian
%
%   Input
%   x - X coordinates of fixations
%   y - Y coordinates of fixations
%
%   Output
%   fit - Fit optimal s0/s1, likelihood, BIC, etc

params = fminsearch(@resMVN,[1 1],[],data); % start at the null hypothesis

fit.params = params;
fit.likelihood = resMVN(params,data);
fit.BIC = -2*log(fit.likelihood) + log(size(data,1));

function res = resMVN(params,data)

likelihood = mvnpdf(data,[0 0],params);
res = -sum(log(likelihood));