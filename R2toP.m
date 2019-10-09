% Obtain p-value from R2 of a fitting against the constant-mean model
% R2 can be array
% K = number of parameters (e.g., K=2 for Y~1+X)
% F = F-statistics 
% See github/study/Statistics/21 Continuous - Linear regression.mlx

function [p, F] = R2toP(R2, n, K)

F = R2./(1-R2)*(n-K)/(K-1);
p = 1 - fcdf(F,K-1,n-K);




