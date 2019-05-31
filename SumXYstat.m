
% see OneNote > Study > Statistics > Z = aX+bY
% covr = covariance 

function [mu, std] = SumXYstat(mu1, std1, mu2, std2, n, covr, a, b)
if nargin < 8
    b = 1;
end
if nargin < 7
    a = 1;
end

    mu = a * mu1 + b * mu2;
    
    se1 = std1/sqrt(n);
    se2 = std2/sqrt(n);
    
    se = sqrt( a^2*se1^2 + b^2*se2^2 + 2*a*b*covr );
    std = se.*sqrt(n);
    
