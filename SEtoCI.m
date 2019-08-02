function ci = SEtoCI(se, n, bSTD, conf)

if nargin < 4
    conf = 0.95;
end
if nargin < 3
    bSTD = false;
end

if bSTD
    se = se / sqrt(n);
end

ci = tinv(conf,n-1)*se;