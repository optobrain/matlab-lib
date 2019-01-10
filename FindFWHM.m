%% find FWHM of a convex data in terms of point number

% bYminZero = false : min of Y can be considered to be zero?

function ret = FindFWHM(y,bYminZero)

if nargin < 2
    bYminZero = false;
end

    if bYminZero
        yh = max(y)/2;
    else
        yh = mean([max(y) min(y)]);
    end
    y1 = abs(diff(sign(y-yh)));
    ret = abs(diff(find(y1==2)));
    