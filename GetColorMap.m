%% Get [n 3] array for colormap
% mode: 'doppler' 

function ret = GetColorMap(mode, n)

if nargin < 2
    n = 64;
end
if nargin < 1
    mode = 'doppler';
end

    ret = zeros(n,3);
    if strcmp(mode,'doppler')
        % bright blue
        ret((n/8+1):n/2,3) = linspace(1,0,n/2-n/8);
        ret(1:n/8,3) = 1;
        ret(1:n/8,1:2) = repmat(linspace(0.5,0,n/8)',[1 2]);
        % bright red
        ret([(n/2+1):(n-n/8)],1) = linspace(0,1,n/2-n/8);
        ret([(n-n/8+1):n],1) = 1;
        ret([(n-n/8+1):n],2:3) = repmat(linspace(0,0.5,n/8)',[1 2]);
    elseif strcmp(mode,'doppler_white')
        ret = ones(n,3);
        ret(1:n/2,1:2) = repmat(linspace(0,1,n/2)',[1 2]);
        ret(n/2+1:end,2:3) = repmat(linspace(1,0,n/2)',[1 2]);
    end
    